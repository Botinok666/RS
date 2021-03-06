#include <intrin.h>
#include "rsssse3.h"

int GetSupportedExtensions()
//No: 0; SSSE3: 1; AVX2: 3
{
    int info[4], nIds, result = 0;
    __cpuidex(info, 0, 0);
    nIds = info[0];
    if (nIds >= 1)
    {
        __cpuidex(info, 1, 0);
        if ((info[2] & ((int)1 << 9)) != 0) //SSSE3
            result |= 1;
    }
    if (nIds >= 7)
    {
        __cpuidex(info, 7, 0);
        if ((info[1] & ((int)1 << 5)) != 0) //AVX2
            result |= 2;
    }
    return result;
}

void InitSSSE3(uint8_t* coefsu, const uint8_t count, uint8_t* lut)
{
    if (count < 2 || count > MAX_T * 2) return;
	uint64_t offset = (uint64_t)lut & 0x1f;
	if (offset) //Array is not aligned to 32 bytes boundary
	{
		offset = 0x20 - offset;
		lut[0] = (uint8_t)offset; //Store offset
		lut += offset; //Offset LUT pointer
	}

    uint8_t* lutSSE = lut + SSE_LUT_SSE_OFFSET;
    for (int i = 0; i < 256; i++)
    {
        //Fill SSE "multiply vector to scalar" table
        uint8_t* usse = lutSSE + i * 32LL, *lsse = usse + 16;
        for (uint8_t j = 0; j < 16; j++)
        {
            *usse++ = GFMul(j << 4, (uint8_t)i);
            *lsse++ = GFMul(j, (uint8_t)i);
        }
    }
    uint16_t* lutLog = (uint16_t*)lut;
    uint8_t* lutExp = lut + ALU_LUT_EXP_OFFSET;
    uint8_t x = 1, index = 0;
    lutExp[0] = x; //0:255 - log, 256:767 - exp, 768:1792 - zero
    for (int i = 0; i < 255; i++) {
        uint8_t y = GFMul(x, 2);
        lutExp[++index] = y;
        lutExp[index + 255] = y;
        x = y;
    }
    lutLog[0] = 511; //Log(0) = inf
    for (int i = 0; i < 255; i++)
        lutLog[lutExp[i]] = i;
    memset(lutExp + 511, 0, 1024);
    
    uint8_t restmp[256];
    memset(restmp, 1, 256); //Set result initially to 1
    uint8_t* lutRoots = lut + SSE_LUT_ROOTS_OFFSET;
    for (int j = 0; j < 255; j++)
    {
        for (int i = 0; i < 256; i++)
            restmp[i] = (uint8_t)lutExp[lutLog[restmp[i]] + i];
        memcpy_s(lutRoots, 256, restmp, 256);
        lutRoots += 256;
    }
    memset(lutRoots, 0, 32); //Clear next 32 bytes so during unaligned load there will be no random data

	//We need to align array at 32 bytes boundary
	offset = 0x20 - ((uint64_t)coefsu & 0x1f);
	coefsu[0] = (uint8_t)offset; //Save offset
	uint8_t* coefs = coefsu + offset - 1; //coefs[0] will be placed outside of aligned array
    memset(coefs, 0, MAX_T * 2 + 1);
	
    coefs[count] = 1;
    coefs[count - 1] = 1;
    for (int k = 2; k <= count; k++)
    {
        coefs[count - k] = 1;
        uint16_t s;
        for (int i = count - k + 2; i <= count; i++)
        {
            s = lutLog[coefs[i - 1]] + k - 1;
            coefs[i - 1] = coefs[i] ^ lutExp[s];
        }
        s = lutLog[coefs[count]] + k - 1;
        coefs[count] = (uint8_t)lutExp[s];
    }
}

static void inline GF_mvs_SSSE3(__m128i* vtt, __m128i* vio, __m128i* vmask, __m128i* lmul, __m128i* umul)
{
    *vtt = _mm_and_si128(*vio, *vmask); //Lower 4 bits
    *vtt = _mm_shuffle_epi8(*lmul, *vtt); //Multiply lower 4 bits
    *vio = _mm_srli_epi64(*vio, 4);
    *vio = _mm_and_si128(*vio, *vmask); //Upper 4 bits
    *vio = _mm_shuffle_epi8(*umul, *vio); //Multiply upper 4 bits
    *vio = _mm_xor_si128(*vio, *vtt); //Result of multiplication
}

int EncodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer)
{
    int length = (n - k - 1) >> 4;
	if (lut[0] != 0xff) //LUT table requires alignment
		lut += lut[0];
    if ((uint64_t)lut & 0xf)
        return -1; //LUT still misaligned? Error must be thrown on upper level
	coefs += coefs[0]; //Align coefficients table
    if ((uint64_t)coefs & 0xf)
        return -1; //Coefficients still misaligned? Error must be thrown on upper level
	
    uint8_t btemp[255 + 16];
    uint8_t* lutSSE = lut + SSE_LUT_SSE_OFFSET, * oj = btemp;
    memcpy_s(btemp, k, buffer, k);
    memset(btemp + k, 0, n - k);

    __m128i vmask = _mm_set1_epi8(0xf);
    //Remainder calculation: long division
    for (int j = 0; j < k; j++)
    {
        int idx = *oj++ * 32;
        __m128i* lutv = (__m128i*)(lutSSE + idx);
        __m128i umul = _mm_load_si128(lutv);
        __m128i lmul = _mm_load_si128(lutv + 1);

        __m128i* ol = (__m128i*)oj;
        __m128i* cl = (__m128i*)coefs;
        for (int l = 0; l <= length; l++)
        {
            __m128i vio = _mm_load_si128(cl++), vtt;
            GF_mvs_SSSE3(&vtt, &vio, &vmask, &lmul, &umul);

            vtt = _mm_loadu_si128(ol);
            vtt = _mm_xor_si128(vtt, vio);
            _mm_storeu_si128(ol++, vtt);
        }
    }
    memcpy_s(buffer + k, n - k, btemp + k, n - k);
    return 0;
}
int DecodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    int scount = n - k, t = scount >> 1, steps = (scount - 1) >> 4;
	if (lut[0] != 0xff) //LUT table requires alignment
		lut += lut[0];
    if ((uint64_t)lut & 0xf)
        return -1; //LUT still misaligned? Error must be thrown on upper level

    uint16_t* lutLog = (uint16_t*)lut;
    uint8_t* lutExp = lut + ALU_LUT_EXP_OFFSET;
    uint8_t* lutSSE = lut + SSE_LUT_SSE_OFFSET;
    __declspec(align(16)) uint8_t lambda[2 * MAX_T], syn[2 * MAX_T];
	__declspec(align(16)) uint8_t b[2 * MAX_T], Lm[2 * MAX_T];

    /*Syndrome calculation*/
	__m128i vz = _mm_setzero_si128();
    memset(syn, buffer[n - 1], MAX_T * 2);
    uint8_t* roots = lut + SSE_LUT_ROOTS_OFFSET;
    __m128i* ls = (__m128i*)syn;
    __m128i vmask = _mm_set1_epi8(0xf);
    for (int j = n - 2; j >= 0; j--)
    {
        int idx = buffer[j] * 32;
        __m128i* lutv = (__m128i*)(lutSSE + idx);
        __m128i umul = _mm_load_si128(lutv);
        __m128i lmul = _mm_load_si128(lutv + 1);

        __m128i* lr = (__m128i*)roots;
        for (int m = 0; m <= steps; m++)
        {
            __m128i vio = _mm_load_si128(lr++), vt;
            //Multiply: {roots}^(n - j) * buffer[j]
            GF_mvs_SSSE3(&vt, &vio, &vmask, &lmul, &umul);
            vt = _mm_load_si128(ls);
            vt = _mm_xor_si128(vt, vio);
            _mm_store_si128(ls++, vt);
        }
        ls = (__m128i*)syn;
        roots += 256; //Hard to explain
    }
    //Check syndromes we've got
    uint16_t syn2[2 * MAX_T];
    uint16_t hasErrors = 0;
    for (int j = 0; j < scount; j++)
    {
        hasErrors |= syn[j];
        syn2[j] = lutLog[syn[j]];
    }
    if (!hasErrors) return 0;

    //Lambda calculation: Berlekamp's method
    //Set initial value of b and lambda to 1, and shift b left
	__m128i* bvec = (__m128i*)b;
	__m128i* lvec = (__m128i*)lambda;
	for (int m = 0; m <= steps; m++)
	{
		_mm_store_si128(bvec++, vz);
		_mm_store_si128(lvec++, vz);
	}
    b[1] = 1;
    lambda[0] = 1;
    int l = 0;
    for (int r = 1; r <= scount; r++)
    {
        uint16_t* si = syn2 + r - 1;
        uint8_t *li = lambda;
        uint8_t delta = lutExp[*si--];
        for (int m = 0; m < l; m++)
        {
            uint16_t lm = lutLog[*++li] + *si--;
            delta ^= lutExp[lm];
        }

        int sx = r <= t ? r - 1 : t;
        int sr = sx >> 4; //Optimization for less copy operations
        if (delta)
        {
			int ri = delta;
			__m128i* lutv = (__m128i*)(lutSSE + ri * 32LL);
			__m128i umul = _mm_load_si128(lutv);
			__m128i lmul = _mm_load_si128(lutv + 1);
			
			bvec = (__m128i*)b;
			lvec = (__m128i*)lambda;
			__m128i* mvec = (__m128i*)Lm;
            for (int m = 0; m <= sr; m++)
            {
				__m128i vx = _mm_load_si128(bvec++), vxl;
				//Multiply: {b} * delta
                GF_mvs_SSSE3(&vxl, &vx, &vmask, &lmul, &umul);

				vxl = _mm_load_si128(lvec);
				_mm_store_si128(mvec++, vxl); //Copy lambda to Lm
				vxl = _mm_xor_si128(vxl, vx);
				_mm_store_si128(lvec++, vxl); //Save result to lambda
            }
            if (2 * l < r)
            {
                l = r - l;
				ri = 255 - lutLog[ri];
				ri = lutExp[ri]; //delta^-1
                lutv = (__m128i*)(lutSSE + ri * 32LL);
				umul = _mm_load_si128(lutv);
				lmul = _mm_load_si128(lutv + 1);
				
				bvec = (__m128i*)b;
				mvec = (__m128i*)Lm;
				for (int m = 0; m <= sr; m++)
				{
					__m128i vx = _mm_load_si128(mvec++), vt;
					//Multiply: {Lm} * delta^-1
                    GF_mvs_SSSE3(&vt, &vx, &vmask, &lmul, &umul);

					_mm_store_si128(bvec++, vx); //Save result to b
				}
            }
        }
		//Shift b left
        sr = (sx + 1) >> 4;
		bvec = (__m128i*)b + sr;
		__m128i vx = _mm_load_si128(bvec);
		vx = _mm_slli_si128(vx, 1);
		_mm_store_si128(bvec--, vx);
		if (!sr) continue;
		uint8_t* bi = (uint8_t*)bvec + 1;
		for (int m = 0; m < sr; m++)
		{
			vx = _mm_load_si128(bvec--);
			_mm_storeu_si128((__m128i*)bi, vx);
			bi -= 16;
		}
		b[0] = 0;
    }
    int nerr = 0;
    for (int m = t; m > 0; m--) //Find deg(lambda)
    {
        if (lambda[m])
		{
            nerr = m;
			break;
		}
    }
    //if (nerr > (scount >> 1))
    //    return -2; //deg(lambda) > t? Uncorrectable error pattern occured

    uint16_t omega[2 * MAX_T];
    /* Omega calculation */
    //Omega must be calculated only up to {nerr} power
    for (int m = 0; m <= nerr; m++)
    {
        uint8_t og = syn[m];
        for (int l = 1; l <= m; l++)
        {
            uint16_t li = lutLog[lambda[l]] + syn2[m - l];
            og ^= lutExp[li];
        }
        omega[m] = lutLog[og];
    }

    /* Chien search
    We need to substitute x[] to lambda(x), where x[] = x^-n+1, x^-n+2, ..., 1 = x^256-n, x^256-n+1, ..., 1 */
    steps = (n - 1) >> 4;
    uint8_t efound = 0, ecorr = 0;

    __declspec(align(16)) uint8_t rtemp[256];
    memset(rtemp, 1, n);
    roots = lut + SSE_LUT_ROOTS_OFFSET + 256 - n;
    for (int j = 0; j < nerr; j++)
    {
        int idx = lambda[j + 1] * 32;
        __m128i* lutv = (__m128i*)(lutSSE + idx);
        __m128i umul = _mm_load_si128(lutv);
        __m128i lmul = _mm_load_si128(lutv + 1);

        __m128i* lr = (__m128i*)roots;
        ls = (__m128i*)rtemp;
        for (int m = 0; m <= steps; m++)
        {
            __m128i vio = _mm_loadu_si128(lr++), vt;
            //Multiply: {roots}^(n - j) * lambda[j + 1]
            GF_mvs_SSSE3(&vt, &vio, &vmask, &lmul, &umul);
            vt = _mm_load_si128(ls);
            vt = _mm_xor_si128(vt, vio);
            _mm_store_si128(ls++, vt);
        }
        roots += 256; //Hard to explain
    }
    //Check what we've calculated
    for (int j = 0; j < n; j++)
    {
        if (rtemp[j]) continue;
        efound++;
        if (j < k) 
        {
            int xIdx = 256 - n + j;
            uint8_t s = 0, y = lutExp[omega[0]];
            for (int l = 1; l <= nerr; l++)
            {
                int ecx = xIdx * l;
                mod255(&ecx);
                y ^= lutExp[ecx + omega[l]]; //Omega(X^-1)
                if (l & 1)
                    s ^= lutExp[ecx + lutLog[lambda[l]]]; //Lambda'(X^-1) * X^-1
            }
            if (!s) return -2;

            s = 255 - (uint8_t)lutLog[s];
            s = lutExp[s + lutLog[y]];
            b[ecorr] = (uint8_t)j;
            Lm[ecorr++] = s;
        }
    }
    
    if (efound != nerr)
        return -2;
    for (int j = 0; j < ecorr; j++)
        buffer[b[j]] ^= Lm[j];

    return nerr;
}
