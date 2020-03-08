#include <intrin.h>
#include "rsssse3.h"
#define ERROR_CHECKING
// Slow multiply, using shifting
// Polynomial x^8 + x^4 + x^3 + x^2 + 1
uint8_t GFMul3(uint8_t a, uint8_t b)
{
    uint8_t r = 0, t;
    while (a)
    {
        if (a & 1)
            r ^= b;
        t = b & 0x80;
        b <<= 1;
        if (t)
            b ^= 0x1d;
        a >>= 1;
    }
    return r;
}

static inline void GFMulSSE2(__m128i *va, __m128i *vb, __m128i *vz, __m128i *vs, __m128i *vgp)
{
    __m128i vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
    *va = _mm_add_epi8(*va, *va);

    vx = _mm_cmpgt_epi8(*vz, *va);
    vx = _mm_and_si128(vx, *vb);
    *vs = _mm_xor_si128(*vs, vx);
    vx = _mm_cmpgt_epi8(*vz, *vb);
    vx = _mm_and_si128(vx, *vgp);
    *vb = _mm_add_epi8(*vb, *vb);
    *vb = _mm_xor_si128(*vb, vx);
}

int GetLatestSupportedExtension()
//No: 0; SSSE3: 1; AVX2: 2
{
    int info[4], nIds;
    __cpuidex(info, 0, 0);
    nIds = info[0];
    if (nIds >= 7)
    {
        __cpuidex(info, 7, 0);
        if ((info[1] & ((int)1 << 5)) != 0) //AVX2
            return 2;
    }
    else if (nIds >= 1)
    {
        __cpuidex(info, 1, 0);
        if ((info[2] & ((int)1 << 9)) != 0) //SSSE3
            return 1;
    }
    return 0;
}

void InitSSSE3(uint8_t* coefsu, const uint8_t count, uint8_t* lut)
{
    if (count < 2) return;
	uint64_t offset = (uint64_t)lut & 0xf;
	if (offset) //Array is not aligned to 16 bytes boundary
	{
		offset = 0x10 - offset;
		lut[0] = (uint8_t)offset; //Store offset
		lut += offset; //Offset LUT pointer
	}
    uint8_t* lutExp = lut + SSE_LUT_EXP_OFFSET, *lutSSE = lut + SSE_LUT_SSE_OFFSET;
    uint8_t x = 1, index = 0;
    lutExp[0] = x; //0:255 - log, 256:767 - exp, 768:8960 - sse
    for (int i = 0; i < 255; i++) {
        uint8_t y = GFMul3(x, 2);
        lutExp[++index] = y;
        lutExp[index + 255] = y;
        x = y;
    }
    lutExp[767] = 0;

    lut[0] = 255; //Log(0) = inf
    for (int i = 0; i < 256; i++)
    {
        lut[lutExp[i]] = i;
        //Fill SSE "multiply vector to scalar" table
        uint8_t* usse = lutSSE + i * 32LL, *lsse = usse + 16;
        for (uint8_t j = 0; j < 16; j++)
        {
            *usse++ = GFMul3(j << 4, (uint8_t)i);
            *lsse++ = GFMul3(j, (uint8_t)i);
        }
    }
    __m128i* lutRev = (__m128i*)(lut + SSE_LUT_REV_OFFSET);
    __m128i* roots = (__m128i*)lutExp;
    __m128i vmask = _mm_set1_epi8(0xf);
    __m128i vlrev = _mm_set_epi8(15, 7, 11, 3, 13, 5, 9, 1, 14, 6, 10, 2, 12, 4, 8, 0);
    __m128i vhrev = _mm_set_epi8(-16, 112, -80, 48, -48, 80, -112, 16, -32, 96, -96, 32, -64, 64, -128, 0);
    //Fill table with roots with reverse ordered bits
    for (int i = 0; i < 256; i += 16)
    {
        __m128i a = _mm_loadu_si128(roots++);
        __m128i c = _mm_and_si128(a, vmask);
        c = _mm_shuffle_epi8(vhrev, c);
        a = _mm_srli_epi64(a, 4);
        a = _mm_and_si128(a, vmask);
        a = _mm_shuffle_epi8(vlrev, a);
        a = _mm_or_si128(a, c);
        _mm_store_si128(lutRev++, a);
    }
    _mm_store_si128(lutRev, _mm_setzero_si128()); //Clear next 16 bytes so during unaligned load there will be no random data
	
    uint8_t* lutLog = lut;
	//We need to align array at 16 bytes boundary
	offset = 0x10 - ((uint64_t)coefsu & 0xf);
	coefsu[0] = (uint8_t)offset; //Save offset
	uint8_t* coefs = coefsu + offset - 1; //coefs[0] will be placed outside of aligned array
	__m128i vz = _mm_setzero_si128();
	__m128i* cptr = (__m128i*)(coefsu + offset);
	for (int i = 0; i < MAX_S_SSE / 16; i++) //Clear coefs
		_mm_store_si128(cptr++, vz);
	
    coefs[count] = 1;
    coefs[count - 1] = 1;
    for (int k = 2; k <= count; k++)
    {
        coefs[count - k] = 1;
        uint16_t s, v = k - 1;
        for (int i = count - k + 2; i <= count; i++)
        {
            s = coefs[i - 1];
            if (!s) continue;
            s = lutLog[s] + v;
            coefs[i - 1] = coefs[i] ^ lutExp[s];
        }
        s = coefs[count];
        if (!s) continue;
        s = lutLog[s] + v;
        coefs[count] = lutExp[s];
    }
}

int EncodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer)
{
    int length = n - k;
	if (lut[0] != 0xff) //LUT table requires alignment
		lut += lut[0];
    if ((uint64_t)lut & 0xf)
        return -3; //LUT still misaligned? Error must be thrown on upper level
	coefs += coefs[0]; //Align coefficients table
    if ((uint64_t)coefs & 0xf)
        return -4; //Coefficients still misaligned? Error must be thrown on upper level
	
    length = (length - 1) >> 4;
    uint8_t btemp[255 + 16];
    uint8_t* lutExp = lut + SSE_LUT_EXP_OFFSET, * lutLog = lut, * lutSSE = lut + SSE_LUT_SSE_OFFSET;
    uint8_t* oj = btemp;
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
            __m128i ucoefs = _mm_load_si128(cl++);
            __m128i lcoefs = _mm_and_si128(ucoefs, vmask); //Lower 4 bits
            lcoefs = _mm_shuffle_epi8(lmul, lcoefs); //Multiply lower 4 bits
            ucoefs = _mm_srli_epi64(ucoefs, 4);
            ucoefs = _mm_and_si128(ucoefs, vmask); //Upper 4 bits
            ucoefs = _mm_shuffle_epi8(umul, ucoefs); //Multiply upper 4 bits
            ucoefs = _mm_xor_si128(ucoefs, lcoefs); //Result of multiplication

            lcoefs = _mm_loadu_si128(ol);
            lcoefs = _mm_xor_si128(lcoefs, ucoefs);
            _mm_storeu_si128(ol++, lcoefs);
        }
    }
    memcpy_s(buffer + k, n - k, btemp + k, n - k);
    return 0;
}
int DecodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    int scount = n - k, steps = (scount - 1) >> 4;
	if (lut[0] != 0xff) //LUT table requires alignment
		lut += lut[0];
    if ((uint64_t)lut & 0xf)
        return -3; //LUT still misaligned? Error must be thrown on upper level

    uint8_t* lutLog = lut, * lutExp = lut + SSE_LUT_EXP_OFFSET, * lutSSE = lut + SSE_LUT_SSE_OFFSET;
    __declspec(align(16)) uint8_t lambda[MAX_S_SSE], omega[MAX_S_SSE], syn[MAX_S_SSE];
	__declspec(align(16)) uint8_t b[MAX_S_SSE], Lm[MAX_S_SSE];

    /*Syndrome calculation*/
    int hasNoErrors = 0xffff;
	__m128i vz = _mm_setzero_si128();
	_mm_store_si128((__m128i*)lambda, vz);
    memset(lambda, 0xff, scount & 0xf); //Temporal storage for root mask
    __m128i vgp = _mm_set1_epi8(0x1d); //Generator polynomial
    __m128i* lrev = (__m128i*)(lut + SSE_LUT_REV_OFFSET);
    __m128i* vsyn = (__m128i*)syn;
    for (int j = 0; j <= steps; j++)
    {
        __m128i vroot = _mm_load_si128(lrev++);
        __m128i vs = _mm_setzero_si128();
        uint8_t* inm = buffer;
        for (int m = 0; m < n - 1; m++)
        {
            __m128i va = vroot;
            __m128i vb = _mm_xor_si128(vs, _mm_set1_epi8(*inm++)); //s ^= in[m]
            //Perform full multiplication: s *= root, or s = b * root
            vs = vz; //Clear s
            GFMulSSE2(&va, &vb, &vz, &vs, &vgp);
        }
        vs = _mm_xor_si128(vs, _mm_set1_epi8(*inm));
        if (j == steps && lambda[0]) //Non-zero lambda[0] means that scount % 0xf != 0
            vs = _mm_and_si128(vs, _mm_load_si128((__m128i*)lambda));
        hasNoErrors &= _mm_movemask_epi8(_mm_cmpeq_epi8(vs, vz));
        _mm_store_si128(vsyn++, vs);
    }
    if (hasNoErrors == 0xffff) return 0;

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
    __m128i vmask = _mm_set1_epi8(0xf);
    for (int r = 1; r <= scount; r++)
    {
        uint8_t* si = syn + r - 1, *li = lambda;
        uint8_t delta = *si;
        for (int m = 0; m < l; m++)
        {
            li++;
            si--;
			if (*li && *si)
			{
				uint16_t idx = lutLog[*li];
				idx += lutLog[*si];
				delta ^= lutExp[idx];
			}
        }
        uint8_t sr = (r > 17) ? ((r - 2) >> 4) : 0; //Optimization for less copy operations
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
				__m128i vx = _mm_load_si128(bvec++);
				//Multiply vector b by scalar delta
				__m128i vxl = _mm_and_si128(vx, vmask); //Lower 4 bits
				vxl = _mm_shuffle_epi8(lmul, vxl); //Multiply lower 4 bits
				vx = _mm_srli_epi64(vx, 4);
				vx = _mm_and_si128(vx, vmask); //Upper 4 bits
				vx = _mm_shuffle_epi8(umul, vx); //Multiply upper 4 bits
				vx = _mm_xor_si128(vx, vxl); //Result of multiplication

				vxl = _mm_load_si128(lvec);
				_mm_store_si128(mvec++, vxl); //Copy lambda to Lm
				vxl = _mm_xor_si128(vxl, vx);
				_mm_store_si128(lvec++, vxl); //Save result to lambda
            }
            if (2 * l <= r - 1)
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
					__m128i vx = _mm_load_si128(mvec++);
					//Multiply vector Lm by scalar delta^-1
					__m128i vxl = _mm_and_si128(vx, vmask); //Lower 4 bits
					vxl = _mm_shuffle_epi8(lmul, vxl); //Multiply lower 4 bits
					vx = _mm_srli_epi64(vx, 4);
					vx = _mm_and_si128(vx, vmask); //Upper 4 bits
					vx = _mm_shuffle_epi8(umul, vx); //Multiply upper 4 bits
					vx = _mm_xor_si128(vx, vxl); //Result of multiplication

					_mm_store_si128(bvec++, vx); //Save result to b
				}
            }
        }
		//Shift b left
		bvec = (__m128i*)b + steps;
		__m128i vx = _mm_load_si128(bvec);
		vx = _mm_slli_si128(vx, 1);
		_mm_store_si128(bvec--, vx);
		if (!steps) continue;
		uint8_t* bi = (uint8_t*)bvec + 1;
		for (int m = 0; m < steps; m++)
		{
			vx = _mm_load_si128(bvec--);
			_mm_storeu_si128((__m128i*)bi, vx);
			bi -= 16;
		}
		b[0] = 0;
    }
    int nerr = 0;
    for (int m = scount - 1; m > 0; m--) //Find deg(lambda)
    {
        if (lambda[m])
		{
            nerr = m;
			break;
		}
    }
    if (nerr > (scount >> 1))
        return -2; //deg(lambda) > t? Uncorrectable error pattern occured

    /* Omega calculation */
    //Omega must be calculated only up to {nerr} power
    for (int m = 0; m <= nerr; m++)
    {
        uint8_t og = syn[m];
        uint8_t* li = lambda, * si = syn + m;
        for (int l = 0; l < m; l++)
        {
            li++;
            si--;
            if (*li && *si)
            {
                uint16_t idx = lutLog[*li];
                idx += lutLog[*si];
                og ^= lutExp[idx];
            }
        }
        omega[m] = og;
    }

    /* Chien search
    We need to substitute x[] to labmda(x), where x[] = x^-n+1, x^-n+2, ..., 1 = x^256-n, x^256-n+1, ..., 1 */
#ifdef ERROR_CHECKING
    steps = (n - 1) >> 4;
    uint8_t efound = 0, ecorr = 0;
#else
    steps = (k - 1) >> 4;
#endif // ERROR_CHECKING

    lrev = (__m128i*)(lut + SSE_LUT_REV_OFFSET + 256 - n);
    for (int j = 0; j <= steps; j++)
    {
        uint8_t* li = lambda + nerr;
        __m128i vs = _mm_set1_epi8(*li--); //s = lambda[nerr]
        __m128i vroot = _mm_loadu_si128(lrev++);
        for (int m = 0; m < nerr; m++)
        {
            __m128i va = vroot;
            __m128i vb = vs;
            //Perform full multiplication: s *= root, or s = b * root
            vs = vz;
            GFMulSSE2(&va, &vb, &vz, &vs, &vgp);
            va = _mm_set1_epi8(*li--);
            vs = _mm_xor_si128(vs, va); //s ^= lambda[nerr - m]
        }
        vs = _mm_cmpeq_epi8(vs, vz);
        int mask = _mm_movemask_epi8(vs), kv = j << 4;
        while (mask)
        {
            if (mask & 1)
            {
#ifdef ERROR_CHECKING
                if (kv < k) {
#endif // ERROR_CHECKING
                    int xIdx = 256 - n + kv;
                    uint8_t s = 0, y = omega[0];
                    for (int l = 1; l <= nerr; l++)
                    {
                        int ecx = xIdx * l;
                        ecx = (ecx >> 8) + (ecx & 0xff);
                        ecx = (ecx >> 8) + (ecx & 0xff);
                        if (omega[l])
                            y ^= lutExp[ecx + lutLog[omega[l]]]; //Omega(X^-1)
                        if ((l & 1) && lambda[l])
                            s ^= lutExp[ecx + lutLog[lambda[l]]]; //Lambda'(X^-1) * X^-1
                    }
                    if (!s) return -4;

                    s = 255 - lutLog[s];
                    s = lutExp[s + lutLog[y]];
#ifdef ERROR_CHECKING
                    b[ecorr] = (uint8_t)kv;
                    Lm[ecorr++] = s;
                }
                efound++;
#else
                buffer[kv] ^= s;
#endif // ERROR_CHECKING
            }
            mask >>= 1;
            kv++;
        }
    }
#ifdef ERROR_CHECKING
    if (efound != nerr)
        return -3;
    for (int j = 0; j < ecorr; j++)
        buffer[b[j]] ^= Lm[j];
#endif // ERROR_CHECKING

    return nerr;
}
