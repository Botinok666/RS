#include <intrin.h>
#include "rsavx2.h"

static void inline GF_mvs_AVX2(__m256i* vtt, __m256i* vio, __m256i* vmask, __m256i* lmul, __m256i* umul)
{
    *vtt = _mm256_and_si256(*vio, *vmask); //Lower 4 bits
    *vtt = _mm256_shuffle_epi8(*lmul, *vtt); //Multiply lower 4 bits
    *vio = _mm256_srli_epi64(*vio, 4);
    *vio = _mm256_and_si256(*vio, *vmask); //Upper 4 bits
    *vio = _mm256_shuffle_epi8(*umul, *vio); //Multiply upper 4 bits
    *vio = _mm256_xor_si256(*vio, *vtt); //Result of multiplication
}

int EncodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer)
{
    int length = (n - k - 1) >> 5;
    if (lut[0] != 0xff) //LUT table requires alignment
        lut += lut[0];
    if ((uint64_t)lut & 0x1f)
        return -1; //LUT still misaligned? Error must be thrown on upper level
    coefs += coefs[0]; //Align coefficients table
    if ((uint64_t)coefs & 0x1f)
        return -1; //Coefficients still misaligned? Error must be thrown on upper level

    uint8_t btemp[255 + 32];
    uint8_t* lutSSE = lut + SSE_LUT_SSE_OFFSET, * oj = btemp;
    memcpy_s(btemp, k, buffer, k);
    memset(btemp + k, 0, n - k);

    _mm256_zeroall();
    __m256i vmask = _mm256_set1_epi8(0xf);
    //Remainder calculation: long division
    for (int j = 0; j < k; j++)
    {
        int idx = *oj++ * 32;
        __m256i* lutv = (__m256i*)(lutSSE + idx);
        __m256i lmul = _mm256_load_si256(lutv);
        __m256i umul = _mm256_permute2x128_si256(lmul, lmul, 0);
        lmul = _mm256_permute2x128_si256(lmul, lmul, 0x11);

        __m256i* ol = (__m256i*)oj;
        __m256i* cl = (__m256i*)coefs;
        for (int l = 0; l <= length; l++)
        {
            __m256i vio = _mm256_load_si256(cl++), vtt;
            GF_mvs_AVX2(&vtt, &vio, &vmask, &lmul, &umul);

            vtt = _mm256_loadu_si256(ol);
            vtt = _mm256_xor_si256(vtt, vio);
            _mm256_storeu_si256(ol++, vtt);
        }
    }
    memcpy_s(buffer + k, n - k, btemp + k, n - k);
    return 0;
}

int DecodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    int scount = n - k, t = scount >> 1, steps = (scount - 1) >> 5;
    if (lut[0] != 0xff) //LUT table requires alignment
        lut += lut[0];
    if ((uint64_t)lut & 0x1f)
        return -1; //LUT still misaligned? Error must be thrown on upper level

    uint16_t* lutLog = (uint16_t*)lut;
    uint8_t* lutExp = lut + ALU_LUT_EXP_OFFSET;
    uint8_t* lutSSE = lut + SSE_LUT_SSE_OFFSET; 
    //uint8_t* lutLog = lut, * lutExp = lut + SSE_LUT_EXP_OFFSET, * lutSSE = lut + SSE_LUT_SSE_OFFSET;
    __declspec(align(32)) uint8_t lambda[2 * MAX_T], syn[2 * MAX_T];
    __declspec(align(32)) uint8_t b[2 * MAX_T], Lm[2 * MAX_T];

    /*Syndrome calculation*/
    _mm256_zeroall();
    __m256i vz = _mm256_setzero_si256();
    memset(syn, buffer[n - 1], MAX_T * 2);
    uint8_t* roots = lut + SSE_LUT_ROOTS_OFFSET;
    __m256i* ls = (__m256i*)syn;
    __m256i vmask = _mm256_set1_epi8(0xf);
    for (int j = n - 2; j >= 0; j--)
    {
        int idx = buffer[j] * 32;
        __m256i* lutv = (__m256i*)(lutSSE + idx);
        __m256i lmul = _mm256_load_si256(lutv);
        __m256i umul = _mm256_permute2x128_si256(lmul, lmul, 0);
        lmul = _mm256_permute2x128_si256(lmul, lmul, 0x11);

        __m256i* lr = (__m256i*)roots;
        for (int m = 0; m <= steps; m++)
        {
            __m256i vio = _mm256_load_si256(lr++), vtt;
            //Multiply: {roots}^(n - j) * buffer[j]
            GF_mvs_AVX2(&vtt, &vio, &vmask, &lmul, &umul);
            vtt = _mm256_load_si256(ls);
            vtt = _mm256_xor_si256(vtt, vio);
            _mm256_store_si256(ls++, vtt);
        }
        ls = (__m256i*)syn;
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
    __m256i* bvec = (__m256i*)b;
    __m256i* lvec = (__m256i*)lambda;
    for (int m = 0; m <= steps; m++)
    {
        _mm256_store_si256(bvec++, vz);
        _mm256_store_si256(lvec++, vz);
    }
    b[1] = 1;
    lambda[0] = 1;
    int l = 0;
    for (int r = 1; r <= scount; r++)
    {
        uint16_t* si = syn2 + r - 1;
        uint8_t* li = lambda;
        uint8_t delta = lutExp[*si--];
        for (int m = 0; m < l; m++)
        {
            uint16_t lm = lutLog[*++li] + *si--;
            delta ^= lutExp[lm];
        }        
        
        int sx = r <= t ? r - 1 : t;
        int sr = sx >> 5; //Optimization for less copy operations
        if (delta)
        {
            int ri = delta;
            __m256i* lutv = (__m256i*)(lutSSE + ri * 32LL);
            __m256i lmul = _mm256_load_si256(lutv);
            __m256i umul = _mm256_permute2x128_si256(lmul, lmul, 0);
            lmul = _mm256_permute2x128_si256(lmul, lmul, 0x11);

            bvec = (__m256i*)b;
            lvec = (__m256i*)lambda;
            __m256i* mvec = (__m256i*)Lm;
            for (int m = 0; m <= sr; m++)
            {
                __m256i vio = _mm256_load_si256(bvec++), vtt;
                //Multiply vector b by scalar delta
                GF_mvs_AVX2(&vtt, &vio, &vmask, &lmul, &umul);

                vtt = _mm256_load_si256(lvec);
                _mm256_store_si256(mvec++, vtt); //Copy lambda to Lm
                vtt = _mm256_xor_si256(vtt, vio);
                _mm256_store_si256(lvec++, vtt); //Save result to lambda
            }
            if (2 * l < r)
            {
                l = r - l;
                ri = 255 - lutLog[ri];
                ri = lutExp[ri]; //delta^-1
                lutv = (__m256i*)(lutSSE + ri * 32LL);
                lmul = _mm256_load_si256(lutv);
                umul = _mm256_permute2x128_si256(lmul, lmul, 0);
                lmul = _mm256_permute2x128_si256(lmul, lmul, 0x11);

                bvec = (__m256i*)b;
                mvec = (__m256i*)Lm;
                for (int m = 0; m <= sr; m++)
                {
                    __m256i vio = _mm256_load_si256(mvec++), vtt;
                    //Multiply vector b by scalar delta
                    GF_mvs_AVX2(&vtt, &vio, &vmask, &lmul, &umul);

                    _mm256_store_si256(bvec++, vio); //Save result to b
                }
            }
        }
        //Shift b left
        sr = (sx + 1) >> 5;
        bvec = (__m256i*)b + sr;
        __m256i vx = _mm256_load_si256(bvec);
        uint8_t tmp8 = ((uint8_t*)bvec)[15];
        vx = _mm256_slli_si256(vx, 1);
        _mm256_store_si256(bvec, vx);
        ((uint8_t*)bvec)[16] = tmp8;
        if (!sr) continue;
        uint8_t* bi = (uint8_t*)--bvec + 1;
        for (int m = 0; m < sr; m++)
        {
            vx = _mm256_load_si256(bvec--);
            _mm256_storeu_si256((__m256i*)bi, vx);
            bi -= 32;
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

    /* Omega calculation */
    //Omega must be calculated only up to {nerr} power
    uint16_t omega[2 * MAX_T];
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
    We need to substitute x[] to labmda(x), where x[] = x^-n+1, x^-n+2, ..., 1 = x^256-n, x^256-n+1, ..., 1 */
    steps = (n - 1) >> 5;
    uint8_t efound = 0, ecorr = 0;
    __declspec(align(32)) uint8_t rtemp[256];
    memset(rtemp, 1, n);
    roots = lut + SSE_LUT_ROOTS_OFFSET + 256 - n;
    for (int j = 0; j < nerr; j++)
    {
        int idx = lambda[j + 1] * 32;
        __m256i* lutv = (__m256i*)(lutSSE + idx);
        __m256i lmul = _mm256_load_si256(lutv);
        __m256i umul = _mm256_permute2x128_si256(lmul, lmul, 0);
        lmul = _mm256_permute2x128_si256(lmul, lmul, 0x11);

        __m256i* lr = (__m256i*)roots;
        ls = (__m256i*)rtemp;
        for (int m = 0; m <= steps; m++)
        {
            __m256i vio = _mm256_loadu_si256(lr++), vtt;
            //Multiply: {roots}^(n - j) * buffer[j]
            GF_mvs_AVX2(&vtt, &vio, &vmask, &lmul, &umul);
            vtt = _mm256_load_si256(ls);
            vtt = _mm256_xor_si256(vtt, vio);
            _mm256_store_si256(ls++, vtt);
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
