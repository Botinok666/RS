#include <intrin.h>
#include "rsavx2.h"
#define ERROR_CHECKING
// Slow multiply, using shifting
// Polynomial x^8 + x^4 + x^3 + x^2 + 1
uint8_t GFMul4(uint8_t a, uint8_t b)
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

static inline void GFMulAVX2(__m256i* va, __m256i* vb, __m256i* vz, __m256i* vs, __m256i* vgp)
{
    __m256i vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
    vx = _mm256_cmpgt_epi8(*vz, *vb);
    vx = _mm256_and_si256(vx, *vgp);
    *vb = _mm256_add_epi8(*vb, *vb);
    *vb = _mm256_xor_si256(*vb, vx);
    *va = _mm256_add_epi8(*va, *va);

    vx = _mm256_cmpgt_epi8(*vz, *va);
    vx = _mm256_and_si256(vx, *vb);
    *vs = _mm256_xor_si256(*vs, vx);
}

void InitAVX2(uint8_t* coefsu, const uint8_t count, uint8_t* lut)
{
    if (count < 2) return;
    uint64_t offset = (uint64_t)lut & 0x1f;
    if (offset) //Array is not aligned to 32 bytes boundary
    {
        offset = 0x20 - offset;
        lut[0] = (uint8_t)offset; //Store offset
        lut += offset; //Offset LUT pointer
    }
    uint8_t* lutExp = lut + SSE_LUT_EXP_OFFSET, * lutSSE = lut + SSE_LUT_SSE_OFFSET;
    uint8_t x = 1, index = 0;
    lutExp[0] = x; //0:255 - log, 256:767 - exp, 768:8960 - sse
    for (int i = 0; i < 255; i++) {
        uint8_t y = GFMul4(x, 2);
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
        uint8_t* usse = lutSSE + i * 32LL, * lsse = usse + 16;
        for (uint8_t j = 0; j < 16; j++)
        {
            *usse++ = GFMul4(j << 4, (uint8_t)i);
            *lsse++ = GFMul4(j, (uint8_t)i);
        }
    }
    __m256i* lutRev = (__m256i*)(lut + SSE_LUT_REV_OFFSET);
    __m256i* roots = (__m256i*)lutExp;
    __m256i vmask = _mm256_set1_epi8(0xf);
    __m256i vlrev = _mm256_broadcastsi128_si256(
        _mm_set_epi8(15, 7, 11, 3, 13, 5, 9, 1, 14, 6, 10, 2, 12, 4, 8, 0));
    __m256i vhrev = _mm256_broadcastsi128_si256(
        _mm_set_epi8(-16, 112, -80, 48, -48, 80, -112, 16, -32, 96, -96, 32, -64, 64, -128, 0));
    //Fill table with roots with reverse ordered bits
    for (int i = 0; i < 256; i += 32)
    {
        __m256i a = _mm256_load_si256(roots++);
        __m256i c = _mm256_and_si256(a, vmask);
        c = _mm256_shuffle_epi8(vhrev, c);
        a = _mm256_srli_epi64(a, 4);
        a = _mm256_and_si256(a, vmask);
        a = _mm256_shuffle_epi8(vlrev, a);
        a = _mm256_or_si256(a, c);
        _mm256_store_si256(lutRev++, a);
    }
    _mm256_store_si256(lutRev, _mm256_setzero_si256()); //Clear next 32 bytes so during unaligned load there will be no random data

    uint8_t* lutLog = lut;
    //We need to align array at 16 bytes boundary
    offset = 0x20 - ((uint64_t)coefsu & 0x1f);
    coefsu[0] = (uint8_t)offset; //Save offset
    uint8_t* coefs = coefsu + offset - 1; //coefs[0] will be placed outside of aligned array
    __m128i vz = _mm_setzero_si128();
    __m128i* cptr = (__m128i*)(coefsu + offset);
    for (int i = 0; i < AVX_COEFS_SIZE / 16; i++) //Clear coefs
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

void inline GF_mvs_AVX2(__m256i* vio, __m256i* vt, __m256i* vmask, __m256i* lmul, __m256i* umul)
{
    *vio = _mm256_and_si256(*vt, *vmask); //Lower 4 bits
    *vio = _mm256_shuffle_epi8(*lmul, *vio); //Multiply lower 4 bits
    *vt = _mm256_srli_epi64(*vt, 4);
    *vt = _mm256_and_si256(*vt, *vmask); //Upper 4 bits
    *vt = _mm256_shuffle_epi8(*umul, *vt); //Multiply upper 4 bits
    *vt = _mm256_xor_si256(*vt, *vio); //Result of multiplication
}

int EncodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer)
{
    int length = n - k;
    if (lut[0] != 0xff) //LUT table requires alignment
        lut += lut[0];
    if ((uint64_t)lut & 0x1f)
        return -3; //LUT still misaligned? Error must be thrown on upper level
    coefs += coefs[0]; //Align coefficients table
    if ((uint64_t)coefs & 0x1f)
        return -4; //Coefficients still misaligned? Error must be thrown on upper level

    length = (length - 1) >> 5;
    uint8_t btemp[255 + 32];
    uint8_t* lutExp = lut + SSE_LUT_EXP_OFFSET, * lutLog = lut, * lutSSE = lut + SSE_LUT_SSE_OFFSET;
    uint8_t* oj = btemp;
    memcpy_s(btemp, k, buffer, k);
    memset(btemp + k, 0, n - k);

    __m256i vmask = _mm256_set1_epi8(0xf);
    //Remainder calculation: long division
    for (int j = 0; j < k; j++)
    {
        int idx = *oj++ * 32;
        __m128i* lutv = (__m128i*)(lutSSE + idx);
        __m256i umul = _mm256_broadcastsi128_si256(_mm_load_si128(lutv));
        __m256i lmul = _mm256_broadcastsi128_si256(_mm_load_si128(lutv + 1));

        __m256i* ol = (__m256i*)oj;
        __m256i* cl = (__m256i*)coefs;
        for (int l = 0; l <= length; l++)
        {
            __m256i ucoefs = _mm256_load_si256(cl++), lcoefs;            
            /*__m256i lcoefs = _mm256_and_si256(ucoefs, vmask); //Lower 4 bits
            lcoefs = _mm256_shuffle_epi8(lmul, lcoefs); //Multiply lower 4 bits
            ucoefs = _mm256_srli_epi64(ucoefs, 4);
            ucoefs = _mm256_and_si256(ucoefs, vmask); //Upper 4 bits
            ucoefs = _mm256_shuffle_epi8(umul, ucoefs); //Multiply upper 4 bits
            ucoefs = _mm256_xor_si256(ucoefs, lcoefs); //Result of multiplication*/
            GF_mvs_AVX2(&ucoefs, &lcoefs, &vmask, &lmul, &umul);

            lcoefs = _mm256_loadu_si256(ol);
            lcoefs = _mm256_xor_si256(lcoefs, ucoefs);
            _mm256_storeu_si256(ol++, lcoefs);
        }
    }
    memcpy_s(buffer + k, n - k, btemp + k, n - k);
    return 0;
}

int DecodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    int scount = n - k, steps = (scount - 1) >> 5;
    if (lut[0] != 0xff) //LUT table requires alignment
        lut += lut[0];
    if ((uint64_t)lut & 0x1f)
        return -3; //LUT still misaligned? Error must be thrown on upper level

    uint8_t* lutLog = lut, * lutExp = lut + SSE_LUT_EXP_OFFSET, * lutSSE = lut + SSE_LUT_SSE_OFFSET;
    __declspec(align(32)) uint8_t lambda[2 * MAX_T], omega[2 * MAX_T], syn[2 * MAX_T];
    __declspec(align(32)) uint8_t b[2 * MAX_T], Lm[2 * MAX_T];

    /*Syndrome calculation*/
    int hasNoErrors = -1;
    __m256i vz = _mm256_setzero_si256();
    _mm256_store_si256((__m256i*)lambda, vz);
    memset(lambda, 0xff, scount & 0x1f); //Temporal storage for root mask
    __m256i vgp = _mm256_set1_epi8(0x1d); //Generator polynomial
    __m256i* lrev = (__m256i*)(lut + SSE_LUT_REV_OFFSET);
    __m256i* vsyn = (__m256i*)syn;
    for (int j = 0; j <= steps; j++)
    {
        __m256i vroot = _mm256_load_si256(lrev++);
        __m256i vs = vz;
        uint8_t* inm = buffer;
        for (int m = 0; m < n - 1; m++)
        {
            __m256i va = vroot;
            __m256i vb = _mm256_xor_si256(vs, _mm256_set1_epi8(*inm++)); //s ^= in[m]
            //Perform full multiplication: s *= root, or s = b * root
            vs = vz; //Clear s
            GFMulAVX2(&va, &vb, &vz, &vs, &vgp);
        }
        vs = _mm256_xor_si256(vs, _mm256_set1_epi8(*inm));
        if (j == steps && lambda[0]) //Non-zero lambda[0] means that scount % 0xf != 0
            vs = _mm256_and_si256(vs, _mm256_load_si256((__m256i*)lambda));
        hasNoErrors &= _mm256_movemask_epi8(_mm256_cmpeq_epi8(vs, vz));
        _mm256_store_si256(vsyn++, vs);
    }
    if (hasNoErrors == -1) return 0;

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
    __m256i vmask = _mm256_set1_epi8(0xf);
    for (int r = 1; r <= scount; r++)
    {
        uint8_t* si = syn + r - 1, * li = lambda;
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
        uint8_t sr = (r > 33) ? ((r - 2) >> 5) : 0; //Optimization for less copy operations
        if (delta)
        {
            int ri = delta;
            __m128i* lutv = (__m128i*)(lutSSE + ri * 32LL);
            __m256i umul = _mm256_broadcastsi128_si256(_mm_load_si128(lutv));
            __m256i lmul = _mm256_broadcastsi128_si256(_mm_load_si128(lutv + 1));

            bvec = (__m256i*)b;
            lvec = (__m256i*)lambda;
            __m256i* mvec = (__m256i*)Lm;
            for (int m = 0; m <= sr; m++)
            {
                __m256i vx = _mm256_load_si256(bvec++), vxl;
                //Multiply vector b by scalar delta
                GF_mvs_AVX2(&vx, &vxl, &vmask, &lmul, &umul);

                vxl = _mm256_load_si256(lvec);
                _mm256_store_si256(mvec++, vxl); //Copy lambda to Lm
                vxl = _mm256_xor_si256(vxl, vx);
                _mm256_store_si256(lvec++, vxl); //Save result to lambda
            }
            if (2 * l <= r - 1)
            {
                l = r - l;
                ri = 255 - lutLog[ri];
                ri = lutExp[ri]; //delta^-1
                lutv = (__m128i*)(lutSSE + ri * 32LL);
                umul = _mm256_broadcastsi128_si256(_mm_load_si128(lutv));
                lmul = _mm256_broadcastsi128_si256(_mm_load_si128(lutv + 1));

                bvec = (__m256i*)b;
                mvec = (__m256i*)Lm;
                for (int m = 0; m <= sr; m++)
                {
                    __m256i vx = _mm256_load_si256(mvec++), vxl;
                    //Multiply vector b by scalar delta
                    GF_mvs_AVX2(&vx, &vxl, &vmask, &lmul, &umul);

                    _mm256_store_si256(bvec++, vx); //Save result to b
                }
            }
        }
        //Shift b left
        bvec = (__m256i*)b + steps;
        __m256i vx = _mm256_load_si256(bvec);
        vx = _mm256_slli_si256(vx, 1);
        _mm256_store_si256(bvec--, vx);
        if (!steps) continue;
        uint8_t* bi = (uint8_t*)bvec + 1;
        for (int m = 0; m < steps; m++)
        {
            vx = _mm256_load_si256(bvec--);
            _mm256_storeu_si256((__m256i*)bi, vx);
            bi -= 32;
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
        __m256i vs = _mm256_set1_epi8(*li--); //s = lambda[nerr]
        __m256i vroot = _mm256_loadu_si256(lrev++);
        for (int m = 0; m < nerr; m++)
        {
            __m256i va = vroot;
            __m256i vb = vs;
            //Perform full multiplication: s *= root, or s = b * root
            vs = vz;
            GFMulAVX2(&va, &vb, &vz, &vs, &vgp);
            va = _mm256_set1_epi8(*li--);
            vs = _mm256_xor_si256(vs, va); //s ^= lambda[nerr - m]
        }
        vs = _mm256_cmpeq_epi8(vs, vz);
        int mask = _mm256_movemask_epi8(vs), kv = j << 5;
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
