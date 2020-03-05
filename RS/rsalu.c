//#define SSE_TEST
#define ERROR_CHECKING

#include "rsalu.h"
//#include <stdio.h>
#ifdef SSE_TEST
#include <intrin.h>
#endif // SSE_TEST

// Slow multiply, using shifting
// Polynomial x^8 + x^4 + x^3 + x^2 + 1
uint8_t GFMul(uint8_t a, uint8_t b)
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

// Fill the Exp/Log table
void FillRLUT(void* LUT) {
    uint8_t x = 1, index = 0;
    uint16_t* lutExp = (uint16_t*)LUT + 256, * lutLog = (uint16_t*)LUT;
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
    memset(lutExp + 511, 0, 1025);
}
void FillCoefficents(void* Coefs, const uint8_t count, void* lut)
{
    if (count < 2) return;

    uint16_t* lutExp = (uint16_t*)lut + 256, * lutLog = (uint16_t*)lut, * coefs = (uint16_t*)Coefs;
    coefs[count] = lutExp[0];
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
        coefs[count] = lutExp[s];
    }
    for (int j = 0; j <= count; j++)
        coefs[j] = lutLog[coefs[j]];
}
int __cdecl Encode(const uint8_t n, const uint8_t k, void* LUT, void* Coefs, uint8_t* buffer)
{
    //Required size of scratch buffer is n
    int length = n - k + 1;
    //if (length < 3 || n < 4)
    //    return -1;
    uint16_t* lutExp = (uint16_t*)LUT + 256, *lutLog = (uint16_t*)LUT;
    uint8_t btemp[255 + 16];
    memcpy_s(btemp, k, buffer, k);
    memset(btemp + k, 0, n - k);
    //Remainder calculation: long division
    uint8_t* oj = btemp;
    uint16_t* coefs = (uint16_t*)Coefs + 1;
    for (int j = 0; j < k; j++)
    {
        int mul = lutLog[*oj++];
        uint8_t* ol = oj;
        uint16_t* cl = coefs;
        for (int l = 1; l < length; l++)
            *ol++ ^= (uint8_t)lutExp[*cl++ + mul];
    }
    memcpy_s(buffer + k, n - k, btemp + k, n - k);
    return 0;
}
int __cdecl Decode(const uint8_t n, const uint8_t k, void* lut, uint8_t* buffer)
{
    int scount = n - k;
    //if (scount < 2 || n < 4)
    //    return -1;
    uint16_t* lutExp = (uint16_t*)lut + 256, * lutLog = (uint16_t*)lut;

#ifndef SSE_TEST
    uint16_t hasErrors = 0;
    uint16_t lambda[MAX_T_ALU], syn[MAX_T_ALU];
    /*Syndrome calculation: Horner's method*/   
    {
        uint16_t stemp[255];
        for (int j = 0; j < n; j++)
            stemp[j] = lutLog[buffer[j]];
        for (int j = 0; j < scount; j++)
        {
            uint16_t s = 0;
            int v = n - 1;
            for (int m = 0; m <= v; m++)
            {
                int ecx = (v - m) * j;
                ecx = (ecx >> 8) + (ecx & 0xff);
                ecx = (ecx >> 8) + (ecx & 0xff);
                s ^= lutExp[stemp[m] + ecx];
            }
            hasErrors |= s;
            syn[j] = lutLog[s];
        }
    }
    if (!hasErrors) return 0;
#else
    __declspec(align(16)) uint16_t lambda[MAX_T_ALU], syn[MAX_T_ALU];
    int steps = (scount - 1) >> 3;
    int hasNoErrors = 0xffff;
    /*Syndrome calculation: Horner's method*/
    {
        uint16_t stemp[256];
        __declspec(align(16)) uint16_t vbuf[8];
        for (int j = 0; j < n; j++)
            stemp[j] = lutLog[buffer[j]];

        int v = n - 1;
        __m128i* vsyn = (__m128i*)syn;
        __m128i* vstmp = (__m128i*)stemp;
        __m128i vj = _mm_set_epi16(7, 6, 5, 4, 3, 2, 1, 0);
        __m128i vz = _mm_setzero_si128();
        __m128i v1 = _mm_set1_epi16(1);
        __m128i v8 = _mm_set1_epi16(8);
        __m128i vmask = _mm_set1_epi16(0xff);

        _mm_store_si128((__m128i*)lambda, vz);
        memset(lambda, 0xff, (long long)(scount & 0x7) << 1); //Temporal storage for root mask

        for (int j = 0; j <= steps; j++)
        {
            __m128i vs = vz;
            __m128i vv = _mm_set1_epi16(n);

            for (int m = 0; m <= v; m++)
            {
                vv = _mm_sub_epi16(vv, v1);
                __m128i ve = _mm_mullo_epi16(vv, vj); //(v-m)*j
                __m128i vx = _mm_and_si128(ve, vmask); //e & 0xff
                ve = _mm_srli_epi16(ve, 8); //e >> 8
                ve = _mm_add_epi16(ve, vx);
                vx = _mm_and_si128(ve, vmask); //e & 0xff
                ve = _mm_srli_epi16(ve, 8); //e >> 8
                ve = _mm_add_epi16(ve, vx);

                vx = _mm_set1_epi16(stemp[m]);
                vx = _mm_add_epi16(vx, ve);

                _mm_store_si128((__m128i*)vbuf, vx);
                vbuf[0] = lutExp[vbuf[0]];
                vbuf[1] = lutExp[vbuf[1]];
                vbuf[2] = lutExp[vbuf[2]];
                vbuf[3] = lutExp[vbuf[3]];
                vbuf[4] = lutExp[vbuf[4]];
                vbuf[5] = lutExp[vbuf[5]];
                vbuf[6] = lutExp[vbuf[6]];
                vbuf[7] = lutExp[vbuf[7]];
                vx = _mm_load_si128((__m128i*)vbuf);

                vs = _mm_xor_si128(vs, vx);
            }
            if (j == steps && lambda[0]) //Non-zero lambda[0] means that scount % 0x7 != 0
                vs = _mm_and_si128(vs, _mm_load_si128((__m128i*)lambda));
            hasNoErrors &= _mm_movemask_epi8(_mm_cmpeq_epi16(vs, vz));
            _mm_store_si128(vsyn++, vs);
            vj = _mm_add_epi16(vj, v8);
        }

        for (int j = 0; j < scount; j++)
            syn[j] = lutLog[syn[j]];
    }
    if (hasNoErrors == 0xffff) return 0;
#endif // !SSE_TEST

    //printf_s("Syndromes: ");
    //for (int j = 0; j < scount; j++)
    //    printf_s("%d ", lutExp[syn[j]]);
    //printf_s("\n");

    //Lambda calculation: Berlekamp's method
    uint16_t omega[MAX_T_ALU], b[MAX_T_ALU], Lm[MAX_T_ALU];
    //Set initial value of b and lambda to 1
    memset(b, 0, scount * sizeof(uint16_t));
    b[1] = 1;
    memset(lambda, 0, scount * sizeof(uint16_t));
    lambda[0] = 1;
    int l = 0;
    /* First, calculate delta = S_r-1 ^ lambda_1 * S_r-2 ^ ... ^ lambda_l * S_r-l-1
    If delta == 0, continue
    Else save lambda to Lm, then calculate new lambda:
    lambda[] = lambda_1 ^ (b_0 * delta), lambda_2 ^ (b_1 * delta), ...
    Now b array needs an update
    If 2*l <= r-1:
        delta = 1 / delta; b[] = Lm_0 * delta, Lm_1 * delta, ...
    else:
        shift b[] to the left, b_0 = 0
    */
    for (int r = 1; r <= scount; r++)
    {
        uint16_t* si = syn + r - 1, * li = lambda;
        uint16_t delta = lutExp[*si--];
        for (int m = 0; m < l; m++)
        {
            uint16_t lm = lutLog[*++li];
            delta ^= lutExp[lm + *si--];
        }
        int c = (r < scount) ? r + 1 : r;
        if (delta)
        {
            li = lambda;
            uint16_t dlg = lutLog[delta], blg;
            Lm[0] = *li++;
            for (int m = 1; m < c; m++)
            {
                blg = lutLog[b[m]];
                Lm[m] = *li;
                *li++ ^= (uint8_t)lutExp[blg + dlg];
            }
            if (2 * l <= r - 1)
            {
                l = r - l;
                dlg = 255 - dlg; //Inverse
                for (int m = 0; m < c; m++)
                {
                    blg = lutLog[Lm[m]];
                    b[m] = (uint8_t)lutExp[blg + dlg];
                }
            }
        }
        for (int m = c - 1; m > 0; m--) //Shift
            b[m] = b[m - 1];
        b[0] = 0;
    }

    int nerr = 0;
    for (int m = 0; m < scount; m++) //Find deg(lambda)
    {
        uint16_t lg = lambda[m];
        if (lg)
            nerr = m;
        lambda[m] = lutLog[lg]; //Convert to indexes
    }
    //printf_s("\nLambda: ");
    //for (int j = 0; j <= nerr; j++)
    //    printf_s("%d ", lutExp[lambda[j]]);
    //printf_s("\n");
    if (nerr > (scount >> 1))
        return -2; //deg(lambda) > t? Uncorrectable error pattern occured

    /* Omega calculation */
    //Omega must be calculated only up to t power
    for (int m = 0; m <= nerr; m++)
    {
        uint16_t og = lutExp[syn[m]];
        for (int l = 1; l <= m; l++)
        {
            uint16_t li = lambda[l];
            og ^= lutExp[li + syn[m - l]];
        }
        omega[m] = lutLog[og]; //Convert to indexes
    }
    //printf_s("Omega: ");
    //for (int j = 0; j <= nerr; j++)
    //    printf_s("%d ", lutExp[omega[j]]);
    //printf_s("\n");

    uint16_t* xterm = syn;
    /* Chien search
    We need to substitute x[] to labmda(x), where x[] = x^-n+1, x^-n+2, ..., 1 = x^256-n, x^256-n+1, ..., 1
    Because lambda[] = lambda_0 * x^l-1, lambda_1 * x^l-2, ... lambda_l
    We can precompute indexes in xterm[]: xterm_j = (log(lambda_j+1) + (256 - n) * (j + 1)) mod 255
    Then for each x: B_j = 1 ^ exp(xterm_0) ^ exp(xterm_1) ^ ... ^ exp(xterm_m-1)
    And change xterm[]: xterm_m = (xterm_m + m + 1) mod 255
    */
    int x = 256 - n;
    for (int j = 1; j <= nerr; j++)
    {
        if (lambda[j] == 511)
            xterm[j - 1] = 511;
        else
        {
            int ecx = lambda[j] + j * x;
            ecx = (ecx >> 8) + (ecx & 0xff);
            ecx = (ecx >> 8) + (ecx & 0xff);
            xterm[j - 1] = (uint16_t)ecx;
        }
    }
#ifdef ERROR_CHECKING
    int efound = 0, ecorr = 0;
    for (int j = 0; j < n; j++)
#else
    for (int j = 0; j < k; j++)
#endif // ERROR_CHECKING
    {
        uint16_t p = 1;
        for (int m = 0; m < nerr; m++)
        {
            int ecx = j * (m + 1);
            ecx = (ecx >> 8) + (ecx & 0xff);
            ecx = (ecx >> 8) + (ecx & 0xff);
            p ^= lutExp[xterm[m] + ecx];
        }
        if (!p)
        {
#ifdef ERROR_CHECKING
            if (j < k) {
#endif // ERROR_CHECKING
                int xIdx = 256 - n + j;
                uint16_t s = 0, y = lutExp[omega[0]];
                for (int l = 1; l <= nerr; l++)
                {
                    int ecx = xIdx * l;
                    ecx = (ecx >> 8) + (ecx & 0xff);
                    ecx = (ecx >> 8) + (ecx & 0xff);
                    y ^= lutExp[ecx + omega[l]]; //Omega(X^-1)
                    if (l & 1)
                        s ^= lutExp[ecx + lambda[l]]; //Lambda'(X^-1) * X^-1
                }
                if (!s) return -4;

                s = 255 - lutLog[s];
                s = lutExp[s + lutLog[y]];
#ifdef ERROR_CHECKING
                b[ecorr] = (uint8_t)j;
                Lm[ecorr++] = s;
            }
            efound++;
#else
            buffer[j] ^= s;
#endif // ERROR_CHECKING
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
