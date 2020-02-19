#include "rssse4.h"
#include <intrin.h>
//#include <stdio.h>
//#define ERROR_CHECKING

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
// Fill the Exp/Log/SSE tables
void FillSSELUT(uint8_t* lut) {
    uint8_t x = 1, index = 0;
    uint8_t* lutExp = lut + 256, *lutSSE = lut + 768;
    lutExp[0] = x; //0:255 - log, 256:767 - exp, 768:8960 - sse
    for (int i = 0; i < 255; i++) {
        uint8_t y = GFMul3(x, 2);
        lutExp[++index] = y;
        lutExp[index + 255] = y;
        x = y;
    }
    lutExp[767] = 0;

    lut[0] = 255; //Log(0) = inf
    for (int i = 0; i < 255; i++)
    {
        lut[lutExp[i]] = i;

        uint8_t* usse = lutSSE + i * 32LL, *lsse = usse + 16;
        for (uint8_t j = 0; j < 16; j++)
        {
            *usse++ = GFMul3(j << 4, (uint8_t)i);
            *lsse++ = GFMul3(j, (uint8_t)i);
        }
    }
}
void FillCoefficentsSSE(uint8_t* coefs, uint8_t count, uint8_t* lut)
{
    if (count < 2) return;

    uint8_t* lutExp = lut + 256, * lutLog = lut;
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
    memset(coefs + count + 1, 0, COEFS_SIZE_SSE - count - 1);
}
int __cdecl EncodeSSE4(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* in, uint8_t* out)
{
    //Required size of scratch buffer is n
    uint8_t length = n - k;
    if (length < 3 || n < 4)
        return -1;
    length >>= 4;
    uint8_t buffer[255 + 16];
    uint8_t* lutExp = lut + 256, * lutLog = lut, * lutSSE = lut + 768;
    uint8_t* oj = buffer;
    __m128i vmask = _mm_set1_epi8(0xf);
    memcpy_s(buffer, k, in, k);
    memcpy_s(out, k, in, k);
    memset(buffer + k, 0, n - k);
    coefs++; //We don't need coefs[0]

    //Remainder calculation: long division
    for (int j = 0; j < k; j++)
    {
        int idx = *oj++ * 32;
        __m128i* lutv = (__m128i*)(lutSSE + idx);
        __m128i umul = _mm_loadu_si128(lutv);
        __m128i lmul = _mm_loadu_si128(lutv + 1);

        __m128i* ol = (__m128i*)oj;
        __m128i* cl = (__m128i*)coefs;
        for (uint8_t l = 0; l <= length; l++)
        {
            __m128i ucoefs = _mm_loadu_si128(cl++);
            __m128i lcoefs = _mm_and_si128(ucoefs, vmask); //Lower 4 bits
            lcoefs = _mm_shuffle_epi8(lmul, lcoefs); //Multiply lower 4 bits
            ucoefs = _mm_srli_epi64(ucoefs, 4);
            ucoefs = _mm_and_si128(ucoefs, vmask); //Upper 4 bits
            ucoefs = _mm_shuffle_epi8(umul, ucoefs); //Multiply upper 4 bits
            ucoefs = _mm_xor_si128(ucoefs, lcoefs); //Result of multiplication

            lcoefs = _mm_loadu_si128(ol);
            lcoefs = _mm_xor_si128(lcoefs, ucoefs);
            _mm_store_si128(ol++, lcoefs);
        }
    }
    memcpy_s(out + k, n - k, buffer + k, n - k);
    return 0;
}
int __cdecl DecodeSSE4(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* scratch, uint8_t* in, uint8_t* out)
{
    //Required size of scratch buffer is 5*scount
    uint8_t scount = n - k;
    if (scount < 2 || n < 4)
        return -1;
    uint8_t* lutExp = lut + 256, * lutLog = lut;
    uint8_t t = scount >> 1; //Error capacity

    uint16_t lambda[MAX_T_SSE], omega[MAX_T_SSE / 2 + 1], syn[MAX_T_SSE];

    /*Syndrome calculation: Horner's method
    Buffer view: B_0, B_1, ...
    Polynomial view: B[] = B_0*x^n-1, B_1*x^n-2, ...
    Roots of generator polynomial: R_0 = 1, R_1 = 2, ...

    Stepper method. S_0(R_y) = 0, S_1(R_y) = exp(y + log(B_0 ^ S_0)), ... S(R_y) = B_n-1 ^ S_n-1
    */
    uint16_t hasErrors = 0;
    for (int j = 0; j < scount; j++)
    {
        uint16_t s = 0;
        uint8_t* inm = in;
        for (int m = 0; m < n - 1; m++)
        {
            s = lutLog[*inm++ ^ s];
            s = lutExp[s + j];
        }
        s ^= *inm;
        hasErrors |= s;
        syn[j] = lutLog[s];
    }

    if (!hasErrors) return 0;

    //Lambda calculation: Berlekamp's method
    uint16_t b[MAX_T_SSE], Lm[MAX_T_SSE];
    //Set initial value of b and lambda to 1
    memset(b, 0, scount * sizeof(uint16_t));
    b[0] = 1;
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
        uint8_t* si = syn + r - 1;
        uint16_t delta = lutExp[*si--];
        for (int m = 1; m <= l; m++)
        {
            uint16_t li = lutLog[lambda[m]];
            delta ^= lutExp[li + *si--];
        }

        if (delta)
        {
            memcpy_s(Lm, scount * sizeof(uint16_t), lambda, scount * sizeof(uint16_t));
            uint16_t dlg = lutLog[delta], blg;
            for (int m = 1; m < scount; m++)
            {
                blg = lutLog[b[m - 1]];
                lambda[m] ^= (uint8_t)lutExp[blg + dlg];
            }
            if (2 * l <= r - 1)
            {
                l = r - l;
                dlg = 255 - dlg; //Inverse
                for (int m = 0; m < scount; m++)
                {
                    blg = lutLog[Lm[m]];
                    b[m] = (uint8_t)lutExp[blg + dlg];
                }
            }
            else
            {
                for (int m = scount - 1; m > 0; m--) //Shift
                    b[m] = b[m - 1];
                b[0] = 0;
            }
        }
        else
        {
            for (int m = scount - 1; m > 0; m--) //Shift
                b[m] = b[m - 1];
            b[0] = 0;
        }
    }
    uint8_t nerr = 0;
    for (int m = 0; m < scount; m++) //Find deg(lambda)
    {
        uint16_t lg = lambda[m];
        if (lg)
            nerr = m;
        lambda[m] = lutLog[lg]; //Convert to indexes
    }
#ifdef ERROR_CHECKING
    if (nerr > t)
        return -2; //deg(lambda) > t? Uncorrectable error pattern occured
#endif // ERROR_CHECKING
    //scratch: lambda, omega, syn

    /* Omega calculation
    O_y = S_y ^ (lambda_1 * S_y-1) ^ ... ^ (lambda_l * S_y-l)
    */
    //Omega must be calculated only up to t power
    for (int m = 0; m <= t; m++)
    {
        uint16_t og = lutExp[syn[m]];
        for (int l = 1; l <= m; l++)
        {
            uint16_t li = lambda[l];
            og ^= lutExp[li + syn[m - l]];
        }
        omega[m] = lutLog[og]; //Convert to indexes
    }
    //scratch: lambda, omega

    uint8_t* xterm = syn;
    /* Chien search
    We need to substitute x[] to labmda(x), where x[] = x^-n+1, x^-n+2, ..., 1 = x^256-n, x^256-n+1, ..., 1
    Because lambda[] = lambda_0 * x^l-1, lambda_1 * x^l-2, ... lambda_l
    We can precompute indexes in xterm[]: xterm_j = (log(lambda_j+1) + (256 - n) * (j + 1)) mod 255
    Then for each x: B_j = 1 ^ exp(xterm_0) ^ exp(xterm_1) ^ ... ^ exp(xterm_m-1)
    And change xterm[]: xterm_m = (xterm_m + m + 1) mod 255
    */
    int x = 256 - n;
    for (int j = 0; j < nerr; j++)
    {
        if (lambda[j + 1] == 511)
        {
            xterm[j] = 511;
            continue;
        }
        xterm[j] = lambda[j + 1] + x * (j + 1);
        xterm[j] %= 255;
    }
#ifdef ERROR_CHECKING
    uint8_t efound = 0;
    for (int j = 0; j < n; j++)
#else
    for (int j = 0; j < k; j++)
#endif // ERROR_CHECKING
    {
        uint8_t p = 1;
        for (int m = 0; m < nerr; m++)
        {
            uint16_t e = xterm[m];
            if (e > 255) continue;
            p ^= lutExp[e];
            e += m + 1;
            if (e >= 255)
                e -= 255;
            xterm[m] = e;
        }
        if (!p)
        {
#ifdef ERROR_CHECKING
            efound++;
#endif // ERROR_CHECKING
            uint16_t xIdx = 256 - n + (uint8_t)j, xp = xIdx;
            uint8_t s = 0, y = (uint8_t)lutExp[omega[0]];
            for (int l = 1; l <= nerr; l++)
            {
                uint16_t lg = omega[l];
                y ^= lutExp[xIdx + lg]; //Omega(X^-1)
                if (l & 1)
                {
                    lg = lambda[l];
                    s ^= lutExp[xIdx + lg]; //Lambda'(X^-1) * X^-1
                }
                xIdx += xp;
                if (xIdx >= 255)
                    xIdx -= 255;
            }
            if (s)
            {
                xp = 255 - lutLog[s];
                s = (uint8_t)lutExp[xp + lutLog[y]];
            }

            in[j] ^= s;
        }
    }
#ifdef ERROR_CHECKING
    if (efound != nerr)
        return -2;
#endif // ERROR_CHECKING

    memcpy_s(out, k, in, k);
    return nerr;
}
