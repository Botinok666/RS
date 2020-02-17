#include "rsalurlut.h"
//#include <stdio.h>
//#define ERROR_CHECKING

// Fast multiply using table lookup
uint8_t inline GFMulFast(uint8_t a, uint8_t b, uint16_t* lutLog, uint16_t* lutExp)
{
    int t = lutLog[a] + lutLog[b];
    return (uint8_t)lutExp[t];
}
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
uint8_t inline GFInvFast(uint8_t a, uint16_t* lutLog, uint16_t* lutExp)
{
    if (!a) return 0;
    int t = 255 - lutLog[a];
    return (uint8_t)lutExp[t];
}
// Fill the Exp/Log table
void FillRLUT(uint16_t* lut) {
    uint8_t x = 1, index = 0;
    uint16_t* lutExp = lut + 256;
    lutExp[0] = x; //0:255 - log, 256:767 - exp, 768:1792 - zero
    for (int i = 0; i < 255; i++) {
        uint8_t y = GFMul(x, 2);
        lutExp[++index] = y;
        lutExp[index + 255] = y;
        x = y;
    }
    lut[0] = 511; //Log(0) = inf
    for (int i = 0; i < 255; i++)
        lut[lutExp[i]] = i;
    memset(lutExp + 511, 0, 1024);
}
void FillCoefficents(uint16_t* coefs, uint8_t count, uint16_t* lut)
{
    if (count < 2) return;

    uint16_t* lutExp = lut + 256, * lutLog = lut;
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
int __cdecl Encode(uint8_t n, uint8_t k, uint16_t* lut, uint16_t* coefs, uint8_t* in, uint8_t* out)
{
    //Required size of scratch buffer is n
    int length = n - k + 1;
    if (length < 3 || n < 4)
        return -1;
    uint16_t* lutExp = lut + 256, *lutLog = lut;
    memcpy_s(out, k, in, k);
    memset(out + k, 0, n - k);
    //Remainder calculation: long division
    for (int j = 0; j < k; j++)
    {
        uint16_t mul = lutLog[out[j]];
        for (int l = 1; l < length; l++)
            out[j + l] ^= (uint8_t)lutExp[coefs[l] + mul];
    }
    memcpy_s(out, k, in, k);
    return 0;
}
int __cdecl Decode(uint8_t n, uint8_t k, uint16_t* lut, uint16_t* scratch, uint8_t* in, uint8_t* out)
{
    //Required size of scratch buffer is 5*scount
    uint8_t scount = n - k;
    if (scount < 2 || n < 4)
        return -1;
    uint16_t* lutExp = lut + 256, * lutLog = lut;
    uint8_t t = scount >> 1; //Error capacity

    uint16_t* lambda = scratch; // new uint8_t[scount]
    uint16_t* omega = lambda + scount; //new uint8_t[t + 1]
    uint16_t* syn = omega + t + 1; //new uint8_t[scount]

    /*Syndrome calculation: substitution method
    Buffer view: B_0, B_1, ...
    Polynomial view: B(x) = B_0*x^n-1, B_1*x^n-2, ...
    Roots of generator polynomial: R_0 = 1, R_1 = 2, ...

    Idea is to compute indexes C[](B_x, R_0) = log(B_0) + log(R_0)*(n-1), log(B_1) + log(R_0)*(n-2), ...
    Obviously, C[](B_x, R_0) = log(B_0), log(B_1), ...
    C[](B_x, R_y) = log(B_0) + log(R_y)*(n-1), ... = [log(R_1) = log(R_0) + 1] = C_0 + n - 1, C_1 + n - 2, ...
    Then, S(R_y) = XOR(C[](B_x, R_y)) = exp(C_0) ^ exp(C_1) ^ ...
    */
    //uint16_t stemp[16]; //2n-k + 2n
    //for (int j = 0; j < n; j++)
    //{
    //    if (!btemp[j]) continue;
    //    stemp[j] = lutLog[btemp[j]];
    //}
    //for (int j = 0; j < scount; j++)
    //{
    //    syn[j] = 0;
    //    for (int m = 0; m < n; m++)
    //    {
    //        if (!btemp[m]) continue;
    //        syn[j] ^= lutExp[stemp[m]];
    //        stemp[m] += n - m - 1;
    //        if (stemp[m] >= 255)
    //            stemp[m] -= 255;
    //    }
    //}
    //printf_s("Syndromes alt: ");
    //for (int j = 0; j < scount; j++)
    //    printf_s("%d ", syn[j]);
    //printf_s("\n");
    
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
        for (int m = 0; m < n - 1; m++)
        {
            s = lutLog[in[m] ^ s];
            s = lutExp[s + j];
        }
        s ^= in[n - 1];
        hasErrors |= s;
        syn[j] = lutLog[s];
    }

    //printf_s("Syndromes: ");
    //for (int j = 0; j < scount; j++)
    //    printf_s("%d ", syn[j]);
    //printf_s("\n");

    if (!hasErrors) return 0;

    //Lambda calculation: Berlekamp's method
    uint16_t* b = syn + scount, *Lm = b + scount; //5 * scount
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
        int idx = r - 1;
        uint16_t delta = lutExp[syn[idx--]];
        for (int m = 1; m <= l; m++)
        {
            uint16_t li = lutLog[lambda[m]];
            delta ^= lutExp[li + syn[idx--]];
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
    for (int m = scount - 1; m > 0; m--) //Find deg(lambda)
    {
        if (lambda[m])
        {
            nerr = m;
            break;
        }
    }
    if (nerr > t)
        return -2; //deg(lambda) > t? Uncorrectable error pattern occured
    //scratch: lambda, omega, syn

    //printf_s("Lambda: ");
    //for (int j = 0; j < lcount; j++)
    //    printf_s("%d ", lambda[j]);
    //printf_s("\n");

    //char err = 0;
    //for (int j = t; j < scount; j++)
    //{
    //    int ls = 0;
    //    for (int m = 1; m <= t; m++)
    //        ls ^= GFMulFast(lambda[m], syn[j - m], lutLog, lutExp);
    //    //printf_s("S=%d, L=%d\n", syn[j], ls);
    //    if (ls != syn[j]) return -2;
    //}
    //if (err)
    //    printf_s("Lambda check failed\n");
    //else
    //    printf_s("Lambda check OK\n");

    /* Omega calculation
    O_y = S_y ^ (lambda_1 * S_y-1) ^ ... ^ (lambda_l * S_y-l)
    */
    //Omega must be calculated only up to t power
    for (int m = 0; m <= t; m++)
    {
        omega[m] = lutExp[syn[m]];
        for (int l = 1; l <= m; l++)
        {
            uint16_t li = lutLog[lambda[l]];
            omega[m] ^= lutExp[li + syn[m - l]];
        }
    }
    //scratch: lambda, omega

    uint16_t* xterm = (uint16_t*)(omega + t + 1);// new uint16_t[t]
    /* Chien search
    We need to substitute x[] to labmda(x), where x[] = x^-n+1, x^-n+2, ..., 1 = x^256-n, x^256-n+1, ..., 1
    Because lambda[] = lambda_0 * x^l-1, lambda_1 * x^l-2, ... lambda_l
    We can precompute indexes in xterm[]: xterm_j = (log(lambda_j+1) + (256 - n) * (j + 1)) mod 255
    Then for each x: B_j = 1 ^ exp(xterm_0) ^ exp(xterm_1) ^ ... ^ exp(xterm_m-1)
    And change xterm[]: xterm_m = (xterm_m + m + 1) mod 255
    */
    int x = 256 - n;
    for (int j = 0; j < t; j++)
    {
        if (!lambda[j + 1])
        {
            xterm[j] = 256;
            continue;
        }
        xterm[j] = lutLog[lambda[j + 1]] + x * (j + 1);
        xterm[j] %= 255;
    }
    uint8_t efound = 0;
#ifdef ERROR_CHECKING
    for (int j = 0; j < n; j++)
#else
    for (int j = 0; j < k; j++)
#endif // ERROR_CHECKING
    {
        uint8_t p = 1;
        for (int m = 0; m < t; m++)
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
            efound++;
            uint16_t xIdx = 256 - n + (uint8_t)j, xp = xIdx;
            uint8_t s = 0, y = (uint8_t)omega[0];
            for (int l = 1; l <= t; l++)
            {
                uint16_t lg = lutLog[omega[l]];
                y ^= lutExp[xIdx + lg]; //Omega(X^-1)
                if (l & 1)
                {
                    lg = lutLog[lambda[l]];
                    s ^= lutExp[xIdx + lg]; //Lambda'(X^-1) * X^-1
                }
                xIdx += xp;
                if (xIdx >= 255)
                    xIdx -= 255;
            }
            s = GFInvFast(s, lutLog, lutExp);
            s = GFMulFast(y, s, lutLog, lutExp);

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
