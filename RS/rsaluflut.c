#include "rsaluflut.h"
//#include <stdio.h>
union u16u8
{
    uint16_t u16;
    struct u8s
    {
        uint8_t a, b;
    } u8;
};
// Fast multiply using table lookup
uint8_t inline GFMulFast2(uint8_t a, uint8_t b, uint8_t* lut)
{
    uint16_t t = (uint16_t)a | ((uint16_t)b << 8);
    return lut[t];
}
// Slow multiply, using shifting
// Polynomial x^8 + x^4 + x^3 + x^2 + 1
uint8_t GFMul2(uint8_t a, uint8_t b)
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
uint8_t inline GFInvFast2(uint8_t a, uint8_t* lut)
{
    int t = (int)a | FL_INV_OFFSET;
    return lut[t];
}
uint8_t inline GFIdxFast2(uint8_t a, uint8_t* lut)
{
    int t = (int)a | FL_GFV_OFFSET;
    return lut[t];
}
// Fill the Exp/Log table
void FillFLUT(uint8_t* lut) {
    lut[1 << 16] = 0; //Inverse of 0 is 0
    uint8_t mul = 1;
    for (int i = 0; i < 256; i++) 
    {
        for (int j = i; j < 256; j++)
        {
            uint8_t y = GFMul2((uint8_t)i, (uint8_t)j);
            lut[i | (j << 8)] = y;
            lut[j | (i << 8)] = y;
            if (y == 1)
            {
                lut[FL_INV_OFFSET | j] = i;
                lut[FL_INV_OFFSET | i] = j;
            }
        }
        lut[FL_GFV_OFFSET | i] = mul;
        mul = GFMul2(2, mul);
    }
}
void FillCoefficentsFL(uint8_t* coefs, uint8_t count, uint8_t* lut)
{
    if (count < 2) return;

    coefs[count] = 1;
    coefs[count - 1] = 1;
    union u16u8 idx;
    for (int k = 2; k <= count; k++)
    {
        coefs[count - k] = 1;
        idx.u8.b = GFIdxFast2(k - 1, lut);
        for (int i = count - k + 2; i <= count; i++)
        {
            idx.u8.a = coefs[i - 1];
            coefs[i - 1] = coefs[i] ^ lut[idx.u16];
        }
        idx.u8.a = coefs[count];
        coefs[count] = lut[idx.u16];
    }
}
int __cdecl EncodeFL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* in, uint8_t* out)
{
    int length = n - k + 1;
    if (length < 3 || n < 4)
        return -1;
    union u16u8 idx;
    memcpy_s(out, k, in, k);
    memset(out + k, 0, n - k);
    //Remainder calculation: long division
    uint8_t* oj = out;
    for (int j = 0; j < k; j++)
    {
        idx.u8.b = *oj++; //mul
        uint8_t* ol = oj, * cl = coefs + 1;
        for (int l = 1; l < length; l++)
        {
            idx.u8.a = *cl++;
            *ol++ ^= lut[idx.u16];
        }
    }
    memcpy_s(out, k, in, k);
    return 0;
}
int __cdecl DecodeFL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* scratch, uint8_t* in, uint8_t* out)
{
    uint8_t scount = n - k;
    if (scount < 2 || n < 4)
        return -1;
    union u16u8 idx;
    uint8_t t = scount >> 1; //Error capacity

    uint8_t* lambda = scratch;
    uint8_t* omega = lambda + scount;
    uint8_t* syn = omega + t + 1;

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

    uint8_t hasErrors = 0;
    for (int j = 0; j < scount; j++)
    {
        uint8_t s = 0;
        idx.u8.b = GFIdxFast2(j, lut);
        for (int m = 0; m < n - 1; m++)
        {
            idx.u8.a = in[m] ^ s;
            s = lut[idx.u16];
        }
        s ^= in[n - 1];
        hasErrors |= s;
        syn[j] = s;
    }
    if (!hasErrors) return 0;

    //Lambda calculation: Berlekamp's method
    uint8_t* b = syn + scount, *Lm = b + scount; //4n-3k+2 + 2n-2k+2
    //Set initial value of b and lambda to 1
    memset(b + 1, 0, (size_t)scount - 1);
    b[0] = 1;
    memset(lambda + 1, 0, (size_t)scount - 1);
    lambda[0] = 1;
    int l = 0;

    for (int r = 1; r <= scount; r++)
    {
        uint8_t* si = syn + r - 1;
        uint8_t delta = *si--;
        for (int m = 1; m <= l; m++)
        {
            idx.u8.a = lambda[m];
            idx.u8.b = *si--;
            delta ^= lut[idx.u16];
        }

        if (delta)
        {
            memcpy_s(Lm, scount, lambda, scount);
            idx.u8.b = delta;
            for (int m = 1; m < scount; m++)
            {
                idx.u8.a = b[m - 1];
                lambda[m] ^= lut[idx.u16];
            }

            if (2 * l <= r - 1)
            {
                l = r - l;
                idx.u8.b = GFInvFast2(delta, lut);
                for (int m = 0; m < scount; m++)
                {
                    idx.u8.a = Lm[m];
                    b[m] = lut[idx.u16];
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

    for (int m = 0; m <= t; m++)
    {
        omega[m] = syn[m];
        for (int l = 1; l <= m; l++)
        {
            idx.u8.a = lambda[l];
            idx.u8.b = syn[m - l];
            omega[m] ^= lut[idx.u16];
        }
    }
    //scratch: lambda, omega

    uint8_t* xpos = omega + t + 1;
    uint8_t* xterm = xpos + t, *xmul = xterm + t;

    uint8_t x = GFIdxFast2(256 - n, lut);
    xterm[0] = x;
    xmul[0] = 2;
    union u16u8 idm;
    idm.u8.b = 2;
    idx.u8.b = x;
    for (int j = 1; j < nerr; j++)
    {
        idm.u8.a = xmul[j - 1];
        xmul[j] = lut[idm.u16]; //x, x^2, x^3, ...
        idx.u8.a = xterm[j - 1];
        xterm[j] = lut[idx.u16]; //x^p, (x^p)^2, (x^p)^3, ...
    }
    for (int j = 0; j < nerr; j++)
    {
        idx.u8.a = lambda[j + 1];
        idx.u8.b = xterm[j];
        xterm[j] = lut[idx.u16]; //l_1*x^p, l_2*(x^p)^2, ...
    }

    for (int j = 0; j < k; j++)
    {
        uint8_t p = 1;
        for (int m = 0; m < nerr; m++)
        {
            idm.u8.a = xterm[m];
            p ^= idm.u8.a;
            idm.u8.b = xmul[m];
            xterm[m] = lut[idm.u16]; //l_1*x^(p+1), l_2*(x^(p+1))^2, ...
        }
        if (!p)
        {
            x = GFIdxFast2(256 - n + j, lut);
            uint8_t y = omega[0], s = 0;
            idx.u8.b = x;
            for (int l = 1; l <= nerr; l++)
            {
                idx.u8.a = omega[l];
                y ^= lut[idx.u16]; //Omega(X^-1)
                if (l & 1)
                {
                    idx.u8.a = lambda[l];
                    s ^= lut[idx.u16]; //Lambda'(X^-1) * X^-1
                }
                idx.u8.a = x;
                idx.u8.b = lut[idx.u16];
            }
            idx.u8.a = y;
            idx.u8.b = GFInvFast2(s, lut);
            s = lut[idx.u16];

            in[j] ^= s;
        }
    }
    memcpy_s(out, k, in, k);
    return nerr;
}
