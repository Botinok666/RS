#include "rsaluflut.h"
#include <stdio.h>
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
    //Required size of scratch buffer is n
    int length = n - k + 1;
    if (length < 3 || n < 4)
        return -1;
    //uint8_t* btemp = (uint8_t*)(coefs + length);
    union u16u8 idx;
    memcpy_s(out, k, in, k);
    memset(out + k, 0, n - k);
    //Remainder calculation: long division
    for (int j = 0; j < k; j++)
    {
        idx.u8.b = out[j]; //mul
        for (int l = 1; l < length; l++)
        {
            idx.u8.a = coefs[l];
            out[j + l] ^= lut[idx.u16];
        }
    }
    memcpy_s(out, k, in, k);
    return 0;
}
int __cdecl DecodeFL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* scratch, uint8_t* in, uint8_t* out)
{
    //Required size of scratch buffer is 6n-5k+4
    int scount = n - k, lcount = scount + 1;
    if (scount < 2 || n < 4)
        return -1;

    union u16u8 idx;
    //uint8_t* btemp = scratch;
    uint8_t* lambda = scratch; // btemp + n;
    uint8_t* omega = lambda + lcount;
    uint8_t* syn = omega + scount;
    memcpy_s(out, k, in, k);

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

    for (int j = 0; j < scount; j++)
    {
        uint8_t s = 0;
        idx.u8.b = GFIdxFast2(j, lut);
        for (int m = 0; m < n - 1; m++)
        {
            idx.u8.a = in[m] ^ s;
            s = lut[idx.u16];
        }
        syn[j] = s ^ in[n - 1];
    }

    //printf_s("Syndromes: ");
    //for (int j = 0; j < scount; j++)
    //    printf_s("%d ", syn[j]);
    //printf_s("\n");

    uint8_t hasErrors = 0;
    for (int j = 0; j < scount; j++)
        hasErrors |= syn[j];
    if (!hasErrors) return 0;

    //Lambda calculation: Berlekamp's method
    uint8_t* b = syn + scount, * Lm = b + lcount; //4n-3k+2 + 2n-2k+2
    //Set initial value of b and lambda to 1
    memset(b + 1, 0, (size_t)lcount - 1);
    b[0] = 1;
    memset(lambda + 1, 0, (size_t)lcount - 1);
    lambda[0] = 1;
    int l = 0;

    for (int r = 1; r <= scount; r++)
    {
        int sidx = r - 1;
        uint8_t delta = syn[sidx--];
        for (int m = 1; m <= l; m++)
            delta ^= GFMulFast2(lambda[m], syn[sidx--], lut);

        if (delta)
        {
            memcpy_s(Lm, lcount, lambda, lcount);
            idx.u8.b = delta;
            for (int m = 1; m < lcount; m++)
            {
                idx.u8.a = b[m - 1];
                lambda[m] ^= lut[idx.u16];
            }

            if (2 * l <= r - 1)
            {
                l = r - l;
                idx.u8.b = GFInvFast2(delta, lut);
                for (int m = 0; m < lcount; m++)
                {
                    idx.u8.a = Lm[m];
                    b[m] = lut[idx.u16];
                }
            }
            else
            {
                for (int m = lcount - 1; m > 0; m--) //Shift
                    b[m] = b[m - 1];
                b[0] = 0;
            }
        }
        else
        {
            for (int m = lcount - 1; m > 0; m--) //Shift
                b[m] = b[m - 1];
            b[0] = 0;
        }
    }
    int nerr = 0, t = scount / 2;
    for (int m = lcount - 1; m > 0; m--) //Find deg(lambda)
    {
        if (lambda[m])
        {
            nerr = m;
            break;
        }
    }
    if (nerr > t)
        return -nerr; //deg(lambda) > t? Uncorrectable error pattern occured
    //scratch: btemp, lambda, omega, syn

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
    //scratch: btemp, lambda, omega

    int xlength = t; //We need to analyze only t errors
    uint8_t* xterm = omega + lcount, *xmul = xterm + xlength;

    uint8_t x = GFIdxFast2(256 - n, lut);
    xterm[0] = x;
    xmul[0] = 2;
    union u16u8 idm;
    idm.u8.b = 2;
    idx.u8.b = x;
    for (int j = 1; j < xlength; j++)
    {
        idm.u8.a = xmul[j - 1];
        xmul[j] = lut[idm.u16]; //x, x^2, x^3, ...
        idx.u8.a = xterm[j - 1];
        xterm[j] = lut[idx.u16]; //x^p, (x^p)^2, (x^p)^3, ...
    }
    for (int j = 0; j < xlength; j++)
        xterm[j] = GFMulFast2(lambda[j + 1], xterm[j], lut); //l_1*x^p, l_2*(x^p)^2, ...

    for (int j = 0; j < n; j++)
    {
        in[j] = 1;
        for (int m = 0; m < xlength; m++)
        {
            idm.u8.a = xterm[m];
            in[j] ^= idm.u8.a;
            idm.u8.b = xmul[m];
            xterm[m] = lut[idm.u16]; //l_1*x^(p+1), l_2*(x^(p+1))^2, ...
        }
    }

    //printf_s("Chien search: \n");
    //for (int j = 0; j < n; j++)
    //    printf_s("%d ", btemp[j]);
    //printf_s("\n");

    for (int m = 0; m < k; m++)
    {
        if (in[m]) continue;
        x = GFIdxFast2(256 - n + m, lut);
        uint8_t y = omega[0], s = 0;
        idx.u8.b = x;
        for (int l = 1; l <= t; l++)
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

        out[m] ^= s;
    }
    
    return nerr;
}
