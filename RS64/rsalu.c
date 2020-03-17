#include "rsalu.h"
//#include <stdio.h>

void InitALU(uint8_t* Coefs, const uint8_t count, uint8_t* lut)
{
    if (count < 2) return;

    uint16_t* lutExp = (uint16_t*)lut + ALU_LUT_EXP_OFFSET, * lutLog = (uint16_t*)lut, * coefs = (uint16_t*)Coefs;
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
    memset(lutExp + 511, 0, 2048);
    //Calculate coefficients
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
int EncodeALU(const uint8_t n, const uint8_t k, uint8_t* LUT, uint8_t* Coefs, uint8_t* buffer)
{
    int length = n - k + 1;
    uint16_t* lutExp = (uint16_t*)LUT + ALU_LUT_EXP_OFFSET, *lutLog = (uint16_t*)LUT;
    uint16_t* coefs = (uint16_t*)Coefs + 1;
    uint8_t btemp[255];
    memcpy_s(btemp, k, buffer, k);
    memset(btemp + k, 0, n - k);
    //Remainder calculation: long division
    uint8_t* oj = btemp;
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
int DecodeALU(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    int scount = n - k;
    uint16_t* lutExp = (uint16_t*)lut + ALU_LUT_EXP_OFFSET, * lutLog = (uint16_t*)lut;

    uint16_t hasErrors = 0;
    uint16_t lambda[2 * MAX_T], syn[2 * MAX_T];
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
                ecx = (ecx >> 8) + (ecx & 0xff); //Close analog to (ecx %= 255), but faster
                ecx = (ecx >> 8) + (ecx & 0xff); //Result [255] is possible, that's acceptable here
                s ^= lutExp[stemp[m] + ecx];
            }
            hasErrors |= s;
            syn[j] = lutLog[s];
        }
    }
    //printf_s("\nSyndromes: ");
    //for (int j = 0; j < scount; j++)
    //    printf_s("%d ", lutExp[syn[j]]);
    //printf_s("\n");
    if (!hasErrors) return 0;

    //Lambda calculation: Berlekamp's method
    uint16_t omega[2 * MAX_T], b[2 * MAX_T], Lm[2 * MAX_T];
    //Set initial value of b and lambda to 1
    memset(b, 0, scount * sizeof(uint16_t));
    b[1] = 1;
    memset(lambda, 0, scount * sizeof(uint16_t));
    lambda[0] = 1;
    int l = 0;
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
            if (2 * l < r)
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
    if (nerr > (scount >> 1))
        return -2; //deg(lambda) > t? Uncorrectable error pattern occured

    /* Omega calculation */
    //Omega must be calculated only up to {nerr} power
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

    uint16_t* xterm = syn;
    /* Chien search */
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
    int efound = 0, ecorr = 0;
    for (int j = 0; j < n; j++)
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
            if (j < k) {
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
                b[ecorr] = (uint8_t)j;
                Lm[ecorr++] = s;
            }
            efound++;
        }
    }
    if (efound != nerr)
        return -3;
    for (int j = 0; j < ecorr; j++)
        buffer[b[j]] ^= Lm[j];

    return nerr;
}
