#include <stdint.h>
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
