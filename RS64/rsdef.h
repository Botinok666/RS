#pragma once
#define MAX_T 48

#define ALU_LUT_EXP_OFFSET	512LL
#define ALU_LUT_SIZE 2048LL
#define ALU_COEFS_SIZE (2 * (2 * MAX_T + 1))

#define SSE_LUT_SSE_OFFSET 2048LL
#define SSE_LUT_ROOTS_OFFSET 10240LL
#define SSE_LUT_SIZE (SSE_LUT_ROOTS_OFFSET + 256 * 255 + 64)
#define SSE_COEFS_SIZE (2 * MAX_T + 32)

#ifdef __cplusplus
extern "C" {
#endif 
	uint8_t GFMul(uint8_t a, uint8_t b);
#ifdef __cplusplus
}
#endif

static void inline mod255(int* x)
{
	*x = (*x >> 8) + (*x & 0xff); //Close analog to (x %= 255), but faster
	*x = (*x >> 8) + (*x & 0xff); //Result [255] is possible, that's acceptable here
}
