#pragma once
#define MAX_T 48

#define ALU_LUT_EXP_OFFSET	256LL
#define ALU_LUT_SIZE 3584LL
#define ALU_COEFS_SIZE (2 * (2 * MAX_T + 1))

#define SSE_LUT_EXP_OFFSET 256LL
#define SSE_LUT_SSE_OFFSET 768LL
#define SSE_LUT_ROOTS_OFFSET 8960LL
#define SSE_LUT_SIZE (8960 + 256 * 255 + 64)
#define SSE_LUT_REV_OFFSET 8960LL
#define SSE_COEFS_SIZE (2 * MAX_T + 32)

#ifdef __cplusplus
extern "C" {
#endif 
	uint8_t GFMul(uint8_t a, uint8_t b);
#ifdef __cplusplus
}
#endif
