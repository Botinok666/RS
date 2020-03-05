#pragma once
#include <stdint.h>
#include <string.h>
#define SSE_LUT_SIZE 9248LL
#define MAX_T_SSE 96LL
#define SSE_LUT_EXP_OFFSET 256LL
#define SSE_LUT_SSE_OFFSET 768LL
#define SSE_LUT_REV_OFFSET 8960LL
#define COEFS_SIZE_SSE (MAX_T_SSE + 16)
#ifdef __cplusplus
extern "C" {
#endif 
	int EncodeSSSE3(uint8_t n, uint8_t k, void* lut, void* coefs, uint8_t* buffer);
	int DecodeSSSE3(uint8_t n, uint8_t k, void* lut, uint8_t* buffer);
	void FillSSELUT(uint8_t* lut);
	void FillCoefficentsSSE(void* coefs, uint8_t count, void* lut);
#ifdef __cplusplus
}
#endif
