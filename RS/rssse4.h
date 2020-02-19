#pragma once
#include <stdint.h>
#include <string.h>
#define SSE_LUT_SIZE 8960LL
#define MAX_T_SSE 96LL
#define COEFS_SIZE_SSE (MAX_T_SSE + 1LL)
#define SCRATCH_SIZE_SSE(n, k) (5 * (n - k))
#ifdef __cplusplus
extern "C" {
#endif 
	int EncodeSSE4(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* in, uint8_t* out);
	int DecodeSSE4(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* scratch, uint8_t* in, uint8_t* out);
	void FillSSELUT(uint8_t* lut);
	void FillCoefficentsSSE(uint8_t* coefs, uint8_t count, uint8_t* lut);
#ifdef __cplusplus
}
#endif
