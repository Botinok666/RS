#pragma once
#include <stdint.h>
#include <string.h>

#define SSE_LUT_EXP_OFFSET 256LL
#define SSE_LUT_SSE_OFFSET 768LL
#define SSE_LUT_REV_OFFSET 8960LL
#define MAX_S_SSE 96LL
#define SSE_LUT_SIZE 9248LL
#define SSE_COEFS_SIZE (MAX_S_SSE + 16)
#ifdef __cplusplus
extern "C" {
#endif 
	__declspec(dllexport) int GetLatestSupportedExtension();
	__declspec(dllexport) void InitSSSE3(uint8_t* coefsu, const uint8_t count, uint8_t* lut);
	__declspec(dllexport) int EncodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
	__declspec(dllexport) int DecodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
#ifdef __cplusplus
}
#endif
