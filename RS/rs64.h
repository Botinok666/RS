#pragma once
#include <stdint.h>
#define MAX_T 48
#define ALU_LUT_SIZE 3584LL
#define ALU_COEFS_SIZE (2 * (2 * MAX_T + 1))
#define SSE_LUT_SIZE (8960 + 256 * 255 + 16)// 9248LL
#define SSE_COEFS_SIZE (2 * MAX_T + 16)
#ifdef __cplusplus
extern "C" {
#endif 
	__declspec(dllimport) void InitALU(uint8_t* Coefs, const uint8_t count, uint8_t* lut);
	__declspec(dllimport) int EncodeALU(const uint8_t n, const uint8_t k, uint8_t* LUT, uint8_t* Coefs, uint8_t* buffer);
	__declspec(dllimport) int DecodeALU(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
	__declspec(dllimport) int GetLatestSupportedExtension();
	__declspec(dllimport) void InitSSSE3(uint8_t* coefsu, const uint8_t count, uint8_t* lut);
	__declspec(dllimport) int EncodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
	__declspec(dllimport) int DecodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
#ifdef __cplusplus
}
#endif