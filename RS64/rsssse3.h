#pragma once
#include <stdint.h>
#include <string.h>
#include "rsdef.h"

#ifdef __cplusplus
extern "C" {
#endif 
	__declspec(dllexport) int GetSupportedExtensions();
	__declspec(dllexport) void InitSSSE3(uint8_t* coefsu, const uint8_t count, uint8_t* lut);
	__declspec(dllexport) int EncodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
	__declspec(dllexport) int DecodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
#ifdef __cplusplus
}
#endif
