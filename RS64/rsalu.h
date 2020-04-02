#pragma once
#include <stdint.h>
#include <string.h>
#include "rsdef.h"
#ifdef __cplusplus
extern "C" {
#endif 
	__declspec(dllexport) void InitALU(uint8_t* Coefs, const uint8_t count, uint8_t* lut);
	__declspec(dllexport) int EncodeALU(const uint8_t n, const uint8_t k, uint8_t* LUT, uint8_t* Coefs, uint8_t* buffer);
	__declspec(dllexport) int DecodeALU(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
#ifdef __cplusplus
}
#endif
