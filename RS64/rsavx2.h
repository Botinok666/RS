#pragma once
#include <stdint.h>
#include <string.h>
#include "rsdef.h"

#ifdef __cplusplus
extern "C" {
#endif 
	__declspec(dllexport) int EncodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
	__declspec(dllexport) int DecodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
#ifdef __cplusplus
}
#endif
