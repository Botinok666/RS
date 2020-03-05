#pragma once
#include <stdint.h>
#include <string.h>
#define FULL_LUT_SIZE (256 * 256 + 512)
#define MAX_T_ALUF 96LL
#define FL_INV_OFFSET (1 << 16)
#define FL_GFV_OFFSET ((1 << 16) + 256)
#define COEFS_SIZE_FLUT (MAX_T_ALUF + 1)
#ifdef __cplusplus
extern "C" {
#endif 
	int EncodeFL(uint8_t n, uint8_t k, void* lut, void* coefs, uint8_t* buffer);
	int DecodeFL(uint8_t n, uint8_t k, void* lut, uint8_t* buffer);
	void FillCoefficentsFL(void* coefs, uint8_t count, void* lut); 
	void FillFLUT(uint8_t* lut);
#ifdef __cplusplus
}
#endif