#pragma once
#include <stdint.h>
#include <string.h>
#define FULL_LUT_SIZE (256 * 256 + 512)
#define FL_INV_OFFSET (1 << 16)
#define FL_GFV_OFFSET ((1 << 16) + 256)
#define COEFS_SIZE_FLUT(n, k) (n - k + 2 + n)
#define SCRATCH_SIZE_FLUT(n, k) (6 * n - 5 * k + 4)
#ifdef __cplusplus
extern "C" {
#endif 
	int EncodeFL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* in, uint8_t* out);
	int DecodeFL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* scratch, uint8_t* in, uint8_t* out);
	void FillCoefficentsFL(uint8_t* coefs, uint8_t count, uint8_t* lut); 
	void FillFLUT(uint8_t* lut);
#ifdef __cplusplus
}
#endif