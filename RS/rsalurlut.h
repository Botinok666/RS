#pragma once
#include <stdint.h>
#include <string.h>
#define REDUCED_LUT_SIZE 1792
#define COEFS_SIZE_RLUT(n, k) (n - k + 1)
#define SCRATCH_SIZE_RLUT(n, k) (5 * (n - k))
#ifdef __cplusplus
extern "C" {
#endif 
int Encode(uint8_t n, uint8_t k, uint16_t* lut, uint16_t* coefs, uint8_t* in, uint8_t* out);
int Decode(uint8_t n, uint8_t k, uint16_t* lut, uint16_t* scratch, uint8_t* in, uint8_t* out);
void FillRLUT(uint16_t* lut);
void FillCoefficents(uint16_t* coefs, uint8_t count, uint16_t* lut);
#ifdef __cplusplus
}
#endif