#pragma once
#include <stdint.h>
#include <string.h>
#define REDUCED_LUT_SIZE 1792LL
#define MAX_T_ALU 96LL
#define COEFS_SIZE_RLUT (MAX_T_ALU + 1LL)
#define SCRATCH_SIZE_RLUT(n, k) (5LL * (n - k))
#ifdef __cplusplus
extern "C" {
#endif 
int Encode(const uint8_t n, const uint8_t k, uint16_t* lut, uint16_t* coefs, uint8_t* in, uint8_t* out);
int Decode(const uint8_t n, const uint8_t k, uint16_t* lut, uint16_t* scratch, uint8_t* in, uint8_t* out);
void FillRLUT(uint16_t* lut);
void FillCoefficents(uint16_t* coefs, const uint8_t count, uint16_t* lut);
#ifdef __cplusplus
}
#endif