#pragma once
#include <stdint.h>
#include <string.h>
#define REDUCED_LUT_SIZE 3584LL
#define MAX_T_ALU 96LL
#define COEFS_SIZE_RLUT (2 * (MAX_T_ALU + 1LL))
#ifdef __cplusplus
extern "C" {
#endif 
int Encode(const uint8_t n, const uint8_t k, void* lut, void* coefs, uint8_t* buffer);
int Decode(const uint8_t n, const uint8_t k, void* lut, uint8_t* buffer);
void FillRLUT(void* lut);
void FillCoefficents(void* coefs, const uint8_t count, void* lut);
#ifdef __cplusplus
}
#endif