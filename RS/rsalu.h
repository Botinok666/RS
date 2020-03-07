#pragma once
#include <stdint.h>
#include <string.h>
#define ALU_LUT_SIZE 3584LL
#define MAX_S_ALU 96LL
#define ALU_LUT_EXP_OFFSET	256LL
#define COEFS_SIZE_ALU (2 * (MAX_S_ALU + 1LL))
#ifdef __cplusplus
extern "C" {
#endif 
int Encode(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
int Decode(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
void FillRLUT(uint8_t* lut);
void FillCoefficents(uint8_t* coefs, const uint8_t count, uint8_t* lut);
#ifdef __cplusplus
}
#endif