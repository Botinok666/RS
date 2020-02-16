#pragma once
#include <stdint.h>
uint16_t NextRand();
void IntroduceSingleError(uint8_t n, uint8_t* buffer);
void IntroduceDistributedErrors(uint8_t errorsCount, uint8_t n, uint8_t* buffer);
void IntroduceSequencedErrors(uint8_t errorsCount, uint8_t n, uint8_t* buffer);
