#pragma once 
#include "rs_crc_decl.h"
void init_rs(int n_rs, int k_rs);
int encode_rs(dtype *CodingData, dtype *bb, int n_rs, int k_rs);
int dec_rs(dtype *CodedData, dtype *ParityBuf, int n_rs, int k_rs);
