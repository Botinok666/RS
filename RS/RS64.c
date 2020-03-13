#include "rs64.h"
#include <string.h>

static uint8_t N = 0, K = 0;
uint8_t lut[SSE_LUT_SIZE], coefs[SSE_COEFS_SIZE];
static int (*EncodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*, uint8_t*);
static int (*DecodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*);

int InitRS(int n, int k)
{
	if (((n - k) >> 1) > MAX_T || 4 > n - k || n > 255 || k > 251)
		return -1;
	InitSSSE3(coefs, n - k, lut);
	if (GetSupportedExtensions() == AVX2_SUPPORTED)
	{
		EncodeData = &EncodeAVX2;
		DecodeData = &DecodeAVX2;
	}
	else
	{
		EncodeData = &EncodeSSSE3;
		DecodeData = &DecodeSSSE3;
	}
	N = (uint8_t)n, K = (uint8_t)k;
	return 0;
}

int EncodeRS(uint8_t* CodedData, uint8_t* ParityBuf, int n, int k)
{
	if (ParityBuf == NULL)
		return EncodeData(N, K, lut, coefs, CodedData);
	uint8_t buffer[255];
	memcpy_s(buffer, K, CodedData, K);
	int result = EncodeData(N, K, lut, coefs, buffer);
	memcpy_s(CodedData, K, buffer, K);
	memcpy_s(ParityBuf, N - K, buffer + K, N - K);
	return result;
}

int DecodeRS(uint8_t* CodedData, uint8_t* ParityBuf, int n, int k)
{
	if (ParityBuf == NULL)
		return DecodeData(N, K, lut, CodedData);
	uint8_t buffer[255];
	memcpy_s(buffer, K, CodedData, K);
	memcpy_s(buffer + K, N - K, ParityBuf, N - K);
	int result = DecodeData(N, K, lut, buffer);
	memcpy_s(CodedData, K, buffer, K);
	return result;
}
