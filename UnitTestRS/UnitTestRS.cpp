#include "pch.h"
#include "CppUnitTest.h"
#include "../RS/rsalu.h"
#include "../RS/rsalu.c"
#include <vector>
//#include "../RS/rsaluflut.h"
//#include "../RS/rsaluflut.c"
#include "../RS/rsssse3.h"
#include "../RS/rsssse3.c"
#include <random>
#include <set>

uint16_t NextRand()
{
	static std::random_device rd;
	static std::mt19937 rgen(rd());
	return rgen();
}
void IntroduceSingleError(uint8_t n, uint8_t* buffer)
{
	uint16_t randn = NextRand();
	int idx = randn % n;
	randn >>= 8; //Use upper byte for error value
	if (!randn)
		randn = 1; //We need to flip at least one bit
	buffer[idx] ^= randn;
}
void IntroduceDistributedErrors(uint8_t errorsCount, uint8_t n, uint8_t* buffer)
{
	std::set<uint8_t> st; //Set required to create exact number of errors as specified
	uint16_t randn;
	int idx;
	while (st.size() < errorsCount)
	{
		randn = NextRand();
		idx = randn % n;
		randn >>= 8; //Use upper byte for error value
		if (!randn)
			randn = 1; //We need to flip at least one bit
		if (st.insert(idx).second)
			buffer[idx] ^= (uint8_t)randn; //If idx is new in current sequence
	}
}
void IntroduceSequencedErrors(uint8_t errorsCount, uint8_t n, uint8_t* buffer)
{
	uint16_t randn = NextRand();
	int idx = randn % n;
	if (idx + errorsCount >= n)
		idx = n - errorsCount - 1;
	for (int j = 0; j < errorsCount; j++)
	{
		randn = NextRand() >> 8; //Use upper byte for error value
		if (!randn)
			randn = 1; //We need to flip at least one bit
		buffer[idx + j] ^= randn;
	}
}

using namespace Microsoft::VisualStudio::CppUnitTestFramework;
//#define MSG(msg) [&]{ std::wstringstream _s; _s << msg; return _s.str(); }().c_str()

namespace Microsoft {
	namespace VisualStudio {
		namespace CppUnitTestFramework
		{
			template<> inline std::wstring ToString<unsigned short>(const unsigned short& t) { RETURN_WIDE_STRING(t); }
		}
	}
}

namespace UnitTestRS
{
	TEST_CLASS(UnitTestRS)
	{
	public:
		
		const std::vector<std::pair<uint8_t, uint8_t>> nk = 
		{ {15, 11}, {15, 7}, {135, 129}, {135, 105}, {135, 67}, {255, 247}, {255, 187}, {255, 159} };
		const std::vector<uint8_t> sdval = { 1, 15, 88, 135, 255 };
		const int repeats = 16;

		void TestAllOnesOrZeros(uint8_t* Coefs, uint8_t* LUT, void (*FillCoeffs)(uint8_t*, uint8_t, uint8_t*),
			int (*EncodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*, uint8_t*),
			int (*DecodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*))
		{
			uint8_t buffer[255], bufferOrig[255];
			wchar_t message[200];
			for (int r = 0; r < repeats; r++)
			{
				for (auto i : nk)
				{
					size_t n = i.first, k = i.second;
					uint8_t t = (uint8_t)((n - k) / 2); //Error capacity
					FillCoeffs(Coefs, (uint8_t)(n - k), LUT);
					std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };
					std::vector<uint8_t> fillval = { 0, 0xff };

					for (auto v : fillval)
					{
						//Test with all zeros
						memset(buffer, v, k);
						int enc = EncodeData((uint8_t)n, (uint8_t)k, LUT, Coefs, buffer);
						Assert::AreEqual(0, enc, L"Encoder returned wrong result");
						memcpy(bufferOrig, buffer, n);

						memcpy(buffer, bufferOrig, n);
						int dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
						Assert::AreEqual(0, dec, L"Decoder returned wrong result");
						Assert::AreEqual(0, memcmp(buffer, bufferOrig, k), L"Data shouldn't have been changed");
						for (auto e : errnom) //Introduce correctable error patterns
						{
							memcpy(buffer, bufferOrig, n);
							IntroduceSingleError((uint8_t)n, buffer);
							dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
							swprintf_s(message, L"Decoder returned wrong result (single error)");
							Assert::AreEqual(1, dec, message);
							int comp = memcmp(buffer, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (single error)");
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceDistributedErrors(e, (uint8_t)n, buffer);
							dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
							swprintf_s(message, L"Decoder returned wrong result (random %d errors)", e);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(buffer, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (random %d errors)", e);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceSequencedErrors(e, (uint8_t)n, buffer);
							dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
							swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors)", e);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(buffer, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (sequenced %d errors)", e);
							Assert::AreEqual(0, comp, message);
						}
					}
				}
			}
		}
		void TestSingleByteData(uint8_t* Coefs, uint8_t* LUT, void (*FillCoeffs)(uint8_t*, uint8_t, uint8_t*),
			int (*EncodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*, uint8_t*),
			int (*DecodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*))
		{
			uint8_t buffer[255], bufferOrig[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity
				FillCoeffs(Coefs, (uint8_t)(n - k), LUT);
				std::vector<uint8_t> sdpos = { 0, (uint8_t)(k / 2), (uint8_t)(k - 1) };
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };

				for (auto v : sdval)
				{
					for (auto p : sdpos)
					{
						memset(buffer, 0, k); //Test with k-1 zeros
						buffer[p] = v; //Different values at different positions
						swprintf_s(positionInfo, L"n %d, k %d, val %d, pos %d", (uint8_t)n, (uint8_t)k, v, p);

						int enc = EncodeData((uint8_t)n, (uint8_t)k, LUT, Coefs, buffer);
						swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, enc, message);
						memcpy(bufferOrig, buffer, n);

						memcpy(buffer, bufferOrig, n);
						int dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
						swprintf_s(message, L"Decoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, dec, message);
						int comp = memcmp(buffer, bufferOrig, k);
						swprintf_s(message, L"Data shouldn't have been changed. %s", positionInfo);
						Assert::AreEqual(0, comp, message);

						for (auto e : errnom) //Introduce correctable error patterns
						{
							memcpy(buffer, bufferOrig, n);
							IntroduceSingleError((uint8_t)n, buffer);
							dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
							swprintf_s(message, L"Decoder returned wrong result (single error). %s", positionInfo);
							Assert::AreEqual(1, dec, message);
							comp = memcmp(buffer, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (single error). %s", positionInfo);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceDistributedErrors(e, (uint8_t)n, buffer);
							dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
							swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(buffer, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (random %d errors). %s", e, positionInfo);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceSequencedErrors(e, (uint8_t)n, buffer);
							dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
							swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(buffer, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual(0, comp, message);
						}
					}
				}
			}
		}
		void TestRandomData(uint8_t* Coefs, uint8_t* LUT, void (*FillCoeffs)(uint8_t*, uint8_t, uint8_t*),
			int (*EncodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*, uint8_t*),
			int (*DecodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*))
		{
			uint8_t buffer[255], bufferOrig[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (int r = 0; r < repeats; r++)
			{
				for (auto i : nk)
				{
					size_t n = i.first, k = i.second;
					uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

					FillCoeffs(Coefs, (uint8_t)(n - k), LUT);
					std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };

					for (int p = 0; p < k; p++)
						buffer[p] = (uint8_t)NextRand();
					swprintf_s(positionInfo, L"n %d, k %d", (uint8_t)n, (uint8_t)k);

					int enc = EncodeData((uint8_t)n, (uint8_t)k, LUT, Coefs, buffer);
					swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
					Assert::AreEqual(0, enc, message);
					memcpy(bufferOrig, buffer, n);

					memcpy(buffer, bufferOrig, n);
					int dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
					swprintf_s(message, L"Decoder returned wrong result. %s", positionInfo);
					Assert::AreEqual(0, dec, message);
					int comp = memcmp(buffer, bufferOrig, k);
					swprintf_s(message, L"Data shouldn't have been changed. %s", positionInfo);
					Assert::AreEqual(0, comp, message);

					for (auto e : errnom) //Introduce correctable error patterns
					{
						memcpy(buffer, bufferOrig, n);
						IntroduceSingleError((uint8_t)n, buffer);
						dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
						swprintf_s(message, L"Decoder returned wrong result (single error). %s", positionInfo);
						Assert::AreEqual(1, dec, message);
						comp = memcmp(buffer, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (single error). %s", positionInfo);
						Assert::AreEqual(0, comp, message);

						memcpy(buffer, bufferOrig, n);
						IntroduceDistributedErrors(e, (uint8_t)n, buffer);
						dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
						swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
						Assert::AreEqual((int)e, dec, message);
						comp = memcmp(buffer, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (random %d errors). %s", e, positionInfo);
						Assert::AreEqual(0, comp, message);

						memcpy(buffer, bufferOrig, n);
						IntroduceSequencedErrors(e, (uint8_t)n, buffer); 
						dec = DecodeData((uint8_t)n, (uint8_t)k, LUT, buffer);
						swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
						Assert::AreEqual((int)e, dec, message);
						comp = memcmp(buffer, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (sequenced %d errors). %s", e, positionInfo);
						Assert::AreEqual(0, comp, message);
					}
				}
			}
		}
		TEST_METHOD(TestAllOnesOrZerosALU)
		{
			uint8_t lut[ALU_LUT_SIZE];
			uint8_t coefs[COEFS_SIZE_ALU];
			FillRLUT(lut);
			TestAllOnesOrZeros(coefs, lut, &FillCoefficents, &Encode, &Decode);
		}

		TEST_METHOD(TestSingleByteDataALU)
		{
			uint8_t lut[ALU_LUT_SIZE];
			uint8_t coefs[COEFS_SIZE_ALU];
			FillRLUT(lut);
			TestSingleByteData(coefs, lut, &FillCoefficents, &Encode, &Decode);
		}		

		TEST_METHOD(TestRandomDataALU)
		{
			uint8_t lut[ALU_LUT_SIZE];
			uint8_t coefs[COEFS_SIZE_ALU];
			FillRLUT(lut);
			TestRandomData(coefs, lut, &FillCoefficents, &Encode, &Decode);
		}

		//TEST_METHOD(TestAllOnesOrZerosALUFL)
		//{
		//	uint8_t* flut = new uint8_t[FULL_LUT_SIZE];
		//	FillFLUT(flut);
		//	uint8_t coefs[COEFS_SIZE_FLUT];
		//	TestAllOnesOrZeros(coefs, flut, &FillCoefficentsFL, &EncodeFL, &DecodeFL);
		//	delete[] flut;
		//}

		//TEST_METHOD(TestSingleByteDataALUFL)
		//{
		//	uint8_t* flut = new uint8_t[FULL_LUT_SIZE];
		//	FillFLUT(flut);
		//	uint8_t coefs[COEFS_SIZE_FLUT];
		//	TestSingleByteData(coefs, flut, &FillCoefficentsFL, &EncodeFL, &DecodeFL);
		//	delete[] flut;
		//}
		//
		//TEST_METHOD(TestRandomDataALUFL)
		//{
		//	uint8_t* flut = new uint8_t[FULL_LUT_SIZE];
		//	FillFLUT(flut);
		//	uint8_t coefs[COEFS_SIZE_FLUT];
		//	TestRandomData(coefs, flut, &FillCoefficentsFL, &EncodeFL, &DecodeFL);
		//	delete[] flut;
		//}
		TEST_METHOD(TestAllOnesOrZerosSSSE3)
		{
			uint8_t* luts = new uint8_t[SSE_LUT_SIZE];
			FillSSELUT(luts);
			uint8_t coefs[COEFS_SIZE_SSE];
			TestAllOnesOrZeros(coefs, luts, &FillCoefficentsSSE, &EncodeSSSE3, &DecodeSSSE3);
			delete[] luts;
		}

		TEST_METHOD(TestSingleByteDataSSSE3)
		{
			uint8_t* luts = new uint8_t[SSE_LUT_SIZE];
			FillSSELUT(luts);
			uint8_t coefs[COEFS_SIZE_SSE];
			TestSingleByteData(coefs, luts, &FillCoefficentsSSE, &EncodeSSSE3, &DecodeSSSE3);
			delete[] luts;
		}
		
		TEST_METHOD(TestRandomDataSSSE3)
		{
			uint8_t* luts = new uint8_t[SSE_LUT_SIZE];
			FillSSELUT(luts);
			uint8_t coefs[COEFS_SIZE_SSE];
			TestRandomData(coefs, luts, &FillCoefficentsSSE, &EncodeSSSE3, &DecodeSSSE3);
			delete[] luts;
		}
	};
}
