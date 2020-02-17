#include "pch.h"
#include "CppUnitTest.h"
#include "../RS/rsalurlut.h"
#include "../RS/rsalurlut.c"
#include <vector>
#include "../RS/rsaluflut.h"
#include "../RS/rsaluflut.c"
#include "common.h"
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
		
		std::vector<std::pair<uint8_t, uint8_t>> nk = 
		{ {15, 11}, {15, 7}, {135, 131}, {135, 67}, {255, 247}, {255, 187}, {255, 127} };
		std::vector<uint8_t> sdval = { 1, 15, 88, 135, 255 };

		TEST_METHOD(TestAllOnesOrZerosALU)
		{
			uint16_t lut[REDUCED_LUT_SIZE];
			FillRLUT(lut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];

			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity
				uint16_t* coefs = new uint16_t[COEFS_SIZE_RLUT(n, k)]; //Minimum size
				FillCoefficents(coefs, (uint8_t)(n - k), lut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_RLUT(n, k)];
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };
				std::vector<uint8_t> fillval = { 0, 0xff };

				for (auto v : fillval)
				{
					//Test with all zeros
					memset(buffer, v, k);
					int enc = Encode((uint8_t)n, (uint8_t)k, lut, coefs, buffer, out);
					Assert::AreEqual(0, enc, L"Encoder returned wrong result");
					memcpy(bufferOrig, out, n);

					memcpy(buffer, bufferOrig, n);
					int dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
					Assert::AreEqual(0, dec, L"Decoder returned wrong result");
					Assert::AreEqual(0, memcmp(out, bufferOrig, k), L"Data shouldn't have been changed");
					for (auto e : errnom) //Introduce correctable error patterns
					{
						memcpy(buffer, bufferOrig, n);
						IntroduceSingleError((uint8_t)n, buffer);
						dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result (single error)");
						Assert::AreEqual(1, dec, message);
						int comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (single error)");
						Assert::AreEqual(0, comp, message);

						memcpy(buffer, bufferOrig, n);
						IntroduceDistributedErrors(e, (uint8_t)n, buffer);
						dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result (random %d errors)", e);
						Assert::AreEqual((int)e, dec, message);
						comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (random %d errors)", e);
						Assert::AreEqual(0, comp, message);

						memcpy(buffer, bufferOrig, n);
						IntroduceSequencedErrors(e, (uint8_t)n, buffer);
						dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors)", e);
						Assert::AreEqual((int)e, dec, message);
						comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (sequenced %d errors)", e);
						Assert::AreEqual(0, comp, message);
					}
				}

				delete[] coefs, scratch;
			}
		}

		TEST_METHOD(TestSingleByteDataALU)
		{
			uint16_t lut[REDUCED_LUT_SIZE];
			FillRLUT(lut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

				uint16_t* coefs = new uint16_t[COEFS_SIZE_RLUT(n, k)]; //Minimum size
				FillCoefficents(coefs, (uint8_t)(n - k), lut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_RLUT(n, k)];
				std::vector<uint8_t> sdpos = { 0, (uint8_t)(k / 2), (uint8_t)(k - 1) };
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };

				for (auto v : sdval)
				{
					for (auto p : sdpos)
					{
						memset(buffer, 0, k); //Test with k-1 zeros
						buffer[p] = v; //Different values at different positions
						swprintf_s(positionInfo, L"n %d, k %d, val %d, pos %d", (uint8_t)n, (uint8_t)k, v, p);

						int enc = Encode((uint8_t)n, (uint8_t)k, lut, coefs, buffer, out);
						swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, enc, message);
						memcpy(bufferOrig, out, n);

						memcpy(buffer, bufferOrig, n);
						int dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, dec, message);
						int comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data shouldn't have been changed. %s", positionInfo);
						Assert::AreEqual(0, comp, message);

						for (auto e : errnom) //Introduce correctable error patterns
						{
							memcpy(buffer, bufferOrig, n);
							IntroduceSingleError((uint8_t)n, buffer);
							dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (single error). %s", positionInfo);
							Assert::AreEqual(1, dec, message);
							comp = memcmp(out, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (single error). %s", positionInfo);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceDistributedErrors(e, (uint8_t)n, buffer);
							dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(out, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (random %d errors). %s", e, positionInfo);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceSequencedErrors(e, (uint8_t)n, buffer);
							dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(out, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual(0, comp, message);
						}
					}
				}
				delete[] coefs, scratch;
			}
		}		
		TEST_METHOD(NegativeTestSingleByteDataALU)
		{
			uint16_t lut[REDUCED_LUT_SIZE];
			FillRLUT(lut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

				uint16_t* coefs = new uint16_t[COEFS_SIZE_RLUT(n, k)]; //Minimum size
				FillCoefficents(coefs, (uint8_t)(n - k), lut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_RLUT(n, k)];
				std::vector<uint8_t> sdpos = { 0, (uint8_t)(k / 2), (uint8_t)(k - 1) };
				std::vector<uint8_t> errbad = { (uint8_t)(t + 1), (uint8_t)(t + t / 2) };

				for (auto v : sdval)
				{
					for (auto p : sdpos)
					{
						memset(buffer, 0, k); //Test with k-1 zeros
						buffer[p] = v; //Different values at different positions
						swprintf_s(positionInfo, L"n %d, k %d, val %d, pos %d", (uint8_t)n, (uint8_t)k, v, p);

						int enc = Encode((uint8_t)n, (uint8_t)k, lut, coefs, buffer, out);
						swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, enc, message);
						memcpy(bufferOrig, out, n);

						for (auto e : errbad) //Introduce uncorrectable error patterns
						{
							memcpy(buffer, bufferOrig, n);
							IntroduceDistributedErrors(e, (uint8_t)n, buffer);
							int dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
							Assert::AreNotEqual(-2, dec, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceSequencedErrors(e, (uint8_t)n, buffer);
							dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
							Assert::AreNotEqual(-2, dec, message);
						}
					}
				}
				delete[] coefs, scratch;
			}
		}

		TEST_METHOD(TestRandomDataALU)
		{
			uint16_t lut[REDUCED_LUT_SIZE];
			FillRLUT(lut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

				uint16_t* coefs = new uint16_t[COEFS_SIZE_RLUT(n, k)]; //Minimum size
				FillCoefficents(coefs, (uint8_t)(n - k), lut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_RLUT(n, k)];
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };

				for (int p = 0; p < k; p++)
					buffer[p] = (uint8_t)NextRand();
				swprintf_s(positionInfo, L"n %d, k %d", (uint8_t)n, (uint8_t)k);

				int enc = Encode((uint8_t)n, (uint8_t)k, lut, coefs, buffer, out);
				swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
				Assert::AreEqual(0, enc, message);
				memcpy(bufferOrig, out, n);

				memcpy(buffer, bufferOrig, n);
				int dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
				swprintf_s(message, L"Decoder returned wrong result. %s", positionInfo);
				Assert::AreEqual(0, dec, message);
				int comp = memcmp(out, bufferOrig, k);
				swprintf_s(message, L"Data shouldn't have been changed. %s", positionInfo);
				Assert::AreEqual(0, comp, message);

				for (auto e : errnom) //Introduce correctable error patterns
				{
					memcpy(buffer, bufferOrig, n);
					IntroduceSingleError((uint8_t)n, buffer);
					dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
					swprintf_s(message, L"Decoder returned wrong result (single error). %s", positionInfo);
					Assert::AreEqual(1, dec, message);
					comp = memcmp(out, bufferOrig, k);
					swprintf_s(message, L"Data was not corrected (single error). %s", positionInfo);
					Assert::AreEqual(0, comp, message);

					memcpy(buffer, bufferOrig, n);
					IntroduceDistributedErrors(e, (uint8_t)n, buffer);
					dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
					swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
					Assert::AreEqual((int)e, dec, message);
					comp = memcmp(out, bufferOrig, k);
					swprintf_s(message, L"Data was not corrected (random %d errors). %s", e, positionInfo);
					Assert::AreEqual(0, comp, message);

					memcpy(buffer, bufferOrig, n);
					IntroduceSequencedErrors(e, (uint8_t)n, buffer);
					dec = Decode((uint8_t)n, (uint8_t)k, lut, scratch, buffer, out);
					swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
					Assert::AreEqual((int)e, dec, message);
					comp = memcmp(out, bufferOrig, k);
					swprintf_s(message, L"Data was not corrected (sequenced %d errors). %s", e, positionInfo);
					Assert::AreEqual(0, comp, message);
				}
					
				delete[] coefs, scratch;
			}
		}

		uint8_t flut[FULL_LUT_SIZE];

		TEST_METHOD(TestAllOnesOrZerosALUFL)
		{
			FillFLUT(flut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];

			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity
				uint8_t* coefs = new uint8_t[COEFS_SIZE_FLUT(n, k)]; //Minimum size
				FillCoefficentsFL(coefs, (uint8_t)(n - k), flut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_FLUT(n, k)];
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };
				std::vector<uint8_t> fillval = { 0, 0xff };

				for (auto v : fillval)
				{
					//Test with all zeros
					memset(buffer, v, k);
					int enc = EncodeFL((uint8_t)n, (uint8_t)k, flut, coefs, buffer, out);
					Assert::AreEqual(0, enc, L"Encoder returned wrong result");
					memcpy(bufferOrig, out, n);

					memcpy(buffer, bufferOrig, n);
					int dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
					Assert::AreEqual(0, dec, L"Decoder returned wrong result");
					Assert::AreEqual(0, memcmp(out, bufferOrig, k), L"Data shouldn't have been changed");
					for (auto e : errnom) //Introduce correctable error patterns
					{
						memcpy(buffer, bufferOrig, n);
						IntroduceSingleError((uint8_t)n, buffer);
						dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result (single error)");
						Assert::AreEqual(1, dec, message);
						int comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (single error)");
						Assert::AreEqual(0, comp, message);

						memcpy(buffer, bufferOrig, n);
						IntroduceDistributedErrors(e, (uint8_t)n, buffer);
						dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result (random %d errors)", e);
						Assert::AreEqual((int)e, dec, message);
						comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (random %d errors)", e);
						Assert::AreEqual(0, comp, message);

						memcpy(buffer, bufferOrig, n);
						IntroduceSequencedErrors(e, (uint8_t)n, buffer);
						dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors)", e);
						Assert::AreEqual((int)e, dec, message);
						comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data was not corrected (sequenced %d errors)", e);
						Assert::AreEqual(0, comp, message);
					}
				}

				delete[] coefs, scratch;
			}
		}

		TEST_METHOD(TestSingleByteDataALUFL)
		{
			FillFLUT(flut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

				uint8_t* coefs = new uint8_t[COEFS_SIZE_FLUT(n, k)]; //Minimum size
				FillCoefficentsFL(coefs, (uint8_t)(n - k), flut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_FLUT(n, k)];
				std::vector<uint8_t> sdpos = { 0, (uint8_t)(k / 2), (uint8_t)(k - 1) };
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };

				for (auto v : sdval)
				{
					for (auto p : sdpos)
					{
						memset(buffer, 0, k); //Test with k-1 zeros
						buffer[p] = v; //Different values at different positions
						swprintf_s(positionInfo, L"n %d, k %d, val %d, pos %d", (uint8_t)n, (uint8_t)k, v, p);

						int enc = EncodeFL((uint8_t)n, (uint8_t)k, flut, coefs, buffer, out);
						swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, enc, message);
						memcpy(bufferOrig, out, n);

						memcpy(buffer, bufferOrig, n);
						int dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
						swprintf_s(message, L"Decoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, dec, message);
						int comp = memcmp(out, bufferOrig, k);
						swprintf_s(message, L"Data shouldn't have been changed. %s", positionInfo);
						Assert::AreEqual(0, comp, message);

						for (auto e : errnom) //Introduce correctable error patterns
						{
							memcpy(buffer, bufferOrig, n);
							IntroduceSingleError((uint8_t)n, buffer);
							dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (single error). %s", positionInfo);
							Assert::AreEqual(1, dec, message);
							comp = memcmp(out, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (single error). %s", positionInfo);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceDistributedErrors(e, (uint8_t)n, buffer);
							dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(out, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (random %d errors). %s", e, positionInfo);
							Assert::AreEqual(0, comp, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceSequencedErrors(e, (uint8_t)n, buffer);
							dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual((int)e, dec, message);
							comp = memcmp(out, bufferOrig, k);
							swprintf_s(message, L"Data was not corrected (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual(0, comp, message);
						}
					}
				}
				delete[] coefs, scratch;
			}
		}
		TEST_METHOD(NegativeTestSingleByteDataALUFL)
		{
			FillFLUT(flut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

				uint8_t* coefs = new uint8_t[COEFS_SIZE_FLUT(n, k)]; //Minimum size
				FillCoefficentsFL(coefs, (uint8_t)(n - k), flut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_FLUT(n, k)];
				std::vector<uint8_t> sdpos = { 0, (uint8_t)(k / 2), (uint8_t)(k - 1) };
				std::vector<uint8_t> errbad = { (uint8_t)(t + 1), (uint8_t)(t + t / 2) };

				for (auto v : sdval)
				{
					for (auto p : sdpos)
					{
						memset(buffer, 0, k); //Test with k-1 zeros
						buffer[p] = v; //Different values at different positions
						swprintf_s(positionInfo, L"n %d, k %d, val %d, pos %d", (uint8_t)n, (uint8_t)k, v, p);

						int enc = EncodeFL((uint8_t)n, (uint8_t)k, flut, coefs, buffer, out);
						swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
						Assert::AreEqual(0, enc, message);
						memcpy(bufferOrig, out, n);

						for (auto e : errbad) //Introduce uncorrectable error patterns
						{
							memcpy(buffer, bufferOrig, n);
							IntroduceDistributedErrors(e, (uint8_t)n, buffer);
							int dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
							Assert::AreEqual(-(int)e, dec, message);

							memcpy(buffer, bufferOrig, n);
							IntroduceSequencedErrors(e, (uint8_t)n, buffer);
							dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
							swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
							Assert::AreEqual(-(int)e, dec, message);
						}
					}
				}
				delete[] coefs, scratch;
			}
		}

		TEST_METHOD(TestRandomDataALUFL)
		{
			FillFLUT(flut);
			uint8_t buffer[255], bufferOrig[255], out[255];
			wchar_t message[200];
			wchar_t positionInfo[50];
			for (auto i : nk)
			{
				size_t n = i.first, k = i.second;
				uint8_t t = (uint8_t)((n - k) / 2); //Error capacity

				uint8_t* coefs = new uint8_t[COEFS_SIZE_FLUT(n, k)]; //Minimum size
				FillCoefficentsFL(coefs, (uint8_t)(n - k), flut);
				uint8_t* scratch = new uint8_t[SCRATCH_SIZE_FLUT(n, k)];
				std::vector<uint8_t> errnom = { (uint8_t)(t / 2), (uint8_t)(t - 1), t };

				for (int p = 0; p < k; p++)
					buffer[p] = (uint8_t)NextRand();
				swprintf_s(positionInfo, L"n %d, k %d", (uint8_t)n, (uint8_t)k);

				int enc = EncodeFL((uint8_t)n, (uint8_t)k, flut, coefs, buffer, out);
				swprintf_s(message, L"Encoder returned wrong result. %s", positionInfo);
				Assert::AreEqual(0, enc, message);
				memcpy(bufferOrig, out, n);

				memcpy(buffer, bufferOrig, n);
				int dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
				swprintf_s(message, L"Decoder returned wrong result. %s", positionInfo);
				Assert::AreEqual(0, dec, message);
				int comp = memcmp(out, bufferOrig, k);
				swprintf_s(message, L"Data shouldn't have been changed. %s", positionInfo);
				Assert::AreEqual(0, comp, message);

				for (auto e : errnom) //Introduce correctable error patterns
				{
					memcpy(buffer, bufferOrig, n);
					IntroduceSingleError((uint8_t)n, buffer);
					dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
					swprintf_s(message, L"Decoder returned wrong result (single error). %s", positionInfo);
					Assert::AreEqual(1, dec, message);
					comp = memcmp(out, bufferOrig, k);
					swprintf_s(message, L"Data was not corrected (single error). %s", positionInfo);
					Assert::AreEqual(0, comp, message);

					memcpy(buffer, bufferOrig, n);
					IntroduceDistributedErrors(e, (uint8_t)n, buffer);
					dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
					swprintf_s(message, L"Decoder returned wrong result (random %d errors). %s", e, positionInfo);
					Assert::AreEqual((int)e, dec, message);
					comp = memcmp(out, bufferOrig, k);
					swprintf_s(message, L"Data was not corrected (random %d errors). %s", e, positionInfo);
					Assert::AreEqual(0, comp, message);

					memcpy(buffer, bufferOrig, n);
					IntroduceSequencedErrors(e, (uint8_t)n, buffer);
					dec = DecodeFL((uint8_t)n, (uint8_t)k, flut, scratch, buffer, out);
					swprintf_s(message, L"Decoder returned wrong result (sequenced %d errors). %s", e, positionInfo);
					Assert::AreEqual((int)e, dec, message);
					comp = memcmp(out, bufferOrig, k);
					swprintf_s(message, L"Data was not corrected (sequenced %d errors). %s", e, positionInfo);
					Assert::AreEqual(0, comp, message);
				}

				delete[] coefs, scratch;
			}
		}
	};
}
