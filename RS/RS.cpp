// RS.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <windows.h>
#include <set>
#include "rs64.h"

#include "rs-jp.h"
#include <intrin.h>

//Wrappers for rs (jpwl) library
int EncodeJPWL(uint8_t n, uint8_t k, uint8_t* coefs, uint8_t* lut, uint8_t* buffer)
{
    return encode_rs(buffer, buffer + k);
}
int DecodeJPWL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    return eras_dec_rs(buffer, nullptr, 0);
}

void PrintLUT(uint16_t* lut)
{
    std::cout << "Log table: \n";
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
            std::cout << *lut++ << ' ';
        std::cout << std::endl;
    }
    std::cout << "Exp table: \n";
    for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
            std::cout << *lut++ << ' ';
        std::cout << std::endl;
    }
}

void PrintStats(_LARGE_INTEGER start, _LARGE_INTEGER end, _LARGE_INTEGER freq, size_t size, const char* testType)
{
    _LARGE_INTEGER elapsed;
    elapsed.QuadPart = end.QuadPart - start.QuadPart;
    elapsed.QuadPart *= 1000;
    elapsed.QuadPart /= freq.QuadPart;
    int msecs = elapsed.LowPart;
    std::cout << testType << " finished in " << msecs / 1000.0f << 's';
    std::cout << ", speed: " << ((size / 1024) * 1000) / msecs << " Kb/s\n";
}
void BenchmarkTests(uint8_t n, uint8_t k, size_t size, uint8_t* memblock, uint8_t* outblock, uint8_t* coefs, uint8_t* lut,
    int (*EncodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*, uint8_t*),
    int (*DecodeData)(uint8_t, uint8_t, uint8_t*, uint8_t*))
{
    std::random_device rd;
    std::mt19937 rgen(rd());
    _LARGE_INTEGER StartingTime, EndingTime, Frequency;
    long long blkCount = size / k + 1;
    uint8_t t = (n - k) >> 1;
    int err = 0;

    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        memcpy_s(&outblock[j * n], k, &memblock[j * k], k);
        EncodeData(n, k, lut, coefs, &outblock[j * n]); //Data expands here
    }
    QueryPerformanceCounter(&EndingTime);
    PrintStats(StartingTime, EndingTime, Frequency, size, "Encoding");
    memset(memblock, 0, size);
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        if (DecodeData(n, k, lut, &outblock[j * n]) != 0)
            err++;
        memcpy_s(&memblock[j * k], k, &outblock[j * n], k); //Data shrinks here
    }
    QueryPerformanceCounter(&EndingTime);
    PrintStats(StartingTime, EndingTime, Frequency, size, "Clear decoding");
    std::cout << "Errors: " << err << std::endl;

    int errorSum = 0;
    //Introduce some errors
    for (int j = 0; j < blkCount; j++)
    {
        uint8_t errors = rgen() % (t + 1);
        for (uint8_t m = 0; m < errors; m++)
        {
            uint16_t randn = rgen();
            int idx = j * n + randn % n; //Random position within a block
            outblock[idx] = randn >> 8; //Upper byte as a random value
        }
        errorSum += errors;
    }
    std::cout << "Introduced approx. " << errorSum << " errors\n";

    err = 0;
    memset(memblock, 0, size);
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        if (DecodeData(n, k, lut, &outblock[j * n]) < 0)
            err++;
        memcpy_s(&memblock[j * k], k, &outblock[j * n], k); //Data shrinks here
    }
    QueryPerformanceCounter(&EndingTime);
    PrintStats(StartingTime, EndingTime, Frequency, size, "Decoding");
    std::cout << "Errors: " << err << std::endl;
}

void RunBenchmark(uint8_t n, uint8_t k, const char* filename)
{
    uint8_t t = (n - k) >> 1;
    if (t < 4 || n < 4) return;

    uint8_t* memblock = nullptr;
    long long blkCount;
    std::streampos size;
    std::ifstream file(filename, std::ios::in | std::ios::binary | std::ios::ate);
    if (file.is_open())
    {
        size = file.tellg();
        blkCount = size / k + 1;
        memblock = new uint8_t[blkCount * n]; //This array will be used for decoding
        file.seekg(0, std::ios::beg);
        file.read((char*)memblock, size);
        file.close();
    }
    else return;
    if (memblock == nullptr)
        return;

    uint8_t* outblock = new uint8_t[blkCount * n];
    uint8_t* origblock = new uint8_t[blkCount * n];
    memcpy_s(origblock, size, memblock, size); //Save original data
    std::cout << "File opened, size " << size << std::endl;

    uint8_t lut[ALU_LUT_SIZE];
    uint8_t coefs[ALU_COEFS_SIZE];
    InitALU(coefs, (uint8_t)(n - k), lut);

    uint8_t luts[SSE_LUT_SIZE];
    uint8_t coefsSSE[SSE_COEFS_SIZE];
    InitSSSE3(coefsSSE, (uint8_t)(n - k), luts);

    std::cout << "\nTesting ALU\n";
    std::cout << std::setprecision(4);
    BenchmarkTests(n, k, size, memblock, outblock, coefs, lut, &EncodeALU, &DecodeALU);
    int check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    std::cout << "\nTesting JPWL\n";
    init_rs(k);
    memcpy_s(memblock, size, origblock, size); //Restore original data
    memset(outblock, 0, size);
    BenchmarkTests(n, k, size, memblock, outblock, nullptr, nullptr, &EncodeJPWL, &DecodeJPWL);
    check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    std::cout << "\nTesting SSSE3\n";
    memcpy_s(memblock, size, origblock, size); //Restore original data
    memset(outblock, 0, size);
    BenchmarkTests(n, k, size, memblock, outblock, coefsSSE, luts, &EncodeSSSE3, &DecodeSSSE3);
    check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    delete[] memblock, outblock, origblock;
}

void TestWithManyErrors(uint8_t n, uint8_t k, uint8_t ta, int rounds)
{
    uint8_t t = (n - k) >> 1;
    t += ta;
    std::random_device rd;
    std::mt19937 rgen(rd());

    uint8_t buffer[255], orig[255], origErr[255];

    uint8_t luts[SSE_LUT_SIZE];
    uint8_t coefsSSE[SSE_COEFS_SIZE];
    InitSSSE3(coefsSSE, (uint8_t)(n - k), luts);

    int ecount[53]; //For returned values in range [-4;48]
    memset(ecount, 0, 53 * sizeof(int));
    for (int m = 0; m < rounds; m++)
    {
        for (int j = 0; j < k; j++) //Generate random data
            buffer[j] = (uint8_t)rgen();
        EncodeSSSE3(n, k, luts, coefsSSE, buffer);
        memcpy_s(orig, n, buffer, n);

        std::set<uint8_t> st; //Set required to create exact number of errors as specified
        while (st.size() < t)
        {
            uint16_t randn = rgen();
            int idx = randn % n;
            randn >>= 8; //Use upper byte for error value
            if (!randn)
                randn = 1; //We need to flip at least one bit
            if (st.insert(idx).second)
                buffer[idx] ^= (uint8_t)randn; //If idx is new in current sequence
        }
        memcpy_s(origErr, n, buffer, n);

        for (auto pos : st)
        {
            for (int i = 0; i < 256; i++)
            {
                if (i == orig[pos]) continue; //Skip original value to not alter the number of errors
                buffer[pos] = (uint8_t)i;
                int result = DecodeSSSE3(n, k, luts, buffer);
                ecount[result + 4]++;
                memcpy_s(buffer, n, origErr, n);
            }
        }
    }
    int totalRounds = 0;
    for (int j = 0; j < 53; j++)
        totalRounds += ecount[j];
    std::cout << "Rounds " << totalRounds << "\nStatistics:\n" << std::setprecision(4);
    double erdet = 0;
    for (int j = 0; j < 53; j++)
    {
        if (!ecount[j]) continue;
        double er = ecount[j] * 100.0 / totalRounds;
        std::cout << (j - 4) << ": " << er << '%' << '\n';
        if (j < 4) 
            erdet += er;
    }
    std::cout << "Total detected errors: " << erdet << '%' << '\n';
}

uint8_t lut[ALU_LUT_SIZE];
uint8_t luts[SSE_LUT_SIZE];

int main()
{
    std::cout << "Hello World!\n";
    int ext = GetLatestSupportedExtension();
    if (ext == 2)
        std::cout << "AVX2 is supported\n";
    else if (ext == 1)
        std::cout << "Only SSSE3 is supported\n";
    else
        std::cout << "Shouldn't be here\n";

    //for (int xt = 8; xt <= 32; xt += 8)
    //{
    //    for (int xj = 1; xj <= xt / 2; xj *= 2)
    //    {
    //        std::cout << "\nt = " << xt << ", e = " << (xj + xt) << '\n';
    //        TestWithManyErrors(255, 255 - (xt * 2), xj, 16);
    //    }
    //}
    //return 0;

    std::cout << "t = 16\n";
    RunBenchmark(255, 223, "C:\\Intel\\meatloaf.jpg"); //t=16
    std::cout << "\nt = 32\n";
    RunBenchmark(255, 191, "C:\\Intel\\meatloaf.jpg"); //t=32
    std::cout << "\nt = 48\n";
    RunBenchmark(255, 159, "C:\\Intel\\meatloaf.jpg"); //t=48
    return 0;
    
    uint8_t n = 255, k = 249;
    //uint8_t buffer[15] = { 11, 12, 13, 14, 15, 16, 17, 0, 0, 0, 21, 0, 0, 0, 0 }; //15, 11
    uint8_t buffer[255];
    memset(buffer, 0, n);
    buffer[0] = 127;
    //uint8_t buffer2[15] = { 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t coefs[ALU_COEFS_SIZE];
    InitALU(coefs, n - k, lut);

    uint8_t coefsSL[SSE_COEFS_SIZE];
    InitSSSE3(coefsSL, n - k, luts);

    init_rs(k);

    uint8_t orig[255];
    EncodeSSSE3(n, k, luts, coefsSL, buffer);
    memcpy_s(orig, n, buffer, n);
    std::cout << "Original buffer: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)orig[j] << ' ';
    std::cout << std::endl;

    memcpy_s(buffer, n, orig, n);
    std::cout << "Check remainder: " << (DecodeALU(n, k, lut, buffer) ? "fail" : "OK") << std::endl;
    int dec;

    memset(buffer, 0, n);
    buffer[0] = 127;
    encode_rs(buffer, buffer + k);
    std::cout << "rs-jp buffer: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;
    //
    //memcpy_s(orig, n, buffer, n);
    //for (int i = 1; i < 256; i++)
    //{
    //    buffer[9] = 105; //112 - S1=0, 167 - S2=0, 81 - S3=0
    //    buffer[17] = 112;
    //    buffer[27] = (uint8_t)i;
    //    buffer[250] = 81;
    //    //std::cout << "Buffer with errors I: \n";
    //    //for (int j = 0; j < n; j++)
    //    //    std::cout << (int)buffer[j] << ' ';
    //    //std::cout << std::endl;

    //    dec = eras_dec_rs(buffer, nullptr, 0);
    //    //dec = DecodeFL(n, k, lutf, buffer);
    //    if (dec >= 0)
    //        std::cout << "i/res: " << i << '/' << dec << ' ';
    //    //std::cout << "Found errors I: " << dec << std::endl;
    //    memcpy_s(buffer, n, orig, n);
    //}        

    buffer[9] = 105; //112 - S1=0, 167 - S2=0, 81 - S3=0
    buffer[17] = 112;
    buffer[27] = 4;
    buffer[250] = 81;
    //std::cout << "Buffer with errors I: \n";
    //for (int j = 0; j < n; j++)
    //    std::cout << (int)buffer[j] << ' ';
    //std::cout << std::endl;

    dec = eras_dec_rs(buffer, nullptr, 0);
    //dec = DecodeFL(n, k, lutf, buffer);
    std::cout << "Found errors I: " << dec << std::endl;

    std::cout << "Corrected I: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    dec = eras_dec_rs(buffer, nullptr, 0);
    //dec = DecodeFL(n, k, lutf, buffer);
    std::cout << "Found errors II: " << dec << std::endl;

    std::cout << "Corrected II: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    memcpy_s(buffer, n, orig, n);
    EncodeSSSE3(n, k, luts, coefsSL, buffer);
    std::cout << "Original buffer V: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    std::cout << "Check remainder: " << (DecodeSSSE3(n, k, luts, buffer) ? "fail" : "OK") << std::endl;

    //buffer[0] = 14;
    buffer[9] = 105; //105
    std::cout << "Buffer with errors V: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    dec = DecodeSSSE3(n, k, luts, buffer);
    std::cout << "Found errors V: " << dec << std::endl;

    std::cout << "Corrected V: \n";
    for (int j = 0; j < k; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;
    return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
