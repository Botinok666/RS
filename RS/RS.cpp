// RS.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <windows.h>
#include <set>
#include "rs64.h"

#include "rs_crc_import.h"

std::ofstream results;

//Wrappers for rs (jpwl) library
int EncodeJPWL(uint8_t n, uint8_t k, uint8_t* coefs, uint8_t* lut, uint8_t* buffer)
{
    return encode_rs(buffer, buffer + k, n, k);
}
int DecodeJPWL(uint8_t n, uint8_t k, uint8_t* lut, uint8_t* buffer)
{
    return dec_rs(buffer, buffer + k, n, k);
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
    if (results.is_open())
        results << '\t' << msecs;
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
            randn >>= 8; //Upper byte as a random value
            outblock[idx] ^= randn ? randn : 1;
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
    if (t < 2 || n < 4) return;

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
    if (results.is_open())
    {
        results << "\nRS(" << (int)n << ", " << (int)k << "), file size " << (size >> 10) << "Kb\n";
        results << "Coder\tEncode, ms\tClear decode, ms\tDecode, ms";
    }

    uint8_t lut[ALU_LUT_SIZE];
    uint8_t coefs[ALU_COEFS_SIZE];
    InitALU(coefs, (uint8_t)(n - k), lut);

    uint8_t* luts = new uint8_t[SSE_LUT_SIZE];
    uint8_t coefsSSE[SSE_COEFS_SIZE];
    InitSSSE3(coefsSSE, (uint8_t)(n - k), luts);

    init_rs(n, k);
    int check;
    std::cout << std::setprecision(4);
    //if (results.is_open())
    //    results << "\nJPWL";
    //std::cout << "\nTesting JPWL\n";
    //memcpy_s(memblock, size, origblock, size); //Restore original data
    //memset(outblock, 0, size);
    //BenchmarkTests(n, k, size, memblock, outblock, nullptr, nullptr, &EncodeJPWL, &DecodeJPWL);
    //check = memcmp(memblock, origblock, size);
    //std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    //if (results.is_open())
    //    results << "\nALU";
    //std::cout << "\nTesting ALU\n";
    //BenchmarkTests(n, k, size, memblock, outblock, coefs, lut, &EncodeALU, &DecodeALU);
    //memcpy_s(memblock, size, origblock, size); //Restore original data
    //memset(outblock, 0, size);
    //check = memcmp(memblock, origblock, size);
    //std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    if (results.is_open())
        results << "\nSSSE3";
    std::cout << "\nTesting SSSE3\n";
    memcpy_s(memblock, size, origblock, size); //Restore original data
    memset(outblock, 0, size);
    BenchmarkTests(n, k, size, memblock, outblock, coefsSSE, luts, &EncodeSSSE3, &DecodeSSSE3);
    check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    if (results.is_open())
        results << "\nAVX2";
    std::cout << "\nTesting AVX2\n";
    memcpy_s(memblock, size, origblock, size); //Restore original data
    memset(outblock, 0, size);
    BenchmarkTests(n, k, size, memblock, outblock, coefsSSE, luts, &EncodeAVX2, &DecodeAVX2);
    check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    delete[] memblock, outblock, origblock, luts;
}

void TestWithManyErrors(uint8_t n, uint8_t k, uint8_t ta, int rounds)
{
    uint8_t t = (n - k) >> 1;
    t += ta;
    std::random_device rd;
    std::mt19937 rgen(rd());

    uint8_t buffer[255], orig[255], origErr[255];

    uint8_t* luts = new uint8_t[SSE_LUT_SIZE];
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
    delete[] luts;
}

uint8_t lut[ALU_LUT_SIZE];
uint8_t luts[SSE_LUT_SIZE];

int main()
{
    std::cout << "Hello World!\n";
    int ext = GetSupportedExtensions();
    if (ext == 3)
        std::cout << "AVX2 is supported\n";
    else if (ext == 1)
        std::cout << "Only SSSE3 is supported\n";
    else
        std::cout << "Shouldn't be here\n";

    //for (int xt = 2; xt <= 8; xt += 2)
    //{
    //    int n = 32;
    //    for (int xj = 1; xj <= xt / 2; xj *= 2)
    //    {
    //        std::cout << "\nt = " << xt << ", e = " << (xj + xt) << '\n';
    //        TestWithManyErrors(n, n - (xt * 2), xj, 32);
    //    }
    //}
    //return 0;

    results.open("C:\\Intel\\bench.txt", std::ios::out | std::ios::trunc);
    std::set<uint8_t> testN = { 37, 43, 45, 51, 53, 75, 85, 96 };// { 48, 64, 96, 128 }; // {40, 48, 56, 64, 80, 96, 112, 128};
    for (auto tn : testN) 
    {
        std::cout << "n = " << (int)tn << ", t = " << (((int)tn - 32) >> 1) << '\n';
        RunBenchmark(tn, 32, "C:\\Intel\\meatloaf.jpg");
    }

    //std::cout << "\nt = 32\n";
    //RunBenchmark(255, 191, "C:\\Intel\\meatloaf.jpg"); //t=32
    //std::cout << "\nt = 48\n";
    //RunBenchmark(255, 159, "C:\\Intel\\meatloaf.jpg"); //t=48
    return 0;
    
    uint8_t n = 64, k = 32;
    uint8_t orig[64], buffer[64] = {
        236,111,118,204,64,22,9,96,31,196,236,13,245,4,254,191,
        199,175,190,155,2,63,167,203,30,138,40,0,150,211,18,34,
        26,71,202,31,250,92,44,85,65,198,119,140,51,162,18,47,
        239,24,225,59,200,43,59,74,218,135,76,243,107,212,180,24};// = { 10,11,12,13,14,15,16,17,18,19,20,0,0,0,0 };
    /*memset(buffer, 0, n);
    buffer[0] = 127;*/


    uint8_t coefs[ALU_COEFS_SIZE];
    InitALU(coefs, n - k, lut);

    uint8_t coefsSL[SSE_COEFS_SIZE];
    InitSSSE3(coefsSL, n - k, luts);

    init_rs(n, k);

    InitRS(n, k);
    //EncodeRS(buffer, NULL, 0, 0);

    //EncodeSSSE3(n, k, luts, coefsSL, buffer);
    memcpy_s(orig, n, buffer, n);
    std::cout << "Original buffer: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)orig[j] << ' ';
    std::cout << std::endl;

    memcpy_s(buffer, n, orig, n);
    DecodeSSSE3(n, k, luts, buffer);
    memcpy_s(buffer, n, orig, n);
    DecodeAVX2(n, k, luts, buffer);

    std::cout << "Check remainder: " << (DecodeRS(buffer, NULL, n, k) ? "fail" : "OK") << std::endl;
    int dec;
    buffer[3] = 112;
    buffer[13] = 85;
    //buffer[27] = 167;
    //buffer[37] = 81;
    std::cout << "Buffer with errors: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    DecodeRS(buffer, NULL, 0, 0);

    std::cout << "\nCorrected: \n";
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
    //std::cout << "Buffer with errors I: \n";
    //for (int j = 0; j < n; j++)
    //    std::cout << (int)buffer[j] << ' ';
    //std::cout << std::endl;

    dec = dec_rs(buffer, buffer + k, n, k);
    //dec = DecodeFL(n, k, lutf, buffer);
    std::cout << "Found errors I: " << dec << std::endl;

    std::cout << "Corrected I: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    dec = dec_rs(buffer, buffer + k, n, k);
    //dec = DecodeFL(n, k, lutf, buffer);
    std::cout << "Found errors II: " << dec << std::endl;

    std::cout << "Corrected II: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    memcpy_s(buffer, n, orig, n);
    EncodeSSSE3(n, k, luts, coefsSL, buffer);
    std::cout << "Original buffer SSE: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    std::cout << "Check remainder: " << (DecodeAVX2(n, k, luts, buffer) ? "fail" : "OK") << std::endl;

    buffer[3] = 112;
    buffer[13] = 85;
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
