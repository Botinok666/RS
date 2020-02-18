// RS.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <random>
#include <windows.h>
#include "rsalurlut.h"
#include "rsaluflut.h"

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

// testMul: go through all possible products of two bytes
int TestMul(uint8_t* fullLut, uint16_t* reducedLut)
{
    uint16_t* lutExp = reducedLut + 256, *lutLog = reducedLut;
    for (int i = 0; i < 256; i++) 
    {
        for (int j = 0; j < 256; j++) 
        {
            uint8_t x = fullLut[i | (j << 8)];
            uint8_t y = (uint8_t)lutExp[lutLog[i] + lutLog[j]];
            if (x != y)
            {
                std::cout << "a: " << i << ", b: " << j << std::endl;
                return -1;
            }
        }
        if (i)
        {
            uint8_t inv = fullLut[(1 << 16) | i];
            if (fullLut[inv | (i << 8)] != 1)
                std::cout << "Inv a: " << i << std::endl;
        }
    }
    return 0;
}

void RunBenchmark(uint8_t n, uint8_t k, const char* filename)
{
    uint8_t t = (n - k) >> 1;
    if (t < 4 || n < 4) return;

    uint8_t* memblock = nullptr;
    int blkCount;
    std::streampos size;
    std::ifstream file(filename, std::ios::in | std::ios::binary | std::ios::ate);
    if (file.is_open())
    {
        size = file.tellg();
        blkCount = (int)(size / k) + 1;
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

    uint16_t* lut = new uint16_t[REDUCED_LUT_SIZE];
    uint8_t* flut = new uint8_t[FULL_LUT_SIZE];
    FillRLUT(lut);
    FillFLUT(flut);
    uint16_t* coefs = new uint16_t[COEFS_SIZE_RLUT(n, k)]; //Minimum size
    FillCoefficents((uint16_t*)coefs, (uint8_t)(n - k), lut);
    uint8_t* coefsFL = new uint8_t[COEFS_SIZE_FLUT(n, k)];
    FillCoefficentsFL(coefsFL, (uint8_t)(n - k), flut);
    uint16_t* scratch = new uint16_t[SCRATCH_SIZE_RLUT(n, k)];
    std::random_device rd;
    std::mt19937 rgen(rd()); 
    LARGE_INTEGER StartingTime, EndingTime, ElapsedMilliseconds;
    LARGE_INTEGER Frequency;

    std::cout << "\nTesting ALU with reduced LUT\n";
    std::cout << std::setprecision(4);
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        Encode(n, k, lut, coefs, &memblock[j * k], &outblock[j * n]); //Data expands here
    }
    QueryPerformanceCounter(&EndingTime);
    ElapsedMilliseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMilliseconds.QuadPart *= 1000;
    ElapsedMilliseconds.QuadPart /= Frequency.QuadPart;
    int msecs = ElapsedMilliseconds.LowPart;
    std::cout << "Encoding finished in " << msecs / 1000.0f << 's';
    std::cout << ", speed: " << ((size / 1024) * 1000) / msecs << " Kb/s\n";

    int errorSum = 0;
    //Introduce some errors
    for (int j = 0; j < blkCount; j++)
    {
        uint8_t errors = rgen() % t;
        for (uint8_t m = 0; m < errors; m++)
        {
            uint16_t randn = rgen();
            int idx = j * n + randn % n; //Random position within a block
            outblock[idx] = randn >> 8; //Upper byte as a random value
        }
        errorSum += errors;
    }
    std::cout << "Introduced approx. " << errorSum << " errors\n";

    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        Decode(n, k, lut, scratch, &outblock[j * n], &memblock[j * k]); //Data shrinks here
    }    
    QueryPerformanceCounter(&EndingTime);
    ElapsedMilliseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMilliseconds.QuadPart *= 1000;
    ElapsedMilliseconds.QuadPart /= Frequency.QuadPart;
    msecs = ElapsedMilliseconds.LowPart;
    std::cout << "Decoding finished in " << msecs / 1000.0f << 's';
    std::cout << ", speed: " << ((size / 1024) * 1000) / msecs << " Kb/s\n";

    int check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    std::cout << "\nTesting ALU with full LUT\n";
    memcpy_s(memblock, size, origblock, size); //Restore original data
    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        EncodeFL(n, k, flut, coefsFL, &memblock[j * k], &outblock[j * n]); //Data expands here
    }
    QueryPerformanceCounter(&EndingTime);
    ElapsedMilliseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMilliseconds.QuadPart *= 1000;
    ElapsedMilliseconds.QuadPart /= Frequency.QuadPart;
    msecs = ElapsedMilliseconds.LowPart;
    std::cout << "Encoding finished in " << msecs / 1000.0f << 's';
    std::cout << ", speed: " << ((size / 1024) * 1000) / msecs << " Kb/s\n";

    errorSum = 0;
    //Introduce some errors
    for (int j = 0; j < blkCount; j++)
    {
        uint8_t errors = rgen() % t;
        for (uint8_t m = 0; m < errors; m++)
        {
            uint16_t randn = rgen();
            int idx = j * n + randn % n; //Random position within a block
            outblock[idx] = randn >> 8; //Upper byte as a random value
        }
        errorSum += errors;
    }
    std::cout << "Introduced approx. " << errorSum << " errors\n";

    QueryPerformanceFrequency(&Frequency);
    QueryPerformanceCounter(&StartingTime);
    for (int j = 0; j < blkCount; j++)
    {
        DecodeFL(n, k, flut, (uint8_t*)scratch, &outblock[j * n], &memblock[j * k]); //Data shrinks here
    }    
    QueryPerformanceCounter(&EndingTime);
    ElapsedMilliseconds.QuadPart = EndingTime.QuadPart - StartingTime.QuadPart;
    ElapsedMilliseconds.QuadPart *= 1000;
    ElapsedMilliseconds.QuadPart /= Frequency.QuadPart;
    msecs = ElapsedMilliseconds.LowPart;
    std::cout << "Decoding finished in " << msecs / 1000.0f << 's';
    std::cout << ", speed: " << ((size / 1024) * 1000) / msecs << " Kb/s\n";

    check = memcmp(memblock, origblock, size);
    std::cout << "Data integrity check " << (check ? "failed\n" : "OK\n");

    delete[] lut, flut, coefs, scratch, memblock, outblock;
}

uint16_t lut[REDUCED_LUT_SIZE];
uint8_t flut[FULL_LUT_SIZE];

int main()
{
    std::cout << "Hello World!\n";
    RunBenchmark(255, 223, "C:\\Intel\\alisa.jpg"); //t=16
    RunBenchmark(255, 191, "C:\\Intel\\alisa.jpg"); //t=32
    RunBenchmark(255, 159, "C:\\Intel\\alisa.jpg"); //t=48
    return 0;

    FillRLUT(lut);
    FillFLUT(flut);
    
    //if (TestMul(flut, lut))
    //    std::cout << "Multiplication test failed\n";
    //else
    //    std::cout << "Multiplication test OK\n";
    uint8_t n = 15, k = 11;
    //uint8_t buffer[15] = { 11, 12, 13, 14, 15, 16, 17, 0, 0, 0, 0, 0, 0, 0, 0 }; //15, 7
    uint8_t buffer[15] = { 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint8_t buffer2[15] = { 0, 0, 0, 0, 0, 255, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    uint16_t coefs[COEFS_SIZE_RLUT(15, 11)];
    FillCoefficents((uint16_t*)coefs, n - k, lut);
    uint16_t scratch[SCRATCH_SIZE_RLUT(15, 11)];
    //memcpy_s(scratch, COEFS_SIZE_RLUT(15, 11), coefs, COEFS_SIZE_RLUT(15, 11));

    uint8_t coefsFL[COEFS_SIZE_FLUT(15, 11)];
    FillCoefficentsFL(coefsFL, n - k, flut);

    uint8_t out[15], orig[15];
    EncodeFL(n, k, flut, coefsFL, buffer, out);
    memcpy_s(orig, 15, out, 15);
    std::cout << "Original buffer: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)orig[j] << ' ';
    std::cout << std::endl;

    memcpy_s(buffer, 15, orig, 15);
    //std::cout << "Check remainder: " << (Decode(n, k, lut, scratch, buffer) ? "fail" : "OK") << std::endl;
    std::cout << "Check remainder: " << (DecodeFL(n, k, flut, (uint8_t*)scratch, buffer, out) ? "fail" : "OK") << std::endl;

    memcpy_s(buffer, 15, orig, 15);
    buffer[0] = 14;
    buffer[3] = 14; //112 - S1=0, 167 - S2=0, 81 - S3=0
    std::cout << "Buffer with errors I: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)buffer[j] << ' ';
    std::cout << std::endl;

    //int dec = Decode(n, k, lut, scratch, buffer);
    int dec = DecodeFL(n, k, flut, (uint8_t*)scratch, buffer, out);
    std::cout << "Found errors I: " << dec << std::endl;

    std::cout << "Corrected I: \n";
    for (int j = 0; j < k; j++)
        std::cout << (int)out[j] << ' ';
    std::cout << std::endl;

    //buffer[2] = ~buffer[2];
    //buffer[4] = ~buffer[4];
    //buffer[12] = ~buffer[12];
    Encode(n, k, lut, coefs, buffer2, out);
    std::cout << "Original buffer V: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)out[j] << ' ';
    std::cout << std::endl;
    out[4] = 42;
    //out[6] = ~out[6];
    out[12] = 105; //105
    std::cout << "Buffer with errors V: \n";
    for (int j = 0; j < n; j++)
        std::cout << (int)out[j] << ' ';
    std::cout << std::endl;

    memcpy(buffer2, out, n);
    dec = Decode(n, k, lut, scratch, buffer2, out);
    std::cout << "Found errors V: " << dec << std::endl;

    std::cout << "Corrected V: \n";
    for (int j = 0; j < k; j++)
        std::cout << (int)out[j] << ' ';
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
