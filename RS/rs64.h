#pragma once
#include <stdint.h>
//Константы зависят от заданных в rsdef.h, изменение только в одном месте приведёт к ошибкам
#define MAX_T 48
#define AVX2_SUPPORTED 3
#define SSSE3_SUPPORTED 1
#define ALU_LUT_SIZE 3584LL
#define ALU_COEFS_SIZE (2 * (2 * MAX_T + 1))
#define SSE_LUT_SIZE (11776 + 256 * 255 + 64)
#define SSE_COEFS_SIZE (2 * MAX_T + 32)
#ifdef __cplusplus
extern "C" {
#endif 
	//Функции из DLL
	/* ALU версия "для порядка", т.к. поддержка инструкций SSSE3 начинается с:
	 * Intel Core 2 Duo (2006), AMD Bulldozer (2011)
	 * Использует свои таблицы, перед использованием они должны быть заполнены InitALU */
	__declspec(dllimport) void InitALU(uint8_t* Coefs, const uint8_t count, uint8_t* lut);
	__declspec(dllimport) int EncodeALU(const uint8_t n, const uint8_t k, uint8_t* LUT, uint8_t* Coefs, uint8_t* buffer);
	__declspec(dllimport) int DecodeALU(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
	/* Определение поддержки расширений
	 * Возвращает 1 для SSSE3 и 3 для AVX2 */
	__declspec(dllimport) int GetSupportedExtensions();
	/* SSSE3 версия использует расширенную таблицу LUT (~74Кб) для применения умножения типа "вектор * скаляр"
	 * При инициализации производится выравнивание данных в таблицах */
	__declspec(dllimport) void InitSSSE3(uint8_t* coefsu, const uint8_t count, uint8_t* lut);
	__declspec(dllimport) int EncodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
	__declspec(dllimport) int DecodeSSSE3(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);
	/* AVX2 версия использует те же таблицы, что и SSSE3. Поддержка начинается с:
	 * Intel Haswell (2013), AMD Excavator (2015) */
	__declspec(dllimport) int EncodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* coefs, uint8_t* buffer);
	__declspec(dllimport) int DecodeAVX2(const uint8_t n, const uint8_t k, uint8_t* lut, uint8_t* buffer);

	//Функции из RS64.c
	/* InitRS заполняет таблицы и выбирает версию кодера/декодера
	 * Однако, некоторые процессоры (например, AMD Zen) эмулируют 256-бит инструкции на 128-бит блоках
	 * Поэтому на самом деле разницы по скорости между SSSE3 и AVX2 может и не быть
	 * Функция возвращает 0, если n и k подходят, и -1 в противном случае */
	int InitRS(int n, int k);
	/* EncodeRS должна вызываться только после InitRS, иначе будет исключение NPE 
	 * Желательно хранить данные и контрольные байты в одном массиве размером n
	 * Тогда указатель на такой массив передавать первым, вторым передавать NULL
	 * Это позволит избежать бестолковых копирований данных туда-сюда
	 * Функция возвращает:
	 *	0, если всё в порядке
	 *	-1, если обнаружена ошибка выравнивания таблиц
	 * Передаваемые числа n и k игнорируются (используются сохранённые при инициализации) */
	int EncodeRS(uint8_t* CodedData, uint8_t* ParityBuf, int n, int k);
	/* DecodeRS должна вызываться только после InitRS, иначе будет исключение NPE
	 * Желательно хранить данные и контрольные байты в одном массиве размером n
	 * Тогда указатель на такой массив передавать первым, вторым передавать NULL
	 * Это позволит избежать бестолковых копирований данных туда-сюда
	 * Функция возвращает:
	 *	[1;(n-k)/2] - количество найденных ошибок
	 *	0, если в кодовом слове нет ошибок
	 *	-1, если обнаружена ошибка выравнивания таблиц
	 *	-2, если степень лямбда-функции получилась выше, чем разрешающая способность кода
	 *	-3, если количество найденных ошибок не соответствует степени лямбда-функции
	 *	-4, если значение какой-либо ошибки получается равным 0
	 * Ошибки [-2;-4] сообщают о превышении разрешающей способности кода, данные в этом случае никак не меняются
	 * Передаваемые числа n и k игнорируются (используются сохранённые при инициализации) */
	int DecodeRS(uint8_t* CodedData, uint8_t* ParityBuf, int n, int k);
	/* Ошибка выравнивания таблиц не должна возникать при нормальной работе
	 * Она может теоретически возникнуть, если использовать библиотеку с managed языком C#
	 * Пока это никак не проверялось */
#ifdef __cplusplus
}
#endif