using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;

namespace Demo
{
	class RSCS
	{
		private const int MAX_T = 48;
		private const int SSSE3_SUPPORTED = 1;
		private const int AVX2_SUPPORTED = 3;
		private const long ALU_LUT_SIZE = 2048L;
		private const long ALU_COEFS_SIZE = 2 * (2 * MAX_T + 1);
		private const long SSE_LUT_SIZE = 10240 + 256 * 255 + 64;
		private const long SSE_COEFS_SIZE = 2 * MAX_T + 32;

		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int GetSupportedExtensions();
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern void InitALU(ref byte coefs, byte count, ref byte lut);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int EncodeALU(byte n, byte k, ref byte lut, ref byte Coefs, ref byte buffer);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int DecodeALU(byte n, byte k, ref byte lut, ref byte buffer);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern void InitSSSE3(ref byte coefs, byte count, ref byte lut);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int EncodeSSSE3(byte n, byte k, ref byte lut, ref byte Coefs, ref byte buffer);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int DecodeSSSE3(byte n, byte k, ref byte lut, ref byte buffer);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int EncodeAVX2(byte n, byte k, ref byte lut, ref byte Coefs, ref byte buffer);
		[DllImport("RS64.dll", CallingConvention = CallingConvention.Cdecl)]
		private static extern int DecodeAVX2(byte n, byte k, ref byte lut, ref byte buffer);

		private readonly byte n;
		private readonly byte k;
		private readonly byte[] lut, coefs;
		private delegate void Init(ref byte coefs, byte count, ref byte lut);
		private delegate int Encode(byte n, byte k, ref byte lut, ref byte Coefs, ref byte buffer);
		private delegate int Decode(byte n, byte k, ref byte lut, ref byte buffer);
		private readonly Init InitTables;
		private readonly Encode EncodeData;
		private readonly Decode DecodeData;
		public enum Coders { ALU, SSSE3, AVX2 }
		public RSCS(byte n, byte k)
		{
			int t = (n - k) >> 1;
			if (t > MAX_T || 2 > t)
				throw new ArgumentOutOfRangeException("t=(n-k)/2 must be in range [2;" + MAX_T.ToString() + "]");
			int ext = GetSupportedExtensions();
			if (ext == SSSE3_SUPPORTED || ext == AVX2_SUPPORTED)
			{
				coefs = new byte[SSE_COEFS_SIZE];
				lut = new byte[SSE_LUT_SIZE];
				InitTables = InitSSSE3;
			}
			else
			{
				coefs = new byte[ALU_COEFS_SIZE];
				lut = new byte[ALU_LUT_SIZE];
				InitTables = InitALU;
			}
			InitTables(ref coefs[0], (byte)(n - k), ref lut[0]);
			if (ext == AVX2_SUPPORTED && t >= 12)
			{
				EncodeData = EncodeAVX2;
				DecodeData = DecodeAVX2;
				Version = Coders.AVX2.ToString();
			}
			else if (ext == SSSE3_SUPPORTED || ext == AVX2_SUPPORTED)
			{
				EncodeData = EncodeSSSE3;
				DecodeData = DecodeSSSE3;
				Version = Coders.SSSE3.ToString();
			}
			else
			{
				EncodeData = EncodeALU;
				DecodeData = DecodeALU;
				Version = Coders.ALU.ToString();
			}
			this.n = n;
			this.k = k;
		}
		public RSCS(byte n, byte k, Coders forceCoder)
		{
			int t = (n - k) >> 1;
			if (t > MAX_T || 2 > t)
				throw new ArgumentOutOfRangeException("t=(n-k)/2 must be in range [2;" + MAX_T.ToString() + "]");
			int ext = GetSupportedExtensions();
			if ((forceCoder == Coders.AVX2 && ext != AVX2_SUPPORTED)
				|| (forceCoder == Coders.SSSE3 && ext == 0))
				throw new PlatformNotSupportedException("Your CPU must support these instructions");
			if (forceCoder == Coders.SSSE3 || forceCoder == Coders.AVX2)
			{
				coefs = new byte[SSE_COEFS_SIZE];
				lut = new byte[SSE_LUT_SIZE];
				InitTables = InitSSSE3;
			}
			else
			{
				coefs = new byte[ALU_COEFS_SIZE];
				lut = new byte[ALU_LUT_SIZE];
				InitTables = InitALU;
			}
			InitTables(ref coefs[0], (byte)(n - k), ref lut[0]);
			if (forceCoder == Coders.AVX2)
			{
				EncodeData = EncodeAVX2;
				DecodeData = DecodeAVX2;
			}
			else if (forceCoder == Coders.SSSE3)
			{
				EncodeData = EncodeSSSE3;
				DecodeData = DecodeSSSE3;
			}
			else
			{
				EncodeData = EncodeALU;
				DecodeData = DecodeALU;
			}
			Version = forceCoder.ToString();
			this.n = n;
			this.k = k;
		}
		public int EncodeBuffer(ref byte buffer)
		{
			int result = EncodeData(n, k, ref lut[0], ref coefs[0], ref buffer);
			if (result == -1)
			{
				InitTables(ref coefs[0], (byte)(n - k), ref lut[0]);
				result = EncodeData(n, k, ref lut[0], ref coefs[0], ref buffer);
			}
			return result;
		}
		public int DecodeBuffer(ref byte buffer)
		{
			int result = DecodeData(n, k, ref lut[0], ref buffer);
			if (result == -1)
			{
				InitTables(ref coefs[0], (byte)(n - k), ref lut[0]);
				result = DecodeData(n, k, ref lut[0], ref buffer);
			}
			return result;
		}
		public string Version { get; private set; }
	}
}
