#pragma once
#include"DEFINE.h"

class BLOCK {
public:
	double* Y;
	double* Cb;
	double* Cr;
	BLOCK(int _x, int _y) {
		Y = new double[_x * _y];
		Cb = new double[_x * _y];
		Cr = new double[_x * _y];
	}
};

typedef struct
{
	uint8_t raw;
	uint8_t code_length;
	uint16_t code_word;
} coefficients;

extern std::vector<coefficients> luma_dc_huffman_table;
extern std::vector<coefficients> chroma_dc_huffman_table;
extern std::vector<coefficients> luma_ac_huffman_table;
extern std::vector<coefficients> chroma_ac_huffman_table;
extern const uint8_t DcLuminanceCodesPerBitsize[16];
extern const uint8_t DcLuminanceValues[12];
extern const uint8_t AcLuminanceCodesPerBitsize[16];
extern const uint8_t AcLuminanceValues[162];
extern const uint8_t DcChrominanceCodesPerBitsize[16];
extern const uint8_t DcChrominanceValues[12];
extern const uint8_t AcChrominanceCodesPerBitsize[16];
extern const uint8_t AcChrominanceValues[162];

void padding_pixels(BYTE*& pixels, int &x, int &y);
void encoding(std::ofstream& out, BYTE* Data, int width, int height);
void JPEG_header_encode(std::ofstream &out, int width, int height);
void JPEG_ending_encode(std::ofstream &out);

class JPEG_Encoder
{
public:
	int width, height;
	JPEG_Encoder(std::ofstream &out, BYTE* pixels, int x, int y) {
		width = x;
		height = y;
		JPEG_header_encode(out, width, height);
		padding_pixels(pixels, width, height);
		encoding(out, pixels, width, height);
		JPEG_ending_encode(out);
	}
};