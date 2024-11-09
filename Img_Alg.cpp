#include "Img_Alg.h"
#include "DEFINE.h"
#include <cmath>
#include <cstring>
#include<iostream>

//std::ofstream out_log_dyx("dyx_log.txt", std::ios::binary | std::ios::out);

void padding_pixels(BYTE*& pixels, int &x, int &y) {
	int new_x = x;
	int new_y = y;
	if (x % 8 != 0) new_x = (x / 8 + 1) * 8;
	if (y % 8 != 0) new_y = (y / 8 + 1) * 8;
	BYTE* new_pixels = new BYTE[new_x * new_y * 3];
	memset(new_pixels, 0, new_x * new_y * 3);
	for (int i = 0; i < new_y; i++) {
		for (int j = 0; j < new_x; j++) {
            if (i < y && j < x) {
                new_pixels[(i * new_x + j) * 3] = pixels[(i * x + j) * 3];
                new_pixels[(i * new_x + j) * 3 + 1] = pixels[(i * x + j) * 3 + 1];
                new_pixels[(i * new_x + j) * 3 + 2] = pixels[(i * x + j) * 3 + 2];
            }
            else if (i < y) {
                new_pixels[(i * new_x + j) * 3] = pixels[(i * x + x - 1) * 3];
				new_pixels[(i * new_x + j) * 3 + 1] = pixels[(i * x + x - 1) * 3 + 1];
				new_pixels[(i * new_x + j) * 3 + 2] = pixels[(i * x + x - 1) * 3 + 2];
            }else if (j < x) {
                new_pixels[(i * new_x + j) * 3] = pixels[((y - 1) * x + j) * 3];
                new_pixels[(i * new_x + j) * 3 + 1] = pixels[((y - 1) * x + j) * 3 + 1];
                new_pixels[(i * new_x + j) * 3 + 2] = pixels[((y - 1) * x + j) * 3 + 2];
            }
            else {
    			new_pixels[(i * new_x + j) * 3] = pixels[((y - 1) * x + x - 1) * 3];
				new_pixels[(i * new_x + j) * 3 + 1] = pixels[((y - 1) * x + x - 1) * 3 + 1];
				new_pixels[(i * new_x + j) * 3 + 2] = pixels[((y - 1) * x + x - 1) * 3 + 2];
            }
        }
	}
    x = new_x;
    y = new_y;
	delete[] pixels;
	pixels = new_pixels;
}

double rgb2Y(BYTE r, BYTE g, BYTE b) {
	return 0.299 * r + 0.587 * g + 0.114 * b - 128;
}
double rgb2Cb(BYTE r, BYTE g, BYTE b) {
	return -0.16874 * r - 0.33126 * g + 0.5 * b;
}
double rgb2Cr(BYTE r, BYTE g, BYTE b) {
	return 0.5 * r - 0.41869 * g - 0.08131 * b;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define PI M_PI
#define DCT_size 8

// Standard Quantization Table for Y component
std::vector<int> default_luma_table =
{
    16, 11, 10, 16,  24,  40,  51,  61,
    12, 12, 14, 19,  26,  58,  60,  55,
    14, 13, 16, 24,  40,  57,  69,  56,
    14, 17, 22, 29,  51,  87,  80,  62,
    18, 22, 37, 56,  68, 109, 103,  77,
    24, 35, 55, 64,  81, 104, 113,  92,
    49, 64, 78, 87, 103, 121, 120, 101,
    72, 92, 95, 98, 112, 100, 103,  99,
};

// Standard Quantization Table for Cb and Cr components
std::vector<int> default_chroma_table =
{
    17, 18, 24, 47, 99, 99, 99, 99,
    18, 21, 26, 66, 99, 99, 99, 99,
    24, 26, 56, 99, 99, 99, 99, 99,
    47, 66, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
    99, 99, 99, 99, 99, 99, 99, 99,
};

#define qt 50
int get_cur(std::vector<int>& default_table, int k) {
    double ret = default_table[k];
    ret = ret * ((qt < 50) ? (50.0 / qt) : (2 - qt / 50.0));
    ret = (ret < 1) ? 1 : ret;
    ret = (ret > 255) ? 255 : ret;
    return (int)ret;
}

BYTE luma_table[DCT_size * DCT_size];
BYTE chroma_table[DCT_size * DCT_size];
std::vector<int> zigzag = {
    0, 1, 5, 6, 14, 15, 27, 28,
    2, 4, 7, 13, 16, 26, 29, 42,
    3, 8, 12, 17, 25, 30, 41, 43,
    9, 11, 18, 24, 31, 40, 44, 53,
    10, 19, 23, 32, 39, 45, 52, 54,
    20, 22, 33, 38, 46, 51, 55, 60,
    21, 34, 37, 47, 50, 56, 59, 61,
    35, 36, 48, 49, 57, 58, 62, 63
};
std::vector<int> zigzag_inv = {
    0, 1, 8, 16, 9, 2, 3, 10,
    17, 24, 32, 25, 18, 11, 4, 5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13, 6, 7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63

};

double ck(int k){
    return (k == 0) ? sqrt(1.0 / DCT_size) : sqrt(2.0 / DCT_size);
}

void DCT(double* in)
{
    double sum = 0.0;
    double out[DCT_size * DCT_size];
    double cos_table[DCT_size][DCT_size];

    for (int i = 0; i < DCT_size; i++) {
        for (int j = 0; j < DCT_size; j++) {
            cos_table[i][j] = cos((2 * i + 1) * j * M_PI / (2 * DCT_size));
        }
    }

    for (int m = 0; m < DCT_size; m++) {
        for (int n = 0; n < DCT_size; n++) {
            for (int i = 0; i < DCT_size; i++) {
                for (int j = 0; j < DCT_size; j++) {
                    sum += in[i * DCT_size + j] * cos_table[j][n] * cos_table[i][m];
                }
            }
            out[m * DCT_size + n] = sum * ck(m) * ck(n);
            sum = 0.0;
        }
    }

    for (int i = 0; i < DCT_size * DCT_size; i++) {
        in[i] = out[i];
    }
}
int write_file(std::ofstream& out, std::vector<uint8_t>& data) {
    out.write(reinterpret_cast<const char*>(data.data()), data.size());
    return (int)data.size();
}
int write_u8(std::ofstream& out, uint8_t data) {
    out.put(data);
    return 1;
}
int write_bits(std::ofstream& out, uint32_t data, int len, int flush) {
    static uint32_t bit_ptr = 0;
    static uint32_t bitbuf = 0x00000000;
    uint8_t w = 0x00;

    bitbuf |= data << (32 - bit_ptr - len);
    bit_ptr += len;

    while (bit_ptr >= 8)
    {
        w = (uint8_t)((bitbuf & 0xFF000000) >> 24);
        write_u8(out, w);
        if (w == 0xFF)
        {
            write_u8(out, 0x00);
        }
        bitbuf <<= 8;
        bit_ptr -= 8;
    }

    if (flush)
    {
        w = (uint8_t)((bitbuf & 0xFF000000) >> 24);
        write_u8(out, w);
    }
    return 0;
}
//flush = 1, flush the buffer
void get_bit_code(int val, int& bit, int& code) {
    int abs_val = val;
    if (val < 0) {
        abs_val = -val;
        val -= 1;
    }
    bit = 1;
    while (abs_val >>= 1) bit++;
    code = (uint16_t)(val & ((1 << bit) - 1));
    //out_log_dyx << val << " " << bit << " " << code << std::endl;
}


void dpcm_rle(std::ofstream& out, double* in, int& lst_dc, std::vector<coefficients>* dc_co, std::vector <coefficients>* ac_co) {
    int dc_data = in[0] - lst_dc;
    lst_dc = in[0];
    if (dc_data == 0) {
        //out_log_dyx << "DC: " << dc_data << std::endl;
        write_bits(out, (*dc_co)[0].code_word, (*dc_co)[0].code_length, 0);
	}
	else {
		int bit, code;
		get_bit_code(dc_data, bit, code);
        //out_log_dyx << "DC: " << dc_data << std::endl;
        write_bits(out, (*dc_co)[bit].code_word, (*dc_co)[bit].code_length, 0);
        write_bits(out, code, bit, 0);
    }

    int count = 0;
    int last = 0;
    for (int i = 1; i <= DCT_size * DCT_size - 1; i++) {
        //out_log_dyx << "AC: " << in[i] << std::endl;
    }
    for (int i = DCT_size * DCT_size - 1; i >= 1; i--) {
		if (in[i] != 0) {
			last = i;
			break;
		}
	}

	for (int i = 1; i <= last; i++) {
        //out_log_dyx << "count: " << count << " in[i]: " << in[i] << std::endl;
        while (in[i] == 0) {
            count++;
            if (count > 15) {
                //out_log_dyx << "AC: fuck " << 0xF0 << " " 
                //    //<< (int)(*ac_co)[0xF0].code_word << " " << (int)(*ac_co)[0xF0].code_length 
                //    << std::endl;
                write_bits(out, (*ac_co)[0xF0].code_word, (*ac_co)[0xF0].code_length, 0);
                count = 0;
            }
            i++;
        }
        int bit, code;
        get_bit_code(in[i], bit, code);
        int run_size = count << 4 | bit;
        //out_log_dyx << "AC: " << in[i] << " " << count << " " << code << " " << bit << std::endl;
        //out_log_dyx << (int)(*ac_co)[run_size].code_word << " " << (int)(*ac_co)[run_size].code_length << std::endl;
        write_bits(out, (*ac_co)[run_size].code_word, (*ac_co)[run_size].code_length, 0);
        write_bits(out, code, bit, 0);
        count = 0;
	}
	if (last != DCT_size * DCT_size - 1)  write_bits(out, (*ac_co)[0x00].code_word, (*ac_co)[0x00].code_length, 0);
    //write_bits(out, 0, 0, 1);
    //out_log_dyx << "end_block" << std::endl;
}

int lst_dc_Y = 0;
int lst_dc_Cb = 0;
int lst_dc_Cr = 0;
void block_encoding(std::ofstream& out, BLOCK* cur) {
    //DCT
    DCT(cur->Y);
    DCT(cur->Cb);
    DCT(cur->Cr);
    //Quantization

    double* Y = new double[DCT_size * DCT_size];
    double* Cb = new double[DCT_size * DCT_size];
    double* Cr = new double[DCT_size * DCT_size];
        
    
    for(int i = 0; i < DCT_size * DCT_size; i++) {
        cur->Y[i] = cur->Y[i] / get_cur(default_luma_table, i);
        cur->Cb[i] = cur->Cb[i] / get_cur(default_chroma_table, i);
        cur->Cr[i] = cur->Cr[i] / get_cur(default_chroma_table, i);
        //out_log_dyx << cur->Y[i] << " ";
	}

    for (int i = 0; i < DCT_size * DCT_size; i++) {
        Y[i] = cur->Y[zigzag_inv[i]];
        Y[i] = (int)(Y[i] + (Y[i] >= 0? 0.5: -0.5));
        Cb[i] = cur->Cb[zigzag_inv[i]];
        Cb[i] = (int)(Cb[i] + (Cb[i] >= 0 ? 0.5 : -0.5));
        Cr[i] = cur->Cr[zigzag_inv[i]];
        Cr[i] = (int)(Cr[i] + (Cr[i] >= 0 ? 0.5 : -0.5));
    }
    //out_log_dyx << std::endl <<  "???" << std::endl;
    //DPCM_RLE
    //out_log_dyx << "\nY: ";
    //for (int i = 0; i < DCT_size * DCT_size; i++) out_log_dyx << (int)Y[i] << " ";
    //out_log_dyx << "\nCb: ";
    //for (int i = 0; i < DCT_size * DCT_size; i++) out_log_dyx << (int)Cb[i] << " ";
    //out_log_dyx << "\nCr: ";
    //for (int i = 0; i < DCT_size * DCT_size; i++) out_log_dyx << (int)Cr[i] << " ";
    dpcm_rle(out, Y, lst_dc_Y, &luma_dc_huffman_table, &luma_ac_huffman_table);
    dpcm_rle(out, Cb, lst_dc_Cb, &chroma_dc_huffman_table, &chroma_ac_huffman_table);
    dpcm_rle(out, Cr, lst_dc_Cr, &chroma_dc_huffman_table, &chroma_ac_huffman_table);
    //free Y, Cb, Cr
    delete[] Y;
    delete[] Cb;
    delete[] Cr;
}

void encoding(std::ofstream& out, BYTE* Data, int width, int height) {
    for (int i = 0; i < height; i += 8) {
        for (int j = 0; j < width; j += 8) {
            BLOCK* Cur_Block = new BLOCK(8, 8);
            for (int u = 0; u < 8; u++) {
                for (int v = 0; v < 8; v++) {
                    BYTE r = Data[((i + u) * width + j + v) * 3];
                    BYTE g = Data[((i + u) * width + j + v) * 3 + 1];
                    BYTE b = Data[((i + u) * width + j + v) * 3 + 2];
                    //out_log_dyx << "r: " << (int)r << " g: " << (int)g << " b: " << (int)b << std::endl;
                    Cur_Block->Y[u * 8 + v] = rgb2Y(r, g, b);
                    Cur_Block->Cb[u * 8 + v] = rgb2Cb(r, g, b);
                    Cur_Block->Cr[u * 8 + v] = rgb2Cr(r, g, b);
                    //out_log_dyx << "rgbtest" << rgb2Y(220, 127, 121);
                    //out_log_dyx << "Y: " << (int)Cur_Block->Y[u * 8 + v] << " Cb: " << (int)Cur_Block->Cb[u * 8 + v] << " Cr: " << (int)Cur_Block->Cr[u * 8 + v] << std::endl;
                    //out_log_dyx << "r: " << (int)r << " g: " << (int)g << " b: " << (int)b << std::endl;
                    //if (g != 0 && b != 0) out_log_dyx << "speical case: i: " << i << " j: " << j << " u: " << u << " v: " << v << " r: " << (int)r << " g: " << (int)g << " b: " << (int)b << 
                    //    " pixelpos:" << (i + u) * width + j + v << std::endl;
                }
            }
            
            //out_log_dyx << "i: " << i << " j: " << j << std::endl;
            block_encoding(out, Cur_Block);
        }
    }
    write_bits(out, 0, 0, 1);//flush the buffer
}
void JPEG_header_encode(std::ofstream &out, int width, int height){
    //write SOI
    out.put(0xFF);
    out.put(0xD8);//marker

    //write APP0
    out.put(0xFF);
    out.put(0xE0);//marker
    out.put(0x00);
    out.put(0x10);//length
    out.put(0x4A);
    out.put(0x46);
    out.put(0x49);
    out.put(0x46);
    out.put(0x00);//JFIF
    out.put(0x01);
    out.put(0x01);//version
    out.put(0x00);//units
    out.put(0x00);
    out.put(0x01);//Xdensity
    out.put(0x00);
    out.put(0x01);//Ydensity
    out.put(0x00);//Xthumbnail
    out.put(0x00);//Ythumbnail

    //write DQT
    //write Y_DQT
    out.put(0xFF);
    out.put(0xDB);//marker
    out.put(0x00);
    out.put(0x84);//length
    out.put(0x00);//QTY
    for (int i = 0; i < DCT_size * DCT_size; i++) luma_table[i] = get_cur(default_luma_table, zigzag_inv[i]);
    for(int i =0 ; i < DCT_size * DCT_size; i++)  write_bits(out, luma_table[i], 8, 0);
    out.put(0x01);//QTC
    for (int i = 0; i < DCT_size * DCT_size; i++)  chroma_table[i] = get_cur(default_chroma_table, zigzag_inv[i]);
    for (int i = 0; i < DCT_size * DCT_size; i++) write_bits(out, chroma_table[i], 8, 0);

	//write SOF0
	out.put(0xFF);
	out.put(0xC0);//marker
	out.put(0x00);
	out.put(0x11);//length
	out.put(0x08);//precision

    write_bits(out, (height >> 8), 8, 0);
    write_bits(out, (height & 0xFF), 8, 0);
    write_bits(out, (width >> 8), 8, 0);
    write_bits(out, (width & 0xFF), 8, 0);
	out.put(0x03);//components
	out.put(0x01);//Y
	out.put(0x11);//sampling
    out.put(0x00);
	out.put(0x02);//Cb
	out.put(0x11);//sampling
    out.put(0x01);
	out.put(0x03);//Cr
	out.put(0x11);//sampling
    out.put(0x01);

    //write DHT
    out.put(0xFF);
    out.put(0xC4);//marker
    out.put(0x01);
    out.put(0xA2);//length
    out.put(0x00);//Y_DC
    for (int i = 0; i < 16; i++) write_bits(out, DcLuminanceCodesPerBitsize[i], 8, 0);
    for (int i = 0; i < 12; i++) write_bits(out, DcLuminanceValues[i], 8, 0);
    out.put(0x10);//Y_AC
    for (int i = 0; i < 16; i++) write_bits(out, AcLuminanceCodesPerBitsize[i], 8, 0);
    for (int i = 0; i < 162; i++) write_bits(out, AcLuminanceValues[i], 8, 0);
    out.put(0x01);//CbCr_DC
    for (int i = 0; i < 16; i++) write_bits(out, DcChrominanceCodesPerBitsize[i], 8, 0);
    for (int i = 0; i < 12; i++) write_bits(out, DcChrominanceValues[i], 8, 0);
    out.put(0x11);//CbCr_AC
    for (int i = 0; i < 16; i++) write_bits(out, AcChrominanceCodesPerBitsize[i], 8, 0);
    for (int i = 0; i < 162; i++) write_bits(out, AcChrominanceValues[i], 8, 0);

//    //write SOS
    out.put(0xFF);
    out.put(0xDA);//marker
    out.put(0x00);
    out.put(0x0C);//length
    out.put(0x03);//components
    out.put(0x01);//Y
    out.put(0x00);//Y_DC_AC
    out.put(0x02);//Cb
    out.put(0x11);//Cb_DC_AC
    out.put(0x03);//Cr
    out.put(0x11);//Cr_DC_AC
    out.put(0x00);//SS
    out.put(0x3F);//SE
    out.put(0x00);//AH_AL
}

void JPEG_ending_encode(std::ofstream &out) {
    //write EOI
	out.put(0xFF);
	out.put(0xD9);//marker
}