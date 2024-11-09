#include "PicReader.h"
#include "Img_Alg.h"
#include "define.h"
#define _CRT_SECURE_NO_WARNINGS

FILE* file;

void myOutput(unsigned char oneByte) { fputc(oneByte, file); }

void JPEG(const std::string &input, const std::string &output) {
    PicReader Img_Reader;
    BYTE* data = nullptr;
    UINT x, y;

    const char* input_cstr = input.c_str(); // Convert std::string to const char*
    Img_Reader.readPic(input_cstr);
	Img_Reader.getData(data, x, y);
	auto pixels = new unsigned char[x * y * 3];
	for (int i = 0; i < x * y; i++) {
		pixels[i * 3 + 0] = data[i * 4 + 0];
		pixels[i * 3 + 1] = data[i * 4 + 1];
		pixels[i * 3 + 2] = data[i * 4 + 2];
	}

	std::ofstream out(output, std::ios::out | std::ios::binary);
	
	JPEG_Encoder Img_Data(out, pixels, (int)x, (int)y);
}


void IMG_COMPRESS(const std::string &address) {
	std::string input =  address;
    std::string inputWithoutExtension = input.substr(0, input.find_last_of("."));
	std::string output = inputWithoutExtension + ".jpeg";
	JPEG(input, output);
	return;
}

void IMG_READ(const std::string &address) {
    std::string input = address;
    std::string inputWithoutExtension = input.substr(0, input.find_last_of("."));
    std::string output = inputWithoutExtension + ".jpeg";
	system(output.c_str());

    return;
}

int main(int argc, char* argv[]) {
	std::string command = argv[1];
	std::string address = argv[2];
	if (command == "-compress") IMG_COMPRESS(address);
	else if (command == "-read") IMG_READ(address);
	else std::cerr << "Invalid command" << std::endl;
}