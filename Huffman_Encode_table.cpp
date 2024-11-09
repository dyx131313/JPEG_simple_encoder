#include "Img_Alg.h"

std::vector<coefficients> luma_dc_huffman_table = {
    {0, 2, 0b00},
    {1, 3, 0b010},
    {2, 3, 0b011},
    {3, 3, 0b100},
    {4, 3, 0b101},
    {5, 3, 0b110},
    {6, 4, 0b1110},
    {7, 5, 0b11110},
    {8, 6, 0b111110},
    {9, 7, 0b1111110},
    {10, 8, 0b11111110},
    {11, 9, 0b111111110},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0}
};

std::vector<coefficients> chroma_dc_huffman_table = {
    {0, 2, 0b00},
    {1, 2, 0b01},
    {2, 2, 0b10},
    {3, 3, 0b110},
    {4, 4, 0b1110},
    {5, 5, 0b11110},
    {6, 6, 0b111110},
    {7, 7, 0b1111110},
    {8, 8, 0b11111110},
    {9, 9, 0b111111110},
    {10, 10, 0b1111111110},
    {11, 11, 0b11111111110},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0}
};

std::vector<coefficients> luma_ac_huffman_table = {
    {0, 4, 0b1010},
    {1, 2, 0b00},
    {2, 2, 0b01},
    {3, 3, 0b100},
    {4, 4, 0b1011},
    {5, 5, 0b11010},
    {6, 7, 0b1111000},
    {7, 8, 0b11111000},
    {8, 10, 0b1111110110},
    {9, 16, 0b1111111110000010},
    {10, 16, 0b1111111110000011},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {17, 4, 0b1100},
    {18, 5, 0b11011},
    {19, 7, 0b1111001},
    {20, 9, 0b111110110},
    {21, 11, 0b11111110110},
    {22, 16, 0b1111111110000100},
    {23, 16, 0b1111111110000101},
    {24, 16, 0b1111111110000110},
    {25, 16, 0b1111111110000111},
    {26, 16, 0b1111111110001000},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {33, 5, 0b11100},
    {34, 8, 0b11111001},
    {35, 10, 0b1111110111},
    {36, 12, 0b111111110100},
    {37, 16, 0b1111111110001001},
    {38, 16, 0b1111111110001010},
    {39, 16, 0b1111111110001011},
    {40, 16, 0b1111111110001100},
    {41, 16, 0b1111111110001101},
    {42, 16, 0b1111111110001110},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {49, 6, 0b111010},
    {50, 9, 0b111110111},
    {51, 12, 0b111111110101},
    {52, 16, 0b1111111110001111},
    {53, 16, 0b1111111110010000},
    {54, 16, 0b1111111110010001},
    {55, 16, 0b1111111110010010},
    {56, 16, 0b1111111110010011},
    {57, 16, 0b1111111110010100},
    {58, 16, 0b1111111110010101},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {65, 6, 0b111011},
    {66, 10, 0b1111111000},
    {67, 16, 0b1111111110010110},
    {68, 16, 0b1111111110010111},
    {69, 16, 0b1111111110011000},
    {70, 16, 0b1111111110011001},
    {71, 16, 0b1111111110011010},
    {72, 16, 0b1111111110011011},
    {73, 16, 0b1111111110011100},
    {74, 16, 0b1111111110011101},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {81, 7, 0b1111010},
    {82, 11, 0b11111110111},
    {83, 16, 0b1111111110011110},
    {84, 16, 0b1111111110011111},
    {85, 16, 0b1111111110100000},
    {86, 16, 0b1111111110100001},
    {87, 16, 0b1111111110100010},
    {88, 16, 0b1111111110100011},
    {89, 16, 0b1111111110100100},
    {90, 16, 0b1111111110100101},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {97, 7, 0b1111011},
    {98, 12, 0b111111110110},
    {99, 16, 0b1111111110100110},
    {100, 16, 0b1111111110100111},
    {101, 16, 0b1111111110101000},
    {102, 16, 0b1111111110101001},
    {103, 16, 0b1111111110101010},
    {104, 16, 0b1111111110101011},
    {105, 16, 0b1111111110101100},
    {106, 16, 0b1111111110101101},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {113, 8, 0b11111010},
    {114, 12, 0b111111110111},
    {115, 16, 0b1111111110101110},
    {116, 16, 0b1111111110101111},
    {117, 16, 0b1111111110110000},
    {118, 16, 0b1111111110110001},
    {119, 16, 0b1111111110110010},
    {120, 16, 0b1111111110110011},
    {121, 16, 0b1111111110110100},
    {122, 16, 0b1111111110110101},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {129, 9, 0b111111000},
    {130, 15, 0b111111111000000},
    {131, 16, 0b1111111110110110},
    {132, 16, 0b1111111110110111},
    {133, 16, 0b1111111110111000},
    {134, 16, 0b1111111110111001},
    {135, 16, 0b1111111110111010},
    {136, 16, 0b1111111110111011},
    {137, 16, 0b1111111110111100},
    {138, 16, 0b1111111110111101},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {145, 9, 0b111111001},
    {146, 16, 0b1111111110111110},
    {147, 16, 0b1111111110111111},
    {148, 16, 0b1111111111000000},
    {149, 16, 0b1111111111000001},
    {150, 16, 0b1111111111000010},
    {151, 16, 0b1111111111000011},
    {152, 16, 0b1111111111000100},
    {153, 16, 0b1111111111000101},
    {154, 16, 0b1111111111000110},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {161, 9, 0b111111010},
    {162, 16, 0b1111111111000111},
    {163, 16, 0b1111111111001000},
    {164, 16, 0b1111111111001001},
    {165, 16, 0b1111111111001010},
    {166, 16, 0b1111111111001011},
    {167, 16, 0b1111111111001100},
    {168, 16, 0b1111111111001101},
    {169, 16, 0b1111111111001110},
    {170, 16, 0b1111111111001111},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {177, 10, 0b1111111001},
    {178, 16, 0b1111111111010000},
    {179, 16, 0b1111111111010001},
    {180, 16, 0b1111111111010010},
    {181, 16, 0b1111111111010011},
    {182, 16, 0b1111111111010100},
    {183, 16, 0b1111111111010101},
    {184, 16, 0b1111111111010110},
    {185, 16, 0b1111111111010111},
    {186, 16, 0b1111111111011000},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {193, 10, 0b1111111010},
    {194, 16, 0b1111111111011001},
    {195, 16, 0b1111111111011010},
    {196, 16, 0b1111111111011011},
    {197, 16, 0b1111111111011100},
    {198, 16, 0b1111111111011101},
    {199, 16, 0b1111111111011110},
    {200, 16, 0b1111111111011111},
    {201, 16, 0b1111111111100000},
    {202, 16, 0b1111111111100001},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {209, 11, 0b11111111000},
    {210, 16, 0b1111111111100010},
    {211, 16, 0b1111111111100011},
    {212, 16, 0b1111111111100100},
    {213, 16, 0b1111111111100101},
    {214, 16, 0b1111111111100110},
    {215, 16, 0b1111111111100111},
    {216, 16, 0b1111111111101000},
    {217, 16, 0b1111111111101001},
    {218, 16, 0b1111111111101010},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {225, 16, 0b1111111111101011},
    {226, 16, 0b1111111111101100},
    {227, 16, 0b1111111111101101},
    {228, 16, 0b1111111111101110},
    {229, 16, 0b1111111111101111},
    {230, 16, 0b1111111111110000},
    {231, 16, 0b1111111111110001},
    {232, 16, 0b1111111111110010},
    {233, 16, 0b1111111111110011},
    {234, 16, 0b1111111111110100},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {240, 11, 0b11111111001},
    {241, 16, 0b1111111111110101},
    {242, 16, 0b1111111111110110},
    {243, 16, 0b1111111111110111},
    {244, 16, 0b1111111111111000},
    {245, 16, 0b1111111111111001},
    {246, 16, 0b1111111111111010},
    {247, 16, 0b1111111111111011},
    {248, 16, 0b1111111111111100},
    {249, 16, 0b1111111111111101},
    {250, 16, 0b1111111111111110}
};

std::vector<coefficients> chroma_ac_huffman_table = {
    {0, 2, 0b00}, // 0/0 (EOB)
    {1, 2, 0b01}, // 0/1
    {2, 3, 0b100}, // 0/2
    {3, 4, 0b1010}, // 0/3
    {4, 5, 0b11000}, // 0/4
    {5, 5, 0b11001}, // 0/5
    {6, 6, 0b111000}, // 0/6
    {7, 7, 0b1111000}, // 0/7
    {8, 9, 0b111110100}, // 0/8
    {9, 10, 0b1111110110}, // 0/9
    {10, 12, 0b111111110100}, // 0/A
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {17, 4, 0b1011}, // 1/1
    {18, 6, 0b111001}, // 1/2
    {19, 8, 0b11110110}, // 1/3
    {20, 9, 0b111110101}, // 1/4
    {21, 11, 0b11111110110}, // 1/5
    {22, 12, 0b111111110101}, // 1/6
    {23, 16, 0b1111111110001000}, // 1/7
    {24, 16, 0b1111111110001001}, // 1/8
    {25, 16, 0b1111111110001010}, // 1/9
    {26, 16, 0b1111111110001011}, // 1/A
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {33, 5, 0b11010}, // 2/1
    {34, 8, 0b11110111}, // 2/2
    {35, 10, 0b1111110111}, // 2/3
    {36, 12, 0b111111110110}, // 2/4
    {37, 15, 0b111111111000010}, // 2/5
    {38, 16, 0b1111111110001100}, // 2/6
    {39, 16, 0b1111111110001101}, // 2/7
    {40, 16, 0b1111111110001110}, // 2/8
    {41, 16, 0b1111111110001111}, // 2/9
    {42, 16, 0b1111111110010000}, // 2/A
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {49, 5, 0b11011}, // 3/1
    {50, 8, 0b11111000}, // 3/2
    {51, 10, 0b1111111000}, // 3/3
    {52, 12, 0b111111110111}, // 3/4
    {53, 16, 0b1111111110010001}, // 3/5
    {54, 16, 0b1111111110010010}, // 3/6
    {55, 16, 0b1111111110010011}, // 3/7
    {56, 16, 0b1111111110010100}, // 3/8
    {57, 16, 0b1111111110010101}, // 3/9
    {58, 16, 0b1111111110010110}, // 3/A
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {65, 6, 0b111010}, // 4/1
    {66, 9, 0b111110110}, // 4/2
    {67, 16, 0b1111111110010111}, // 4/3
    {68, 16, 0b1111111110011000}, // 4/4
    {69, 16, 0b1111111110011001}, // 4/5
    {70, 16, 0b1111111110011010}, // 4/6
    {71, 16, 0b1111111110011011}, // 4/7
    {72, 16, 0b1111111110011100}, // 4/8
    {73, 16, 0b1111111110011101}, // 4/9
    {74, 16, 0b1111111110011110}, // 4/A
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {81, 6, 0b111011}, // 5/1
    {82, 10, 0b1111111001}, // 5/2
    {83, 16, 0b1111111110011111}, // 5/3
    {84, 16, 0b1111111110100000}, // 5/4
    {85, 16, 0b1111111110100001}, // 5/5
    {86, 16, 0b1111111110100010}, // 5/6
    {87, 16, 0b1111111110100011}, // 5/7
    {88, 16, 0b1111111110100100}, // 5/8
    {89, 16, 0b1111111110100101}, // 5/9
    {90, 16, 0b1111111110100110}, // 5/A
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {0, 0, 0},
    {97, 7, 0b1111001}, // 6/1
    {98, 11, 0b11111110111}, // 6/2
    {99, 16, 0b1111111110100111}, // 6/3
    {100, 16, 0b1111111110101000}, // 6/4
    {101, 16, 0b1111111110101001}, // 6/5
    {102, 16, 0b1111111110101010}, // 6/6
    {103, 16, 0b1111111110101011}, // 6/7
    {104, 16, 0b1111111110101100}, // 6/8
    {105, 16, 0b1111111110101101}, // 6/9
    {106, 16, 0b1111111110101110}, // 6/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {113, 7, 0b1111010}, // 7/1
    {114, 11, 0b11111111000}, // 7/2
    {115, 16, 0b1111111110101111}, // 7/3
    {116, 16, 0b1111111110110000}, // 7/4
    {117, 16, 0b1111111110110001}, // 7/5
    {118, 16, 0b1111111110110010}, // 7/6
    {119, 16, 0b1111111110110011}, // 7/7
    {120, 16, 0b1111111110110100}, // 7/8
    {121, 16, 0b1111111110110101}, // 7/9
    {122, 16, 0b1111111110110110}, // 7/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {129, 8, 0b11111001}, // 8/1
    {130, 16, 0b1111111110110111}, // 8/2
    {131, 16, 0b1111111110111000}, // 8/3
    {132, 16, 0b1111111110111001}, // 8/4
    {133, 16, 0b1111111110111010}, // 8/5
    {134, 16, 0b1111111110111011}, // 8/6
    {135, 16, 0b1111111110111100}, // 8/7
    {136, 16, 0b1111111110111101}, // 8/8
    {137, 16, 0b1111111110111110}, // 8/9
    {138, 16, 0b1111111110111111}, // 8/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {145, 9, 0b111110111}, // 9/1
    {146, 16, 0b1111111111000000}, // 9/2
    {147, 16, 0b1111111111000001}, // 9/3
    {148, 16, 0b1111111111000010}, // 9/4
    {149, 16, 0b1111111111000011}, // 9/5
    {150, 16, 0b1111111111000100}, // 9/6
    {151, 16, 0b1111111111000101}, // 9/7
    {152, 16, 0b1111111111000110}, // 9/8
    {153, 16, 0b1111111111000111}, // 9/9
    {154, 16, 0b1111111111001000}, // 9/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {161, 9, 0b111111000}, // A/1
    {162, 16, 0b1111111111001001}, // A/2
    {163, 16, 0b1111111111001010}, // A/3
    {164, 16, 0b1111111111001011}, // A/4
    {165, 16, 0b1111111111001100}, // A/5
    {166, 16, 0b1111111111001101}, // A/6
    {167, 16, 0b1111111111001110}, // A/7
    {168, 16, 0b1111111111001111}, // A/8
    {169, 16, 0b1111111111010000}, // A/9
    {170, 16, 0b1111111111010001}, // A/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {177, 9, 0b111111001}, // B/1
    {178, 16, 0b1111111111010010}, // B/2
    {179, 16, 0b1111111111010011}, // B/3
    {180, 16, 0b1111111111010100}, // B/4
    {181, 16, 0b1111111111010101}, // B/5
    {182, 16, 0b1111111111010110}, // B/6
    {183, 16, 0b1111111111010111}, // B/7
    {184, 16, 0b1111111111011000}, // B/8
    {185, 16, 0b1111111111011001}, // B/9
    {186, 16, 0b1111111111011010}, // B/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {193, 9, 0b111111010}, // C/1
    {194, 16, 0b1111111111011011}, // C/2
    {195, 16, 0b1111111111011100}, // C/3
    {196, 16, 0b1111111111011101}, // C/4
    {197, 16, 0b1111111111011110}, // C/5
    {198, 16, 0b1111111111011111}, // C/6
    {199, 16, 0b1111111111100000}, // C/7
    {200, 16, 0b1111111111100001}, // C/8
    {201, 16, 0b1111111111100010}, // C/9
    {202, 16, 0b1111111111100011}, // C/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {209, 11, 0b11111111001}, // D/1
    {210, 16, 0b1111111111100100}, // D/2
    {211, 16, 0b1111111111100101}, // D/3
    {212, 16, 0b1111111111100110}, // D/4
    {213, 16, 0b1111111111100111}, // D/5
    {214, 16, 0b1111111111101000}, // D/6
    {215, 16, 0b1111111111101001}, // D/7
    {216, 16, 0b1111111111101010}, // D/8
    {217, 16, 0b1111111111101011}, // D/9
    {218, 16, 0b1111111111101100}, // D/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {225, 14, 0b11111111100000}, // E/1
    {226, 16, 0b1111111111101101}, // E/2
    {227, 16, 0b1111111111101110}, // E/3
    {228, 16, 0b1111111111101111}, // E/4
    {229, 16, 0b1111111111110000}, // E/5
    {230, 16, 0b1111111111110001}, // E/6
    {231, 16, 0b1111111111110010}, // E/7
    {232, 16, 0b1111111111110011}, // E/8
    {233, 16, 0b1111111111110100}, // E/9
    {234, 16, 0b1111111111110101}, // E/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    {240, 10, 0b1111111010}, // F/0 (ZRL)
    {241, 15, 0b111111111000011}, // F/1
    {242, 16, 0b1111111111110110}, // F/2
    {243, 16, 0b1111111111110111}, // F/3
    {244, 16, 0b1111111111111000}, // F/4
    {245, 16, 0b1111111111111001}, // F/5
    {246, 16, 0b1111111111111010}, // F/6
    {247, 16, 0b1111111111111011}, // F/7
    {248, 16, 0b1111111111111100}, // F/8
    {249, 16, 0b1111111111111101}, // F/9
    {250, 16, 0b1111111111111110}, // F/A
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
    { 0, 0, 0 },
};


const uint8_t DcLuminanceCodesPerBitsize[16] = { 0,1,5,1,1,1,1,1,1,0,0,0,0,0,0,0 };   // sum = 12
const uint8_t DcLuminanceValues[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 };         // => 12 codes
const uint8_t AcLuminanceCodesPerBitsize[16] = { 0,2,1,3,3,2,4,3,5,5,4,4,0,0,1,125 }; // sum = 162
const uint8_t AcLuminanceValues[162] =                                        // => 162 codes
{ 0x01,0x02,0x03,0x00,0x04,0x11,0x05,0x12,0x21,0x31,0x41,0x06,0x13,0x51,0x61,0x07,0x22,0x71,0x14,0x32,0x81,0x91,0xA1,0x08, // 16*10+2 symbols because
  0x23,0x42,0xB1,0xC1,0x15,0x52,0xD1,0xF0,0x24,0x33,0x62,0x72,0x82,0x09,0x0A,0x16,0x17,0x18,0x19,0x1A,0x25,0x26,0x27,0x28, // upper 4 bits can be 0..F
  0x29,0x2A,0x34,0x35,0x36,0x37,0x38,0x39,0x3A,0x43,0x44,0x45,0x46,0x47,0x48,0x49,0x4A,0x53,0x54,0x55,0x56,0x57,0x58,0x59, // while lower 4 bits can be 1..A
  0x5A,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6A,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x83,0x84,0x85,0x86,0x87,0x88,0x89, // plus two special codes 0x00 and 0xF0
  0x8A,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9A,0xA2,0xA3,0xA4,0xA5,0xA6,0xA7,0xA8,0xA9,0xAA,0xB2,0xB3,0xB4,0xB5,0xB6, // order of these symbols was determined empirically by JPEG committee
  0xB7,0xB8,0xB9,0xBA,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,0xC8,0xC9,0xCA,0xD2,0xD3,0xD4,0xD5,0xD6,0xD7,0xD8,0xD9,0xDA,0xE1,0xE2,
  0xE3,0xE4,0xE5,0xE6,0xE7,0xE8,0xE9,0xEA,0xF1,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,0xF8,0xF9,0xFA };
// Huffman definitions for second DC/AC tables (chrominance / Cb and Cr channels)
const uint8_t DcChrominanceCodesPerBitsize[16] = { 0,3,1,1,1,1,1,1,1,1,1,0,0,0,0,0 };   // sum = 12
const uint8_t DcChrominanceValues[12] = { 0,1,2,3,4,5,6,7,8,9,10,11 };         // => 12 codes (identical to DcLuminanceValues)
const uint8_t AcChrominanceCodesPerBitsize[16] = { 0,2,1,2,4,4,3,4,7,5,4,4,0,1,2,119 }; // sum = 162
const uint8_t AcChrominanceValues[162] =                                        // => 162 codes
{ 0x00,0x01,0x02,0x03,0x11,0x04,0x05,0x21,0x31,0x06,0x12,0x41,0x51,0x07,0x61,0x71,0x13,0x22,0x32,0x81,0x08,0x14,0x42,0x91, // same number of symbol, just different order
  0xA1,0xB1,0xC1,0x09,0x23,0x33,0x52,0xF0,0x15,0x62,0x72,0xD1,0x0A,0x16,0x24,0x34,0xE1,0x25,0xF1,0x17,0x18,0x19,0x1A,0x26, // (which is more efficient for AC coding)
  0x27,0x28,0x29,0x2A,0x35,0x36,0x37,0x38,0x39,0x3A,0x43,0x44,0x45,0x46,0x47,0x48,0x49,0x4A,0x53,0x54,0x55,0x56,0x57,0x58,
  0x59,0x5A,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6A,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7A,0x82,0x83,0x84,0x85,0x86,0x87,
  0x88,0x89,0x8A,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9A,0xA2,0xA3,0xA4,0xA5,0xA6,0xA7,0xA8,0xA9,0xAA,0xB2,0xB3,0xB4,
  0xB5,0xB6,0xB7,0xB8,0xB9,0xBA,0xC2,0xC3,0xC4,0xC5,0xC6,0xC7,0xC8,0xC9,0xCA,0xD2,0xD3,0xD4,0xD5,0xD6,0xD7,0xD8,0xD9,0xDA,
  0xE2,0xE3,0xE4,0xE5,0xE6,0xE7,0xE8,0xE9,0xEA,0xF2,0xF3,0xF4,0xF5,0xF6,0xF7,0xF8,0xF9,0xFA };