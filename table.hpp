#include "TwoDArray.hpp"

class QuantizationTable
{
public:
    int luminance_quantization_table[8][8] ={
      {16,  11,  10,  16,  24,  40,  51,  61},
      {12,  12,  14,  19,  26,  58,  60,  55},
      {14,  13,  16,  24,  40,  57,  69,  56},
      {14,  17,  22,  29,  51,  87,  80,  62},
      {18,  22,  37,  56,  68, 109, 103,  77},
      {24,  35,  55,  64,  81, 104, 113,  92},
      {49,  64,  78,  87, 103, 121, 120, 101},
      {72,  92,  95,  98, 112, 100, 103,  99}
        };
    int chrominance_quantization_table[8][8] ={
      {17,  18,  24,  47,  99,  99,  99,  99},
      {18,  21,  26,  66,  99,  99,  99,  99},
      {24,  26,  56,  99,  99,  99,  99,  99},
      {47,  66,  99,  99,  99,  99,  99,  99},
      {99,  99,  99,  99,  99,  99,  99,  99},
      {99,  99,  99,  99,  99,  99,  99,  99},
      {99,  99,  99,  99,  99,  99,  99,  99},
      {99,  99,  99,  99,  99,  99,  99,  99},
    };

    std::vector<TwoDArray<unsigned int>> DQT[4];

    void UpSampling(TwoDArray<unsigned int> &arr, TwoDArray<unsigned int> &arr2, int size) {
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                int a = arr(i/2, j/2); 
                int b = arr((i-1)/2, j/2);
                int c = arr(i/2, (j-1)/2);
                int d = arr((i-1)/2, (j-1)/2);
                arr2(i, j) = a + b + c + d;
            }
        }
    }

    QuantizationTable() 
    {
	puts("123");
        TwoDArray<unsigned int> c(8, 8);
        TwoDArray<unsigned int> l(8, 8);
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++) {
                c(i, j) = chrominance_quantization_table[i][j];
                l(i, j) = luminance_quantization_table[i][j];
            }

	puts("4");
        DQT[0].push_back(std::move(l));
        DQT[0].push_back(std::move(c));
	puts("456");
        
        // create larger quantization table
        for (int s = 16, j = 0; s <= 64; s*=2, j++) {
            TwoDArray<unsigned int> l_temp(s,s);
            TwoDArray<unsigned int> c_temp(s,s);
		printf("j: %d\n", j);
            UpSampling(DQT[j][0], l_temp, s);
            UpSampling(DQT[j][1], c_temp, s);
            DQT[j+1].push_back(std::move(l_temp));
            DQT[j+1].push_back(std::move(c_temp));     
        }
    }
};
