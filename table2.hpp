#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <sstream>
#include <cmath>
#include <assert.h>

template <class T>
class TwoDArray
{
public:
	std::unique_ptr<std::unique_ptr<T[]>[]>     smartPtr2D;
	std::unique_ptr<T[]>                        smartPtr1D;
	int width;
	int height;	
	TwoDArray() {}
	TwoDArray(int width, int height): width(width), height(height) 
	{
		
		smartPtr2D = std::make_unique<std::unique_ptr<T[]>[]>(height);
		for (int i = 0; i < height; i++)
		{
			smartPtr1D = std::make_unique<T[]>(width);
			for(int j = 0; j < width; j++)
				smartPtr1D[j] = 0;
		    smartPtr2D[i] = std::move(smartPtr1D);
		}
	}
	T& operator()(int i, int j)
	{
		return smartPtr2D[i][j];
	}
	const T& operator () (int i, int j) const
	{
	    return smartPtr2D[i][j];
	}

	void show()
	{
		for(int i = 0; i < height; i++)
		{
			for (int j = 0; j < width ; j++)
			{
				if (typeid(T).name() == typeid(unsigned char).name())
					std::cout << int(smartPtr2D[i][j]) << " ";
				else
					std::cout << std::fixed << smartPtr2D[i][j] << " ";
			}
			std::cout << std::endl;
		}
	}
};


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

    QuantizationTable() : QuantizationTable(80) {}

    QuantizationTable(int quality) // 1~100
    {
        TwoDArray<unsigned int> c(8, 8);
        TwoDArray<unsigned int> l(8, 8);

	unsigned int qualityScale;
	if (quality < 50)
	    qualityScale = 5000 / qualityScale;
	else
	    qualityScale = 200 - quality * 2;
        for (int i = 0; i < 8; i++)
            for (int j = 0; j < 8; j++) {
                c(i, j) = (chrominance_quantization_table[i][j] * qualityScale + 50) / 100;
                l(i, j) = (luminance_quantization_table[i][j] * qualityScale + 50) / 100;
		if (c(i, j) == 0)
		    c(i, j) = 1;
	        if (c(i, j) > 255)
		    c(i, j) = 255;
		if (l(i, j) == 0)
		    l(i, j) = 1;
	        if (l(i, j) > 255)
		    l(i, j) = 255;
            }

        DQT[0].push_back(std::move(l));
        DQT[0].push_back(std::move(c));
        
        // create larger quantization table
        for (int s = 16, j = 0; s <= 64; s*=2, j++) {
            TwoDArray<unsigned int> l_temp(s,s);
            TwoDArray<unsigned int> c_temp(s,s);
            UpSampling(DQT[j][0], l_temp, s);
            UpSampling(DQT[j][1], c_temp, s);
            DQT[j+1].push_back(std::move(l_temp));
            DQT[j+1].push_back(std::move(c_temp));     
        }
    }
};
