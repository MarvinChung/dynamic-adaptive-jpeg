#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <sstream>
#include <cmath>
#include <assert.h>
#include "zigzagtable.hpp"

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
    
    std::vector<T> zigzag()
    {
        std::vector<T> zz;
        int *zigzagOrder = NULL;
        if(width == 8 && height == 8)
            zigzagOrder = zigzagOrder8x8;
        else if (width == 16 && height == 16)
            zigzagOrder = zigzagOrder16x16;
        else if (width == 32 && height == 32)
            zigzagOrder = zigzagOrder32x32;
        else if (width == 64 && height == 64)
            zigzagOrder = zigzagOrder64x64;
        
        assert(zigzagOrder!=NULL);
        
        for(int i = 0; i < width*height; i++)
        {
            int row_idx = zigzagOrder[i]/width;
            int col_idx = zigzagOrder[i]%width;
            zz.push_back(smartPtr2D[row_idx][col_idx]);
        }
        return zz;
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