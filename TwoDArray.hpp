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

// 	inline TwoDArray<T>& operator/=(const TwoDArray<T>& arr)
//     {
//     	assert(this->width == arr.width);
//     	assert(this->height == arr.height);

//         for(int i = 0; i < height; i++)
//         	for(int j = 0; j < width; j++)
//         	{
//         		smartPtr2D[i][j] /= arr[i][j];
//         	}

//         return *this;
//     }

// 	TwoDArray<T>& operator*=(const TwoDArray<T>& arr)
//     {
//     	assert(this->width == arr.width);
//     	assert(this->height == arr.height);

//         for(int i = 0; i < height; i++)
//         	for(int j = 0; j < width; j++)
//         	{
//         		smartPtr2D[i][j] *= arr[i][j];
//         	}

//         return *this;
//     }

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