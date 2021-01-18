#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <sstream>


class Image
{
public:
	std::unique_ptr<std::unique_ptr<std::unique_ptr<unsigned char[]>[]>[]>     smartPtr3D;
	std::unique_ptr<std::unique_ptr<unsigned char[]>[]>     smartPtr2D;
	std::unique_ptr<unsigned char[]>                        smartPtr1D;		
	int width;
	int height;
	Image() {};
	Image(int width, int height): width(width), height(height) 
	{
		std::cout << "width: " << width << " height: " << height << std::endl;
		smartPtr3D = std::make_unique< std::unique_ptr<std::unique_ptr<unsigned char[]>[]>[] >(height);
		for (int i = 0; i < height; i++)
		{
		    smartPtr2D = std::make_unique<std::unique_ptr<unsigned char[]>[]>(width);
		    for (int j = 0; j < width; j++)
		    {
		    	smartPtr1D = std::make_unique<unsigned char[]>(3);
		    	smartPtr2D[j] = std::move(smartPtr1D);
		    }
		    smartPtr3D[i] = std::move(smartPtr2D);
		}
	}
	unsigned char& operator()(int i, int j, int k)
	{
		return smartPtr3D[i][j][k];
	}
	const unsigned char& operator () (int i, int j, int k) const
	{
	    return smartPtr3D[i][j][k];
	}
	void show()
	{
		for(int i = 0; i < height; i++)
		{
			for (int j = 0; j < width ; j++)
			{
				for(int k = 0; k < 3; k++)
					printf("%d ",smartPtr3D[i][j][k]);
				std::cout << std::endl;
			}
		}
	}
	void printXY(int i , int j)
	{
		for(int k = 0; k < 3; k++)
			printf("%d ",smartPtr3D[i][j][k]);
		std::cout << std::endl;
	}

	void printXYZ(int i , int j, int k)
	{
		printf("%d\n",smartPtr3D[i][j][k]);
	}

	
};

class PPM_Image_Reader{
public:
	Image m_image;
	PPM_Image_Reader() {}
	PPM_Image_Reader(const std::string& input)
	{
		int image_pixel_ct = 0;
		std::ifstream infile(input);
		int iter = 0;
		for( std::string line; std::getline(infile, line); iter++)
		{
			//std::cout << line << "end line ";
		    if(iter == 0 || iter == 2)
		    {
		    	continue;
		    }
		    else if(iter == 1)
		    {	
		    	std::vector<std::string> tmp;
		    	tmp.push_back(line);
		    	std::istringstream iss(tmp[0]);
  				std::string s;
  				for(int ct=0; getline(iss, s, ' '); ct++)
  				{
  					tmp.push_back(s);
  				}
		    	
		    	m_image = Image(std::stoi(tmp[1]), std::stoi(tmp[2]));
		    }
		    else
		    {
		    	std::vector<std::string> tmp;
		    	tmp.push_back(line);
		    	std::istringstream iss(tmp[0]);
  				std::string s;
				int x = image_pixel_ct/m_image.width;
				int y = image_pixel_ct%m_image.width;
  				for(int ct=0; getline(iss, s, ' '); ct++)
  				{
  					tmp.push_back(s);

  					m_image(x, y, ct) = std::stoi(tmp[ct+1]);
  				}
#ifdef DEBUG
  				m_image.printXY(x, y);
#endif 
		    	image_pixel_ct++;
		    }
		}
	}
};

// class JPGEncoder{
// public:
// 	PPM_Image_Reader m_image_reader;
// 	JPGEncoder(PPM_Image_Reader& image_reader) : m_image_reader(image_reader)
// 	{

// 	}
// }

int main()
{
	PPM_Image_Reader("balls.ppm");
}