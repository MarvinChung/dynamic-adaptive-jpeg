#include "table.hpp"
#define MERGE_THERSHOLD 200

class Scanner
{
public:
	FILE *m_fp;
    unsigned char m_buffer[2];
    unsigned char m_bit_buffer;
    unsigned char m_bit_count;
    Scanner(FILE *fp): m_fp(fp) {
        m_bit_count = 0;
    }
    unsigned char readByte()
    {
        //check_reach_EOF();
        fread(&m_buffer[0], 1, 1, m_fp);

        return m_buffer[0];
    }

    unsigned short read2Bytes()
    {
        //check_reach_EOF();
        fread(m_buffer, 1, 2, m_fp);
        return  m_buffer[0] * 256 + m_buffer[1];
    }

    unsigned char readBit() {
        if (m_bit_count == 0) {
            fread(&m_bit_buffer, 1, 1, m_fp);
            if (m_bit_buffer == 0xFF) {
                unsigned char check;
                fread(&check, 1, 1, m_fp);
                if (check != 0x00) {
                    fprintf(stderr, "non 0xFF00!\n");
                }
            }
        }
        bool ret = m_bit_buffer & (1 << (7 - m_bit_count));
        m_bit_count = (m_bit_count == 7 ? 0 : m_bit_count + 1);
        //printf("ret:%d\n",ret);
        return (unsigned char)ret;
    }

    void seek(int length)
    {
        fseek(m_fp, length, SEEK_CUR);
    }

    bool check_reach_EOF()
    {
        if(feof(m_fp))
        {
            puts("reachEOF!");
        }
        return feof(m_fp);
    }
};

template <class T>
class Image
{
public:

	std::unique_ptr<std::unique_ptr<std::unique_ptr<T[]>[]>[]>     smartPtr3D;
	std::unique_ptr<std::unique_ptr<T[]>[]>     smartPtr2D;
	std::unique_ptr<T[]>                        smartPtr1D;		
	int width;
	int height;
	Image() {}
	Image(int width, int height): width(width), height(height) 
	{
		std::cout << "width: " << width << " height: " << height << std::endl;
		smartPtr3D = std::make_unique< std::unique_ptr<std::unique_ptr<T[]>[]>[] >(height);
		for (int i = 0; i < height; i++)
		{
		    smartPtr2D = std::make_unique<std::unique_ptr<T[]>[]>(width);
		    for (int j = 0; j < width; j++)
		    {
		    	smartPtr1D = std::make_unique<T[]>(3);
		    	for(int k = 0; k < 3; k++)
		    		smartPtr1D[k] = 0;
		    	smartPtr2D[j] = std::move(smartPtr1D);
		    }
		    smartPtr3D[i] = std::move(smartPtr2D);
		}
	}
	T& operator()(int i, int j, int k)
	{
		return smartPtr3D[i][j][k];
	}
	const T& operator () (int i, int j, int k) const
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
					if (typeid(T).name() == typeid(unsigned char).name())
						std::cout << int(smartPtr3D[i][j][k]) << " ";
					else
						std::cout << std::fixed << smartPtr3D[i][j][k] << " ";
				std::cout << std::endl;
			}
		}
	}
	void printXY(int i , int j)
	{
		for(int k = 0; k < 3; k++)
			if (typeid(T).name() == typeid(unsigned char).name())
				std::cout << int(smartPtr3D[i][j][k]) << " ";
			else
				std::cout << std::fixed << smartPtr3D[i][j][k] << " ";
		std::cout << std::endl;
	}

	void printXYZ(int i , int j, int k)
	{
		if (typeid(T).name() == typeid(unsigned char).name())
			std::cout << int(smartPtr3D[i][j][k]) << " ";
		else
			std::cout << std::fixed << smartPtr3D[i][j][k] << " ";
	}

	Image<double> RGB2YCbCr()
	{
#ifdef DEBUG
		std::cout << "template_name:" << typeid(T).name() << std::endl;
#endif
		// only RGB (unsigned char) images are allowed to be changed to YCbCr (double)
		assert(typeid(T).name() == typeid(unsigned char).name());
		Image<double> YCbCr_Image(width, height);
		for(int i = 0; i < height; i++)
		{
			for(int j = 0; j < width; j++)
			{
				//Y:0 Cb:1 Cr:2
				YCbCr_Image(i,j,0) = 0.299    * smartPtr3D[i][j][0] + 0.587    * smartPtr3D[i][j][1] + 0.114    * smartPtr3D[i][j][2];
				YCbCr_Image(i,j,1) = 128 - 0.168736 * smartPtr3D[i][j][0] - 0.331264 * smartPtr3D[i][j][1] + 0.5      * smartPtr3D[i][j][2];
				YCbCr_Image(i,j,2) = 128 + 0.5      * smartPtr3D[i][j][0] - 0.418688 * smartPtr3D[i][j][1] - 0.081312 * smartPtr3D[i][j][2];
			}
		}
		return YCbCr_Image;
	}

	Image<unsigned char> YCbCr2RGB()
	{
		assert(typeid(T).name() == typeid(double).name());
		Image<unsigned char> RGB_Image(width, height);

		for(int i = 0; i < height; i++)
		{
			for(int j = 0; j < width; j++)
			{
				// change from double to int 
				double Y = smartPtr3D[i][j][0];
				double Cb = smartPtr3D[i][j][1];
				double Cr = smartPtr3D[i][j][2];


				int r = (int) (Y + 1.40200 * (Cr - 0x80));
				int g = (int) (Y - 0.34414 * (Cb - 0x80) - 0.71414 * (Cr - 0x80));
				int b = (int) (Y + 1.77200 * (Cb - 0x80));

				RGB_Image(i,j,0) = r;
				RGB_Image(i,j,1) = g;
				RGB_Image(i,j,2) = b;

			}
		}
		return RGB_Image;
	}

};


class PPM_Image_Reader{
public:
	Image<unsigned char> m_image;
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
		    	
		    	m_image = Image<unsigned char>(std::stoi(tmp[1]), std::stoi(tmp[2]));
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

class JPGEncoder{
public:
	PPM_Image_Reader m_image_reader;
    QuantizationTable m_q_table;

	JPGEncoder(PPM_Image_Reader& image_reader) 
	{
		m_image_reader = std::move(image_reader);

	}

	/* DCT rows and columns separately
	 *
	 * C(i) = a(i)/2 * sum for x=0 to N-1 of
	 *   s(x) * cos( pi * i * (2x + 1) / 2N )
	 *
	 * if x = 0, a(x) = 1/sqrt(2)
	 *      else a(x) = 1
	 */
	void dct_1d(double *in, double *out, const int count)
	{
		int x, u;

		for (u=0; u<count; u++)
		{
			double z = 0;

			for (x=0; x<count; x++)
			{
				z += in[x] * cos(M_PI * (double)u * (double)(2*x+1)
					/ (double)(2*count));
			}

			if (u == 0) z *= 1.0 / sqrt(2.0);
			out[u] = z/2.0;
		}
	}


	template<class T>
	void dct_2d(int channel, TwoDArray<T>& DCT_Block, const Image<T>& img, const int row_idx, const int col_idx, const int block_size)
	{
		// Code from https://unix4lyfe.org/dct/
		// I use slower code since I can't understand the faster version

		int i,j;
		auto in = std::make_unique<double[]>(block_size);
		auto out = std::make_unique<double[]>(block_size);
		TwoDArray<double> rows(block_size, block_size);
		/* transform rows */
		for (j=0; j<block_size; j++)
		{
			for (i=0; i<block_size; i++){
				if(row_idx+i<img.height && col_idx+j<img.width)
					in[i] = img(row_idx+i, col_idx+j, channel);
				else
					in[i] = 0;
			}
			dct_1d(in.get(), out.get(), block_size);
			for (i=0; i<block_size; i++) rows(j,i) = out[i];
		}

		/* transform columns */
		for (j=0; j<block_size; j++)
		{
			for (i=0; i<block_size; i++)
				in[i] = rows(i,j);
			dct_1d(in.get(), out.get(), block_size);
			for (i=0; i<block_size; i++) DCT_Block(i,j) = out[i];
		}
	}

	#define COEFFS(Cu,Cv,u,v) { \
	if (u == 0) Cu = 1.0 / sqrt(2.0); else Cu = 1.0; \
	if (v == 0) Cv = 1.0 / sqrt(2.0); else Cv = 1.0; \
	}

	template<class T>
	void idct(int channel,const TwoDArray<T>& DCT_Block, Image<T>& img, const int row_idx, const int col_idx, const int block_size)
	{
		int u,v,x,y;
		/* iDCT */
		for (y=0; y<block_size; y++)
			for (x=0; x<block_size; x++)
			{
				double z = 0.0;
				for (v=0; v<block_size; v++)
					for (u=0; u<block_size; u++)
					{
						double S, q;
						double Cu, Cv;
						
						COEFFS(Cu,Cv,u,v);
						S = DCT_Block(v,u);

						q = Cu * Cv * S *
							cos((double)(2*x+1) * (double)u * M_PI/16.0) *
							cos((double)(2*y+1) * (double)v * M_PI/16.0);

						z += q;
					}

				z /= 4.0;
				//When RGB has this constraint, YCbCr is still ok
				if (z > 255.0) z = 255.0;
				if (z < 0) z = 0.0;

				// pixel(tga, x+xpos, y+ypos) = (uint8_t) z;
				if(row_idx+y<img.height && col_idx+x<img.width){
					img(row_idx+y, col_idx+x, channel) = z;
				}
			}
	}

	template<class T>
	void DEBUG_DRAW(Image<T>& debug_YCbCr_img)
	{
		Image<unsigned char> debug_RGB_img = debug_YCbCr_img.YCbCr2RGB();
		std::ofstream ofs("debug.ppm");
		// draw idct image
		// Draw Image
	    std::cout << "=====================" << std::endl << "Start drawing debug image ppm" << std::endl;

	    std::cout << "P3" << std::endl << debug_RGB_img.width << " " << debug_RGB_img.height << std::endl << "255" << std::endl;
	    ofs << "P3" << std::endl << debug_RGB_img.width << " " << debug_RGB_img.height << std::endl << "255" << std::endl;
	    for(int i = 0; i < debug_RGB_img.height; i++)
	        for(int j = 0; j < debug_RGB_img.width; j++)
	        {
	            assert(debug_RGB_img(i,j,0) <= 255);
	            assert(debug_RGB_img(i,j,0) <= 255);
	            assert(debug_RGB_img(i,j,1) <= 255);
	            //std::cout << pic[i][j][0] << " " << pic[i][j][1] << " " << pic[i][j][2] << "\n";
	            ofs << int(debug_RGB_img(i,j,0)) << " " << int(debug_RGB_img(i,j,1)) << " " << int(debug_RGB_img(i,j,2)) << "\n";
	        }
	}

	static double one_norm_distance(const std::vector<double> &A, const std::vector<double>& B){
		assert(A.size() == B.size());
		double sum = 0;
		for(int i=0; i<A.size(); i++){
			sum += abs(A[i] - B[i]);
		}
		return sum;
	}

	template<class T>
	static std::vector<double> left_up_corner(const TwoDArray<T>& blk, double divide_factor){
		// we only fetch first 16 items.
		std::vector<T> v;
		for(int i=0; i<4; i++){
			for(int j=0; j<4; j++){
				v.push_back(blk(i,j) / divide_factor);
			}
		}
		return v;
	}
	template<class T>
	void adaptive_merge(const Image<T>& img, std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks[]){
		// Should divisible by 8. (Padding first)
		assert(img.width % 8 == 0);
		assert(img.height % 8 == 0);
		TwoDArray<char> merge_blk(img.width/8, img.height/8);
		/*
			merge_blk: 
				0 -> 8x8 blocks
				1 -> 16x16 blocks
				2 -> 32x32 blocks
				3 -> 64x64 blocks
				4 -> merged
			Initial: zeros
		*/
		int cul_size[4] = {0};
		for(int i=0; i<img.height; i+=8){
			for(int j=0; j<img.width; j+=8){
				for(int size=3; size>0; size--){
					int blk_size = 8 << size;
					// out of range -> continue
					if(i + blk_size > img.height || j + blk_size > img.width)
						continue;
					// merge the block that merged -> continue
					bool ok_to_merge = true;
					for(int ii=i/8; ii<i/8 + (1<<size) && ok_to_merge; ii++){
						for(int jj=j/8; jj<j/8 + (1<<size) && ok_to_merge; jj++){
							if(merge_blk(ii,jj) != 0){
								ok_to_merge = false;
							}
						}	
					}
					if(!ok_to_merge)
						continue;

					// Calc DCT after merging
					TwoDArray<double> merge_DCT_block(blk_size, blk_size);
					dct_2d(0, merge_DCT_block, img, i, j, blk_size);
					std::vector<double> merge_DCT_info = left_up_corner(merge_DCT_block, (1 << (2*size)));
					double distance = 0;
					for(int ii=i; ii<i+blk_size; ii+=8){
						for(int jj=j; jj<j+blk_size; jj+=8){
							// Calc DCT before merging, and calc the L1 distance.
							TwoDArray<double> single_DCT_block(8, 8);
							dct_2d(0, single_DCT_block, img, ii, jj, 8);
							std::vector<double> single_DCT_info = left_up_corner(single_DCT_block, 1);
							distance += one_norm_distance(merge_DCT_info, single_DCT_info);
						}
					}
					// sum distance -> mean distance
					distance /= (1 << (2*size));
					if(distance < MERGE_THERSHOLD){
						// this block merged.
						for(int ii=i/8; ii<i/8 + (1<<size); ii++){
							for(int jj=j/8; jj<j/8 + (1<<size); jj++){
								// 4: MERGED
								merge_blk(ii,jj) = 4;
							}	
						}
						merge_blk(i/8,j/8) = size;
					}
				}
				if(merge_blk(i/8,j/8) != 4)
					cul_size[merge_blk(i/8,j/8)] ++;
				printf("8: %d\t 16:%d\t 32:%d\t64:%d\t completeness:( %.2f / 100%%) \r", 
					cul_size[0], cul_size[1], cul_size[2], cul_size[3],
					100 * (double)(cul_size[0] + cul_size[1] * 4 + cul_size[2] * 16 + cul_size[3] * 64) / (img.width * img.height / 64));
			}
		}
		puts("Merge Finished");
		// FILE: for visualize
		// FILE *fout = fopen("merge.txt", "w");
		for(int channel=0; channel<3; channel++){
			for(int i=0; i<img.height/8; i++){
				for(int j=0; j<img.width/8; j++){
					// fprintf(fout, "%d", merge_blk(i,j));
					if( merge_blk(i,j) != 4){
						int blk_size = 8 << merge_blk(i,j);
						TwoDArray<double> DCT_Block(blk_size, blk_size);
						dct_2d(channel, DCT_Block, img, i*8, j*8, blk_size);
                        // Test
                        TwoDArray<double> DCT_Block2(8, 8);
                        for(int ii=0; ii<blk_size; ii++){
                            for(int jj=0; jj<blk_size; jj++){
                                DCT_Block(ii,jj) = DCT_Block(ii, jj) / (1 << (2 * merge_blk(i,j)));
                            }
                        }
                        /*
                        TwoDArray<double> DCT_Block3(8, 8);
                        dct_2d(channel, DCT_Block3, img, i*8, j*8, 8);
                        for(int ii=0; ii<8; ii++){
                            for(int jj=0; jj<8; jj++){
                                printf("%.2lf/%.2lf/%.2lf\t", DCT_Block(ii,jj), DCT_Block2(ii,jj), DCT_Block3(ii,jj));
                            }
                        }
                        */
						
						DCT_channel_Blocks[channel].push_back(std::make_tuple(blk_size, i*8, j*8, std::move(DCT_Block)));
                        
                        
                        // Test
//                         DCT_channel_Blocks[channel].push_back(std::make_tuple(blk_size, i*8, j*8, std::move(DCT_Block)));
						// DCT_channel_Blocks[channel].push_back(std::make_tuple(8, i, j, std::move(DCT_Block)));
					
					}
				}
				// fprintf(fout, "\n");
			}
		}	
	}

	template<class T>
	void DCT(const Image<T>& img, std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks[])
	{
		
		Image<T> debug_img(img.width, img.height);

		puts("start DCT");
		// Let's first play with fix size 8, therefore block size not yet useful
		for(int channel = 0; channel < 3; channel++)
		{
			std::cout << "start channel: " << channel << std::endl;
			// first don't consider sampling rate
			// We can only use Y (index 0) then become grayscale image 
			for(int i = 0; i < img.height/8; i++)
			{
				for(int j = 0; j < img.width/8; j++)
				{
					TwoDArray<double> DCT_Block(8, 8);
					dct_2d(channel, DCT_Block, img, i*8, j*8, 8);
                    // for debug
					// idct(channel, DCT_Block, debug_img, i*8, j*8, 8);
					DCT_channel_Blocks[channel].push_back(std::make_tuple(8, i*8, j*8, std::move(DCT_Block)));
				}
			}
		}

		// for debug show
		// DEBUG_DRAW(debug_img);

	}

	void Quantize(std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks[])
	{
        // Y:0 Cb:1 Cr:2
        // Table[0]: Y, Table[1]:Cb,Cr
        int flag = 0;
        for (int i = 0; i < DCT_channel_Blocks[0].size(); i++)
        {
            flag = 0;
            for (int j = 0; j < 4; j++)
            {
                // check 4 table and find the table with correct size
                if( std::get<3>(DCT_channel_Blocks[0][i]).width == m_q_table.DQT[j][0].width && std::get<3>(DCT_channel_Blocks[0][i]).height == m_q_table.DQT[j][0].height)
                {

                    flag = 1;
//                     for(int h = 0; h < std::get<3>(DCT_channel_Blocks[0][i]).height; h++)
//                         for(int w = 0; w < std::get<3>(DCT_channel_Blocks[0][i]).width; w++)
//                             std::get<3>(DCT_channel_Blocks[0][i])(h,w) /= m_q_table.DQT[j][0](h,w); 
                    for( int h = 0; h < 8 ; h++)
                        for(int w = 0; w < 8; w++)
                            std::get<3>(DCT_channel_Blocks[0][i])(h,w) /= m_q_table.DQT[j][0](h,w); 
                    break;
                }
                
                        
           }
           assert(flag);
        }

        

        for (int i = 0; i < DCT_channel_Blocks[1].size(); i++)
        {
            flag = 0;
            for (int j = 0; j < 4; j++)
            {
                if( std::get<3>(DCT_channel_Blocks[1][i]).width == m_q_table.DQT[j][1].width && std::get<3>(DCT_channel_Blocks[1][i]).height == m_q_table.DQT[j][1].height)
                {
                    flag = 1;
//                     for(int h = 0; h < std::get<3>(DCT_channel_Blocks[1][i]).height; h++)
//                         for(int w = 0; w < std::get<3>(DCT_channel_Blocks[1][i]).width; w++)
//                             std::get<3>(DCT_channel_Blocks[1][i])(h,w) /= m_q_table.DQT[j][1](h,w); 
                    for( int h = 0; h < 8 ; h++)
                        for(int w = 0; w < 8; w++)
                            std::get<3>(DCT_channel_Blocks[1][i])(h,w) /= m_q_table.DQT[j][1](h,w); 
                    break;
                }
            }
            assert(flag);
        }

        


        for (int i = 0; i < DCT_channel_Blocks[2].size(); i++)
        {
            flag = 0;
            for (int j = 0; j < 4; j++)
            {
                if( std::get<3>(DCT_channel_Blocks[2][i]).width == m_q_table.DQT[j][1].width && std::get<3>(DCT_channel_Blocks[2][i]).height == m_q_table.DQT[j][1].height)
                {
                    flag = 1;
//                     for(int h = 0; h < std::get<3>(DCT_channel_Blocks[2][i]).height; h++)
//                         for(int w = 0; w < std::get<3>(DCT_channel_Blocks[2][i]).width; w++)
//                             std::get<3>(DCT_channel_Blocks[2][i])(h,w) /= m_q_table.DQT[j][1](h,w); 
                    for( int h = 0; h < 8 ; h++)
                        for(int w = 0; w < 8; w++)
                            std::get<3>(DCT_channel_Blocks[2][i])(h,w) /= m_q_table.DQT[j][1](h,w); 
                    break;
                }
            }
            assert(flag);
        }
        
        
//         Image<double> debug_img(m_image_reader.m_image.width, m_image_reader.m_image.height);
        
//         for ( int channel = 0; channel < 3; channel++)
//         {
//             for (int i = 0; i < DCT_channel_Blocks[channel].size(); i++)
//             {
//                 int block_size = std::get<0>(DCT_channel_Blocks[channel][i]);
//                 TwoDArray<double> DCT_Block(block_size, block_size);
//                 int row_idx = std::get<1>(DCT_channel_Blocks[channel][i]);
//                 int col_idx = std::get<2>(DCT_channel_Blocks[channel][i]);

//                 // copy block
//                 for(int h = 0; h < block_size; h++)
//                     for(int w = 0; w < block_size; w++){                      
//                         DCT_Block(h,w) = std::get<3>(DCT_channel_Blocks[channel][i])(h,w); 
//                     }

//                 idct(channel, DCT_Block, debug_img, row_idx, col_idx, block_size);
//             }
//         }
//         DEBUG_DRAW(debug_img);

	}

	void inv_Quantize(std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks[])
	{
        // Y:0 Cb:1 Cr:2
        // Table[0]: Y, Table[1]:Cb,Cr
		std::cout << "start inv quantize" << std::endl;
        int flag = 0;
        for (int i = 0; i < DCT_channel_Blocks[0].size(); i++)
        {
            flag = 0;
            for (int j = 0; j < 4; j++)
            {
                if( std::get<3>(DCT_channel_Blocks[0][i]).width == m_q_table.DQT[j][0].width && std::get<3>(DCT_channel_Blocks[0][i]).height == m_q_table.DQT[j][0].height)
                {

                    flag = 1;
//                     for(int h = 0; h < std::get<3>(DCT_channel_Blocks[0][i]).height; h++)
//                         for(int w = 0; w < std::get<3>(DCT_channel_Blocks[0][i]).width; w++)
//                             std::get<3>(DCT_channel_Blocks[0][i])(h,w) *= m_q_table.DQT[j][0](h,w); 
                    for(int h = 0; h < 8; h++)
                        for(int w = 0; w < 8; w++)
                            std::get<3>(DCT_channel_Blocks[0][i])(h,w) *= m_q_table.DQT[j][0](h,w); 
                    break;
                }
           }
           assert(flag);
        }
        
        
        for (int i = 0; i < DCT_channel_Blocks[1].size(); i++)
        {
            flag = 0;
            for (int j = 0; j < 4; j++)
            {
                if( std::get<3>(DCT_channel_Blocks[1][i]).width == m_q_table.DQT[j][1].width && std::get<3>(DCT_channel_Blocks[1][i]).height == m_q_table.DQT[j][1].height)
                {
                    flag = 1;
//                     for(int h = 0; h < std::get<3>(DCT_channel_Blocks[1][i]).height; h++)
//                         for(int w = 0; w < std::get<3>(DCT_channel_Blocks[1][i]).width; w++)
//                             std::get<3>(DCT_channel_Blocks[1][i])(h,w) *= m_q_table.DQT[j][1](h,w); 
                    for(int h = 0; h < 8; h++)
                        for(int w = 0; w < 8; w++)
                            std::get<3>(DCT_channel_Blocks[1][i])(h,w) *= m_q_table.DQT[j][1](h,w); 
                    break;
                }
            }
            assert(flag);
        }   
        
        

        for (int i = 0; i < DCT_channel_Blocks[2].size(); i++)
        {
            flag = 0;
            for (int j = 0; j < 4; j++)
            {
                if( std::get<3>(DCT_channel_Blocks[2][i]).width == m_q_table.DQT[j][1].width && std::get<3>(DCT_channel_Blocks[2][i]).height == m_q_table.DQT[j][1].height)
                {
                    flag = 1;
//                     for(int h = 0; h < std::get<3>(DCT_channel_Blocks[2][i]).height; h++)
//                         for(int w = 0; w < std::get<3>(DCT_channel_Blocks[2][i]).width; w++)
//                             std::get<3>(DCT_channel_Blocks[2][i])(h,w) *= m_q_table.DQT[j][1](h,w);
                    for(int h = 0; h < 8; h++)
                        for(int w = 0; w < 8; w++)
                            std::get<3>(DCT_channel_Blocks[2][i])(h,w) *= m_q_table.DQT[j][1](h,w); 
                    break;
                }
            }
            assert(flag);
        }
        
        // debug check
        puts("DEBUG DRAW inv Quantize");
        Image<double> debug_img(m_image_reader.m_image.width, m_image_reader.m_image.height);
        
        for ( int channel = 0; channel < 3; channel++)
        {
            for (int i = 0; i < DCT_channel_Blocks[channel].size(); i++)
            {
                int block_size = std::get<0>(DCT_channel_Blocks[channel][i]);
                TwoDArray<double> DCT_Block(block_size, block_size);
                int row_idx = std::get<1>(DCT_channel_Blocks[channel][i]);
                int col_idx = std::get<2>(DCT_channel_Blocks[channel][i]);

                // copy block
                for(int h = 0; h < block_size; h++)
                    for(int w = 0; w < block_size; w++){                      
                        DCT_Block(h,w) = std::get<3>(DCT_channel_Blocks[channel][i])(h,w); 
                    }

                idct(channel, DCT_Block, debug_img, row_idx, col_idx, block_size);
            }
        }
        DEBUG_DRAW(debug_img);
	}
    
    int check_AC_nextZeros(std::vector<double>& ac_vector, int start_idx)
    {
        int ct = 0;
        for(int i = start_idx; i < ac_vector.size(); i++)
        {
            if(int(ac_vector[i]) != 0){
                break;
            }
            ct++;
        }
        return ct;
    }
    
    void DataStream(std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks[])
    {
        std::cout << "Start DataStream" << std::endl;
        std::vector<bool> bitstream[3];
        std::ofstream m_ofs("output.txt", std::ios_base::out);
        for(int channel = 0; channel < 3; channel++)
            for(int i = 0; i < DCT_channel_Blocks[channel].size(); i++)
            {
                // 00 -> 8x8
                // 01 -> 16x16
                // 10 -> 32x32
                // 11 -> 64x64
                if(std::get<0>(DCT_channel_Blocks[channel][i]) == 8)
                {
                    bitstream[channel].push_back(false);
                    bitstream[channel].push_back(false);
                }
                else if(std::get<0>(DCT_channel_Blocks[channel][i]) == 16)
                {
                    bitstream[channel].push_back(false);
                    bitstream[channel].push_back(true);
                }
                else if(std::get<0>(DCT_channel_Blocks[channel][i]) == 32)
                {
                    bitstream[channel].push_back(true);
                    bitstream[channel].push_back(false);
                }
                else if(std::get<0>(DCT_channel_Blocks[channel][i]) == 64)
                {
                    bitstream[channel].push_back(true);
                    bitstream[channel].push_back(true);
                }
                else{
                   puts("holy assert");
                    assert(1);
                }
            }
        std::cout << "Start DataStream2" << std::endl;
        std::vector<bool> dc_stream[3];
        std::vector<bool> ac_stream[3];
        for(int channel = 0; channel < 3; channel++)
        {
            double previous;
            // get chrominance/luminance huffman table 
            std::vector<bool> *huffman_table_dc = (channel==0)?DCLuminanceToCode:DCChrominanceToCode;
            
            std::map<unsigned char, std::vector<bool>> huffman_table_ac = 
                     (channel==0)?ACLuminanceToCode:ACChrominanceToCode;
            
			int number;
			unsigned current;
            // DPCM (Differential Pulse Coded Modulation)
            for(int i = 0; i < DCT_channel_Blocks[channel].size(); i++)
            {
                std::vector<double> temp = std::get<3>(DCT_channel_Blocks[channel][i]).zigzag();
                
                // Deal with DC
                // std::cout << "[Deal with DC] channel: " << channel  <<"i: " << i << std::endl;
				
				// First one
                if(i == 0)
                {
                //std::cout << "[Deal with DC] channel: aaaaaa" << std::endl;
					number   = (int(temp[0]));
					current  = std::abs(int(temp[0]));
                    previous = temp[0];
                }
                else // differential part
                {
                // std::cout << "[Deal with DC] channel: bbbbbbbbb" << std::endl;
                    number  = (int(temp[i] - previous));
                    current = std::abs(int(temp[i] - previous));
                    previous = temp[i];
					//std::cout << "channel: " << channel << ", cuurent: " << number << std::endl;
				}

                int bit_size = 0;
				// find number of bits
				while(current != 0)
				{
                //std::cout << "[Deal with DC] channel: bb " << current << std::endl;
					bit_size++;
					current = current >> 1;
				}
				// write to dc bitstream
				for(int tt = 0; tt < huffman_table_dc[bit_size].size(); tt++)
					dc_stream[channel].push_back(huffman_table_dc[bit_size][tt]);
				// fixed length encoding
				if(number > 0)
				{
                        //	std::cout << "[Deal with AC] 11 " << std::endl;
					for(int tt = bit_size-1; tt >=0; tt--)
						dc_stream[channel].push_back((bool)(number>>tt)&1);
				}
				else
				{
                        //	std::cout << "[Deal with AC] 10 " << std::endl;
					for(int tt = bit_size-1; tt >=0; tt--)
						dc_stream[channel].push_back((bool)((~number)>>tt)&1);
				}
                
                // Deal with AC
                
                // if all zero
                int blockSize = std::get<0>(DCT_channel_Blocks[channel][i]);
                int numberOfZeros;
                      
                int cursor = 1;
                while (cursor < blockSize*blockSize) {
                    numberOfZeros = check_AC_nextZeros(temp, cursor);
                    //std::cout << "[Deal with AC] cursor: " << cursor << std::endl;
                    //std::cout << "[Deal with AC] num 0" << numberOfZeros << std::endl;
                    if (numberOfZeros == blockSize*blockSize - cursor) {
                        //std::cout << "[Deal with AC] All zero" << std::endl;
                        std::vector<bool> ac_ret = huffman_table_ac[0x00];
                        for (int s = 0; s < ac_ret.size(); s++)
                            ac_stream[channel].push_back( ac_ret[s] );
                        cursor += numberOfZeros;
                    } else if (numberOfZeros >= 16) {
                        //std::cout << "[Deal with AC] >= 16" << std::endl;
                        std::vector<bool> ac_ret = huffman_table_ac[0xf0];
                        for (int s = 0; s < ac_ret.size(); s++)
                            ac_stream[channel].push_back( ac_ret[s] );
                        cursor += 16;
                    } else {
                        // upper 4 bits representing number of zeros
                        //std::cout << "[Deal with AC] test1" << std::endl;
                        cursor += numberOfZeros;

                        int current = std::abs(int(temp[cursor]));
                        int bit_size = 0;
                        // find number of bits
                        while(current != 0)
                        {
                            bit_size++;
                            current = current >> 1;
                        }
                        //std::cout << "[Deal with AC] test2" << std::endl;
                        // write to ac bitstream
                        unsigned char ac_index = (numberOfZeros << 4) + bit_size;
                        for(int tt = 0; tt < huffman_table_ac[ac_index].size(); tt++)
                            ac_stream[channel].push_back(huffman_table_ac[ac_index][tt]);
                        
                        // fixed length encoding
                        int number = int(temp[cursor]);
                        //std::cout << "[Deal with AC] test3: " << bit_size << std::endl;
                        if(number > 0)
                        {
                        	//std::cout << "[Deal with AC] 1 " << std::endl;
                            for(int tt = bit_size-1; tt >=0; tt--) {
                                ac_stream[channel].push_back((bool)(number>>tt)&1);
                        	//std::cout << "[Deal with AC] 2 " << std::endl;
                                
                            }
                        }
                        else
                        {
                        	//std::cout << "[Deal with AC] 3 " << std::endl;
                            for(int tt = bit_size-1; tt >=0; tt--) {
                                ac_stream[channel].push_back((bool)((~number)>>tt)&1);
                        	//std::cout << "[Deal with AC] 4 " << std::endl;
                                
                            }
                        }

                        cursor += 1;
                    }
                    
                }
                //std::cout << "[Deal with AC] End cursor: " << cursor << std::endl;
             
                
            }
        }
        //std::cout << "[Deal with AC] out" << std::endl;
        
        for(int channel = 0; channel < 3; channel++)
        {
        //std::cout << "[Deal with AC] zzz" << std::endl;
            for(int i = 0; i < bitstream[channel].size(); i++)
                m_ofs << bitstream[channel][i]?"1":"0";
            for(int i = 0; i < dc_stream[channel].size(); i++)
                m_ofs << dc_stream[channel][i]?"1":"0";
            for(int i = 0; i < ac_stream[channel].size(); i++)
                m_ofs << ac_stream[channel][i]?"1":"0";
            m_ofs << "\n\n=========================================================\n\n";
        }
        m_ofs << std::endl;
         
    }

	FILE *m_fp;
	unsigned char m_buffer[2];
	unsigned char m_bit_buffer;
	unsigned char m_bit_count;
	int which_one[4] = {8, 16, 32, 64};

	template<class T>
	void decode(const Image<T>& img,
			std::vector< std::tuple<int, int, int, TwoDArray<double> > > tuples[3]) {
		int height = img.height;
		int width = img.width;
		m_fp = fopen("output.txt", "r");
		Scanner scanner(m_fp);

		TwoDArray<bool> bit_Map(img.width, img.height);
		for (int channel = 0; channel < 3; channel++) {
			for (int i = 0; i < img.height; i++)
				for (int j = 0; j < img.width; j++)
					bit_Map(i, j) = 0;
			int filled = 0;
			
			// reconstruct map tuples of (size, i, j)
			int current_cursor = 0;
			int current_width = current_cursor % width;
			int current_height = current_cursor / width;
			std::cout << "decode: 1" << std::endl;
			while (filled < width*height) {
				bool bit1 = scanner.readBit();
				bool bit2 = scanner.readBit();
				int code = bit1 << 1 + bit2;
				std::cout << "decode: 2, " << code << std::endl;
				for (; current_cursor < width*height; current_cursor++) {
					current_width = current_cursor % width;
					current_height = current_cursor / width;
					// std::cout << "decode: 2.5, " << current_width << ", " << current_height << std::endl;
					if (bit_Map(current_height, current_width)) {
						current_cursor++;
					}
					else {
						break;
					}
				}
				int currentSize = which_one[code];
				for (int i = 0; i < currentSize; i++)
					for (int j = 0; j < currentSize; j++)
						bit_Map(current_height+i, current_width+j) = 1;
				TwoDArray<double> tt(currentSize, currentSize);
				tuples[channel].push_back(std::make_tuple(currentSize, current_height, current_width, std::move(tt)));
				filled += currentSize*currentSize;
				std::cout << "decode: 4, " << filled << std::endl;
			}
			
			std::cout << "tuple size: " << tuples[channel].size() << std::endl;
			// scan DC
			// get DC value, then read fixed-length encoded value
			// if first bit is zero --> negative
			for (int i = 0; i < tuples[channel].size(); i++) {
				;
			}
			std::cout << "tuple size1: " << tuples[channel].size() << std::endl;

			// scan AC
			// get AC, get blockSize, if 0x00, 0xf0: special case
			//                      else like DC, 
			for (int i = 0; i < tuples[channel].size(); i++) {
				;
			}
			std::cout << "tuple size2: " << tuples[channel].size() << std::endl;
		}
	}
    
	void run(int mode)
	{
		Image<double> YCbCr_Image = m_image_reader.m_image.RGB2YCbCr();
		puts("finsh RGB2YCbCr");
        std::cout << "original image width:" << m_image_reader.m_image.width << " height: " << m_image_reader.m_image.height << std::endl;
        std::cout << "YCbCr_Image image width:" << YCbCr_Image.width << " height: " << YCbCr_Image.height << std::endl;

		// tuple (block size, row_idx, col_idx, DCT_Block)
		std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks[3]; 
		if (mode)
			adaptive_merge(YCbCr_Image, DCT_channel_Blocks);
		else
			DCT(YCbCr_Image, DCT_channel_Blocks);
        
        Image<double> debug_img(m_image_reader.m_image.width, m_image_reader.m_image.height);
        
        /**for ( int channel = 0; channel < 3; channel++)
        {
            for (int i = 0; i < DCT_channel_Blocks[channel].size(); i++)
            {
                int block_size = std::get<0>(DCT_channel_Blocks[channel][i]);
                TwoDArray<double> DCT_Block(block_size, block_size);
                int row_idx = std::get<1>(DCT_channel_Blocks[channel][i]);
                int col_idx = std::get<2>(DCT_channel_Blocks[channel][i]);

                // copy block
                for(int h = 0; h < block_size; h++)
                    for(int w = 0; w < block_size; w++){                      
                        DCT_Block(h,w) = std::get<3>(DCT_channel_Blocks[channel][i])(h,w); 
                    }

                // idct(channel, DCT_Block, debug_img, row_idx, col_idx, block_size);
            }
        }
        DEBUG_DRAW(debug_img);**/

        
// 		std::get<3>(DCT_channel_Blocks[0][0]).show();
	    Quantize(DCT_channel_Blocks);
        DataStream(DCT_channel_Blocks);
        inv_Quantize(DCT_channel_Blocks);

#ifdef DEBUG
        std::vector<double> hello = std::get<3>(DCT_channel_Blocks[0][1]).zigzag();
        for(int i = 0; i < hello.size(); i++)
        {
            std::cout << hello[i] << " ";
        }
        std::cout << std::endl;
#endif
		std::vector< std::tuple<int, int, int, TwoDArray<double> > > DCT_channel_Blocks_2[3];
		// decode(YCbCr_Image, DCT_channel_Blocks_2);
	}

};

int main(int argc, char **argv)
{
	PPM_Image_Reader reader(argv[1]);
	JPGEncoder jpg_encoder(reader);
	jpg_encoder.run(atoi(argv[2]));

}
