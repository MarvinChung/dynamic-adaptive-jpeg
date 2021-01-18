#include <fstream> 
#include <assert.h>
#include <map>
#include <iostream>
#include <math.h>
#include<vector>
#include <algorithm>

#define DC_SYMBOL 0
#define AC_SYMBOL 1
// reference: https://github.com/MROS/jpeg_decoder/blob/master/main.cpp  

void print_byte(unsigned char read_char)
{
    std::cout << std::hex << (int)read_char << std::endl;
}

struct RGB {
    unsigned char R, G, B;
};

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

typedef struct {
    unsigned char color_id;
    unsigned char h_sample_rate;
    unsigned char v_sample_rate;
    unsigned char quant_id;
}Channel;

typedef struct {
    int width;
    int height;
}Image;

class acCode  {
public:
    unsigned char len;
    unsigned char zeros;
    int value;
    acCode(unsigned char len, unsigned char zeros, int value): len(len), zeros(zeros), value(value) {}
};

typedef double BLOCK[8][8];      

class JPGDecoder
{
public:
    Scanner m_scanner;
    std::ofstream& m_ofs;
    // map store node_index, code -> value
    // Example: 7, 11010 -> 0x31
    // 4 huffman table DC[2] AC[2]
    std::map<std::pair<unsigned char, unsigned int>, unsigned char> m_huffTable[2][2];
    double m_cos_cache[200];
    // At most 4 DQT section, may have 64 or 128 elements.
    int m_quantTable[4][128];
    // Y:1 Cb:2 Cr:3, 0 is redundancy
    Channel m_channel[4];
    Image m_image;
    BLOCK mcu[4][2][2];
    int m_maxHeight;
    int m_maxWidth;
    // Y:1 Cb:2 Cr:3, 0 is redundancy
    int m_DC_values[4];

    // Use for changing huffman 
    // DC id0  DC id1
    std::map<int ,int > m_DC_statistics_table[2];

    JPGDecoder(Scanner& scanner, std::ofstream& ofs): m_scanner(scanner), m_ofs(ofs) {
        init_cos_cache();
        m_maxWidth = 0;
        m_maxWidth = 0;
        m_DC_values[0] = 0;
        m_DC_values[1] = 0;
        m_DC_values[2] = 0;
        m_DC_values[3] = 0;
    }
    void read_jpg();

    void read2end();
    void readAPP();
    void readCOM();
    void readSOF();
    void readDQT();
    void readSOS();
    void readDHT();

    void readBlock(int channel_id, int h, int w);
    void readMCU();

    void readCompressedData();
    int readDC(int channel_id);
    acCode readAC(int channel_id);
    unsigned char findHuffman(int DCorAC, int table_id);


    void init_cos_cache();
    void quantify(); 
    void zigzag();
    void idct();

    RGB **toRGB() {
        RGB **ret = (RGB **)malloc(sizeof(RGB **) * m_maxHeight * 8);
        for (int i = 0; i < m_maxHeight * 8; i++) {
            ret[i] = (RGB *)malloc(sizeof(RGB *) * m_maxWidth * 8);
        }
        for (int i = 0; i < m_maxHeight * 8; i++) {
            for (int j = 0; j < m_maxWidth * 8; j++) {
                double Y = trans(1, i, j);
                double Cb = trans(2, i, j);
                double Cr = trans(3, i, j);
                ret[i][j].R = chomp(Y + 1.402*Cr + 128);
                ret[i][j].G = chomp(Y - 0.34414*Cb - 0.71414*Cr + 128);
                ret[i][j].B = chomp(Y + 1.772*Cb + 128);
            }
        }
        return ret;
    }
    double cc(int i, int j) {
        if (i == 0 && j == 0) {
            return 1.0/2.0;
        } else if (i == 0 || j == 0) {
            return 1.0/sqrt(2.0);
        } else {
            return 1.0;
        }
    }
    double c(int i) {
        static double x = 1.0/sqrt(2.0);
        if (i == 0) {
            return x;
        } else {
            return 1.0;
        }
    }
    unsigned char chomp(double x) {
        if (x > 255.0) {
            return 255;
        } else if (x < 0) {
            return 0;
        } else {
            return (unsigned char) x;
        }
    }
    double trans(int id, int h, int w) {
        int vh = h * m_channel[id].v_sample_rate / m_maxHeight;
        int vw = w * m_channel[id].h_sample_rate / m_maxWidth;
        return mcu[id][vh / 8][vw / 8][vh % 8][vw % 8];
    }

    void show() {
        printf("*************** mcu show ***********************\n");
        for (int id = 1; id <= 3; id++) {
            for (int h = 0; h < m_channel[id].v_sample_rate; h++) {
                for (int w = 0; w < m_channel[id].h_sample_rate; w++) {
                    printf("mcu id: %d, %d %d\n", id, h, w);
                    for (int i = 0; i < 8; i++) {
                        for (int j = 0; j < 8; j++) {
                            printf("%lf ", mcu[id][h][w][i][j]);
                        }
                        printf("\n");
                    }
                }
            }
        }
    };
    
};


void JPGDecoder::read2end()
{
    int length = m_scanner.read2Bytes();
    unsigned char read_char;
    printf("length: %d\n", length);
    length -= 2;
    
    while(length > 0)
    {
        read_char = m_scanner.readByte();
        //print_byte(read_char);
        length -= 1;
    }
    // unsigned read_char;
    // int ct = 0;
    // while(read_char != 0xFF)
    // {
    //  ct++;
    //  read_char = m_scanner.readByte();
    //  print_byte(read_char);
    // }
    // std::cout << ct << std::endl;
    // exit(0);
}

void JPGDecoder::readAPP()
{

    int length = m_scanner.read2Bytes();
    unsigned char read_char;
    int version[2];
    int x_density[2];
    int y_density[2];
    
    char identifier[5];
    for(int i = 0; i < 5; i++)
    {
        read_char = m_scanner.readByte();
        identifier[i] = read_char;
    }

    version[0] = m_scanner.readByte();
    version[1] = m_scanner.readByte();
    m_scanner.seek(1);
    x_density[0] = m_scanner.readByte();
    x_density[1] = m_scanner.readByte();
    y_density[0] = m_scanner.readByte();
    y_density[1] = m_scanner.readByte();

    printf("length: %d\n", length);
    printf("identifier: %s\n", identifier);
    printf("version: %d.%d\n", version[0], version[1]);
    printf("x_density: %d\n", x_density[0] * 16 + x_density[1]);
    printf("y_density: %d\n", y_density[0] * 16 + y_density[1]);
    m_scanner.seek(length-14);
    
}


void JPGDecoder::readCOM()
{
    read2end();
}

void JPGDecoder::readSOF()
{
    //read2end();
    
    int length = m_scanner.read2Bytes();
    unsigned char precision = m_scanner.readByte();
    int height    = m_scanner.read2Bytes();
    int width     = m_scanner.read2Bytes();
    int channel   = m_scanner.readByte();
    std::cout << "image height:" << height << std::endl << "image width:" << width << std::endl << "image channel(YCbCr):" << channel << std::endl;
    m_image.height = height;
    m_image.width = width;
    for(int i = 0; i < 3; i++)
    {
        // color id: 
        //     Y:1, Cb:2, Cr:3
        unsigned char color_id = m_scanner.readByte();
        unsigned char sample_rate = m_scanner.readByte();
        // upper 4 bits : horizontal
        // lower 4 bits : vertical
        unsigned char h_sample_rate = sample_rate  >> 4;
        unsigned char v_sample_rate = sample_rate  & 0x0F;
        unsigned char quant_id      = m_scanner.readByte();
        printf("color_id: %d\n", color_id);
        printf("h_sample_rate: %d\n", h_sample_rate);
        printf("v_sample_rate: %d\n", v_sample_rate);
        printf("quant_id: %d\n\n", quant_id);

        m_channel[color_id].color_id = color_id;
        m_channel[color_id].h_sample_rate = h_sample_rate;
        m_channel[color_id].v_sample_rate = v_sample_rate;
        m_channel[color_id].quant_id      = quant_id;
        m_maxHeight = (m_maxHeight > h_sample_rate ? m_maxHeight : h_sample_rate);
        m_maxWidth = (m_maxWidth > v_sample_rate ? m_maxWidth : v_sample_rate);
    }
}

void JPGDecoder::readDQT()
{
    //read2end();
    int length = m_scanner.read2Bytes();
    std::cout << "length: " << length << std::endl;
    unsigned char read_char;
    length -= 2;
    while (length > 0) {
        // A DQT section may have several DQT tables
        puts("** DQT table:");
        read_char = m_scanner.readByte();
        length--;

        // upper 4 bits : the size of each quantize value
        // lower 4 bits : ID
        unsigned precision = ((read_char>>4) == 0) ? 8: 16; 
        unsigned char id = read_char & 0x0F;
        printf("precision：%d\n", precision);
        printf("Quant ID: %d\n", id);
        // bits to byte
        precision /= 8;
        // remain byte = 64 * precision
        for (int i = 0; i < 64; i++) {
            unsigned char t = 0;
            for (int p = 0; p < precision; p++) {
                read_char = m_scanner.readByte();
                t = t << 8;
                t += read_char;
            }
            m_quantTable[id][i] = t;
        }
        for (int i = 0; i < 64; i++) {
            if (i % 8 == 0) {
                printf("\n");
            }
            printf("%2d ", m_quantTable[id][i]);
        }
        printf("\n");
        length -= (precision*64);
    }
}

void JPGDecoder::readSOS()
{
    //read2end();

    int length = m_scanner.read2Bytes();
    printf("length: %d\n", length);
    // color_nums should be fix to 3
    unsigned char color_nums = m_scanner.readByte();
    assert(color_nums == 3);
    for (int i = 0; i < 3; i++)
    {
        unsigned char color_id = m_scanner.readByte();
        unsigned char huffman_id = m_scanner.readByte();
        unsigned char DC       = huffman_id >> 4;
        unsigned char AC       = huffman_id & 0x0F; 
        printf("color_id: %d\n", color_id);
        printf("DC: %d\n", DC);
        printf("AC: %d\n\n", AC);
    }
    // jpg baseline standard won't use this
    m_scanner.seek(3);
}

void JPGDecoder::readDHT()
{
    // Assume zero extension
    int length = m_scanner.read2Bytes();
    printf("length: %d\n", length);
    length -= 2;
    while (length > 0) {
        unsigned char read_char = m_scanner.readByte();
        int type = (read_char >> 4);
        printf(type == 0 ? "DC\n" : "AC\n");
        unsigned char id = read_char & 0x0F;
        printf("ID: %d\n", id);

        int leaf_numbers = 0;
        int node_height[1024] = {0};
        int ct = 0;
        for (int i = 0; i < 16; i++) {
            read_char = m_scanner.readByte();
            leaf_numbers += read_char;
            //print_byte(read_char);
            while(read_char > 0)
            {   
                node_height[ct++] = i+1;
                read_char--;
            }
        }

        // ct can be less than 16
        assert(ct == leaf_numbers);
        int huffCode = 0;
        for (int i = 0; i < leaf_numbers; i++) {
            read_char = m_scanner.readByte();
            if (i != 0)
                huffCode = (huffCode + 1) << (node_height[i] - node_height[i-1]);
            m_huffTable[type][id][std::make_pair<unsigned char, unsigned int>(node_height[i], huffCode)] = read_char;
            printf("node_id:%4d height:%4d huffCode:%4d: value:%4d\n", i+1, node_height[i], huffCode, read_char);
        }

        length -= (1 + 16 + leaf_numbers);
    }
}

unsigned char JPGDecoder::findHuffman(int DCorAC, int table_id)
{
    //printf("DCorAC:%d\n", DCorAC);
    // DC:0 AC:1
    // huffman tree height max : 16
    int read_bit = 0; 
    unsigned char codeLen;
    for(int count = 1; count <= 16; count++)
    {
        read_bit = read_bit << 1;
        read_bit += (unsigned int)m_scanner.readBit();
        // printf("read_bit:%d\n",read_bit);
        if (m_huffTable[DCorAC][table_id].find(std::make_pair(count, read_bit)) != m_huffTable[DCorAC][table_id].end()) {
            codeLen = m_huffTable[DCorAC][table_id][std::make_pair(count, read_bit)];
            return codeLen;
        }
    } 
    fprintf(stderr, "key not found\n");
    exit(1);          
}

int JPGDecoder::readDC(int channel_id)
{
    // SSSS + DIFF
    // Y: 0, Cb:1, Cr:2
    // huffman table:
    // Y      -> table_id 0
    // Cb, Cr -> table_id 1
    int codeLen = findHuffman(DC_SYMBOL, channel_id/2);

    if (codeLen == 0)
        return 0;
    unsigned char read_bit;
    int ret = 1;

    // Use one's complement for DC
    unsigned char sign_bit = m_scanner.readBit();
    for (int i = 1; i < codeLen; i++)
    {
        read_bit = m_scanner.readBit();
        ret = ret << 1;
        // if sign_bit = 0 -> negative
        ret += sign_bit ? read_bit : !read_bit;
    }
    ret = sign_bit ? ret: -ret;

    // calc how many times ret has occur to see if the huffman table is good
    // auto iter = m_DC_statistics_table[channel_id/2].find(ret);
    // //std::cout << iter << std::endl;
    // if(iter != m_DC_statistics_table[channel_id/2].end())
    //     iter->second += 1;
    // else
    //     m_DC_statistics_table[channel_id/2].insert(std::pair<int, int>(ret, 1));
    
    return ret;
}

acCode JPGDecoder::readAC(int channel_id)
{
    // ZZ + DIFF
    // Y: 0, Cb:1, Cr:2
    // huffman table:
    // Y      -> table_id 0
    // Cb, Cr -> table_id 1
    unsigned char read_byte = findHuffman(AC_SYMBOL, channel_id/2);
    // upper 4 bits : the numbers of contiguous zeros next
    // lower 4 bits : how many bits should be read after zeros
    
    int zeros = read_byte >> 4;
    int codeLen = read_byte & 0x0F;
    if (read_byte == 0) {
        // When read 0x00 -> All zeros
        return acCode(0,0,0);
    } else if (read_byte == 0xF0) {
        // When read 0xF0 -> 16 zeros
        return acCode(0, 16, 0);
    }
    // Use one's complement for AC
    unsigned char sign_bit = m_scanner.readBit();
    unsigned char read_bit;
    int ret = 1;
    for (int i = 1; i < codeLen; i++)
    {
        read_bit = m_scanner.readBit();
        ret = ret << 1;
        // if sign_bit = 0 -> negative
        ret += sign_bit ? read_bit : !read_bit;
    }
    ret = sign_bit ? ret: -ret;
    return acCode(codeLen, zeros, ret);
}

void JPGDecoder::readBlock(int channel_id, int h, int w)
{
    m_DC_values[channel_id] = readDC(channel_id) + m_DC_values[channel_id];
    mcu[channel_id][h][w][0][0] = m_DC_values[channel_id];
    unsigned int count = 1;
    while (count < 64) {
        acCode ac = readAC(channel_id);
        if (ac.len == 0 && ac.zeros == 16) {
            for (int j = 0; j < ac.zeros; j++) {
                mcu[channel_id][h][w][count/8][count%8] = 0;
                count++;
            }
        } else if (ac.len == 0) {
            break;
        } else {
            for (int j = 0; j < ac.zeros; j++) {
                mcu[channel_id][h][w][count/8][count%8] = 0;
                count++;
            }
            mcu[channel_id][h][w][count/8][count%8] = ac.value;
            count++;
        }
    }
    while (count < 64) {
        mcu[channel_id][h][w][count/8][count%8] = 0;
        count++;
    }
}

void JPGDecoder::readMCU()
{
    // Y, Cb, Cr
    for(int i = 1; i <= 3; i++)
    {
        for (int h = 0; h < m_channel[i].v_sample_rate; h++) 
        {
            for (int w = 0; w < m_channel[i].h_sample_rate ; w++) 
            {
                readBlock(i, h, w);
            }
        }
    }
}

void JPGDecoder::readCompressedData()
{
    int pic[500][500][3];

    // Sequentially draw image by reading each MCU
    int h = (m_image.height-1)/(8 * m_maxHeight) + 1;
    int w = (m_image.width-1)/(8 * m_maxWidth) + 1;
    for (int i = 0; i < h; i++)
    {
        for ( int j = 0; j < w; j++)
        {
            readMCU();
            quantify();
            zigzag();
            idct();
            RGB **b = toRGB();
            // Every block color
            for (int y = i*8*m_maxHeight; y < (i+1)*8*m_maxHeight; y++) {
                for (int x = j*8*m_maxWidth; x < (j+1)*8*m_maxWidth; x++) {
                    int by = y - i*8*m_maxHeight;
                    int bx = x - j*8*m_maxWidth;
                    // printf("mm_y:%d   mm_x:%d\n",(i+1)*8*m_maxHeight,(j+1)*8*m_maxWidth);
                    // printf("R:%d G:%d B:%d\n",b[by][bx].R, b[by][bx].G, b[by][bx].B);
                    // printf("y:%d, x:%d m_image.height:%d, m_image.width:%d\n",y,x,m_image.height,m_image.width);
                    pic[y][x][0] = b[by][bx].R;
                    pic[y][x][1] = b[by][bx].G; 
                    pic[y][x][2] = b[by][bx].B;

                }
            }
            //show();
            //exit(0);
        }
    }

    //Draw Image
    std::cout << "=====================" << std::endl << "Start drawing image ppm" << std::endl;

    std::cout << "P3" << std::endl << m_image.width << " " << m_image.height << std::endl << "255" << std::endl;
    m_ofs << "P3" << std::endl << m_image.width << " " << m_image.height << std::endl << "255" << std::endl;
    for(int i = 0; i < m_image.height; i++)
        for(int j = 0; j < m_image.width; j++)
        {
            assert(pic[i][j][0] <= 255);
            assert(pic[i][j][1] <= 255);
            assert(pic[i][j][2] <= 255);
            //std::cout << pic[i][j][0] << " " << pic[i][j][1] << " " << pic[i][j][2] << "\n";
            m_ofs << pic[i][j][0] << " " << pic[i][j][1] << " " << pic[i][j][2] << "\n";
        }
    // printf("m_maxHeight:%d m_maxWidth:%d\n", m_maxHeight, m_maxWidth);

}   


void JPGDecoder::read_jpg()
{
    unsigned char read_char;
    read_char = m_scanner.readByte();
    while(read_char == 0xFF and !m_scanner.check_reach_EOF())
    {
        read_char = m_scanner.readByte();
        switch (read_char) {
            case 0xD8:
                puts("===================");
                puts("readSOI, Start of Image");
                break;
            case 0xE0:
            case 0xE1:
            case 0xE2:
            case 0xE3:
            case 0xE4:
            case 0xE5:
            case 0xE6:
            case 0xE7:
            case 0xE8:
            case 0xE9:
            case 0xEA:
            case 0xEB:
            case 0xEC:
            case 0xED:
            case 0xEE:
            case 0xEF:
                puts("===================");
                puts("readAPP, application 0"); 
                readAPP();
                break;
            case 0xFE:
                puts("===================");
                puts("readCOM");
                readCOM();
                break;
            case 0xDB:
                puts("===================");
                puts("readDQT");
                readDQT();
                break;
            case 0xC0:
                puts("===================");
                puts("readSOF");
                readSOF();
                break;
            case 0xC4:
                puts("===================");
                puts("readDHT");
                readDHT();
                break;
            case 0xDA:
                puts("===================");
                puts("readSOS");
                readSOS();
                // Last step, reads all MCU
                readCompressedData();
                break;
            case 0xD9:
                break;
            default:
                puts("not");
        
        }
        //print_byte(read_char);
        read_char = m_scanner.readByte();
        //print_byte(read_char);
    }
}

void JPGDecoder::quantify() {
    for (int id = 1; id <= 3; id++) {
        for (int h = 0; h < m_channel[id].v_sample_rate; h++) {
            for (int w = 0; w < m_channel[id].h_sample_rate; w++) {
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        mcu[id][h][w][i][j] *= m_quantTable[m_channel[id].quant_id][i*8 + j];
                    }
                }
            }
        }
    }
}

void JPGDecoder::zigzag() {
    for (int id = 1; id <= 3; id++) {
        for (int h = 0; h < m_channel[id].v_sample_rate; h++) {
            for (int w = 0; w < m_channel[id].h_sample_rate; w++) {
                int zz[8][8] = {
                        { 0,  1,  5,  6, 14, 15, 27, 28},
                        { 2,  4,  7, 13, 16, 26, 29, 42},
                        { 3,  8, 12, 17, 25, 30, 41, 43},
                        { 9, 11, 18, 24, 31, 40, 44, 53},
                        {10, 19, 23, 32, 39, 45, 52, 54},
                        {20, 22, 33, 38, 46, 51, 55, 60},
                        {21, 34, 37, 47, 50, 56, 59, 61},
                        {35, 36, 48, 49, 57, 58, 62, 63}
                };
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        zz[i][j] = mcu[id][h][w][zz[i][j] / 8][zz[i][j] % 8];
                    }
                }
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        mcu[id][h][w][i][j] = zz[i][j];
                    }
                }
            }
        }
    }
};


void JPGDecoder::idct() {
    for (int id = 1; id <= 3; id++) {
        for (int h = 0; h < m_channel[id].v_sample_rate; h++) {
            for (int w = 0; w < m_channel[id].h_sample_rate; w++) {
                double tmp[8][8] = {0};
                double s[8][8] = {};
                for (int j = 0; j < 8; j++) {
                    for (int x = 0; x < 8; x++) {
                        for (int y = 0; y < 8; y++) {
                            s[j][x] += c (y) * mcu[id][h][w][x][y] * m_cos_cache[(j + j + 1) * y];
                        }
                        s[j][x] = s[j][x] / 2.0;
                    }
                }
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        for (int x = 0; x < 8; x++) {
                            tmp[i][j] += c(x) * s[j][x] * m_cos_cache[(i + i + 1) * x];
                        }
                        tmp[i][j] = tmp[i][j] / 2.0;
                    }
                }
                for (int i = 0; i < 8; i++) {
                    for (int j = 0; j < 8; j++) {
                        mcu[id][h][w][i][j] = tmp[i][j];
                    }
                }
            }
        }
    }
}



void JPGDecoder::init_cos_cache() {
    for (int i = 0; i < 200; i++) {
        m_cos_cache[i] = cos(i * M_PI / 16.0);
    }
}

// class HuffmanChanger{
// public:
//     JPGDecoder &m_decoder;
//     std::vector<std::pair<int, int>> m_sorted_list[2];
//     std::vector<std::pair<std::pair<unsigned char, unsigned int>, unsigned char>> m_huffman_list[2];
//     std::map<std::pair<unsigned char, unsigned int>, unsigned char> m_new_huffTable[2][2];
//     std::map<int, int> value_mapper[2];

//     HuffmanChanger(JPGDecoder &decoder) : m_decoder(decoder) {}

//     void change_DC_table();
//     //static bool my_compare1(const std::pair<int, int> &p1, const std::pair<int, int> &p2);
//     //static bool my_compare2(const std::pair<std::pair<unsigned char, unsigned int>, unsigned char> &p1, const std::pair<std::pair<unsigned char, unsigned int>, unsigned char> &p2);
// };


// bool my_compare1(const std::pair<int, int> &p1, const std::pair<int, int> &p2){
//     return p1.second > p2.second;
// }

// bool my_compare2(const std::pair<std::pair<unsigned char, unsigned int>, unsigned char> &p1, const std::pair<std::pair<unsigned char, unsigned int>, unsigned char> &p2){
//     return p1.first.first > p2.first.first;
// }

// void HuffmanChanger::change_DC_table()
// {   
//     for(int i = 0; i < 2; i++)
//     {
//         for(auto itr = m_decoder.m_DC_statistics_table[i].begin(); itr != m_decoder.m_DC_statistics_table[i].end(); ++itr){
//             m_sorted_list[i].push_back(std::make_pair(itr->first, itr->second));
//         }
//         std::sort(m_sorted_list[i].begin(), m_sorted_list[i].end(), my_compare1);
//     }

//     for(int i = 0; i < 2; i++)
//     {
//         for(auto itr = m_decoder.m_huffTable[DC_SYMBOL][i].begin(); itr != m_decoder.m_huffTable[DC_SYMBOL][i].end(); ++itr){
//             //node_height, huffCode, value
//             m_huffman_list[i].push_back(std::make_pair(std::make_pair(itr->first.first, itr->first.second), itr->second));
//         }
//         std::sort(m_huffman_list[i].begin(), m_huffman_list[i].end(), my_compare2);
//     }

//     for(int i = 0; i < 2; i++)
//     {

//     }
    
// }

int main(int argc, char *argv[])
{
    std::ofstream m_ofs("output.ppm", std::ios_base::out);
    puts("please ensure ./a.out img_name");
    FILE *fp;
    fp = fopen(argv[1], "r");
    Scanner scanner(fp);
    JPGDecoder decoder(scanner, m_ofs);
    decoder.read_jpg();
    std::cout << "finish program" <<std::endl <<  "output.ppm is ready" <<std::endl;
    // HuffmanChanger changer(decoder);
    return 0;
}