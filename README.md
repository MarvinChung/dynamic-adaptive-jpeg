# Dynamic Adaptive JPEG

> Team C:
> 
> R09922028 施承志
> 
> R09922071 劉俊緯
> 
> R09922091 張庭與
> 
> R09944021 鍾起鳴


# Introduction

JPEG codec divides image into 8x8 blocks for encoding. However, real-world images often constains large blocks with few details. We would like to design a new variant of JPEC codec that divides image based on its local detail.

# Methods

## 1. Block Division

### 1.1. [Canny edge detection algorithm](https://en.wikipedia.org/wiki/Canny_edge_detector)

Canny sometimes fails to detect smoother edges such as the pillar on the left of the first sample. However, this allows us to not mark the moderately less detailed blocks for further merging.

| Sample | Detected Edge | Initial Block Division |
|:-:|:-:|:-:|
| ![Sample 1](https://i.imgur.com/ACEn8eu.png) | ![Detected edge from sample 1](https://i.imgur.com/LUwrYdU.png) | ![Block division of sample 1](https://i.imgur.com/FVBTEgI.png) |
| ![Sample 2](https://i.imgur.com/nSmpV30.jpg) | ![Detected edge from sample 2](https://i.imgur.com/jQO1Tyg.png) | ![Block division of sample 2](https://i.imgur.com/UBMbuID.jpg) |

- Red stripes: 8x8 blocks with edge
- Blue squares: 8x8 blocks without any edge


### 1.2. Block Merging

Contiguous mergeable 8x8 blocks are merged into larger blocks (16x16, 32x32 or 64x64) if allowed.

- Legend:
    - Red stripes: 8x8 blocks that are not merged
    - Orange squares: 16x16 merged blocks
    - Blue squares: 32x32 merged blocks
    - Green squares 64x64 merged blocks

![](https://i.imgur.com/NtznotY.jpg)


### 1.3. Block Encoding

- 4 different block size are represented with 2 bits
- Blocks are represented with a bitstream, in a top-to-bottom, left-to-right manner.
- Already-coded blocks are skipped for cleaness and unambiguity.

![](https://imgur.com/KD9rD5j.png)![](https://i.imgur.com/kt2Nx9R)![](https://i.imgur.com/kZiOY26)

- Encoded bitstream: `00 10 01 00 00 00 00 00 00 00 01 00 01 10 00 ...`


## 2. DCT

DCT for blocks with different size is trivial as the argument of DCT algorithm is modifiable.


## 3. Quantisation

### 3.1. Table Scaling

<table>
<tr>
<th style="text-align:center">Standard luminance table</th>
<th style="text-align:center">Standard chrominance table</th>
</tr>
<tr>
<td>

```
16  11  10  16  24  40  51  61
12  12  14  19  26  58  60  55
14  13  16  24  40  57  69  56
14  17  22  29  51  87  80  62
18  22  37  56  68 109 103  77
24  35  55  64  81 104 113  92
49  64  78  87 103 121 120 101
72  92  95  98 112 100 103  99
```
</td>
<td>

```
17  18  24  47  99  99  99  99
18  21  26  66  99  99  99  99
24  26  56  99  99  99  99  99
47  66  99  99  99  99  99  99
99  99  99  99  99  99  99  99
99  99  99  99  99  99  99  99
99  99  99  99  99  99  99  99
99  99  99  99  99  99  99  99
```
</td>



</table>

Quantisation tables for larger blocks are upsampled from the standard tables.
```C++
upscaled_table = ((std_table * qualityScale + 50) / 100).clip(1,255)
quantized_value = DCT(block) / upscaled_table
```
<img src="https://i.imgur.com/hlBIeIB.png" height="200px">


The range of DC value from larger blocks are greater than that from smaller blocks. Hence, nomalisation of DC values from blocks with different sizes is performed to keep the values in the same scale.


### 3.2. [Zigzagging](https://stackoverflow.com/questions/3025595/code-golf-zigzag-pattern-scanning)

Line-based encoding introduces the potential 'great hop' of DC value between lines and reduces the efficiency of dDC table. We introduce Zigzagging with the hope to reduce 'great hop' from happening.

Zigzagging index tables for the 4 different block sizes are pre-generated and coded in.

## 4. Intermediate Bitstream

The intermediate bitstream format in our codec is very similiar with that in the JPEG codec.
```
BlockSizes(Y)  | DC(Y)  | AC(Y)  |
BlockSizes(Cb) | DC(Cb) | AC(Cb) |
BlockSizes(Cr) | DC(Cr) | AC(Cr) | 
```

### 4.1. Encoding (Implemented)

- DCT, DQT, zigzag
- Block sizes
    - Written in the format mentioned in section 1.3.
- DC
    - Kept for the first block, and is encoded with DPCM (Differential Pulse Coded Modulation) for the rest of the blocks. 
- AC
    - Scanned with the zigzag pattern, as mentioned in section 2.2.
    - Values:
        - `0xXY`: Next X values are 0s, and the value after the 0s are Y bits long.
        - `0x00`: All remaining values are 0s
        - `0xF0`: Next 16 values are 0s

### 4.2. Decoding (Not Implemented)

- Read blocks sizes
    - Reconstructed from in a top-to-bottom, left-to-right manner (using a bitmap to check whether the block is placed)
    - Skipping already-recognized spaces
- Read DC
    - Read bitstream until it matches one of the value of DC table
    - Then, read number of bit according to the index of matched code
    - Use 1's complement representation to reconstruct the origin DC or DPCM value
- Read AC
    - Read bitstream until it matches one of the value of AC table
    - According to the code, we can get the bytes representing it
    - We parse the byte in this manner
        - For the upper 4 bits, we can reconstruct number of zeros of AC values
        - For the lower 4 bits, we will read number of bits according to it, and reconstruct the origin AC value
- Inverse zigzag, IDQT, then IDCT


## 5. Comparison

> GIMP JPEG spec:
> 
> Quality 90,
> 
> Subsample: 4:4:4
> 
> DCT: Integer
> 
> Optimized, Progressive

- Comparison of different methods (Variation with quality 80)

| File Name  | Original PPM<br>(Bytes) | Variation-8x8<br>(Bytes) | Variation-Adaptive<br>(Bytes) | GIMP JPEG<br>(Bytes) |
| ---------- | ------------ | ------------- | ------------------ | --------- |
| [balls.ppm](https://imgur.com/wS45aXQ.png)  | 10880560  (100%) | 200338  (1.8%)<br>[Image link](https://imgur.com/lTnFqoc.png) | 260329  (2.4%)<br>[Image link](https://imgur.com/IhiK4GU.png) | 260170  (2.4%)<br>[Image link](https://imgur.com/Bz4pZIj.jpg) |
| [test_1.ppm](https://imgur.com/1Nl7lSm.png) | 90426514  (100%) | 1814081  (2.0%)<br>[Image link](https://imgur.com/z0h6iWj.png) | 1095085  (1.2%) | 1708280  (1.9%) |
| [test_2.ppm](https://imgur.com/TsrvnmX.png) | 96039186  (100%) | 2377799  (2.5%)<br>[Image link](https://imgur.com/JquVfLq.png) | 2140897  (2.2%) | 3022289  (3.1%) |
| [test_3.ppm](https://imgur.com/x0fjol7.png) | 89758153  (100%) | 1569522  (1.7%)<br>[Image link](https://imgur.com/zGJuQ0Z.png) | 752873  (0.8%) | 1215101  (1.3%) |


The compress version of all test cases with 8x8 blocks and quality 80 does not deviate from the original ones by over 1 in every RGB channel in every pixel.

In the case of `balls.ppm`, the adaptive compressed image is larger than the 8x8 compressed one. On the other hand, all other test cases shows a decline in file size from the 8x8 compression to the adaptive one.

<table>
<tr>
<th style="text-align:center">Original ppm</th>
<th style="text-align:center">Variation-8x8</th>
</tr>
<tr>
<td><img src="https://imgur.com/wS45aXQ.png"/></td>
<td><img src="https://imgur.com/lTnFqoc.png"/></td>
</tr>
<tr>
<th style="text-align:center">Variation-Adaptive</th>
<th style="text-align:center">GIMP JPEG</th>
</tr>
<tr>
<td><img src="https://imgur.com/IhiK4GU.png"/></td>
<td><img src="https://imgur.com/Bz4pZIj.jpg"/></td>
</tr>
</table>

- Comparison of different quality using 8x8 (quantization table coefficients)

| File Name  | original | 100 | 80 | 60 |
| ---------- | -------- | ------------ | ------------- | ------------------ | --------- |
| balls.ppm | 10880560 (100%) | 1071485 (9.8%) | 200338 (1.8%) | 161873 (1.5%) | 
| test_1.ppm | 90426514 (100%) | 6145913 (6.8%) | 1814081 (2.0%) | 1490622 (1.6%) | 
| test_2.ppm | 96039186 (100%) | 10951132 (11.4%) | 2377799 (2.5%) | 1769368 (1.8%) | 
| test_3.ppm | 89758153 (100%) | 4672494 (5.2%) | 1569522 (1.7%) | 1256616 (1.4%) | 

We support different quality of quantization table. We can see the table list above that the data compression rate increase as we lower the quality. Despite of the fact that we lower the quality, it is hard to discriminate against different quality of images as we show below.

<table>
<tr>
<th>Original Image</th>
<th>Quality 100</th>
</tr>
<tr>
<td><img src="https://imgur.com/wS45aXQ.png"/></td>
<td><img src="https://i.imgur.com/YtLa9Hq.jpg"/></td>
</tr>
<tr>
<th>Quality 80</th>
<th>Quality 60</th>
</tr>
<tr>
<td><img src="https://imgur.com/lTnFqoc.png"/></td>
<td><img src="https://i.imgur.com/TpE8nvb.jpg"/></td>
</tr>
</table>

# Conclusion

The compress version of all test cases with 8x8 blocks and quality 80 does not deviate from the original ones by over 1 in every RGB channel in every pixel, which we consider as a success.

The adaptive compression suffers from severe blocking effect and unnatural artifacts. After inspecting the compressed image, we found that our original thoughts about Canny algorithm is not good enough and that the larger blocks ruins perceived image quality very badly

However, we still learn from the process and become more familiar with the algorithms and maths used in JPEG.


# Appendix A. References

- Class handouts

