# Project Proposal: Dynamic Adaptive JPEG

> Team C:
> R09922028 施承志
>
> R09922071 劉俊緯
>
> R09922091 張庭與
>
> R09944021 鍾起鳴


## Idea Origin

JPEG codec divides image into 8x8 blocks for encoding. However, real-world images often constains large blocks with few details. If we can design a codec which takes image detail into consider, it could improve the compression rate of an real-life image.


## Base Idea

- Modify JPEG so that it can divide image into coding blocks of different sizes based on local image detail.


## Challenge

### Encoding

- An additional table for block division records.
- Larger DC difference between larger blocks


### Encoder

- Image detail recognition
- DCT for different block sizes
- Different quantisation


### Decoder

- DCT for different block sizes


## Expected effect

- Higher compression ratio (i.e. smaller image size) for general images
- Not much more computational costs
- Not much worse in special cases
