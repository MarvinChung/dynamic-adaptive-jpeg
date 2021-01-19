import matplotlib.pyplot as plt
import cv2
import numpy as np
np.set_printoptions(2, suppress=True)

# img = cv2.imread('sample_images/lenna.png')
img = cv2.imread('dynamic-adaptive-jpeg/balls.ppm')
hsv = cv2.cvtColor(img, cv2.COLOR_BGR2YCR_CB)
gray = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY) 
edges = cv2.Canny(gray, 150, 0)

Y = hsv[:, :, 0]
# Y = edges
img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
# plt.imshow(edges, cmap='gray')
# plt.imshow(img)
# plt.show()

# exit()
# plt.imshow(Y, cmap='gray')
plt.imshow(img)
# plt.show()
print(cv2.dct(Y[:16, :16].astype('float')))

exit()

height, width, channel = img.shape
dw8 = np.zeros([height // 8, width // 8, 16])
for i in range(0, height, 8):
    for j in range(0, width, 8):
        # block = np.mean(Y[i:i+8, j:j+8])
        
        img_dct = cv2.dct(Y[i:i+8, j:j+8].astype('float'))[:4, :4]
        img_dct_log = np.log(abs(img_dct))
        img_dct_log = img_dct
        img_dct_log = img_dct_log.reshape([-1])        
        dw8[i//8, j//8] = img_dct_log
        # print("%.1f" % float(img_dct_log[0]) )
        # print(img_dct_log)
        # print(img_dct)
        # exit()
        # plt.annotate("%.1f" % float(img_dct_log[0]) , (i, j))
# for i in range(8, height, 8):
#     plt.plot([i, i], [0, width-1], color='green')
# for i in range(8, width, 8):
#     plt.plot([0, height-1], [i, i], color='green')

for i in range(0, height-1, 16):
    for j in range(0, width-1, 16):
        blocks = []
        big_img = Y[i : i + 16, j : j + 16]
        bigimg_dct = cv2.dct(big_img.astype('float'))[:4, :4].reshape([-1]) / 2
        ds = []
        for di in [0, 1]:
            for dj in [0, 1]:
                # print([i//8+di, j//8+dj])
                blocks.append(dw8[i//8+di, j//8+dj])
                ds.append(np.linalg.norm(dw8[i//8+di, j//8+dj] - bigimg_dct , 1))
                # print(dw8[i//8+di, j//8+dj])
        # print(bigimg_dct)
        # exit()
        # d01 = np.linalg.norm(blocks[0] - blocks[1], 1)
        # d02 = np.linalg.norm(blocks[0] - blocks[2], 1)
        # d13 = np.linalg.norm(blocks[1] - blocks[3], 1)
        # d23 = np.linalg.norm(blocks[2] - blocks[3], 1)
        # ds =  [d01, d02, d13, d23]
        # d = d01 + d02 + d13 + d23
        # print(ds, sum(ds), max(ds))
        # for b in blocks:
        #     print(b)
        # exit()
        if  max(ds) < 200:
            plt.plot([j, j+16, j+16, j, j],[i, i, i+16, i+16, i],  color='blue')
        else:
            plt.plot([j, j+16],[i, i+16], color='red')



plt.show()