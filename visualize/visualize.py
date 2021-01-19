import matplotlib.pyplot as plt
import cv2
import numpy as np
np.set_printoptions(2, suppress=True)

# img = cv2.imread('sample_images/lenna.png')
img = cv2.imread('dynamic-adaptive-jpeg/balls.ppm')
plt.imshow(img)

M = np.zeros(img.shape[:2])
print(M.shape)
with open("dynamic-adaptive-jpeg/merge.txt", 'r') as fin:
    i = -1
    for row in fin:
        i += 1
        for j, cell in enumerate(row.strip()):
            cell = int(cell)
            if cell != 4:
                cell_to_color = ['red', 'orange', 'blue', 'green']
                blk_size = 8 * (2 ** cell)
                ii = i * 8
                jj = j * 8
                # print(i, j, blk_size, [jj, jj+blk_size, jj+blk_size, jj, jj],[ii, ii, ii+blk_size,ii+blk_size, ii])
                if cell != 0:
                    plt.plot([jj+1, jj+blk_size-1, jj+blk_size-1, jj+1, jj+1],[ii+1, ii+1, ii+blk_size-1,ii+blk_size-1, ii+1],  color=cell_to_color[cell])
                else:
                    plt.plot([jj+1, jj+blk_size-1],[ii+1, ii+blk_size-1], color=cell_to_color[cell])

# print(M)
plt.show()