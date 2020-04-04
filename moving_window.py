# Moving window analysis for raster comparison
# Author: Owen Smith

import math
import numpy as np

# i, j = y, x
# a1 = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1]])
#
# a2 = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1],
#                [1, 1, 1, 1, 1, 1, 1, 1, 1]])
# w = 3


def neighbors(im, i, j, d=1):
    n = im[i-d:i+d+1, j-d:j+d+1].flatten()
    return n


def moving_window(array_1, array_2, window):
    for i in range(array_1.shape[0]):
        for j in range(array_1.shape[1]):
            height, width = array_1.shape
            tw = (height * width) / (window * window)
            s = math.fsum([q for q in range(math.ceil(tw))])
            p = math.fsum([array_1[i, j]-array_2[i, j] for i in range(
                     2)])
            fw = (1 / tw) * s * (1 - (p / (2 * window**2)))
            return fw / tw


def multi_mw(array_1, array_2, max_window):
    co = []
    for i in range(1, max_window):
        co.append(moving_window(array_1, array_2, i))
    return co





