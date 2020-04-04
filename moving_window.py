################################################################################
# Title: moving_window.py
# Author: Owen Smith, University of North Georgia
# Purpose: Functions to create a moving window window comarison coeficent for
#          array and raster analysis.
################################################################################

import math
import numpy as np


def neighbors(im, i, j, d=1):
    n = im[i-d:i+d+1, j-d:j+d+1].flatten()
    return n


# TODO: Refactor for raster analysis
def moving_window(array_1, array_2, window):
    for i in range(array_1.shape[0]):
        g = []
        height, width = array_2.shape
        tw = (height * width) / (window * window)
        for j in range(array_1.shape[1]):
            a = neighbors(array_1, i, j)
            b = neighbors(array_2, i, j)
            c = a - b
            d = len(np.nonzero(c)[0]) * 2
            e = math.ceil(2 * (math.sqrt(len(c)) ** 2))
            f = (1 - d / e) if d != 0 else 1 - 0
            g.append(f / tw)
            h = math.fsum(g)

    print('Similarity: ', h)
