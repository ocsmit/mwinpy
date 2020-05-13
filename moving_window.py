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


def moving_window(array_1, array_2):
    """
    Moving window function to take two arrays and apply a 3*3 window
    comparison for each cell and compute both cell similarity percentage and
    similarity percentage for the whole dataset.

    :param array_1: NumPy array 1
    :param array_2: NumPy array 2
    :return: similarity array
    """
    g = []
    w = (3**2)
    height, width = array_2.shape
    tw = (height * width)
    for i in range(array_1.shape[0]):
        for j in range(array_1.shape[1]):
            a = neighbors(array_1, i, j)
            b = neighbors(array_2, i, j)
            c = abs(a - b)
            d = len(np.nonzero(c)[0])
            e = abs((1 - d / w)) if d != 0 else w / w
            print(e)
            g.append(e)
    sim = math.fsum(g) / tw
    sim_arr = np.array(g)
    sim_arr1 = sim_arr.reshape(array_1.shape)
    print('Total similarity: ', sim * 100, '%')
    return sim_arr1
