################################################################################
# Title: moving_window.py
# Author: Owen Smith, University of North Georgia
# Purpose: Functions to create a moving window window comparison coeficent for
#          array and raster analysis.
################################################################################

import math
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt


def neighbors(im, i, j, d=1):
    n = im[i-d:i+d+1, j-d:j+d+1].flatten()
    return n


class MWin:

    def __init__(self, w):
        """
        Class to run a moving window comparison on two raster maps of the same size.

        Parameters
        ----------
            w : int
                Size of the moving window to be used.
                +-----+-------+
                |  1  |  3x3  |
                |  2  |  5x5  |
                |  3  |  7x7  |
                |  4  |  9x9  |
                |  5  | 11x11 |
                |  6  | 13x13 |
                |  7  | 15x15 |
                |  8  | 17x17 |
                |  9  | 19x19 |
                | 10  | 21x21 |
                +-----+-------+

        Attributes
        ----------
            w : int
                Size of the moving window to be used.

        """

        self.w = w

    def load_rasters(self, path_1, path_2):

        path1 = gdal.Open(path_1)
        path2 = gdal.Open(path_2)
        self.spatial = path1
        self.array_1 = path1.GetRasterBand(1).ReadAsArray()
        self.array_2 = path2.GetRasterBand(1).ReadAsArray()

    def load_arrays(self, array_1, array_2):

        self.array_1 = np.array(array_1)
        self.array_2 = np.array(array_2)

    def fit(self):
        vector = []
        w = (((self.w * 2) + 1) ** 2)
        height, width = self.array_2.shape
        tw = (height * width)

        for i in range(self.array_1.shape[0]):
            for j in range(self.array_1.shape[1]):
                a = neighbors(self.array_1, i, j, self.w)
                b = neighbors(self.array_2, i, j, self.w)
                c = abs(a - b)
                d = len(np.nonzero(c)[0])
                e = abs((1 - d / w)) if d != 0 else w / w
                vector.append(e)

        sim = math.fsum(vector) / tw
        self.sim_matrix = np.array(vector).reshape(self.array_1.shape)
        self.vector = np.array(vector)

        print('Total similarity: ', sim * 100, '%')

        self.similarity = sim * 100

    def save_tif(self, out_tif):
        driver = gdal.GetDriverByName('GTiff')
        metadata = driver.GetMetadata()
        shape = self.array_1.shape
        dst_ds = driver.Create(out_tif,
                               shape[1],
                               shape[0],
                               1,
                               gdal.GDT_Float32)
        proj = self.spatial.GetProjection()
        geo = self.spatial.GetGeoTransform()
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(self.out)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.FlushCache()
        dst_ds = None
