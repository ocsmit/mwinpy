################################################################################
# Title: moving_window.py
# Author: Owen Smith, University of North Georgia
# Purpose: Functions to create a moving window window comarison coeficent for
#          array and raster analysis.
################################################################################

import math
import numpy as np
from osgeo import gdal

def neighbors(im, i, j, d=1):
    n = im[i-d:i+d+1, j-d:j+d+1].flatten()
    return n


def moving_window(array_1, array_2, window=1):
    """
    Moving window function to take two arrays and apply a 3*3 window
    comparison for each cell and compute both cell similarity percentage and
    similarity percentage for the whole dataset.
    :param array_1: NumPy array 1
    :param array_2: NumPy array 2
    :param window=1: Window that is applied to the arrays. The default is 1
    which applies a 3x3 window.

            window values:
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

    :return: similarity array
    """
    g = []
    w = (((window * 2) + 1)**2)
    height, width = array_2.shape
    tw = (height * width)
    for i in range(array_1.shape[0]):
        for j in range(array_1.shape[1]):
            a = neighbors(array_1, i, j, window)
            b = neighbors(array_2, i, j, window)
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


def raster_moving_window(raster_1, raster_2, window=1, out_tiff=None):
    """
    Moving window function to for rasters which applies a 3*3 window
    comparison for each cell and compute both cell similarity percentage and
    similarity percentage for the whole dataset. IF out_tif is specified then it
    saves the out put array as a GeoTIFF.
    :param raster_1: path for raster 1
    :param raster_2: path for raster 2
    :param out_tiff: path to save the output similarity GeoTIFF
    :param window=1: Window that is applied to the rasters. The default is 1
    which applies a 3x3 window.

            window values:
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


    :return: similarity array
    """
    path1 = gdal.Open(raster_1)
    path2 = gdal.Open(raster_2)
    arr_1 = path1.GetRasterBand(1).ReadAsArray().astype(np.uint8)
    arr_2 = path2.GetRasterBand(1).ReadAsArray().astype(np.uint8)

    sim = moving_window(arr_1, arr_2, window)
    if out_tiff:
        driver = gdal.GetDriverByName('GTiff')
        metadata = driver.GetMetadata()
        shape = arr_1.shape
        dst_ds = driver.Create(out_tiff,
                               shape[1],
                               shape[0],
                               1,
                               gdal.GDT_Float32)
        proj = path1.GetProjection()
        geo = path2.GetGeoTransform()
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).WriteArray(sim)
        dst_ds.GetRasterBand(1).SetNoDataValue(3)
        dst_ds.FlushCache()
        dst_ds = None
        del(sim)
    if not out_tiff:
        del(sim)
