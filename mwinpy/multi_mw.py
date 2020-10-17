from mwinpy import MWin_multi
import numpy as np
import time
from threading import Thread
import multiprocessing as mp
from functools import partial
from osgeo import gdal
from math import fsum, sqrt
from joblib import Parallel, delayed, dump, load


class new_alg:

    def __init__(self, threads, shape, window):

        self.i, self.j = shape

        self.w = window
        self.threads = threads
        self.n = self.i // threads
        self.rem = self.i % threads
        self.__split_itr = []
        self.vector = []

        itera = 0
        for f in range(threads):
            if f == threads - 1:
                self.__split_itr.append((self.i // threads) + itera +
                                        (self.i % threads))
            else:
                self.__split_itr.append((self.i // threads) + itera)
            itera += (self.i // threads)

    def __split(self):
        iter_range = []
        for j in range(self.threads):
            if j == 0:
                iter_range.append([self.__split_itr[j] -
                                   (self.i // self.threads),
                                   self.__split_itr[j]])
            elif j == self.threads - 1:
                iter_range.append([self.__split_itr[j] -
                                   ((self.i // self.threads) +
                                    (self.i % self.threads)),
                                   self.__split_itr[j]])
            else:
                iter_range.append([self.__split_itr[j] -
                                   (self.i // self.threads),
                                   self.__split_itr[j]])
        return iter_range

    def __neighbors(self, im, i, j, d=1):
        n = im[i - d:i + d + 1, j - d:j + d + 1].flatten()
        return n

    def mw(self, arr1, arr2, sl):
        w = (((self.w * 2) + 1) ** 2)
        arr1_sl, arr2_sl = arr1[sl[0]:sl[1],], arr2[sl[0]:sl[1],]
        vector = []
        for ii in range(arr1_sl.shape[0]):
            for j in range(arr1_sl.shape[1]):
                a = self.__neighbors(arr1_sl, ii, j, self.w)
                b = self.__neighbors(arr2_sl, ii, j, self.w)
                c = abs(a - b)
                d = len(np.nonzero(c)[0])
                e = abs((1 - d / w)) if d != 0 else w / w
                vector.append(e)
        return vector

    def test(self, sl):
        return

    def fit(self, arr, arr2):
        slice = self.__split()
        results = Parallel(n_jobs=4)(delayed(self.mw)(arr, arr2, sl)
                                     for sl in slice)
        #results = Parallel(n_jobs=2)(delayed(sqrt)(i**2) for i in range(10))
        return results


if __name__ == '__main__':
    # arr1 = np.random.randint(2, size=(753, 200))
    # arr2 = np.random.randint(2, size=(753, 200))
    #x = "/home/owen/Data/LC08_L1TP_019035_20191128_20191216_01_T1/LC08_L1TP_019035_20191128_20191216_01_T1_B2.TIF"
    #y = "/home/owen/Data/LC08_L1TP_019035_20191128_20191216_01_T1/LC08_L1TP_019035_20191128_20191216_01_T1_B3.TIF"
    x = "/home/owen/Data/mwpydata/2013.tif"
    y = "/home/owen/Data/mwpydata/2016.tif"
    path1 = gdal.Open(x)
    path2 = gdal.Open(y)
    arr1 = path1.GetRasterBand(1).ReadAsArray()
    arr2 = path2.GetRasterBand(1).ReadAsArray()
    sh = arr1.shape

    mw = new_alg(4, sh, 3)
    test = mw.fit(arr1, arr2)
    #print(Parallel(n_jobs=2)(delayed(sqrt)(i**2) for i in range(10)))
    arr_test = np.array(test[0] + test[1] + test[2] + test[3])
    arr_out = arr_test.reshape(arr1.shape)

