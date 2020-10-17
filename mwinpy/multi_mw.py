from mwinpy import MWin_multi
import matplotlib.pyplot as plt
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
        print(arr1_sl.shape)
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
        print(slice)
        results = Parallel(n_jobs=self.threads)(delayed(self.mw)(arr, arr2, sl)
                                     for sl in slice)
        #results = Parallel(n_jobs=2)(delayed(sqrt)(i**2) for i in range(10))
        sim_list = []
        for i in range(len(results)):
            sim_list += results[i]
        self.sim = fsum(sim_list) / (self.i * self.j)
        return results


if __name__ == '__main__':
    # arr1 = np.random.randint(2, size=(753, 200))
    # arr2 = np.random.randint(2, size=(753, 200))
    x = "/home/owen/Data/mwin/NLCD_2016.tif"
    y = "/home/owen/Data/mwin/NLCD_2013.tif"
    path1 = gdal.Open(x)
    path2 = gdal.Open(y)
    arr1 = path1.GetRasterBand(1).ReadAsArray()
    arr2 = path2.GetRasterBand(1).ReadAsArray()
    sh = arr1.shape
    start = time.time()
    mw = new_alg(10, sh, 3)
    test = mw.fit(arr1, arr2)
    print(mw.sim)
    end = time.time() - start
    print("Seconds: {}".format(end))
    #print(Parallel(n_jobs=2)(delayed(sqrt)(i**2) for i in range(10)))
    test_list = []
    for i in range(len(test)):
        test_list += test[i]
    arr_test = np.array(test_list)
    arr_out = arr_test.reshape(arr1.shape)
    plt.imshow(arr_out)
    plt.show()

