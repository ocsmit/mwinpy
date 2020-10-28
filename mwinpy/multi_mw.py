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


class MWin:

    def __init__(self, threads, window):

        self.w = window
        self.threads = threads
        self.__split_itr = []
        self.vector = []
        self.__d = window // 2

    def __split(self):
        # Create split ranges based of intialist of split numbers
        iter_range = []
        for j in range(self.threads):
            if j == 0:
                iter_range.append([self.__split_itr[j] - (self.cluster),
                                   self.__split_itr[j] + self.__d])
            elif j == self.threads - 1:
                iter_range.append([self.__split_itr[j] - (self.cluster +
                                  self.rem + self.__d), self.__split_itr[j]])
            else:
                iter_range.append([self.__split_itr[j] - (self.cluster +
                                  self.__d), self.__split_itr[j] + self.__d])
        return iter_range

    def __neighbors(self, im, i, j, window):
        d = window // 2
        n = im[i - d:i + d + 1, j - d:j + d + 1].flatten()
        return n

    def __mw(self, arr1, arr2, ii, j, w):
        if arr1[ii][j] and arr2[ii][j] == self.nodata:
            return -1
        else:
            a = self.__neighbors(arr1, ii, j, self.w)
            b = self.__neighbors(arr2, ii, j, self.w)
            c = abs(a - b)
            d = len(np.nonzero(c)[0])
            e = abs((1 - d / w)) if d != 0 else w / w
            return e

    def moving_window(self, arr1, arr2):
        w = (((self.w * 2) + 1) ** 2)
        vector = []
        for ii in range(arr1.shape[0]):
            for j in range(arr1.shape[1]):
                vector.append(self.__mw(arr1, arr2, ii, j, w))

        return vector

    def split_moving_window(self, arr1, arr2, sl):
        w = (((self.w * 2) + 1) ** 2)
        arr1_sl, arr2_sl = arr1[sl[1][0]:sl[1][1],], arr2[sl[1][0]:sl[1][1],]
        vector = []
        if sl[0] == 0:
            for ii in range(arr1_sl.shape[0] - self.__d):
                for j in range(arr1_sl.shape[1]):
                    vector.append(self.__mw(arr1_sl, arr2_sl, ii, j, w))

        elif sl[0] == self.threads - 1:
            for ii in range(self.__d, arr1_sl.shape[0]):
                for j in range(arr1_sl.shape[1]):
                    vector.append(self.__mw(arr1_sl, arr2_sl, ii, j, w))
        else:
            for ii in range(self.__d, arr1_sl.shape[0] - (self.__d)):
                for j in range(arr1_sl.shape[1]):
                    vector.append(self.__mw(arr1_sl, arr2_sl, ii, j, w))
        return vector

    def fit(self, x, y, nodata=None):

        if type(x) and type(y) == str:
            path1 = gdal.Open(x)
            path2 = gdal.Open(y)
            arr1 = path1.GetRasterBand(1).ReadAsArray()
            arr2 = path2.GetRasterBand(1).ReadAsArray()

            self.nodata = path1.GetRasterBand(1).GetNoDataValue()
        elif type(x) and type(y) == np.ndarray:
            arr1 = x
            arr2 = y
            self.nodata = nodata
        else:
            raise TypeError("Inputs x and y must be valid file paths or numpy arrays.")

        self.i, self.j = arr1.shape
        # Create inital list of starts and stops for the split

        self.cluster = self.i // self.threads
        self.rem = self.i % self.threads

        itera = 0
        for f in range(self.threads):
            if f == self.threads - 1:
                self.__split_itr.append((self.cluster) + itera +
                                        (self.rem))
            else:
                self.__split_itr.append((self.cluster) + itera)
            itera += (self.cluster)

        if self.threads == 1:
            results = []
            results.append(self.moving_window(arr1, arr2))
        else:
            slice = self.__split()
            slice_dict = {}
            for i in range(len(slice)):
                slice_dict.update({i: slice[i]})
            results = Parallel(n_jobs=self.threads)(delayed(
                      self.split_moving_window)(arr1, arr2, sl)
                      for sl in slice_dict.items())

        vector = []
        for i in range(len(results)):
            vector += results[i]
        vector_arr = np.array(vector)

        tw = sum(1 for n in vector if n != -1)
        self.sim = np.sum(vector_arr) / (tw)

        mat = vector_arr.reshape([self.i, self.j])
        self.matrix = np.ma.masked_where(mat == -1,
                                    mat,
                                    copy=True)


        return results

    def plot(self, cmap="Greys"):

        plt.imshow(self.matrix, cmap=cmap)
        plt.show()

if __name__ == '__main__':
    # arr1 = np.random.randint(2, size=(753, 200))
    # arr2 = np.random.randint(2, size=(753, 200))
    x = "/home/owen/Data/mwin/nan_2016.tif"
    y = "/home/owen/Data/mwin/nan_2013.tif"
    
    w = [3, 13, 23, 33, 43, 53, 63]
    t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    out_dict = {}
    out_times = []
    
    t = 3
    w = 3  

    # t = int(input("Threads: "))
    # w = int(input("Window: "))
    start = time.time()
    mw = MWin(t, w)
    test = mw.fit(x, y)
    end = time.time() - start
    print(mw.sim)
    mw.plot(cmap="magma")

