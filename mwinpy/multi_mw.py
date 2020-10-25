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

class load_data:

    def __init__(self, raster1, raster2):
        path1 = gdal.Open(x)
        path2 = gdal.Open(y)
        self.arrays = path1.GetRasterBand(1).ReadAsArray(), \
                      path2.GetRasterBand(1).ReadAsArray()
        self.shape = self.arrays[0].shape


class MWin:

    def __init__(self, threads, shape, window):

        self.i, self.j = shape

        self.w = window
        self.threads = threads
        self.cluster = self.i // threads
        self.rem = self.i % threads
        self.__split_itr = []
        self.vector = []
        self.__d = window // 2
        # Create inital list of starts and stops for the split
        itera = 0
        for f in range(threads):
            if f == threads - 1:
                self.__split_itr.append((self.cluster) + itera +
                                        (self.rem))
            else:
                self.__split_itr.append((self.cluster) + itera)
            itera += (self.cluster)

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

    def moving_window(self, arr1, arr2):
        w = (((self.w * 2) + 1) ** 2)
        vector = []
        for ii in range(arr1.shape[0]):
            for j in range(arr1.shape[1]):
                a = self.__neighbors(arr1, ii, j, self.w)
                b = self.__neighbors(arr2, ii, j, self.w)
                c = abs(a - b)
                d = len(np.nonzero(c)[0])
                e = abs((1 - d / w)) if d != 0 else w / w
                vector.append(e)

        return vector


    def split_moving_window(self, arr1, arr2, sl):
        w = (((self.w * 2) + 1) ** 2)
        arr1_sl, arr2_sl = arr1[sl[1][0]:sl[1][1],], arr2[sl[1][0]:sl[1][1],]
        vector = []
        if sl[0] == 0:
            for ii in range(arr1_sl.shape[0] - self.__d):
                for j in range(arr1_sl.shape[1]):
                    a = self.__neighbors(arr1_sl, ii, j, self.w)
                    b = self.__neighbors(arr2_sl, ii, j, self.w)
                    c = abs(a - b)
                    d = len(np.nonzero(c)[0])
                    e = abs((1 - d / w)) if d != 0 else w / w
                    vector.append(e)

        elif sl[0] == self.threads - 1:
            for ii in range(self.__d, arr1_sl.shape[0]):
                for j in range(arr1_sl.shape[1]):
                    a = self.__neighbors(arr1_sl, ii, j, self.w)
                    b = self.__neighbors(arr2_sl, ii, j, self.w)
                    c = abs(a - b)
                    d = len(np.nonzero(c)[0])
                    e = abs((1 - d / w)) if d != 0 else w / w
                    vector.append(e)
        else:
            for ii in range(self.__d, arr1_sl.shape[0] - (self.__d)):
                for j in range(arr1_sl.shape[1]):
                    a = self.__neighbors(arr1_sl, ii, j, self.w)
                    b = self.__neighbors(arr2_sl, ii, j, self.w)
                    c = abs(a - b)
                    d = len(np.nonzero(c)[0])
                    e = abs((1 - d / w)) if d != 0 else w / w
                    vector.append(e)
        return vector

    def fit(self, arr, arr2):
        if self.threads == 1:
            results = []
            results.append(self.moving_window(arr, arr2))
        else:
            slice = self.__split()
            slice_dict = {}
            for i in range(len(slice)):
                slice_dict.update({i: slice[i]})
            results = Parallel(n_jobs=self.threads)(delayed(self.split_moving_window)(arr, arr2, sl)
                                 for sl in slice_dict.items())
        sim_list = []
        for i in range(len(results)):
            sim_list += results[i]
        self.sim = fsum(sim_list) / (self.i * self.j)

        vector = []
        for i in range(len(results)):
            vector += results[i]
        vector_arr = np.array(vector)
        self.matrix = vector_arr.reshape([self.i, self.j])

        return results

    def plot(self, cmap="PiYG"):

        plt.imshow(self.matrix, cmap=cmap)
        plt.show()

if __name__ == '__main__':
    # arr1 = np.random.randint(2, size=(753, 200))
    # arr2 = np.random.randint(2, size=(753, 200))
    x = "/home/owen/Data/mwin/nlcd_L_2016.tif"
    y = "/home/owen/Data/mwin/nlcd_L_2013.tif"

    data = load_data(x, y)
    sh = data.shape
    
    w = [3, 13, 23, 33, 43, 53, 63]
    t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    out_dict = {}
    out_times = []
    
    t = 10
    w = 3

    # t = int(input("Threads: "))
    # w = int(input("Window: "))
    start = time.time()
    mw = MWin(t, sh, w)
    test = mw.fit(data.arrays[0], data.arrays[1])
    end = time.time() - start
    out_times.append(end)
    print(mw.sim)

