import os
import matplotlib.pyplot as plt
import numpy as np
import time
from osgeo import gdal
from joblib import Parallel, delayed

# TODO: Documentation, thorough commenting

class MWin:
    '''
    Moving window comparison algorithm.

    Parameters
    ----------

    threads : int, default=1
        The number of cores to split the algorithm across.

    window : int, default=3
        The window size which will move across the raster map.

    '''
    def __init__(self, window_size=3, n_jobs=1):
        '''
        Initialize moving window

        Parameters
        ----------

        window_size : int, default=3
            Single integer size which indicates the i,j size of the moving
            window. E.g. window_size of 3 is equal to a window size of 3,3.
            Values must be odd numbers due to the way the window size is
            computed where the window_size // 2 is equal to the number of cells
            on each side of a cell is needed for a window.

        n_jobs : int, default=1
            Number of jobs to run in parallel across the CPU. -1 will use all
            available cores.
        '''
        self.__w = window_size
        self.threads = n_jobs
        # Initialize
        self.__split_itr = []
        self.vector = []
        # Get d value which indicates the number of cells to move from the
        # original cell for the window size w.
        self.__d = window_size // 2
        self.verbosity = 0
        self.time = []

    def __verbose(func):
        def wrapper(self, *args, **kwargs):
            if self.verbosity == 1:
                start_time = time.time()
                func(self, *args, **kwargs)
                end_time = time.time() - start_time
                self.time.append(end_time)
                print(f"{end_time} Seconds")
            else:
                func(self, *args, **kwargs)
        return wrapper

    def __check_types(func):
        def wrapper(self, x, y, nodata):
            if type(x) != type(y):
                raise TypeError("Inputs x and y are not the same type")
            func(self, x, y, nodata)
        return wrapper


    def __split(self):
        '''
        Creates list of array position values which are used to split the array
        across the CPU.

        Returns
        -------

        iter_range : list
            List of all split cell ranges with edge artifacts taken into
            account.
        '''
        # Initial list
        iter_range = []
        # Number of splits determined by the number of threads to use.
        for j in range(self.threads):
            # Take into account edge artifacts created by loss of contextual
            # data for windows due to split.
            if j == 0:
                iter_range.append([0, self.__split_itr[j] + self.__d])
            elif j == self.threads - 1:
                iter_range.append([self.__split_itr[j] - (self.cluster +
                                  self.rem + self.__d), self.__split_itr[j]])
            else:
                iter_range.append([self.__split_itr[j] - (self.cluster +
                                  self.__d), self.__split_itr[j] + self.__d])
        return iter_range

    def __neighbors(self, arr, i, j):
        '''
        Computes neighbor vector for window size w with consideration for outer
        edges.

        Parameters
        ----------

        arr : ndarray
            Input matrix to get cell i,j neighbors

        i, j : int
            Cell location i, j

        Returns
        -------

        n : ndarray
            Cell neighborhood vector
        '''

        # Use min and max to account for upper x and y neighbors
        # ranges of cells within the 'd' from the edges
        n = arr[max(i-self.__d,0):min(i+self.__d + 1,arr.shape[0]),
                max(j-self.__d,0):min(j+self.__d + 1,arr.shape[1])].flatten()
        return n

    def __mw(self, arr1, arr2, ii, j):
        '''
        Computes moving window algorithm as described by Costanza 1989.

        Parameters
        ----------

        arr1 : ndarray
            First array for comparison
        arr2 : ndarry
            Second array for comparison
        ii : int
            Cell location i
        j : int
            Cell location j

        Returns
        -------

        similarity : float
            Similarity value for cell i, j with neighborhood of window size
            'w'. Value will range between 0.0 - 1.0.
        '''

        # number of total cells across two neighborhoods with window size w.

        # Assign nodata value of -1 which will be ignored in final coefficient
        if arr1[ii][j] and arr2[ii][j] == self.nodata:
            return -1
        else:
            # Get neighborhood vectors of each cell i, j in arrays 1 and 2.
            a = self.__neighbors(arr1, ii, j)
            b = self.__neighbors(arr2, ii, j)
            # Find number of cells which are different values. A value of
            # 0 indicates it is the same.
            c = abs(a - b)
            d = len(np.nonzero(c)[0])
            # Divide number of cells which are different by the total number of
            # cells in the two neighborhoods. If it is 100% similar assign cell
            # the value of 1.
            w = (((len(c) * 2) + 1) ** 2)
            similarity = (1 - d / w) if d != 0 else 1
            return similarity

    def moving_window(self, arr1, arr2):
        '''
        Moving window implementation where the number of cores utilized is one.

        Parameters
        ----------

        arr1 : ndarray
            First array for comparison
        arr2 : ndarry
            Second array for comparison

        Returns
        -------

        vector : list
            List vector which contains similarity values of each cell

        '''
        vector = []
        for ii in range(self.arr1.shape[0]):
            for j in range(self.arr1.shape[1]):
                vector.append(self.__mw(arr1, arr2, ii, j))
        return vector

    def split_moving_window(self, arr1, arr2, sl):
        '''
        Moving window implementation where the number of cores utilized is
        more than one.

        Parameters
        ----------

        arr1 : ndarray
            First array for comparison
        arr2 : ndarry
            Second array for comparison
        sl : list
            Section of arr1 and arr2 which has assigned to CPU core.

        Returns
        -------

        vector : list
            List vector which contains similarity values of each cell

        '''
        arr1_sl, arr2_sl = arr1[sl[1][0] : sl[1][1], ], \
                           arr2[sl[1][0] : sl[1][1], ]
        vector = []
        if sl[0] == 0:
            for ii in range(0, arr1_sl.shape[0] - self.__d):
                for j in range(arr1_sl.shape[1]):
                    vector.append(self.__mw(arr1_sl, arr2_sl, ii, j))
        elif sl[0] == self.threads - 1:
            for ii in range(self.__d, arr1_sl.shape[0]):
                for j in range(arr1_sl.shape[1]):
                    vector.append(self.__mw(arr1_sl, arr2_sl, ii, j))
        else:
            for ii in range(self.__d, arr1_sl.shape[0] - self.__d):
                for j in range(arr1_sl.shape[1]):
                    vector.append(self.__mw(arr1_sl, arr2_sl, ii, j))
        return vector

    @__check_types
    def __load_data(self, x, y, nodata):

        if type(x) and type(y) == str:
            path1 = gdal.Open(x)
            path2 = gdal.Open(y)
            self.arr1 = path1.GetRasterBand(1).ReadAsArray()
            self.arr2 = path2.GetRasterBand(1).ReadAsArray()
            self.nodata = path1.GetRasterBand(1).GetNoDataValue()
        elif type(x) and type(y) == np.ndarray:
            self.arr1 = x
            self.arr2 = y
            self.nodata = nodata

    @__verbose
    def fit(self, x, y, nodata=None, ):
        '''
        Computes moving window algorithm.

        Parameters
        ----------

        x : ndarray, str
            Must be ndarray or a valid raster file path.
        y : ndarray, str
            Must be ndarray or a valid raster file path.
        nodata : int, default=None
            Only needs to be set if x and y are ndarrays, otherwise the nodata
            value is read from the raster datasets.
        '''
        self.__load_data(x, y, nodata)

        if self.arr1.shape != self.arr2.shape:
            raise Exception("Shape of x and y is not the same")

        self.i, self.j = self.arr1.shape
        # Create initial list of starts and stops for the split

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
            results.append(self.moving_window(self.arr1, self.arr2))
        else:
            slice = self.__split()
            slice_dict = {}
            for i in range(len(slice)):
                slice_dict.update({i: slice[i]})
            results = Parallel(n_jobs=self.threads)(delayed(
                      self.split_moving_window)(self.arr1, self.arr2, sl)
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

    def fit_multi():


        return

    def plot(self, cmap="Greys"):

        plt.imshow(self.matrix, cmap=cmap)
        plt.show()


    def save_tif(self, snap, out_tif):
        if os.path.exists(out_tif):
            os.remove(out_tif)

        snap = gdal.Open(snap)

        driver = gdal.GetDriverByName('GTiff')
        metadata = driver.GetMetadata()
        shape = self.matrix.shape
        dst_ds = driver.Create(out_tif,
                               xsize=shape[1],
                               ysize=shape[0],
                               bands=1,
                               eType=gdal.GDT_Float32)
        proj = snap.GetProjection()
        geo = snap.GetGeoTransform()
        dst_ds.SetGeoTransform(geo)
        dst_ds.SetProjection(proj)
        dst_ds.GetRasterBand(1).SetNoDataValue(-1)
        dst_ds.GetRasterBand(1).WriteArray(self.matrix)
        dst_ds.FlushCache()
        dst_ds = None


if __name__ == '__main__':
    # arr1 = np.random.randint(2, size=(753, 200))
    # arr2 = np.random.randint(2, size=(753, 200))
    x = "/home/owen/Data/mwin/nan_2016.tif"
    y = "/home/owen/Data/mwin/nan_2013.tif"

    w = [3, 13, 23, 33, 43, 53, 63]
    t = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    out_dict = {}
    out_times = []

    t = 4
    w = 3

    # t = int(input("Threads: "))
    # w = int(input("Window: "))
    start = time.time()
    mw = MWin(w, t)
    test = mw.fit(x, y)
    end = time.time() - start
    print(mw.sim)
    #mw.plot(cmap="magma")
    mw.save_tif(x, "/home/owen/tmp/OUT.tif")
