from mwinpy import MWin_multi
import numpy as np
import time
import threading
import multiprocessing as mp
from functools import partial
from osgeo import gdal
from math import fsum

def neighbors(im, i, j, d=1):
    n = im[i-d:i+d+1, j-d:j+d+1].flatten()
    return n

def test(i):
    print(i)
    
def mw(i):
    time.sleep(1)
    path_1 = '/home/owen/Data/Landsat_8/aoi/LC08_L1TP_019036_20190925_20191017_01_T1_B2.TIF'
    path_2 = '/home/owen/Data/Landsat_8/aoi/LC08_L1TP_019036_20190925_20191017_01_T1_B3.TIF'
    path1 = gdal.Open(path_1)
    path2 = gdal.Open(path_2)
    array_1 = path1.GetRasterBand(1).ReadAsArray()
    shape = array_1.shape
    vector = []
    if i - 1 < 0:
        v11 = path1.GetRasterBand(1).ReadAsArray()[i, :]
        v21 = path2.GetRasterBand(1).ReadAsArray()[i, :]

    else:
        v11 = path1.GetRasterBand(1).ReadAsArray()[i - 1, :]
        v21 = path2.GetRasterBand(1).ReadAsArray()[i - 1, :]
    v12 = path1.GetRasterBand(1).ReadAsArray()[i, :]
    v22 = path2.GetRasterBand(1).ReadAsArray()[i, :]
    if i == shape[0] - 1:
        v13 = path1.GetRasterBand(1).ReadAsArray()[i, :]
        v23 = path2.GetRasterBand(1).ReadAsArray()[i, :]
    else:
        v13 = path1.GetRasterBand(1).ReadAsArray()[i + 1, :]
        v23 = path2.GetRasterBand(1).ReadAsArray()[i + 1, :]
    arr1 = np.stack((v11, v12, v13))
    arr2 = np.stack((v21, v22, v23))
    w = 1
    w = (((w * 2) + 1) ** 2)
    for j in range(arr1.shape[1]):
        a = neighbors(arr1, 1, j, w)
        b = neighbors(arr2, 1, j, w)
        c = abs(a - b)
        d = len(np.nonzero(c)[0])
        e = abs((1 - d / w)) if d != 0 else w / w
        vector.append(e)
    print(i, fsum(vector))
    return vector


if __name__ == "__main__":
    x ='/home/owen/Data/Landsat_8/aoi/LC08_L1TP_019036_20190925_20191017_01_T1_B2.TIF'
    y ='/home/owen/Data/Landsat_8/aoi/LC08_L1TP_019036_20190925_20191017_01_T1_B3.TIF'
    vector = []
    w = MWin_multi(1, x, y)
    threads = mp.cpu_count()

    arg = [i for i in range(w.shape[1])]
    it = iter(arg)

    def gen():
        for foo in range(w.shape[0]):
            yield foo
    start = time.time()
    with mp.Pool() as pool:
        L = pool.map(mw, gen())
        print(len(L))
    end = time.time() - start
    print(end)
    #jobs = []
    #for i in range(w.shape[1]):
    #    thread: = mp.Process(target=w.mw, args=(i, ))
    #    jobs.append(thread)
#
    #for j in jobs:
    #    j.start()
#
    #for j in jobs():
    #    j.join()
#
    print('complete')
