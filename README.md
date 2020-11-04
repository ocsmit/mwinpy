# MWinPy
[![Documentation Status](https://readthedocs.org/projects/mwinpy/badge/?version=latest)](https://mwinpy.readthedocs.io/en/latest/?badge=latest)


MWinPy is an implementation of the moving window comparison algorithm designed to work with geospatial data.

Moving window comparisons can be useful for comparing changes between two different data sets while taking into account spatial patterns that a pixel by pixel approach will fail to detect. 

### Usage

```
from mwinpy import MWin

# initialize a 5x5 moving window with 3 cores
mw = MWin(5, 3) 

# Categorical rasters to compare
x, y = "./NLCD_2013.tif", "./NLCD_2016.tif"

# Run 5x5 comparison window over x & y. NoData values read automatically.
mw.fit(x, y)

mw.plot()
```
