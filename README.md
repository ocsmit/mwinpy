# MWinPy
Moving window raster comparison algorithm. 

MWinPy is an implementation of the moving window comparison algorithm designed to work with geospatial data.

Moving window comparisons can be useful for comparing changes between two different data sets while taking into account spatial patterns that a pixel by pixel approach will fail to detect. 

### Usage

```
from mwinpy import MWin

# initialize the moving window
mw = MWin(1) # 1 denotes a 3 x 3 window. 

# Read rasters into arrays and store within MWin class
mw.load_rasters('./raster_1.tif", "./raster_2.tif")

# Run moving window comparing the 2 rasters with a 3 x 3 
# moving window
mw.fit()

# Total Similarity: 82.96 %

# Save tif
mw.save_tif("./out_raster.tif")

# Or assign array to variable to analyze with other packages
sim_arr = mw.out
```

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


### Example output with 2013 & 2016 NLCD data

<img src="https://render.githubusercontent.com/render/math?math=a_1 " width="25" height="25"> | Raster Map 1

<img src="https://user-images.githubusercontent.com/55674113/94374759-61e2f100-00dc-11eb-8f38-0f3019566b04.png" width="465" height="465"/>

<img src="https://render.githubusercontent.com/render/math?math=a_2" width="25" height="25"> | Raster Map 2

<img src="https://user-images.githubusercontent.com/55674113/94374896-855a6b80-00dd-11eb-82af-c6dee9ea547c.png" width="465" height="465"/>

<img src="https://render.githubusercontent.com/render/math?math=F_w" width="25" height="25"> | Similarity

<img src="https://user-images.githubusercontent.com/55674113/94375240-4d085c80-00e0-11eb-950e-7e59aa751342.png" width="484.87" height="465"/>

### Algorithm

<img src="https://user-images.githubusercontent.com/55674113/77957186-1c667800-72a1-11ea-9a5a-408f7372dd69.png"
alt="Algorithm" width="448.7" height="174.7"/>

<img src="https://user-images.githubusercontent.com/55674113/77956975-c1cd1c00-72a0-11ea-99e9-6a41bed1e1fc.png"
width="441" height="515"/>


Citations:

-  [Comparing Raster Map Comparison Algorithms for Spatial Modeling and
 Analysis](https://www.ingentaconnect.com/content/asprs/pers/2005/00000071/00000008/art00008) 
-  [Making more out of pixel-level changeinformation: using a neighbourhood
 approach toimprove land change characterization across largeand heterogeneous areas](https://www.tandfonline.com/doi/full/10.1080/10106049.2018.1458252) 
-  [A computational framework for generalized moving windows and its
 application to landscape pattern analysis](https://www.sciencedirect.com/science/article/abs/pii/S0303243415300337) 
 
