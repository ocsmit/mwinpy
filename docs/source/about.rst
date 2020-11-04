About
=====

MWinPy is a parallelized implementation of the moving window comparison
algorithm ([1]_) designed to ingrated easily with with geocomputational
workflows.

Moving window comuarisons can be useful for comparing changes between two different data sets while taking into account spatial patterns that a pixel by pixel approach will fail to detect.

The index which is implemented and takes into account NoData values is as follows:

.. math::  F_w = \frac{\sum_{s=1}^{t_w}\left[1 - \frac{\sum_{i=1}^{p}|a_{1,i} - a_{2,i}|}{2n_s}\right]}{t_w}

where :math:`F_w` is the correlation between 0 and 1 where 1 indicates complete similarity and 0 indicates none, :math:`w` is the window size, :math:`s` is the index for moving windows, :math:`t_w` is the number of windows with the window size :math:`w`, :math:`a_{1,i}` and :math:`a_{2,i}` represent the numbers of cells with category :math:`i` in rasters 1 and 2, respectively, :math:`n_s` is the number of data cells in the window, and :math:`p` is the number of categories.

The moving window algorithm is further quantified with the weighted multiple resolution index:

.. math:: F_t = \frac{\sum_{w=1}^{n}F_we^{-k(w - 1)}}{\sum_{w=1}^{n}e^{-k(w - 1)}}

where :math:`n` is the total number of window resolutions to iterate, and :math:`k` is the constant weight where :math:`k = 0`  gives all windows the same weight and :math:`k = 1` gives the first several resolutions more weight.

Citation
----------

.. [1] Costanza, R. (1989). Model goodness of fit: a multiple resolution
       procedure. Ecological modelling, 47(3-4), 199-215.
