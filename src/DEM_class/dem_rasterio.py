# pip install rasterio
# pip install numpy
# pip install matplotlib
# pip install scipy
# pip install photutils
# pip install tifffile

import numpy as np
import rasterio as rio
from rasterio.plot import show
import matplotlib.pyplot as plt

from dem_ras_class import Dem_Ras_Class

tiff_file= "../../repos/TopoCurve/DEM_files/Purgatory.tif"

dem_test = Dem_Ras_Class(tiff_file)


dem_test.plot_func(dem_test.fftf_2d([2000,200], "lowpass"))



