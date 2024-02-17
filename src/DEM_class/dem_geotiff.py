#pip install numpy
#pip install geotiff
#pip install Pillow
#pip install scipy
#pip install photutils
#pip install tifffile

from dem_geotiff_class import Dem_Class
from PIL import Image
import numpy as np

tiff_file = "../../repos/TopoCurve/DEM_files/Purgatory.tif"


# Initializing DEM class
dem = Dem_Class(tiff_file)

detrended, plane = dem.detrend()

# Plotting detrended DEM
dem.plot(detrended, "greyscale_dem_detrend.png")
#print(np.mean(detrended))

# Mirroring DEM on all sides
#dimx_f,dimy_f, mirror = dem.mirror_dem()
#dem.plot(mirror, "greyscale_dem_mirror.png")

#Tukey Window
#dem.plot(dem.tukeyWindow(0.5), "tukeyWind.png")

# Padding array
#dem.plot(dem.padding(0.5), "greyscale_dem_padding.png")

# FFT
 #(self, 1/200, "lowpass", 0.5)
#print(dem.FFT([200,2000], "lowpass", 0.5))
dem.plot(dem.FFT([90,100], "lowpass", 0.5), "fftdem.png")

