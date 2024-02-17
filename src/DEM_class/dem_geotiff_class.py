import numpy as np
from scipy import signal
from PIL import Image
import math
from photutils.psf import TukeyWindow
from scipy.fft import fft2, ifft2
from tifffile import TiffFile
from geotiff import GeoTiff

class Dem_Class():
    def __init__(self, tiff_file):
        tif=TiffFile(tiff_file)
        
        
        # Ensure input file is of the right type and contains georeferencing information
        if not tif.is_geotiff:
            raise Exception("Not a geotiff file")
            
        if not tif.geotiff_metadata:
            raise Exception("Metadata missing")
            
        
        # Store projection information               
        self.metadata={'GeogAngularUnitsGeoKey': tif.geotiff_metadata["GeogAngularUnitsGeoKey"],
                       'GeogCitationGeoKey': tif.geotiff_metadata["GeogCitationGeoKey"],
                       'GTCitationGeoKey': tif.geotiff_metadata["GTCitationGeoKey"],
                       'GTModelTypeGeoKey': tif.geotiff_metadata["GTModelTypeGeoKey"],
                       'GTRasterTypeGeoKey': tif.geotiff_metadata["GTRasterTypeGeoKey"],
                       'KeyDirectoryVersion': tif.geotiff_metadata["KeyDirectoryVersion"],
                       'KeyRevision': tif.geotiff_metadata["KeyRevision"],
                       'KeyRevisionMinor': tif.geotiff_metadata["KeyRevisionMinor"],
                       'ModelPixelScale': tif.geotiff_metadata["ModelPixelScale"],
                       'ModelTiepoint': tif.geotiff_metadata["ModelTiepoint"],
                       'ProjectedCSTypeGeoKey': tif.geotiff_metadata["ProjectedCSTypeGeoKey"],
                       'ProjLinearUnitsGeoKey': tif.geotiff_metadata["ProjLinearUnitsGeoKey"],}
        
        crs= tif.geotiff_metadata["ProjectedCSTypeGeoKey"].value
        
        # Pull out array of elevation values and store it as array within the dem class
        gtiff=GeoTiff(tiff_file, crs_code=crs)
        self.z_array = gtiff.read()
        
        # Pull out dimensions of DEM grid
        self.dimx,self.dimy=gtiff.tif_shape
        
        # Assign grid spacing and check to ensure grid spacing is uniform in x and y directions
        self.dx=tif.geotiff_metadata["ModelPixelScale"]
        
        if abs(self.dx[1]-self.dx[0]) <1e-3:
            self.dx_dy=self.dx[0]
        else:
            raise Exception("WARNING: Grid spacing is not uniform in x and y directions!")
        
    def detrend(self):
        self.Z_detrended = signal.detrend(self.z_array)
        self.plane = self.z_array-self.Z_detrended
        return self.Z_detrended, self.plane
    
    def plot(self, input, filename):
        img_array = 255*((input - np.amin(input))/(np.amax(input)- np.amin(input)))
        Image.fromarray(img_array).convert("RGB").save(filename)
    
    def mirror_dem(self):
        detrend, plane = self.detrend()
        mirrored_array = []
        top = np.concatenate((np.rot90(detrend,2),np.flipud(detrend),np.rot90(detrend,2)), axis=1)
        middle = np.concatenate((np.fliplr(detrend),detrend,np.fliplr(detrend)), axis=1)
        bottom = np.concatenate((np.rot90(detrend,2),np.flipud(detrend),np.rot90(detrend,2)), axis=1)
        mirrored_array  = np.concatenate((top, middle, bottom), axis=0)
        self.dimx_ma=len(mirrored_array)
        self.dimy_ma=len(mirrored_array[0])
        return mirrored_array
    
    def tukeyWindow(self, alphaIn):
        mirrored_array = self.mirror_dem()
        taper = TukeyWindow(alpha=alphaIn)
        data = taper((len(mirrored_array), len(mirrored_array[0])))
        tukey_array =np.multiply(data, mirrored_array)

        self.dim_x =  self.z_array.shape[0]
        self.dim_y =  self.z_array.shape[1]
    
        return tukey_array

    def padding(self, alphaIn):
        # Finds next power of two
        tukey_array = self.tukeyWindow(alphaIn)
        if self.dimx_ma>self.dimy_ma:
            a = int(math.log2(self.dimx_ma))
            if 2**a == self.dimx_ma:
                self.powerOfTwo= self.dimx_ma
            self.powerOfTwo= 2**(a+1)
        else:
            a = int(math.log2(self.dimy_ma))
            if 2**a == self.dimy_ma:
                self.powerOfTwo= self.dimy_ma
            self.powerOfTwo= 2**(a+1)

        # Finds difference in dimention of final array and power of two
        self.pad_x_max = math.ceil((self.powerOfTwo -self.dimx_ma)/2)
        self.pad_x_min = math.floor((self.powerOfTwo -self.dimx_ma)/2)
        self.pad_y_max =math.ceil((self.powerOfTwo -self.dimy_ma)/2)
        self.pad_y_min = math.floor((self.powerOfTwo -self.dimy_ma)/2)

        #pads array
        padded_window_array =np.pad(tukey_array, ((self.pad_x_max, self.pad_x_min), (self.pad_y_max, self.pad_y_min)), 'constant', constant_values=(0, 0))
        return padded_window_array

    def FFT(self, filter, filterType, alphaIn):
        padded_window_array= self.padding(alphaIn)
        #Doing fft on the windowed and padded array
        fft_array = fft2(padded_window_array)

        self.dx = np.array(self.dx)
        self.powerOfTwo = np.array(self.powerOfTwo)
        dkx = np.divide(1,(self.dx[0]*self.powerOfTwo))
        dky = np.divide(1,(self.dx[0]*self.powerOfTwo))
       
        xc = self.powerOfTwo/2+1; yc = self.powerOfTwo/2+1 #matrix indices of zero wavenumber

        i, j = np.meshgrid(np.arange(self.powerOfTwo), np.arange(self.powerOfTwo))
        km = np.sqrt(np.square(dky * (i - yc)) + np.square(dkx * (j - xc)))
        
        match filterType:
            case 'lowpass':
                kfilt=np.divide(np.ones_like(filter),filter)
                sigma=abs(kfilt[1]-kfilt[0])/3
                F = np.exp(np.multiply(-1, np.square(km - kfilt[0])) / (2 * sigma**2))
                F=F/np.max(F)
                F[km<kfilt[0]]=1

            case 'highpass':
                kfilt=np.divide(np.ones_like(filter),filter)
                sigma=abs(kfilt[1]-kfilt[0])/3
                F = np.exp(np.multiply(-1, np.square(km - kfilt[0])) / (2 * sigma**2))
                F=F/np.max(F)
                F[km>=kfilt[1]]=1

        ZMWF = np.real(ifft2(np.multiply(fft_array,F)))

        self.Filter = filter
        self.ZFilt = ZMWF[(self.pad_x_max + self.dim_x): -(self.pad_x_min + self.dim_x),(self.pad_y_max + self.dim_y): -(self.pad_y_min + self.dim_y)]+self.plane
        self.ZDiff = self.z_array-self.ZFilt

        return self.ZFilt
    
    