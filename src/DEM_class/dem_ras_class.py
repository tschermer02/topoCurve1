import rasterio as rio
from scipy import signal
import numpy as np
from rasterio.plot import show
import matplotlib.pyplot as plt
import math
from scipy.fft import fft2, fftshift, ifft2
from photutils.psf import TukeyWindow
from tifffile import TiffFile

class Dem_Ras_Class():
    def __init__(self, dem_path):
        
        tif=TiffFile(dem_path)
        
        
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
        
        
        # Pull out array of elevation values and store it as array within the dem class
        dem = rio.open(dem_path)
        self.z_array = dem.read(1).astype('float64')
        
        # Pull out dimensions of DEM grid
        self.dimx=dem.width
        self.dimy=dem.height
        
        # Assign grid spacing and check to ensure grid spacing is uniform in x and y directions
        self.dx = tif.geotiff_metadata["ModelPixelScale"]
        
        if abs(self.dx[1]-self.dx[0]) <1e-3: # This checks if dx and dy are basically equal, T=contn. F=popup error.
            self.dx_dy=self.dx[0]
        else:
            raise Exception("WARNING: Grid spacing is not uniform in x and y directions!")
      

    def detrend(self): #Function that detrends the array and outputs the detrended array and the plane used. 
        
        self.detrended = signal.detrend(self.z_array)
        self.detrended_plane = self.z_array-self.detrended

        return self.detrended, self.detrended_plane
    
    def plot_func(self, input): #Function to plot array

        fig, ax = plt.subplots(1, figsize = (12,12))
        show(input, cmap='Greys_r', ax=ax)
        plt.axis("off")
        plt.show()

    def mirror_array(self): #function that mirrors the array along the sides and rotates at corners.

        detrended, plane = self.detrend() #runs detrened function to get detrended array.

        #Concatenates(merge) top, middle, and bottom arrays horizontally.
        top = np.concatenate((np.rot90(detrended,2), np.flipud(detrended), np.rot90(detrended,2)), axis = 1)
        mid = np.concatenate((np.fliplr(detrended), detrended, np.fliplr(detrended)), axis =1)
        bot = np.concatenate((np.rot90(detrended,2), np.flipud(detrended), np.rot90(detrended,2)), axis =1)

        #Concatenates top, middle, and bottom groups together vertically.
        Zm = np.concatenate((top, mid, bot), axis=0)

        return Zm


    def tukey_window(self): #Function that creates a tukey window and applies it to the mirrored array.

        input = self.mirror_array() #runs mirror function to get imput.

        #Dimensions of mirrored array.
        length = len(input) 
        width = len(input[0])

        taper = TukeyWindow(alpha=0.5) #tukey window with alfa of 0.5.
        data = taper((length, width)) #Creating the tukey window.
        output = data * input #this is the tukey_window * mirrored_array.
        
        # this commented code will cut out the origanal size after the tukey window is applied..
        self.dim_x =  self.z_array.shape[0]
        self.dim_y =  self.z_array.shape[1]
        #output[(self.pad_x_max + self.dim_x): -(self.pad_x_max + self.dim_x),(self.pad_y_max + self.dim_y): -(self.pad_y_max + self.dim_y)]

        return output


    def padding_array(self): #Function that takes the tukey window and pads it to the next power of 2 with 0.

        input = self.tukey_window() #Runs tukey window function for the input.

        #dimentions of the mirrored tukey windowed array.
        x_dim =  len(input)
        y_dim =  len(input[0])
        
        if(x_dim > y_dim): #this uses the dimensions find the next power of 2 for x_dim.
            N = x_dim
            a = int(math.log2(N))
            if 2**a == N:
                self.power_of2 = N
            self.power_of2 = 2**(a + 1)
        
        if(x_dim < y_dim): #this uses the dimensions find the next power of 2 for y_dim.
            N = y_dim
            a = int(math.log2(N))
            if 2**a == N:
                self.power_of2 = N
            self.power_of2 = 2**(a + 1)

        #This is to trim the padding down to fit around the mirrored tukey windowed array, so the x,y dimensions are the power of 2.
        self.pad_x_max = math.ceil((self.power_of2 - x_dim)/2)
        self.pad_x_min = math.floor((self.power_of2 - x_dim)/2)
        self.pad_y_max = math.ceil((self.power_of2 - y_dim)/2)
        self.pad_y_min = math.floor((self.power_of2 - y_dim)/2)

        #This just brings it all together.
        self.array = np.pad(input,((self.pad_x_max, self.pad_x_min), (self.pad_y_max, self.pad_y_min)), 'constant', constant_values= (0,0))
        
        return self.array
        
    
    def fftf_2d(self, filter, filter_type):
        input = self.padding_array()
        self.filter = filter

        # Doing the 2d fourier transformation on the input.
        # https://docs.scipy.org/doc/scipy/tutorial/fft.html#and-n-d-discrete-fourier-transforms
        
        input_fft = fft2(input)
        dkx = np.divide(1,(self.dx * self.power_of2)); dky = np.divide(1,(self.dx * self.power_of2)) # Defining wave number increments.
        
        # Making the radical wave number matrix
        xc = self.power_of2/2+1; yc = self.power_of2/2+1 # 
        [cols, rows] = np.meshgrid(self.power_of2 , self.power_of2) # matrices of column and row indices 
        km = np.square(np.sqrt(dky*(rows - yc))) + np.square(np.sqrt(dkx*(rows - xc))) # km = matrix of radial wavenumbers

        
        match filter_type:

            case 'lowpass':
                kfilt=np.divide(1, filter)
                sigma=np.abs(kfilt[1]-kfilt[0])/3
                F = np.exp(-np.square(km - kfilt[0]) / (2 * sigma**2))
                F[km < kfilt[0]]=1

            case 'highpass':
                kfilt=np.divide(1, filter)
                sigma=np.abs(kfilt[1]-kfilt[0])/3
                F = np.exp(-np.square(km - kfilt[0]) / (2 * sigma**2))
                F[km >= kfilt[1]]=1

        input_fft_real = np.real(ifft2(np.multiply(input_fft,F)))

        '''
        What values to use for filter?
            190M to Hz

        Is the output of fftf_2d is something that goes into a plot function?
            Print out

        # What are this/ what are they doing?

        obj.DEM.ZFilt=ZMWF(r+1:ny-r,c+1:nx-c)+plane
        obj.DEM.ZDiff=obj.DEM.Z-obj.DEM.ZFilt
        obj.Filter=Filter
        '''
        self.z_filt_array = input_fft_real[(self.pad_x_max + self.dim_x): -(self.pad_x_max + self.dim_x),(self.pad_y_max + self.dim_y): -(self.pad_y_max + self.dim_y)]
        self.z_filt_array = self.z_filt_array + self.detrended_plane
        self.z_diff = self.z_array-self.z_filt_array


        return self.z_filt_array





