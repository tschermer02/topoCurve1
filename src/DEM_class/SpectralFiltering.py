import numpy as np
import os
from scipy import signal
from PIL import Image
import math
from photutils.psf import TukeyWindow
from scipy.fft import fft2, ifft2
from tifffile import TiffFile
from geotiff import GeoTiff

class Dem_Class():  # Define a class named Dem_Class
    """
    A class for processing digital elevation models (DEM).
    
    Attributes:
        metadata (dict): Dictionary containing metadata extracted from the GeoTIFF file.
        z_array (numpy.ndarray): Array of elevation values.
        dimx (int): Number of pixels in the x-direction.
        dimy (int): Number of pixels in the y-direction.
        dx (float): Grid spacing.
        dx_dy (float): Grid spacing in both x and y directions.
        dimx_ma (int): Number of pixels in the mirrored DEM in the x-direction.
        dimy_ma (int): Number of pixels in the mirrored DEM in the y-direction.
        dim_x (int): Dimension of the DEM grid in the x-direction.
        dim_y (int): Dimension of the DEM grid in the y-direction.
        powerOfTwo (int): Next power of two after the dimensions of the mirrored array.
        pad_x_max (int): Amount of padding to add to the maximum x-direction.
        pad_x_min (int): Amount of padding to add to the minimum x-direction.
        pad_y_max (int): Amount of padding to add to the maximum y-direction.
        pad_y_min (int): Amount of padding to add to the minimum y-direction.
        Filter (str): Type of filter applied.
        ZFilt (numpy.ndarray): Filtered elevation values.
        ZDiff (numpy.ndarray): Difference between original and filtered elevation values.
    """

    def __init__(self, tiff_file):  # Define the constructor method for the class
        """
        Initializes the Dem_Class object.

        Args:
            tiff_file (str): Path to the GeoTIFF file.
        """
        tif = TiffFile(tiff_file)  # Open the GeoTIFF file

        # Ensure input file is of the right type and contains georeferencing information
        if not tif.is_geotiff:
            raise Exception("Not a geotiff file")  # Raise an exception if it's not a GeoTIFF file
            
        if not tif.geotiff_metadata:
            raise Exception("Metadata missing")  # Raise an exception if metadata is missing
            
        # Store projection information               
        self.metadata = {
            'GeogAngularUnitsGeoKey': tif.geotiff_metadata["GeogAngularUnitsGeoKey"],
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
            'ProjLinearUnitsGeoKey': tif.geotiff_metadata["ProjLinearUnitsGeoKey"],
        }
        
        crs = tif.geotiff_metadata["ProjectedCSTypeGeoKey"].value  # Get the coordinate reference system (CRS)
        
        # Pull out array of elevation values and store it as array within the dem class
        gtiff = GeoTiff(tiff_file, crs_code=crs)  # Create a GeoTiff object
        self.z_array = gtiff.read()  # Read the elevation values
        
        # Pull out dimensions of DEM grid
        self.dimx, self.dimy = gtiff.tif_shape  # Get the dimensions of the DEM grid
        
        # Assign grid spacing and check to ensure grid spacing is uniform in x and y directions
        self.dx = tif.geotiff_metadata["ModelPixelScale"]  # Get the grid spacing
        if abs(self.dx[1] - self.dx[0]) < 1e-3:
            self.dx_dy = self.dx[0]  # If spacing is uniform, assign it to dx_dy
        else:
            raise Exception("WARNING: Grid spacing is not uniform in x and y directions!")  # Raise a warning if grid spacing is not uniform
    
    def detrend(self):  # Define a method to detrend elevation values
        """
        Detrend the elevation values.

        Returns:
            Z_detrended (numpy.ndarray): Detrended elevation values.
            plane (numpy.ndarray): Trend component of the elevation values.
        """
        self.Z_detrended = signal.detrend(self.z_array)  # Detrend the elevation values
        self.plane = self.z_array - self.Z_detrended  # Calculate the trend component
        return self.Z_detrended, self.plane  # Return the detrended values and the trend component
    
    
    def plot(self, input, filename):  # Define a method to plot and save elevation values
        """
        Plot the elevation values and save the plot as an image file.

        Args:
            input (numpy.ndarray): Input elevation values.
            filename (str): Name of the output image file.
        """
        output_dir = '../topoCurve1/src/data/SpectralFilteringOutput/'  # Define the output directory
        os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn't exist
        # Normalize input values to range [0, 255] and convert to uint8
        img_array = 255 * ((input - np.amin(input)) / (np.amax(input) - np.amin(input)))
        # Convert array to RGB image and save it
        Image.fromarray(img_array.astype(np.uint8)).convert("RGB").save(output_dir + filename)

    def mirror_dem(self):  # Define a method to mirror elevation values
        """
        Mirror the elevation values.

        Returns:
            mirrored_array (numpy.ndarray): Mirrored elevation values.
        """
        detrend, plane = self.detrend()  # Detrend the elevation values
        mirrored_array = []  # Initialize an empty array for mirrored values
        # Mirror the detrended array horizontally and vertically
        top = np.concatenate((np.rot90(detrend, 2), np.flipud(detrend), np.rot90(detrend, 2)), axis=1)
        middle = np.concatenate((np.fliplr(detrend), detrend, np.fliplr(detrend)), axis=1)
        bottom = np.concatenate((np.rot90(detrend, 2), np.flipud(detrend), np.rot90(detrend, 2)), axis=1)
        mirrored_array = np.concatenate((top, middle, bottom), axis=0)  # Concatenate the mirrored parts
        self.dimx_ma = len(mirrored_array)  # Store the dimensions of the mirrored array
        self.dimy_ma = len(mirrored_array[0])
        return mirrored_array  # Return the mirrored array
    
    def tukeyWindow(self, alphaIn):  # Define a method to apply a Tukey window
        """
        Apply a Tukey window to the elevation values.

        Args:
            alphaIn (float): Parameter controlling the shape of the Tukey window.

        Returns:
            tukey_array (numpy.ndarray): Elevation values after applying the Tukey window.
        """
        mirrored_array = self.mirror_dem()  # Mirror the elevation values
        taper = TukeyWindow(alpha=alphaIn)  # Create a TukeyWindow object
        data = taper((len(mirrored_array), len(mirrored_array[0])))  # Apply the Tukey window
        tukey_array = np.multiply(data, mirrored_array)  # Multiply the window with the elevation values

        self.dim_x = self.z_array.shape[0]  # Store the dimensions of the original array
        self.dim_y = self.z_array.shape[1]
    
        return tukey_array  # Return the elevation values after applying the Tukey window

    def padding(self, alphaIn):  # Define a method to pad elevation values
        """
        Pad the elevation values.

        Args:
            alphaIn (float): Parameter controlling the shape of the Tukey window.

        Returns:
            padded_window_array (numpy.ndarray): Padded elevation values.
        """
        # Finds next power of two
        tukey_array = self.tukeyWindow(alphaIn)  # Apply Tukey window to elevation values
        if self.dimx_ma > self.dimy_ma:  # Determine the size for padding
            a = int(math.log2(self.dimx_ma))
            if 2**a == self.dimx_ma:
                self.powerOfTwo = self.dimx_ma
            self.powerOfTwo = 2**(a+1)
        else:
            a = int(math.log2(self.dimy_ma))
            if 2**a == self.dimy_ma:
                self.powerOfTwo = self.dimy_ma
            self.powerOfTwo = 2**(a+1)

        # Finds difference in dimension of final array and power of two
        self.pad_x_max = math.ceil((self.powerOfTwo - self.dimx_ma) / 2)
        self.pad_x_min = math.floor((self.powerOfTwo - self.dimx_ma) / 2)
        self.pad_y_max = math.ceil((self.powerOfTwo - self.dimy_ma) / 2)
        self.pad_y_min = math.floor((self.powerOfTwo - self.dimy_ma) / 2)

        # Pads array
        padded_window_array = np.pad(tukey_array, ((self.pad_x_max, self.pad_x_min), (self.pad_y_max, self.pad_y_min)), 'constant', constant_values=(0, 0))
        return padded_window_array  # Return the padded elevation values

    def FFT(self, filter, filterType, alphaIn):
        """
        Apply FFT filtering to the elevation values.

        Args:
            filter (float): Filter parameter.
            filterType (str): Type of filter ('lowpass' or 'highpass').
            alphaIn (float): Parameter controlling the shape of the Tukey window.

        Returns:
            ZFilt (numpy.ndarray): Filtered elevation values.
        """
        padded_window_array = self.padding(alphaIn)  # Pad the elevation values
        # Doing fft on the windowed and padded array
        fft_array = fft2(padded_window_array)  # Compute 2-dimensional FFT of the padded array

        self.dx = np.array(self.dx)  # Convert grid spacing to a NumPy array
        self.powerOfTwo = np.array(self.powerOfTwo)  # Convert the power of two to a NumPy array
        dkx = np.divide(1, (self.dx[0] * self.powerOfTwo))  # Compute the spacing in the frequency domain in x-direction
        dky = np.divide(1, (self.dx[0] * self.powerOfTwo))  # Compute the spacing in the frequency domain in y-direction

        xc = self.powerOfTwo / 2 + 1  # Get the matrix indices of zero wavenumber in x-direction
        yc = self.powerOfTwo / 2 + 1  # Get the matrix indices of zero wavenumber in y-direction

        i, j = np.meshgrid(np.arange(self.powerOfTwo), np.arange(self.powerOfTwo))  # Create a grid of coordinates
        km = np.sqrt(np.square(dky * (i - yc)) + np.square(dkx * (j - xc)))  # Compute the distance from the zero wavenumber

        # Apply filter based on filter type
        match filterType:
            case 'lowpass':
                kfilt = np.divide(np.ones_like(filter), filter)  # Generate the filter kernel
                sigma = abs(kfilt[1] - kfilt[0]) / 3  # Compute sigma for Gaussian filter
                F = np.exp(np.multiply(-1, np.square(km - kfilt[0])) / (2 * sigma**2))  # Compute the filter
                F = F / np.max(F)  # Normalize the filter
                F[km < kfilt[0]] = 1  # Apply lowpass filter
                
            case 'highpass':
                kfilt = np.divide(np.ones_like(filter), filter)  # Generate the filter kernel
                sigma = abs(kfilt[1] - kfilt[0]) / 3  # Compute sigma for Gaussian filter
                F = np.exp(np.multiply(-1, np.square(km - kfilt[0])) / (2 * sigma**2))  # Compute the filter
                F = F / np.max(F)  # Normalize the filter
                F[km >= kfilt[1]] = 1  # Apply highpass filter

        ZMWF = np.real(ifft2(np.multiply(fft_array, F)))  # Apply inverse FFT to get filtered elevation values

        self.Filter = filter  # Store the filter parameter
        # Extract the filtered elevation values and add back the trend component
        self.ZFilt = ZMWF[(self.pad_x_max + self.dim_x): -(self.pad_x_min + self.dim_x), 
                        (self.pad_y_max + self.dim_y): -(self.pad_y_min + self.dim_y)] + self.plane
        self.ZDiff = self.z_array - self.ZFilt  # Compute the difference between original and filtered elevation values

        return self.ZFilt  # Return the filtered elevation values
