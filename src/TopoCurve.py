# topo_curve/topo_curve/topocurve.py
import os
import math
import numpy as np
from scipy.fft import fft2, ifft2
from photutils.psf import TukeyWindow
from PIL import Image
from tifffile import TiffFile
from geotiff import GeoTiff

class TopoCurve ():
    """
    A class for processing digital elevation models (DEM).

    Attributes:
        metadata (dict): Dictionary containing metadata extracted from the GeoTIFF file.
        z_array (numpy.ndarray): Array of elevation values.
        dimx (int): Number of pixels in the x-direction.
        dimy (int): Number of pixels in the y-direction.
        dx (float): Grid spacing.
        dx_dy (float): Grid spacing in both x and y directions.
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
    
    def CurveCalc(self, ZFilt, dx, dy, kt):
        """
        Calculate principal curvatures and curvature features.

        Args:
            ZFilt (numpy.ndarray): Filtered surface data.
            dx (float): Grid spacing in the x direction.
            dy (float): Grid spacing in the y direction.
            kt (float): Threshold value.

        Returns: 
            K1, K2 (tuples): Principal curvatures.
            KM (tuple): Mean curvature.
            KG (tuple): Gaussian curvature.
        
        """

        r = len(ZFilt)  # Number of rows in ZFilt (surface data)
        c = len(ZFilt[0])  # Number of columns in ZFilt
        
        X = np.arange(c) * dx # Generate 1D array of x coordinates
        Y = np.arange(r) * dy  # Generate 1D array of y coordinates
        
        SX, SY = np.meshgrid(X, Y)  # Create 2D arrays of x and y coordinates
        
        du = SX[0, 1] - SX[0, 0]  # Grid spacing in the u direction
        dv = SY[1, 0] - SY[0, 0]  # Grid spacing in the v direction
        
        SZU, SZV = np.gradient(ZFilt, du, dv)  # Compute gradients of ZFilt
        
        SXU = np.ones_like(SX) * dx  # Derivative of x component of surface vector in the u direction
        SXV = np.zeros_like(SX)  # Derivative of x component of surface vector in the v direction
        SYU = np.zeros_like(SY)  # Derivative of y component of surface vector in the u direction
        SYV = np.ones_like(SY) * dy  # Derivative of y component of surface vector in the v direction
        
        SU = np.zeros((r, c, 3))  # Initialize array to store surface vector derivatives
        SV = np.zeros((r, c, 3))  # Initialize array to store surface vector derivatives
        SU[:, :, 0] = SXU  # Store x component of surface vector derivative in SU
        SU[:, :, 1] = SYU  # Store y component of surface vector derivative in SU
        SU[:, :, 2] = SZU  # Store z component of surface vector derivative in SU
        SV[:, :, 0] = SXV  # Store x component of surface vector derivative in SV
        SV[:, :, 1] = SYV  # Store y component of surface vector derivative in SV
        SV[:, :, 2] = SZV  # Store z component of surface vector derivative in SV
        
        E = np.sum(SU * SU, axis=2)  # Compute E tensor component
        F = np.sum(SU * SV, axis=2)  # Compute F tensor component
        G = np.sum(SV * SV, axis=2)  # Compute G tensor component
        
        al = F * G - F ** 2  # Compute intermediate value
        bl = E * G - G * E  # Compute intermediate value
        cl = E * F - F * E  # Compute intermediate value
        
        K1 = (-bl - np.sqrt(np.abs(bl ** 2 - 4 * al * cl))) / (2 * al)  # Compute principal curvature K1
        K2 = (-bl + np.sqrt(np.abs(bl ** 2 - 4 * al * cl))) / (2 * al)  # Compute principal curvature K2
        
        K1[np.abs(K1) <= kt] = 0  # Apply thresholding to K1
        K2[np.abs(K2) <= kt] = 0  # Apply thresholding to K2
        
        KG = K1 * K2  # Compute Gaussian curvature
        KM = 0.5 * (K1 + K2)  # Compute mean curvature
        
        SMAP = np.empty_like(KG) * np.nan  # Initialize surface map
        SDist = [[None] * 9, [None] * 9]  # Initialize distribution of curvature features
        
        # Classify curvature features based on thresholds and store distribution information
        
        CMAP = {  # Create a dictionary to store curvature values
            'KG': KG,  # Gaussian curvature
            'KM': KM,  # Mean curvature
            'K1': K1,  # Principal curvature K1
            'K2': K2   # Principal curvature K2
        }

        # Return principal curvatures and surface map
        return K1, K2, KM, KG

    def plot(self, input, filename):  # Define a method to plot and save elevation values
        """
        Plot the elevation values and save the plot as an image file.

        Args:
            input (numpy.ndarray): Input elevation values.
            filename (str): Name of the output image file.
        """
        output_dir = '../topoCurve1/reports/figures/'  # Define the output directory
        os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn't exist
        # Normalize input values to range [0, 255] and convert to uint8
        img_array = 255 * ((input - np.amin(input)) / (np.amax(input) - np.amin(input)))
        # Convert array to RGB image and save it
        Image.fromarray(img_array.astype(np.uint8)).convert("RGB").save(output_dir + filename)
