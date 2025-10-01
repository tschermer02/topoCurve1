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
        Initializes the Topocurve object.

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
        Calculates curvature attributes for a given surface.

        Args:
            ZFilt (ndarray): 2D array representing the surface.
            dx (float): Grid spacing along the x-axis.
            dy (float): Grid spacing along the y-axis.
            kt (float): Threshold for curvature values.

        Returns:
            tuple: Curvature attributes including K1, K2, KM, KG.
        """

        # Alias ZFilt as SZ for clarity
        SZ = ZFilt
        
        # Get dimensions of ZFilt
        r = len(ZFilt)
        c = len(ZFilt[0])
        
        # Generate X and Y coordinate arrays
        X = np.arange(c) * dx
        Y = np.arange(r) * dy
        
        # Create meshgrid for surface coordinates
        SX, SY = np.meshgrid(X, Y)
        
        # Calculate incremental changes along x and y directions
        du = SX[0, 1] - SX[0, 0]
        dv = SY[1, 0] - SY[0, 0]
        
        # Get shape of SZ
        m, n = SZ.shape
        
        # Initialize arrays for derivatives
        SXU = np.ones_like(SX) * dx
        SXV = np.zeros_like(SX)
        SYU = np.zeros_like(SY)
        SYV = np.ones_like(SY) * dy
        
        # Calculate first-order derivatives of surface
        SZU, SZV = np.gradient(SZ, du, dv)
        
        # Initialize arrays for surface derivatives
        SU = np.zeros((m, n, 3))
        SV = np.zeros((m, n, 3))
        
        # Assign values to surface derivative arrays
        SU[:, :, 0] = SXU
        SU[:, :, 1] = SYU
        SU[:, :, 2] = SZU
        SV[:, :, 0] = SXV
        SV[:, :, 1] = SYV
        SV[:, :, 2] = SZV
        
        # Calculate coefficients for quadratic forms
        E = np.sum(SU * SU, axis=2)
        F = np.sum(SU * SV, axis=2)
        G = np.sum(SV * SV, axis=2)
        
        # Calculate curvature tensor
        CUV = np.cross(SU, SV, axis=2)
        AC = np.sqrt(np.sum(CUV ** 2, axis=2))
        
        # Calculate normal vectors
        NX = CUV[:, :, 0] / AC
        NY = CUV[:, :, 1] / AC
        NZ = CUV[:, :, 2] / AC
        
        # Calculate second-order derivatives of normal vectors
        NXU, NXV = np.gradient(NX, du, dv)
        NYU, NYV = np.gradient(NY, du, dv)
        NZU, NZV = np.gradient(NZ, du, dv)
        
        # Initialize arrays for normal vector derivatives
        NU = np.zeros((m, n, 3))
        NV = np.zeros((m, n, 3))
        
        # Assign values to normal vector derivative arrays
        NU[:, :, 0] = NXU
        NU[:, :, 1] = NYU
        NU[:, :, 2] = NZU
        NV[:, :, 0] = NXV
        NV[:, :, 1] = NYV
        NV[:, :, 2] = NZV
        
        # Calculate coefficients for second fundamental form
        e = -np.sum(NU * SU, axis=2)
        f = -0.5 * (np.sum(NU * SV, axis=2) + np.sum(NV * SU, axis=2))
        g = -np.sum(NV * SV, axis=2)
        
        # Calculate principal curvatures
        K1 = np.zeros_like(SZ)
        K2 = np.zeros_like(SZ)
        K1U = np.zeros_like(SZ)
        K1V = np.zeros_like(SZ)
        K2U = np.zeros_like(SZ)
        K2V = np.zeros_like(SZ)
        
        a = E * G - F ** 2
        b = -(g * E - 2 * f * F + e * G)
        c = e * g - f ** 2

        K1 = (-b + np.sqrt(np.abs(b ** 2 - 4 * a * c))) / (2 * a)
        K2 = (-b - np.sqrt(np.abs(b ** 2 - 4 * a * c))) / (2 * a)
        
        # Apply threshold to curvature values
        K1[np.abs(K1) <= kt] = 0
        K2[np.abs(K2) <= kt] = 0
        
        # Calculate Gaussian and mean curvatures
        KG = K1 * K2
        KM = 0.5 * (K1 + K2)
        
        # Initialize arrays for curvature map and distribution
        SMAP = np.empty_like(KG) * np.nan
        SDist = [[None] * 9, [None] * 9]
        
        # Classify surface regions based on curvature values
        in_ = np.where((KG < 0) & (np.abs(KM) <= 0))
        SMAP[in_] = -4
        SDist[0][0] = ['Perfect Saddles']
        SDist[1][0] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG > 0) & (KM < 0))
        SMAP[in_] = -3
        SDist[0][1] = ['Peaks']
        SDist[1][1] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG < 0) & (KM < 0))
        SMAP[in_] = -2
        SDist[0][2] = ['Antiformal Saddles']
        SDist[1][2] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG == 0) & (KM < 0))
        SMAP[in_] = -1
        SDist[0][3] = ['Antiforms']
        SDist[1][3] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG == 0) & (np.abs(KM) <= 0))
        SMAP[in_] = 0
        SDist[0][4] = ['Planes']
        SDist[1][4] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG == 0) & (KM > 0))
        SMAP[in_] = 1
        SDist[0][5] = ['Synforms']
        SDist[1][5] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG < 0) & (KM > 0))
        SMAP[in_] = 2
        SDist[0][6] = ['Synformal Saddles']
        SDist[1][6] = np.size(in_) / np.size(SMAP)
        
        in_ = np.where((KG > 0) & (KM > 0))
        SMAP[in_] = 3
        SDist[0][7] = ['Basins']
        SDist[1][7] = np.size(in_) / np.size(SMAP)
        
        # Create a dictionary containing curvature attributes
        CMAP = {
            'KG': KG,
            'KM': KM,
            'K1': K1,
            'K2': K2,
            'K1U': K1U,
            'K1V': K1V,
            'K2U': K2U,
            'K2V': K2V
        }

        # Print curvature attributes
        print(K1, K2, KM, KG)

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
