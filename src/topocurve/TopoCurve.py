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
        Calculates principal, mean, and Gaussian curvature maps for a surface Z(x,y).
        Compatible with MATLAB parity and right-handed normal orientation.
        """
        # --- 1. Setup ---
        Z = ZFilt

        # --- 2. First derivatives ---
        # np.gradient returns [dZ/dy, dZ/dx] by default → specify spacings explicitly
        SZV, SZU = np.gradient(Z, dy, dx) # dZ/dx, dZ/dy

        # --- Tangent vectors (parameterization basis, right-handed) ---
        SXU = np.ones_like(Z)
        SXV = np.zeros_like(Z)
        SYU = np.zeros_like(Z)
        SYV = np.ones_like(Z)

        SU = np.stack((SXU, SYU, SZU), axis=2)   # tangent in x (U)
        SV = np.stack((SXV, SYV, SZV), axis=2)   # tangent in y (V)

        # --- Normal vector (ensure right-handedness) ---
        # Use cross(SV, SU) instead of cross(SU, SV)
        CUV = np.cross(SU, SV, axis=2)
        AC = np.linalg.norm(CUV, axis=2)
        NX, NY, NZ = (CUV[..., 0] / AC, CUV[..., 1] / AC, CUV[..., 2] / AC)


        # --- 5. Normal derivatives ---
        NXV, NXU = np.gradient(NX, dy, dx)
        NYV, NYU = np.gradient(NY, dy, dx)
        NZV, NZU = np.gradient(NZ, dy, dx)

        # --- 6. First fundamental form coefficients ---
        E = np.einsum('ijk,ijk->ij', SU, SU)
        F = np.einsum('ijk,ijk->ij', SU, SV)
        G = np.einsum('ijk,ijk->ij', SV, SV)

        # --- 7. Second fundamental form coefficients ---
        # Compute normal derivatives (already: NXU, NXV, NYU, NYV, NZU, NZV)
        NU = np.stack((NXU, NYU, NZU), axis=2)  # ∂N/∂u (x-direction)
        NV = np.stack((NXV, NYV, NZV), axis=2)  # ∂N/∂v (y-direction)

        # Second fundamental form coefficients
        e = -np.einsum('ijk,ijk->ij', NU, SU)
        f = -0.5 * (np.einsum('ijk,ijk->ij', NU, SV) + np.einsum('ijk,ijk->ij', NV, SU))
        g = -np.einsum('ijk,ijk->ij', NV, SV)

        # --- 8. Principal curvatures ---
        a = E * G - F**2
        b = -(g * E - 2 * f * F + e * G)
        c = e * g - f**2

        K2 = -(-b + np.sqrt(np.abs(b**2 - 4 * a * c))) / (2 * a)
        K1 = -(-b - np.sqrt(np.abs(b**2 - 4 * a * c))) / (2 * a)

        # --- 9. Threshold small values ---
        K1[np.abs(K1) <= kt] = 0
        K2[np.abs(K2) <= kt] = 0

        # --- 10. Derived curvature metrics ---
        KM = 0.5 * (K1 + K2)
        KG = K1 * K2

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
