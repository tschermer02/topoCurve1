import numpy as np
import unittest
from TopoCurve import TopoCurve
from SpectralFiltering import SpectralFiltering
class TestTopoCurve(unittest.TestCase):

    def test_Topocurve_initialization(self):
        # Test initialization of TopoCurve class with a valid GeoTIFF file
        tiff_file = "data\DEM_files\Purgatory.tif"
        topocurve = TopoCurve(tiff_file)
        
        # Assert metadata is properly loaded
        self.assertIsInstance(topocurve.metadata, dict)
        
        # Assert elevation array is properly loaded
        self.assertIsInstance(topocurve.z_array, np.ndarray)
        
        # Assert dimensions are properly set
        self.assertIsInstance(topocurve.dimx, int)
        self.assertIsInstance(topocurve.dimy, int)
        
        # Assert grid spacing is properly set
        self.assertIsInstance(topocurve.dx, float)
        self.assertIsInstance(topocurve.dx_dy, float)

    def test_SpectralFiltering_operations(self):
        # Create a mock elevation array
        elevation_array = np.random.rand(100, 100)
        
        # Initialize SpectralFiltering object
        spectral_filtering = SpectralFiltering("data\DEM_files\Purgatory.tif")
        
        # Test detrending operation
        Z_detrended, plane = spectral_filtering.detrend(elevation_array)
        self.assertEqual(Z_detrended.shape, elevation_array.shape)
        self.assertEqual(plane.shape, elevation_array.shape)
        
        # Test mirroring operation
        mirrored_array = spectral_filtering.mirror_dem(elevation_array)
        self.assertEqual(mirrored_array.shape, (300, 300))  # Check dimensions after mirroring
        
        # Test Tukey window application
        tukey_array = spectral_filtering.tukeyWindow(mirrored_array, alphaIn=0.5)
        self.assertEqual(tukey_array.shape, mirrored_array.shape)
        
        # Test padding operation
        padded_array = spectral_filtering.padding(tukey_array, alphaIn=0.5)
        self.assertEqual(padded_array.shape, (512, 512))  # Assuming padding to the next power of two
        
        # Test FFT filtering operation
        dx, dy, ZFilt = spectral_filtering.FFT(padded_array, filter=10, filterType='lowpass', alphaIn=0.5)
        self.assertIsInstance(dx, float)
        self.assertIsInstance(dy, float)
        self.assertIsInstance(ZFilt, np.ndarray)

if __name__ == '__main__':
    unittest.main()
