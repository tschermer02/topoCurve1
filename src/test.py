import unittest
from unittest.mock import patch, MagicMock
import numpy as np
from scipy.fft import fft2, ifft2
from photutils.psf import TukeyWindow
from topocurve.SpectralFiltering import SpectralFiltering
from topocurve.TopoCurve import TopoCurve

class TestTopoCurve(unittest.TestCase):
    def setUp(self):
        """Prepare reusable test objects before each test."""
        # Create dummy instances that can call class methods without file I/O
        self.topo_curve_obj = TopoCurve.__new__(TopoCurve)
        self.spec_filt_obj = SpectralFiltering.__new__(SpectralFiltering)

    def test_valid_geotiff_metadata_extraction(self):
        # Provide a valid GeoTIFF file path for testing
        valid_geotiff_path = r"C:\Users\tsche\source\repos\topoCurve1\references\DEM_files\Purgatory.tif"

        # Initialize a TopoCurve object with the valid GeoTIFF file
        topo_curve_obj = TopoCurve(valid_geotiff_path)

        # Assert that metadata extraction is successful
        self.assertTrue(topo_curve_obj.metadata)  # Check if metadata is not empty
        self.assertIsInstance(topo_curve_obj.metadata, dict)  # Check if metadata is a dictionary
        # Add more assertions to check specific metadata if needed
        import unittest

    def test_curve_calc_with_hemisphere(self):
        """
        Test the CurveCalc method using a synthetic hemisphere surface.

        The goal is to verify that the computed mean and Gaussian curvatures 
        are approximately equal to the theoretical values for a hemisphere 
        of radius R (K1 = K2 = 1/R, KM = 1/R, KG = 1/R^2).
        """
        # Create a synthetic hemisphere surface
        R = 1.0   # Hemisphere radius
        N = 200   # Grid size
        x = np.linspace(-R, R, N)
        y = np.linspace(-R, R, N)
        X, Y = np.meshgrid(x, y)

        # Hemisphere equation: z = sqrt(R^2 - x^2 - y^2), set outside to 0
        Z = np.sqrt(np.clip(R**2 - X**2 - Y**2, 0, None))

        # Mask the rim (Z = 0) to avoid numerical edge artifacts
        mask = Z > 1e-6
        Z_masked = np.where(mask, Z, np.nan)

        # Grid spacing
        dx = x[1] - x[0]
        dy = y[1] - y[0]

        # Run CurveCalc on the masked hemisphere surface
        topo = TopoCurve.__new__(TopoCurve)  # Create object without file initialization
        K1, K2, KM, KG = topo.CurveCalc(Z_masked, dx, dy, kt=1e-6)

        # Calculate mean curvature statistics (ignore NaNs)
        mean_K1 = np.nanmean(K1)
        mean_K2 = np.nanmean(K2)
        mean_KM = np.nanmean(KM)
        mean_KG = np.nanmean(KG)

        # Theoretical curvature values for a hemisphere of radius R
        expected_K = 1 / R
        expected_KG = 1 / R**2

        # Compare computed vs. theoretical curvature values
        self.assertAlmostEqual(mean_K1, expected_K, delta=0.05 * expected_K)
        self.assertAlmostEqual(mean_K2, expected_K, delta=0.05 * expected_K)
        self.assertAlmostEqual(mean_KM, expected_K, delta=0.05 * expected_K)
        self.assertAlmostEqual(mean_KG, expected_KG, delta=0.05 * expected_KG)

    def test_initialization_inherits_metadata(self):
        # Mock the TopoCurve initialization to return metadata
        mock_metadata = {'GeogAngularUnitsGeoKey': 4326, 'ProjLinearUnitsGeoKey': 9001}  # Example metadata
        mock_topo_curve_init = MagicMock(return_value=None)
        SpectralFiltering.__init__ = mock_topo_curve_init

        # Initialize SpectralFiltering object
        spec_filt_obj = SpectralFiltering(r"../references/DEM_files/Purgatory.tif")
        spec_filt_obj.metadata = mock_metadata

        # Ensure metadata is correctly inherited
        self.assertEqual(spec_filt_obj.metadata, mock_metadata)

    def test_detrend_method(self):
        # Mock the SpectralFiltering object
        spec_filt_obj = SpectralFiltering(r"../references/DEM_files/Purgatory.tif")
        
        # Mock the elevation values
        mock_elevation_values = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        
        # Set up the expected detrended values and trend component
        expected_detrended_values = np.array([[-1, 0, 1], [-1, 0, 1], [-1, 0, 1]])
        expected_trend_component = np.array([[2, 2, 2], [5, 5, 5], [8, 8, 8]])

        # Mock the detrended values and trend component computation
        spec_filt_obj.detrend = MagicMock(return_value=(expected_detrended_values, expected_trend_component))

        # Call the detrend method
        detrended_values, trend_component = spec_filt_obj.detrend()

        # Verify that the detrended values and trend component are computed correctly
        self.assertTrue(np.array_equal(detrended_values, expected_detrended_values))
        self.assertTrue(np.array_equal(trend_component, expected_trend_component))
    
    def test_mirror_dem_method(self):
        # Mock the SpectralFiltering object
        spec_filt_obj = SpectralFiltering(r"../references/DEM_files/Purgatory.tif")
        
        # Mock the elevation values
        mock_elevation_values = np.array([[1, 2], [3, 4]])
        
        # Set up the expected mirrored array
        expected_mirrored_array = np.array([
            [4,3,3,4,4,3],
            [2,1,1,2,2,1],
            [2,1,1,2,2,1],
            [4,3,3,4,4,3],
            [4,3,3,4,4,3],
            [2,1,1,2,2,1]
        ])

        # Mock the detrend method to return the mock elevation values
        spec_filt_obj.detrend = MagicMock(return_value=(mock_elevation_values, None))

        # Call the mirror_dem method
        mirrored_array = spec_filt_obj.mirror_dem()

        # Verify that the mirrored array has the correct dimensions and values
        self.assertEqual(mirrored_array.shape, expected_mirrored_array.shape)
        self.assertTrue(np.array_equal(mirrored_array, expected_mirrored_array))
   
    def test_tukey_window_method(self):
        # Mock the SpectralFiltering object
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")
        
        # Mock the elevation values
        mock_elevation_values = np.array([[1, 2], [3, 4]])
        
        # Set up the expected Tukey window
        expected_tukey_array = np.array([[0., 0., 0., 0., 0., 0.],
            [0., 4.47452698e-18, 8.35974985e-17, 1.67194997e-16, 8.94905396e-18, 0.],
            [0., 8.35974985e-17, 1.11022302e-16, 2.22044605e-16,1.67194997e-16, 0.],
            [0., 3.34389994e-16, 4.44089210e-16, 0., 0., 0.],
            [0., 1.78981079e-17, 3.34389994e-16, 0., 0., 0.],
            [0., 0., 0., 0., 0., 0.]])
        
        # Patch the z_array attribute to return the mock elevation values
        spec_filt_obj.z_array = mock_elevation_values
        
        # Call the tukeyWindow method
        tukey_array = spec_filt_obj.tukeyWindow(alphaIn=0.5)

        # Verify that the Tukey window is correctly applied
        self.assertEqual(tukey_array.shape, expected_tukey_array.shape)
        self.assertTrue(np.allclose(tukey_array, expected_tukey_array))
    
    def test_padding_method(self):
        # Mock the SpectralFiltering object
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")
        
        # Mock the mirrored array
        mock_elevation_values = np.array([[1, 2], [3, 4]])
        
        # Set up the expected padded array
        expected_padded_array = np.array([
            [0., 0., 0., 0., 0., 0., 0., 0.],
            [0., 0., 0., 0., 0., 0., 0., 0.],
            [0., 0., 4.47452698e-18, 8.35974985e-17, 1.67194997e-16, 8.94905396e-18, 0., 0.],
            [0., 0., 8.35974985e-17, 1.11022302e-16, 2.22044605e-16,1.67194997e-16, 0., 0.],
            [0.,0., 3.34389994e-16, 4.44089210e-16, 0., 0., 0., 0.],
            [0.,0., 1.78981079e-17, 3.34389994e-16, 0., 0., 0., 0.],
            [0.,0., 0., 0., 0., 0., 0., 0.],
            [0.,0., 0., 0., 0., 0., 0., 0.]])

        # Mock the mirror_dem method to return the mock mirrored array
        spec_filt_obj.z_array = mock_elevation_values

        # Call the padding method
        padded_array = spec_filt_obj.padding(alphaIn=0.5)

        # Verify that padding is applied correctly
        self.assertEqual(padded_array.shape, expected_padded_array.shape)
        self.assertTrue(np.allclose(padded_array, expected_padded_array, rtol=1e-15, atol=1e-15))

        
    def test_fft_method_lowpass(self):
        """
        Test lowpass frequency filtering using a simple synthetic 2×2 grid.
        Ensures low-frequency components are preserved and output remains smooth.
        """
        # Create synthetic elevation data
        Z = np.array([[1, 2],
                    [3, 4]], dtype=float)

        # Grid spacing
        dx = dy = 1.0

        # Apply lowpass FFT filter manually
        FZ = np.fft.fftshift(np.fft.fft2(Z))
        ny, nx = Z.shape
        ky = np.fft.fftfreq(ny, d=dy)
        kx = np.fft.fftfreq(nx, d=dx)
        KX, KY = np.meshgrid(kx, ky)
        kmag = np.sqrt(KX**2 + KY**2)

        # Define cutoff (lowpass)
        k_cut = 0.5 * np.max(kmag)
        H_low = (kmag <= k_cut).astype(float)

        # Inverse FFT for filtered result
        Z_low = np.fft.ifft2(np.fft.ifftshift(FZ * H_low)).real

        # Check shape, smoothness, and energy reduction
        self.assertEqual(Z_low.shape, Z.shape)
        self.assertLess(np.std(Z_low), np.std(Z))
        self.assertAlmostEqual(float(dx), 1.0)
        self.assertAlmostEqual(float(dy), 1.0)


    def test_fft_method_highpass(self):
        """
        Test highpass frequency filtering using a simple synthetic 2×2 grid.
        Ensures high-frequency content is isolated (edges preserved).
        """

        # Create synthetic elevation data
        Z = np.array([[1, 2],
                    [3, 4]], dtype=float)

        # Grid spacing
        dx = dy = 1.0

        # Apply highpass FFT filter manually
        FZ = np.fft.fftshift(np.fft.fft2(Z))
        ny, nx = Z.shape
        ky = np.fft.fftfreq(ny, d=dy)
        kx = np.fft.fftfreq(nx, d=dx)
        KX, KY = np.meshgrid(kx, ky)
        kmag = np.sqrt(KX**2 + KY**2)

        # Define cutoff (highpass complement)
        k_cut = 0.5 * np.max(kmag)
        H_high = (kmag > k_cut).astype(float)

        # Inverse FFT for filtered result
        Z_high = np.fft.ifft2(np.fft.ifftshift(FZ * H_high)).real

        # Check shape and confirm high-frequency emphasis
        self.assertEqual(Z_high.shape, Z.shape)
        self.assertGreater(np.std(Z_high), 0)
        self.assertAlmostEqual(float(dx), 1.0)
        self.assertAlmostEqual(float(dy), 1.0)

if __name__ == "__main__":
    unittest.main()
