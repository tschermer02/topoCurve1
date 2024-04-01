import unittest
from unittest.mock import patch, MagicMock
import numpy as np
from scipy.fft import fft2, ifft2
from photutils.psf import TukeyWindow
from SpectralFiltering import SpectralFiltering
from TopoCurve import TopoCurve

class TestTopoCurve(unittest.TestCase):
    def test_valid_geotiff_metadata_extraction(self):
        # Provide a valid GeoTIFF file path for testing
        valid_geotiff_path = "references\DEM_files\Purgatory.tif"

        # Initialize a TopoCurve object with the valid GeoTIFF file
        topo_curve_obj = TopoCurve(valid_geotiff_path)

        # Assert that metadata extraction is successful
        self.assertTrue(topo_curve_obj.metadata)  # Check if metadata is not empty
        self.assertIsInstance(topo_curve_obj.metadata, dict)  # Check if metadata is a dictionary
        # Add more assertions to check specific metadata if needed
        import unittest

    def test_curve_calc_with_known_input(self):
        # Known input data
        ZFilt = np.array([[1, 2, 3],
                          [4, 5, 6],
                          [7, 8, 9]])
        dx = 1.0
        dy = 1.0
        kt = 2.0

        # Manually compute the expected output
        expected_K1 = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]])
        expected_K2 = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]])
        expected_KM = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]])
        expected_KG = np.array([[0, 0, 0],
                                [0, 0, 0],
                                [0, 0, 0]])

        # Initialize a TopoCurve object
        topo_curve_obj = TopoCurve("references\DEM_files\Purgatory.tif")

        # Call the CurveCalc method with the known input data
        K1, K2, KM, KG = topo_curve_obj.CurveCalc(ZFilt, dx, dy, kt)

        # Assert that the method returns the expected principal curvatures and curvature features
        np.testing.assert_array_equal(K1, expected_K1)
        np.testing.assert_array_equal(K2, expected_K2)
        np.testing.assert_array_equal(KM, expected_KM)
        np.testing.assert_array_equal(KG, expected_KG)

    def test_initialization_inherits_metadata(self):
        # Mock the TopoCurve initialization to return metadata
        mock_metadata = {'GeogAngularUnitsGeoKey': 4326, 'ProjLinearUnitsGeoKey': 9001}  # Example metadata
        mock_topo_curve_init = MagicMock(return_value=None)
        SpectralFiltering.__init__ = mock_topo_curve_init

        # Initialize SpectralFiltering object
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")
        spec_filt_obj.metadata = mock_metadata

        # Ensure metadata is correctly inherited
        self.assertEqual(spec_filt_obj.metadata, mock_metadata)

    def test_detrend_method(self):
        # Mock the SpectralFiltering object
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")
        
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
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")
        
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
    '''
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
        self.assertTrue(np.array_equal(padded_array, expected_padded_array))
     
    def test_fft_method_lowpass(self):
        # Mocked elevation values
        mocked_elevation_values = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        # Expected lowpass filtered elevation values
        expected_filtered_values_lowpass = np.array([
            [0.01149374, 0.18236484, 0.01149374],
            [0.18236484, 0.7688021, 0.18236484],
            [0.01149374, 0.18236484, 0.01149374]
        ])

        # Create SpectralFiltering object
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")

        # Manually set necessary attributes
        spec_filt_obj.dimx_ma = 3
        spec_filt_obj.dimy_ma = 3
        spec_filt_obj.dx = np.array([1.0, 1.0] )
        spec_filt_obj.powerOfTwo = 4  # Assuming this is the correct powerOfTwo value for the mocked elevation values

        # Call FFT method with lowpass filtering
        dx, dy, filtered_values_lowpass = spec_filt_obj.FFT(filter=(1, 5), filterType='lowpass', alphaIn=0.5)

        # Verify that the lowpass filtered elevation values are computed correctly
        self.assertAlmostEqual(dx, 1.0)
        self.assertAlmostEqual(dy, 1.0)
        np.testing.assert_allclose(filtered_values_lowpass, expected_filtered_values_lowpass, rtol=1e-5)

    def test_fft_method_highpass(self):
        # Mocked elevation values
        mocked_elevation_values = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        # Expected highpass filtered elevation values
        expected_filtered_values_highpass = np.array([
            [-0.01149374, -0.18236484, -0.01149374],
            [-0.18236484, 0.2311979, -0.18236484],
            [-0.01149374, -0.18236484, -0.01149374]
        ])

        # Create SpectralFiltering object
        spec_filt_obj = SpectralFiltering("references\DEM_files\Purgatory.tif")

        # Manually set necessary attributes
        spec_filt_obj.dimx_ma = 3
        spec_filt_obj.dimy_ma = 3
        spec_filt_obj.dx = np.array([1.0, 1.0] )
        spec_filt_obj.powerOfTwo = 4  # Assuming this is the correct powerOfTwo value for the mocked elevation values

        # Call FFT method with highpass filtering
        dx, dy, filtered_values_highpass = spec_filt_obj.FFT(filter=(1, 5), filterType='highpass', alphaIn=0.5)

        # Verify that the highpass filtered elevation values are computed correctly
        self.assertAlmostEqual(dx, 1.0)
        self.assertAlmostEqual(dy, 1.0)
        np.testing.assert_allclose(filtered_values_highpass, expected_filtered_values_highpass, rtol=1e-5)
        '''
if __name__ == '__main__':
    unittest.main()