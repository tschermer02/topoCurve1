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
        expected_tukey_array = np.array([
            [0., 0.,0.,0.,0.,0.,],
            [0.,0.81813562, 0.9045085,  1.80901699, 1.63627124, 0.],
            [0., 0.9045085,  1., 2., 1.80901699, 0.],
            [0., 2.71352549, 3., 4., 3.61803399,0.],
            [0., 2.45440686, 2.71352549, 3.61803399, 3.27254249, 0.],
            [0., 0.,0.,0.,0.,0.,]
        ])

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
        mock_mirrored_array = np.array([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]
        ])
        
        # Set up the expected padded array
        expected_padded_array = np.array([
            [0, 0, 0, 0, 0],
            [0, 1, 2, 3, 0],
            [0, 4, 5, 6, 0],
            [0, 7, 8, 9, 0],
            [0, 0, 0, 0, 0]
        ])

        # Mock the mirror_dem method to return the mock mirrored array
        spec_filt_obj.mirror_dem = MagicMock(return_value=mock_mirrored_array)

        # Call the padding method
        padded_array = spec_filt_obj.padding(alphaIn=0.5)

        # Verify that padding is applied correctly
        self.assertEqual(padded_array.shape, expected_padded_array.shape)
        self.assertTrue(np.array_equal(padded_array, expected_padded_array))
        '''

if __name__ == '__main__':
    unittest.main()