from TopoCurve import TopoCurve
from SpectralFiltering import SpectralFiltering

# Instantiate TopoCurve object with a GeoTIFF file
dem = TopoCurve('path/to/your/file.tif')

# Calculate principal curvatures and curvature features
K1, K2, KM, KG = dem.CurveCalc(ZFilt, dx, dy, kt)

# Plot and save elevation values
dem.plot(input_array, 'output_image.png')

# Instantiate SpectralFiltering object with a GeoTIFF file
spectral_filter = SpectralFiltering('path/to/your/file.tif')

# Apply FFT filtering
dx, dy, filtered_elevation = spectral_filter.FFT(filter, filterType, alphaIn)