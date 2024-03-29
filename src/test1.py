from TopoCurve import TopoCurve
from SpectralFiltering import SpectralFiltering

# Instantiate TopoCurve object with a GeoTIFF file
dem = TopoCurve('references\DEM_files\Purgatory.tif')

# Instantiate SpectralFiltering object with a GeoTIFF file
spectral_filter = SpectralFiltering('references\DEM_files\Purgatory.tif')

# Apply FFT filtering
dx, dy, filtered_elevation = spectral_filter.FFT([190, 200], 'lowpass', 0)

# Calculate principal curvatures and curvature features
K1, K2, KM, KG = dem.CurveCalc(filtered_elevation, dx, dy, 0)

# Plot and save elevation values
dem.plot(filtered_elevation, 'output_image.png')
