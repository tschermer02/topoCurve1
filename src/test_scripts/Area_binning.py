#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 11:39:24 2025

@author: ntklema
"""

from topocurve import TopoCurve,SpectralFiltering
import numpy as np
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc

# Define the path to the TIFF file
tiff_file = '/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Umpqua_10m_2.tif'

# Instantiate TopoCurve object
dem = TopoCurve(tiff_file=tiff_file)

# Instantiate SpectralFiltering object
spectral_filter = SpectralFiltering(tiff_file)

#%% Route Flow
from pysheds.grid import Grid
# import seaborn as sns


# Read raw DEM
grid = Grid.from_raster(tiff_file)
dem_ac = grid.read_raster(tiff_file)

# Fill depressions
flooded_dem = grid.fill_depressions(dem_ac)

# Resolve flats
inflated_dem = grid.resolve_flats(flooded_dem)

fdir = grid.flowdir(inflated_dem,routing='dinf')
acc = grid.accumulation(fdir,routing='dinf')*grid.affine._scaling[0]**2

#%% Extract Catchment

fdir2 = grid.flowdir(inflated_dem)
# Specify "pour point"
x,y = 423786,4837890

# Delineate the catchment
catch = grid.catchment(x=x, y=y, fdir=fdir2, xytype='coordinate')
#%%
import matplotlib.pyplot as plt
import matplotlib.colors as colors

fig, ax = plt.subplots(figsize=(8,6))
fig.patch.set_alpha(0)
plt.grid('on', zorder=0)
im = ax.imshow(acc, extent=grid.extent, zorder=2,
               cmap='cubehelix',
               norm=colors.LogNorm(1, acc.max()),
               interpolation='bilinear')
plt.colorbar(im, ax=ax, label='Upstream Area')
plt.title('Flow Accumulation', size=14)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.tight_layout()

#%% Rasterize ROI polygon
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

# Import mask shapefile
mask_path='/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Umpqua_10m_2_Basins.shp'
mask = gpd.read_file(mask_path)

# Open the DEM to use its properties as a template
with rasterio.open(tiff_file) as dem_src:
    dem_transform = dem_src.transform
    dem_crs = dem_src.crs
    dem_shape = dem_src.shape
    dem_meta = dem_src.meta.copy()
    
# Check projection
# if mask.crs != dem_crs:
#     mask = mask.to_crs(dem_crs)    

#      
poly = [i for i in mask.geometry]

BS = rasterize(
    poly,
    out_shape=dem_shape,
    transform=dem_transform,
    fill=0,  # Background value
    all_touched=True, # Include all pixels touched by the polygons
    dtype=rasterio.uint8 # Data type for the output raster
)

#%%

# Apply FFT filtering with a lowpass filter at 190-200
dx, dy, filtered_elevation = spectral_filter.FFT([150, 200], 'lowpass', 0)

# Compute curvature attributes
K = dem.CurveCalc(filtered_elevation, dx, dy, 0)

#%% Print results
KM=K[6]['KM']
KG=K[6]['KG']
LogA=np.log10(acc)

nb = 100

av=np.linspace(2,np.max(LogA[np.isfinite(LogA)]),nb)


km=np.zeros([nb-1])
a=np.zeros([nb-1])
kg=np.zeros([nb-1])



for i in np.arange(nb-1):
    ind=np.where((LogA >= av[i]) & (LogA < av[i+1]))
    a[i]=(av[i]+av[i+1])/2
    km[i]=np.mean(KM[ind])
    kg[i]=np.mean(KG[ind])

    
    







