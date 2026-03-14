#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 11:39:24 2025

@author: ntklema
"""
import os
import topocurve
from topocurve import TopoCurve,SpectralFiltering
import numpy as np
import matplotlib.pyplot as plt
import cmcrameri.cm as cmc

print(TopoCurve.__module__)

package_file_path = topocurve.__file__
print(f"Package file path: {package_file_path}")

# Define the path to the TIFF file
tiff_file = '/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Umpqua_10m_2.tif'

# Instantiate TopoCurve object
dem = TopoCurve(tiff_file=tiff_file)

# Instantiate SpectralFiltering object
spectral_filter = SpectralFiltering(tiff_file)

# Apply FFT filtering with a lowpass filter at 190-200
dx, dy, filtered_elevation = spectral_filter.FFT([150, 200], 'lowpass', 0)

# Compute curvature attributes
K = dem.CurveCalc(filtered_elevation, dx, dy, 0)


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



#%% Print results
KM=K[6]['KM']
KG=K[6]['KG']
LogA=np.log10(acc)
# LogA=LogA*BS
# LogA[np.where(np.isnan(LogA))]=0
SMAP=K[4]

color_list = {
    'B': (0.11, 0.30, 0.95, 1.0),
    'SS': (0.78, 0.95, 0.98, 1.0),
    'S': (0.55, 0.85, 0.90, 1.0),
    'P': (1.00, 1.00, 1.00, 1.0),
    'A': (0.80, 0.20, 0.15, 1.0),
    'AS': (0.95, 0.28, 0.26, 1.0),
    'D': (0.33, 0.00, 0.00, 1.0)
}
#%%
nb = 100

av=np.linspace(1.9,7,nb)


km=np.zeros([nb-1])
a=np.zeros([nb-1])
kg=np.zeros([nb-1])
a_pdf=np.zeros([nb-1])
p_b=np.zeros([nb-1])
p_d=np.zeros([nb-1])
p_as=np.zeros([nb-1])
p_ss=np.zeros([nb-1])



for i in np.arange(nb-1):
    ind=np.where((LogA >= av[i]) & (LogA < av[i+1]) & (BS!=0))
    a[i]=(av[i]+av[i+1])/2
    a_pdf[i]=np.size(ind)/np.size(np.where(BS!=0))
    km[i]=np.mean(KM[ind])
    kg[i]=np.mean(KG[ind])
    
    S=SMAP[ind]

    p_b[i]=np.size(np.where(S==-3))/np.size(S)
    
    p_d[i]=np.size(np.where(S==3))/np.size(S)
    
    p_ss[i]=np.size(np.where(S==-2))/np.size(S)
    
    p_as[i]=np.size(np.where(S==2))/np.size(S)

def ma_conv(data, window_size):
    ma = np.convolve(data, np.ones(window_size) / window_size, mode='same')
    
    return ma
sm_win=3
km = ma_conv(km, sm_win)
kg = ma_conv(kg, sm_win)
a_pdf = ma_conv(a_pdf, sm_win)

#% Shape Class conditional pdfs
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(10**a,p_b,color=color_list['B'],linewidth=2,label='Basins') 
ax.plot(10**a,p_d,color=color_list['D'],linewidth=2,label='Domes')
ax.plot(10**a,p_as,color=color_list['AS'],linewidth=2,label='Antiformal saddles')  
ax.plot(10**a,p_ss,color=color_list['SS'],linewidth=2,label='Synformal saddles')
ax.legend(bbox_to_anchor=(0.55, 0.5), loc='upper left')   
ax.set_xscale('log')   
ax.set_xlabel('Upstream Drainage Area ($m^2$)',fontsize=14)
ax.set_ylabel('P(C|A)',fontsize=14)
# ax.set_ylim(0,0.4)
ax.set_xlim(1e2,1e7)


#%%

fig, ax = plt.subplots(figsize=(10, 6))

ax.fill_between(10**a, a_pdf, color= "b", alpha= 0.8)
ax.set_xscale('log')
ax.set_xlabel('Drainage Area')
ax.set_xlim(10**2,10**7)
ax.set_ylim(0,0.07)

#%%
fig, ax = plt.subplots(figsize=(6.5, 4))

ax.plot(10**a,kg,color=[0.2, 0.4, 0],linewidth=2)
ax.set_xscale('log')
ax.axhline(y=0,color='k',linewidth=0.5)
ax.set_ylabel('Gaussian Curv. (m$^{-2}$)',fontsize=14,color=[0.2, 0.4, 0])
ax.tick_params(axis='y',labelsize='large',colors=[0.2, 0.4, 0])
ax.tick_params(axis='x',labelsize='large')
ax2 = ax.twinx() 
ax2.plot(10**a,km,color=[0.4, 0, 0.6],linewidth=2,ls='--')
ax2.set_ylim(-10e-3,10e-3)
ax2.set_ylabel('Mean Curv. (m$^{-1}$)',fontsize=14,color=[0.4, 0, 0.6])
ax2.tick_params(axis='y',labelsize='large',colors=[0.4, 0, 0.6])
ax2.ticklabel_format(axis='y', style='sci', scilimits=(-3,-3))
ax.set_xlabel('Upstream Drainage Area (m$^2$)',fontsize=14)
ax.set_xlim(10**2,10**7)
ax.set_ylim(-1.5e-5,1.5e-5)

# ax['d'].imshow(K[4])

#%%
f=dem.plot_smap(SMAP, tiff_file, 'title', '/Users/ntklema/Downloads')
f.xlim(-123.975,-123.9)





# Example usage:


