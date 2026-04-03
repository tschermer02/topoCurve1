#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 11:39:24 2025

@author: ntklema
"""

from topocurve import TopoCurve,SpectralFiltering
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.colors as colors
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

# -------------- User inputs for generating TopoCurve object ------------

# Path to the TIFF file (if using geotiff for DEM)
tiff_file = '/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Umpqua_10m_2.tif'

# Define the path to a mask shapefile to select ROI
mask_path='/Users/ntklema/Library/CloudStorage/OneDrive-FortLewisCollege/Research_Projects/Curvature/ESurf_Paper/Data/Umpqua_10m_2_Basins.shp'
mask = gpd.read_file(mask_path)

#------------- Filter DEM and calculate curvatures------------------
# Instantiate TopoCurve object
dem = TopoCurve(tiff_file=tiff_file)

# Instantiate SpectralFiltering object
spectral_filter = SpectralFiltering(tiff_file)

# Apply FFT filtering with a lowpass filter
filter=[150,200] # Low pass filter cutoffs
dx, dy, filtered_elevation = spectral_filter.FFT(filter, 'lowpass', 0)

# Compute curvature attributes
K = dem.CurveCalc(filtered_elevation, dx, dy, 0)

# ----------- Rasterize shapefile with ROI polygon -----------------
# Open the DEM to use its properties as a template
with rasterio.open(tiff_file) as dem_src:
    dem_transform = dem_src.transform
    dem_crs = dem_src.crs
    dem_shape = dem_src.shape
    dem_meta = dem_src.meta.copy()
             
poly = [i for i in mask.geometry]
BS = rasterize(
    poly,
    out_shape=dem_shape,
    transform=dem_transform,
    fill=0,  # Background value
    all_touched=True, # Include all pixels touched by the polygons
    dtype=rasterio.uint8 # Data type for the output raster
)
#%% ---------- Route Flow to generate upstream area grid -----------
from pysheds.grid import Grid

# Read raw DEM
grid = Grid.from_raster(tiff_file)
dem_ac = grid.read_raster(tiff_file)

# Fill depressions
flooded_dem = grid.fill_depressions(dem_ac)

# Resolve flats
inflated_dem = grid.resolve_flats(flooded_dem)

# Route flow (use 'routing' field to adjust flow routing algorithm)
fdir = grid.flowdir(inflated_dem,routing='dinf')
acc = grid.accumulation(fdir,routing='dinf')*grid.affine._scaling[0]**2
LogA=np.log10(acc) # Log transform drainage area data

# fig, ax = plt.subplots(figsize=(8,6))
# fig.patch.set_alpha(0)
# plt.grid('on', zorder=0)
# im = ax.imshow(acc, extent=grid.extent, zorder=2,
#                cmap='cubehelix',
#                norm=colors.LogNorm(1, acc.max()),
#                interpolation='bilinear')
# plt.colorbar(im, ax=ax, label='Upstream Area')
# plt.title('Flow Accumulation', size=14)
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.tight_layout()


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
mn_log_A =1.9
max_log_A=7.3


av=np.linspace(mn_log_A,max_log_A,nb)


a=np.zeros([nb-1])
km=np.zeros([nb-1])
kg=np.zeros([nb-1])
k1=np.zeros([nb-1])
k2=np.zeros([nb-1])
a_pdf=np.zeros([nb-1])
p_b=np.zeros([nb-1])
p_d=np.zeros([nb-1])
p_as=np.zeros([nb-1])
p_ss=np.zeros([nb-1])



for i in np.arange(nb-1):
    ind=np.where((LogA >= av[i]) & (LogA < av[i+1]) & (BS!=0))
    a[i]=(av[i]+av[i+1])/2
    a_pdf[i]=np.size(ind)/np.size(np.where(BS!=0))
    km[i]=np.mean(K[6]['KM'][ind])
    kg[i]=np.mean(K[6]['KG'][ind])
    k1[i]=np.mean(K[6]['K1'][ind])
    k2[i]=np.mean(K[6]['K2'][ind])
    
    S=K[4][ind]
    p_b[i]=np.size(np.where(S==-3))/np.size(S)
    p_d[i]=np.size(np.where(S==3))/np.size(S)
    p_ss[i]=np.size(np.where(S==-2))/np.size(S)
    p_as[i]=np.size(np.where(S==2))/np.size(S)


sm_win=3
a_pdf = np.convolve(a_pdf, np.ones(sm_win) / sm_win, mode='same')
km = np.convolve(km, np.ones(sm_win) / sm_win, mode='same')
kg = np.convolve(kg, np.ones(sm_win) / sm_win, mode='same')
k1 = np.convolve(k1, np.ones(sm_win) / sm_win, mode='same')
k2 = np.convolve(k2, np.ones(sm_win) / sm_win, mode='same')
p_b = np.convolve(p_b, np.ones(sm_win) / sm_win, mode='same')
p_d = np.convolve(p_d, np.ones(sm_win) / sm_win, mode='same')
p_as = np.convolve(p_as, np.ones(sm_win) / sm_win, mode='same')
p_ss = np.convolve(p_ss, np.ones(sm_win) / sm_win, mode='same')


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


#%% Upstream area distribution/slope plot

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

#%%
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

color_list = {
    -3: (0.11, 0.30, 0.95, 1.0),
    -2: (0.78, 0.95, 0.98, 1.0),
    -1: (0.55, 0.85, 0.90, 1.0),
    0: (1.00, 1.00, 1.00, 1.0),
    1: (0.80, 0.20, 0.15, 1.0),
    2: (0.95, 0.28, 0.26, 1.0),
    3: (0.33, 0.00, 0.00, 1.0)
}

lon_lim=[-123.99, -123.91]
lat_lim=[43.71,43.76]

lon_min, lon_max, lat_min, lat_max=dem.get_latlon_extent()
E_min, E_max, N_min, N_max=dem.get_extent()
r,c=K[4].shape
lat=np.linspace(lat_min,lat_max,r)
lon=np.linspace(lon_min,lon_max,c)
N=dem.y
E=dem.x
in_c=np.where((lon>=lon_lim[0]) & (lon<=lon_lim[1]))[0]
in_r=np.where((lat>=lat_lim[0]) & (lat<=lat_lim[1]))[0]

extent=[E[np.min(in_c)],E[np.max(in_c)],N[np.min(in_r)],N[np.max(in_r)]]



#%%
keys = sorted(color_list.keys())
colors = [color_list[k] for k in keys]

cmap = ListedColormap(colors)
bounds = [-3.5, -2.5, -1.5, -0.5, 0.5, 1.5, 2.5, 3.5]
norm = BoundaryNorm(bounds, cmap.N)

Z = dem.z_array
azimuth = 315
altitude = 45
az = np.deg2rad(azimuth)
alt = np.deg2rad(altitude)

dZdy, dZdx = np.gradient(Z, dem.dy, dem.dx)
slope = np.pi / 2 - np.arctan(np.sqrt(dZdx * dZdx + dZdy * dZdy))
aspect = np.arctan2(-dZdx, dZdy)

shaded = (np.sin(alt) * np.sin(slope) +
        np.cos(alt) * np.cos(slope) * np.cos(az - aspect))

hillshade = (shaded - shaded.min()) / (shaded.max() - shaded.min())

plt.figure(figsize=(8, 8))
plt.imshow(hillshade, cmap="gray",extent=dem.get_extent(),origin="lower")
plt.imshow(K[4],cmap=cmap,norm=norm, alpha=0.6, extent=dem.get_extent(),origin="lower")
plt.xlim(E[np.min(in_c)],E[np.max(in_c)])
plt.ylim(N[np.min(in_r)],N[np.max(in_r)])

#%%


in_=np.where(BS==1)
S=K[4][in_]
nb=np.size(np.where(S==-3))/np.size(S)*100
nd=np.size(np.where(S==3))/np.size(S)*100
nss=np.size(np.where(S==-2))/np.size(S)*100
nas=np.size(np.where(S==2))/np.size(S)*100
ns=np.size(np.where(S==-1))/np.size(S)*100
na=np.size(np.where(S==1))/np.size(S)*100
np=np.size(np.where(S==0))/np.size(S)*100



color_list = {
    'a': (0, 0, 1, 1.0),
    'b': (0.4, 0.00, 0.00, 1.0),
    'c': (1, 0, 0, 1.0),
    'd': (0, 1, 1, 1.0)
      
}
keys = sorted(color_list.keys())
colors = [color_list[k] for k in keys]
plt.pie([nb,nd,nas,nss],colors=colors,startangle=5,
        labels=(r'$'+str(np.round(nb,1))+'\%$',
                r'$'+str(np.round(nd,1))+'\%$',
                r'$'+str(np.round(nas,1))+'\%$',
                r'$'+str(np.round(nss,1))+'\%$',))




# Example usage:


