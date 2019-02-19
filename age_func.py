import numpy as np
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import pyresample as pr
import scipy.ndimage as ndimage
from pyproj import Proj, transform



def plot_pcolormesh(lons,lats,var,outname,vmin=None,vmax=None,cmap='jet',label='Variable'):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 66, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    pp = plt.pcolormesh(lons,lats,var,vmin=vmin,vmax=vmax, cmap=cmap, transform=ccrs.PlateCarree())


    # add the colourbar to the bottom of the plot.
    # The first moves the bottom of the map up to 15% of total figure height, 
    # the second makes the new axes for the colourbar, 
    # the third makes the colourbar, and the final adds the label
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(0,1.1,0.1))
    cbar.set_label(label=label,size=14, family='serif')
    
    plt.savefig(outname,bbox_inches='tight')
    
def plot_contour(lons,lats,data,levels=[.15],colors=['purple'],lw=[1.],labels=['Variable'],outname='test.png'):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 66, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    for i in range(len(levels)):   
        cs = plt.contour(lons,lats,data[i],levels=[levels[i]], colors=colors[i], linewidths=lw[i], transform=ccrs.PlateCarree())
        cs.collections[0].set_label(labels[i])
        
    ax1.legend(loc='upper left')
    
    plt.savefig(outname,bbox_inches='tight')
    
def plot_contour_bg(lons,lats,bg,data,levels=[.15],colors=['purple'],lw=[1.],labels=['Variable'],bg_label='Snow_depth',outname='test.png'):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 66, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    #plot 'background'
    pp = plt.pcolormesh(lons,lats,bg,vmin=0,vmax=1, cmap='jet', transform=ccrs.PlateCarree())
    # add the colourbar to the bottom of the plot.
    
    # plot sea ice field
    for i in range(len(levels)):   
        cs = plt.contour(lons,lats,data[i],levels=[levels[i]], colors=colors[i], linewidths=lw[i], transform=ccrs.PlateCarree())
        cs.collections[0].set_label(labels[i])
        
    ax1.legend(loc='upper left')
    
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(0,1.1,0.1))
    cbar.set_label(label=bg_label,size=14, family='serif')

    
    plt.savefig(outname,bbox_inches='tight')
    
def smooth_data(data,lon,lat,coarse_lon,coarse_lat):    
    #smoothen the data for nicer contours
    data = ndimage.gaussian_filter(data, sigma=2, order=0)
    #regrid to equally spaced grid in latlon - otherwise there will be problems with cyclic point in contour plots
    orig_def = pr.geometry.SwathDefinition(lons=lon, lats=lat)
    targ_def = pr.geometry.SwathDefinition(lons=coarse_lon, lats=coarse_lat)
    coarse_def = pr.geometry.SwathDefinition(lons=coarse_lon[::5,::5], lats=coarse_lat[::5,::5])
    data_coarse = pr.kd_tree.resample_nearest(orig_def, data, coarse_def, radius_of_influence=50000, fill_value=0)
    #fill all nans with 0 >> closed contours
    data_coarse = np.nan_to_num(data_coarse)
    data_smooth = pr.kd_tree.resample_gauss(coarse_def, data_coarse, targ_def, radius_of_influence=500000, neighbours=10, sigmas=250000, fill_value=0)
    
    #plot_pcolormesh(lon,lat,myi,'test.png',vmin=0,vmax=1,label='MYI fraction') 
    #plot_pcolormesh(lon_g[::5],lat_g[::5],myi_coarse,'test1.png',vmin=0,vmax=1,label='MYI fraction')    
    #plot_pcolormesh(lon_g,lat_g,myi_smooth,'test2.png',vmin=0,vmax=1,label='MYI fraction')
    #plot_contour(lon_g,lat_g,myi_smooth,'test3.png',levels=[.1], lw=[10], label='MYI extent')

    
    return(data_smooth)


#from pynextsim.projection_info import ProjectionInfo
#mm = ProjectionInfo(f)
###def __init__(self,
        ###ecc    = 0.081816153,
        ###a      = 6378.273e3,
        ###lat_0  = 90.,
        ###lon_0  = -45.,
        ###lat_ts = 60.,
        ###proj='stere',
##print(mm.proj,mm.lat_ts,mm.lat_0,mm.lon_0,mm.a,mm.ecc)
#nextsim_proj = '+proj=%s +lat_ts=%f +lat_0=%f +lon_0=%f +a=%f +e=%f +units=m' %(mm.proj,mm.lat_ts,mm.lat_0,mm.lon_0,mm.a,0.081816153)
#inProj  = Proj(nextsim_proj,preserve_units=True)
#outProj = Proj("+init=EPSG:4326") # WGS84 in degrees
#lonc,latc=transform(inProj,outProj,xc,yc)

