import numpy as np
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy
import cartopy.crs as ccrs
import pandas as pd

inpath='/input_obs_data/CFSR/'
outpath = 'data/outputs/'
outpath_plots = 'plots/'


#number of days (6-h periods) with temperature > 0 (or -5) before January (November-December-possibly also January)


#make the list of all files
fl = sorted(glob(inpath+'cfsr.6h.*.nc', recursive=True))
#print(fl)

date_list = []
warm_day_list = []

#get the land and NP-hole mask (north pole hole got smaller over time!!!)
f = Dataset(fl[0])
temp = f.variables['TMP_L103']
lats = f.variables['lat'][:]
lons = f.variables['lon'][:]

#loop over all the years
#accummulate ice type on maps for each January
years = range(2005,2016)
for yr in years:
    print(yr)
    fl = sorted(glob(inpath+'cfsr.6h.'+str(yr)+'*.nc'))
    
    freq = np.zeros_like(temp[0,:,:])
    times = 0
    
    for i in range(0,len(fl[:])):
        f = fl[i]
        print(f)
        
        #select only Nov-Dec
        mon = f.split('.')[2].split(str(yr))[-1]
        if int(mon) < 11: continue
        print(mon)
        
        f = Dataset(f)
        temp = f.variables['TMP_L103'][:,:,:]-273.15    #convert from K to C
        
        warm = np.where(temp>0,1,0)
        times = times+warm.shape[0]
        freq = freq+np.sum(warm,axis=0)
        
    freq = freq/times
    
    print(times)
    print(freq.shape)
    
    #plot a map
    #%% Make the figure

    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    pp = plt.pcolormesh(lons,lats,freq,vmin=0,vmax=.05, cmap='jet', transform=ccrs.PlateCarree())


    # add the colourbar to the bottom of the plot.
    # The first moves the bottom of the map up to 15% of total figure height, 
    # the second makes the new axes for the colourbar, 
    # the third makes the colourbar, and the final adds the label
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(0,1.1,0.1))
    cbar.set_label(label='warm-6h fraction',size=14, family='serif')
    
    output = 'cfsr_warm_freq_'+str(yr)+'.png'
    plt.savefig(outpath_plots+output,bbox_inches='tight')
    
