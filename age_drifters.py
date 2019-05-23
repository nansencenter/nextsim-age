import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample
import pyproj

from age_func import *

#evaluate neXtSIM drift based on the output in the drifters. Compare to OSI-SAF

#compare April daily means

inpath='/input_obs_data/polona/FRASIL/age_datamor_long/'
inpath='data/drifters/'
outpath = 'data/outputs/'
outpath_plots = 'plots/new/'

icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
drosi_path = '/input_obs_data/data/OSISAF_ice_drift/'

#get OSI-SAF grid
fn = drosi_path+'2007/01/ice-drift_ice_drift_nh_polstere-625_multi-oi_200701011200-200701031200.nc'
f = Dataset(fn)
lat_osi = f.variables['lat'][:]
lon_osi = f.variables['lon'][:]

#get all model drifters
fl = sorted(glob(inpath+'OSISAF_*0415.nc'))
print(fl)

for fn in fl:
    print(fn)
    f = Dataset(fn)
    lats0 = f.variables['latitude'][0,:,0]
    lons0 = f.variables['longitude'][0,:,0]
    lats1 = f.variables['latitude'][1,:,0]
    lons1 = f.variables['longitude'][1,:,0]
    time = f.variables['time'][:]
    #index = f.variables['index'][0,:,0]
    #sic = f.variables['sic'][0,:,0]
    
    #project lat,lon coordinates and calculate displacements
    #use OSI-SAF projection: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
    wgs84=pyproj.Proj("+init=EPSG:4326") 
    nh_stere=pyproj.Proj("+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45")
    x0,y0 = pyproj.transform(wgs84, nh_stere,lons0,lats0)
    x1,y1 = pyproj.transform(wgs84, nh_stere,lons1,lats1)
    dx = x0-x1
    dy = y0-y1
    
    print(dx)
    print(dy)
    
    #put displacements on a regular grid (pyresample) - they should be very close in space and no information lost by interpolation
    swath_def = pyresample.geometry.SwathDefinition(lons=lons0, lats=lats0)
    targ_def = pyresample.geometry.SwathDefinition(lons=lon_osi, lats=lat_osi)
    
    dx_g = pyresample.kd_tree.resample_nearest(swath_def, dx, targ_def, radius_of_influence=62.5)
    dy_g = pyresample.kd_tree.resample_nearest(swath_def, dy, targ_def, radius_of_influence=62.5)

    #get velocities
    tm = (time[0]-time[1])*24*60*60 #this is exactly 2 days anyway
    u = dx_g/tm
    v = dy_g/tm

    speed = np.sqrt(u**2+v**2)
    ##Plot data
    #outname = outpath_plots+'nextsim_drifters_test.png'
    #plot_pcolormesh(lon_osi,lat_osi,speed,outname,cmap='viridis',label='MYI')
    #exit() 
    
    #something is wrong, arrows conmverge at the NP - some plotting of calculation problem!!!
    #quiver plot
    outname = outpath_plots+'nextsim_drifters_test_arrows.png'
    plot_quiver(lon_osi,lat_osi,speed,u,v,outname,cmap='viridis',label='MYI')
    exit()

    

    
    #collect all velocites for the month/winter and correlate with OSI-SAF
    #make correlation maps (for speed and direction)

    exit()
#make PDFs






#make scatter plots

#make correlation maps
