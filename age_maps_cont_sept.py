import numpy as np
from glob import glob
from datetime import datetime
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

from age_func import *

#make April 1 maps
inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
outpath_plots = 'plots/new/'

fl = sorted(glob(inpath+'field*0915T000000Z.bin'))
print(fl)

#make a nice cyclic equally spaced lat-lon grid for plotting contoured fields with cartopy
lat_g = np.arange(90,62,-.1)
lon_g = np.arange(-180,181,1.)   
lon_gm, lat_gm = np.meshgrid(lon_g,lat_g)


snow_bias = []
ridge_bias = []
area_bias = []
dates = []

for f in fl:
    tmp = f.split('_')[-1].split('T')[0]
    date = datetime.strptime(tmp, "%Y%m%d")
    stamp = datetime.strftime(date, "%Y%m%d")+'1200'
    print(date)
    year = str(date.year)
    if int(year)<2005: continue
    
    #neXtSIM grid
    nb = NextsimBin(f)
    x,y = nb.get_xy_grids()
    lon,lat=nb.mapping(x,y,inverse=True)
    
    #neXtSIM data
    ic = nb.get_gridded_vars(['Concentration'],x,y)['Concentration']
    
    #OSI-SAF data
    fosi = icosi_path+year+'_nh_polstere/ice_conc_nh_polstere-100_multi_'+stamp+'.nc'
    f = Dataset(fosi)
    sf_osi = f.variables['status_flag'][0,:,:]
    lat_osi = f.variables['lat'][:]
    lon_osi = f.variables['lon'][:]
    ic_osi =  f.variables['ice_conc'][:]/100
    
    #smoothen the data for nicer contours
    ic_smooth = smooth_data(ic,lon,lat,lon_gm,lat_gm)
    ic_osi_smooth = smooth_data(ic_osi,lon_osi,lat_osi,lon_gm,lat_gm) 
    
    ##constrain with OSI-SAF area
    #landmask = np.where(sf_osi>20,1,0)
    #lm_smooth = smooth_data(landmask,lon_osi,lat_osi,lon_gm,lat_gm)
    #lm_smooth = np.where(lm_smooth>0,1,0)
    #ic_osi_smooth = np.where(lm_smooth,0,ic_osi_smooth)
    #ic_smooth = np.where(lm_smooth,0,ic_smooth)
    
    
    #compare September sea ice extent
    plot_contour(lon_g,lat_g,data=[ic_smooth,ic_osi_smooth],levels=[.15,.15],colors=['red','purple'], lw=[3,3], \
                 labels=['neXtSIM ice extent','OSI-SAF ice extent'],outname=outpath_plots+'extent_sept_'+year+'.png')
    
    #plot_contour_bg(lon_g,lat_g,ic_osi_smooth,data=[ic_smooth,ic_osi_smooth],levels=[.15,.15],colors=['red','purple'], lw=[3,3], \
                 #labels=['neXtSIM ice extent','OSI-SAF ice extent'],bg_label='Ridge ratio',outname=outpath_plots+'extent_sept_'+year+'.png')   
    
    #exit()
