import numpy as np
from glob import glob
from datetime import datetime
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

from age_func import *

#make April 1 maps
inpath = 'data/run04/'
outpath = 'data/outputs/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
tyosi_path = '/input_obs_data/data/OSISAF_ice_type/'
cfsr_path = '/input_obs_data/data/CFSR/'
outpath_plots = 'plots/run04/'

fl = sorted(glob(inpath+'field*0401T000000Z.bin'))
print(fl)

#make a nice cyclic equally spaced lat-lon grid for plotting contoured fields with cartopy
lat_g = np.arange(90,62,-.1)
lon_g = np.arange(-180,181,1.)   
lon_gm, lat_gm = np.meshgrid(lon_g,lat_g)

#read CFSR grid
f = Dataset(cfsr_path+'cfsr.6h.197901.nc')
lat_f = f.variables['lat'][:]
lon_f = f.variables['lon'][:]
lon_fm, lat_fm = np.meshgrid(lon_f,lat_f)

snow_bias = []
ridge_bias = []
warm_bias = []
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
    icthin = nb.get_gridded_vars(['Concentration_thin_ice'],x,y)['Concentration_thin_ice']
    fyi = nb.get_gridded_vars(['Fyi_fraction'],x,y)['Fyi_fraction']
    myi = ic-icthin-fyi
    
    snow = nb.get_gridded_vars(['Snow'],x,y)['Snow']
    rr = nb.get_gridded_vars(['Ridge_ratio'],x,y)['Ridge_ratio']
    
    #OSI-SAF data
    tmp = year+'04301200'
    netcdf_name = 'OSI-SAF_ice_type_cumul_'+tmp+'.nc'
    fosi_cumul = outpath+netcdf_name
    f = Dataset(fosi_cumul)
    myi_osi_cumul = f.variables['myi_freq'][:] 
    lat_osi = f.variables['lat'][:]
    lon_osi = f.variables['lon'][:]

    
    #smoothen the data for nicer contours
    myi_smooth = smooth_data(myi,lon,lat,lon_gm,lat_gm)
    myi_osi_smooth = smooth_data(myi_osi_cumul,lon_osi,lat_osi,lon_gm,lat_gm)

    plot_contour(lon_g,lat_g,data=[myi_smooth,myi_osi_smooth],levels=[.1,.1],colors=['red','purple'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],outname=outpath_plots+'it_contour_'+year+'.png')
    
    
    #plot neXtSIM ridges/snow as background
    snow_smooth = smooth_data(snow,lon,lat,lon_gm,lat_gm)
    rr_smooth = smooth_data(rr,lon,lat,lon_gm,lat_gm)
    
    plot_contour_bg(lon_g,lat_g,rr_smooth,data=[myi_smooth,myi_osi_smooth],levels=[.05,.1],colors=['red','purple'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='Ridge ratio',outname=outpath_plots+'it_contour_bg_'+year+'.png')
    
    plot_contour_bg(lon_g,lat_g,snow_smooth,data=[myi_smooth,myi_osi_smooth],levels=[.05,.1],colors=['red','purple'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='Snow',outname=outpath_plots+'it_contour_bg_snow_'+year+'.png')
    
    #plot CFSR 'fall warm intrusions extent' as background
    yr = int(year)-1
    cfsr = np.load(outpath+'cfsr_warm_freq_'+str(yr))
    cfsr_smooth = smooth_data(cfsr,lon_fm,lat_fm,lon_gm,lat_gm)
    
    plot_contour_bg(lon_g,lat_g,cfsr_smooth,data=[myi_smooth,myi_osi_smooth],levels=[.05,.1],colors=['red','purple'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='warm instrusions freq',outname=outpath_plots+'it_contour_bg_cfsr_'+year+'.png')


    #checking charasteristicns of the difference areas
    #are they more/less ridges, have more/less snow, more/less warm air intrusions?
    nph = lat_gm > 88. #NP hole
    mask = ((myi_smooth>.05) & (myi_osi_smooth<.1) & ~nph) | ((myi_smooth<.05) & (myi_osi_smooth>.1) & ~nph) 
    rr_diff = np.mean(rr_smooth[mask])
    snow_diff = np.mean(snow_smooth[mask])
    warm_diff = np.mean(cfsr_smooth[mask])
    
    ##checking
    #tmp = np.ma.array(rr_smooth,mask=~mask)
    #plot_pcolormesh(lon_g,lat_g,tmp,'test.png')
    #exit()
    
    ridge_bias.append(rr_diff)
    snow_bias.append(snow_diff)
    warm_bias.append(warm_diff)
    dates.append(date)


#save for time series    
np.save(outpath+'snow_bias',np.array(snow_bias))   
np.save(outpath+'ridge_bias',np.array(ridge_bias))
np.save(outpath+'warm_bias',np.array(warm_bias))
np.save(outpath+'dates_bias',np.array(dates))

##load
#dates = np.load(outpath+'dates_bias.npy')
#snow_bias = np.load(outpath+'snow_bias.npy')
#ridge_bias = np.load(outpath+'ridge_bias.npy')
#warm_bias = np.load(outpath+'warm_bias.npy')


#df = pd.DataFrame({ 'mean snow depth' : snow_bias,
                    #'mean ridge ratio' : ridge_bias,
                    #'mean warm e. freq.' : warm_bias}, index=dates)

#fig2 = df.plot().get_figure()
#fig2.savefig(outpath_plots+'bias_test.png')



