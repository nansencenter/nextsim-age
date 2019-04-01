import numpy as np
from glob import glob
from datetime import datetime
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import pandas as pd
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

from age_func import *

#make April average maps
inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
outpath = 'data/outputs/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
tyosi_path = '/input_obs_data/data/OSISAF_ice_type/'
cfsr_path = '/input_obs_data/data/CFSR/'
outpath_plots = 'plots/new/'

#get model grid from moorings for age mask
fn = 'data/Moorings.nc.~7~'
f = Dataset(fn)
lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
age_mask = get_poly_mask(lons,lats)

#make a nice cyclic equally spaced lat-lon grid for plotting contoured fields with cartopy
lat_g = np.arange(90,62,-.1)
lon_g = np.arange(-180,181,1.)   
lon_gm, lat_gm = np.meshgrid(lon_g,lat_g)

#read CFSR grid
f = Dataset(cfsr_path+'cfsr.6h.197901.nc')
lat_f = f.variables['lat'][:]
lon_f = f.variables['lon'][:]
lon_fm, lat_fm = np.meshgrid(lon_f,lat_f)

#OSI-SAF grid
fn = outpath+'OSI-SAF_ice_type_cumul_200504301200.nc'
f = Dataset(fn)
myi_osi_cumul = f.variables['myi_freq'][:] 
lat_osi = f.variables['lat'][:]
lon_osi = f.variables['lon'][:]

#get all daily files
fl = sorted(glob(inpath+'field*0415T000000Z.bin'))

myi_type_list = []
myi_age_list = []
myi_osi_list = []

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
    age = nb.get_gridded_vars(['Age_d'],x,y)['Age_d']
    snow = nb.get_gridded_vars(['Snow'],x,y)['Snow']
    rr = nb.get_gridded_vars(['Ridge_ratio'],x,y)['Ridge_ratio']
    
    #MYI from type
    myi = ic-icthin-fyi
    
    #MYI from 'surface age' (detectable from space)
    age = nb.get_gridded_vars(['Age_d'],x,y)['Age_d']
    if (date.month>9) | ((date.month==9) & (date.day>15)):
        pyr=date.year
    else:
        pyr=date.year-1
    diff = (date-datetime(pyr,9,15)).total_seconds()/60/60/24/356
    print(diff)
    age = age/60/60/24/356
    #fill nans with zeros
    age = np.nan_to_num(age)
    myi_age = np.ma.array(age,mask=age<diff)
    ##Plot myi_age
    #outname = outpath_plots+'age_test.png'
    #plot_pcolormesh(lon,lat,myi_age,outname,vmin=0,vmax=1,cmap='viridis',label='MYI=1')
    ##exit()
      
    #OSI-SAF data
    tmp = year+'04301200'
    netcdf_name = 'OSI-SAF_ice_type_cumul_'+tmp+'.nc'
    fosi_cumul = outpath+netcdf_name
    f = Dataset(fosi_cumul)
    myi_osi_cumul = f.variables['myi_freq'][:]
    myi_osi = np.where(myi_osi_cumul>0,1,0)
    
    #smoothen and regrid the data for nicer contours
    tmp = smooth_data(myi,lon,lat,lon_osi,lat_osi)
    myi_smooth = smooth_data(tmp,lon_osi,lat_osi,lon_gm,lat_gm)
    
    myi_age_smooth = regrid_data(myi_age,lon,lat,lon_gm,lat_gm)
    myi_age_smooth = np.where(myi_age_smooth>.5,1,-1)
    ##Plot myi_age
    #outname = outpath_plots+'age_test_smooth.png'
    #plot_pcolormesh(lon_gm,lat_gm,myi_age_smooth,outname,vmin=0,vmax=1,cmap='viridis',label='MYI=1')
    ##exit()
    
    myi_osi_smooth = regrid_data(myi_osi,lon_osi,lat_osi,lon_gm,lat_gm)
    
    #Plot all age and type data
    plot_contour_bg(lon_g,lat_g,myi_age_smooth,data=[myi_smooth,myi_osi_smooth],levels=[.1,.1],colors=['red','k'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='MYI age',outname=outpath_plots+'it_contour_'+year+'.png')    
    
    
    #plot neXtSIM ridges/snow as background
    snow_smooth = smooth_data(snow,lon,lat,lon_gm,lat_gm)
    rr_smooth = smooth_data(rr,lon,lat,lon_gm,lat_gm)
    
    plot_contour_bg(lon_g,lat_g,rr_smooth,data=[myi_age_smooth,myi_osi_smooth],levels=[.5,.1],colors=['purple','k'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='Ridge ratio',outname=outpath_plots+'it_contour_bg_'+year+'.png')
    
    plot_contour_bg(lon_g,lat_g,snow_smooth,data=[myi_smooth,myi_osi_smooth],levels=[.5,.1],colors=['purple','purple'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='Snow',outname=outpath_plots+'it_contour_bg_snow_'+year+'.png')
    
    #plot CFSR 'fall warm intrusions extent' as background
    yr = int(year)-1
    cfsr = np.load(outpath+'cfsr_warm_freq_'+str(yr))
    cfsr_smooth = smooth_data(cfsr,lon_fm,lat_fm,lon_gm,lat_gm)
    
    plot_contour_bg(lon_g,lat_g,cfsr_smooth,data=[myi_age_smooth,myi_osi_smooth],levels=[.05,.1],colors=['red','k'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='warm instrusions freq',outname=outpath_plots+'it_contour_bg_cfsr_'+year+'.png')


    #mask the areas by age_mask and save for the time series
    #regrid all data to a regulary spaced grid
    age_mask_og = regrid_data(age_mask,lons,lats,lon_osi,lat_osi)
    nph = lat_osi > 88. #NP hole
    mask = (age_mask_og<1) | nph
    
    tmp = np.where(myi_smooth>.1,1,0)
    tmp = regrid_data(tmp,lon_gm,lat_gm,lon_osi,lat_osi)
    myi_type_area = np.sum(np.ma.array(tmp, mask=mask)*62.5*62.5)
    
    tmp = regrid_data(myi_age,lon,lat,lon_osi,lat_osi)
    myi_age_area = np.sum(np.ma.array(tmp, mask=mask)*62.5*62.5)
    
    myi_osi_area = np.sum(np.ma.array(myi_osi_cumul, mask=mask)*62.5*62.5)
    
    #append to lists
    myi_type_list.append(myi_type_area)
    myi_age_list.append(myi_age_area)
    myi_osi_list.append(myi_osi_area)
    
    
    
    

    #checking charasteristicns of the difference areas
    #are there more/less ridges, have more/less snow, more/less warm air intrusions?
    nph = lat_gm > 88. #NP hole
    mask = ((myi_age_smooth>.5) & (myi_osi_smooth<.5) & ~nph) | ((myi_age_smooth<.5) & (myi_osi_smooth>.5) & ~nph) 
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
outfile = outpath+'age_ts' 
np.savez(outfile, dates = np.array(dates),\
    sb = np.array(snow_bias), rb = np.array(ridge_bias), wb = np.array(warm_bias), \
    myit = np.array(myi_type_list), myia = np.array(myi_age_list), myio =np.array(myi_osi_list) )

#load
container = np.load(outpath+'age_ts.npz')
dates = container['dates']
snow_bias = container['sb']
ridge_bias = container['rb']
warm_bias = container['wb']

df = pd.DataFrame({ 'mean snow depth' : snow_bias,
                    'mean ridge ratio' : ridge_bias,
                    'mean warm e. freq.' : warm_bias}, index=dates)

fig2 = df.plot().get_figure()
fig2.savefig(outpath_plots+'bias_test.png')



