import numpy as np
from glob import glob
from datetime import datetime
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

from age_func import *

#make April 1 maps
inpath = 'data/run04/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
tyosi_path = '/input_obs_data/data/OSISAF_ice_type/'
cfsr_path = '/input_obs_data/data/CFSR/'
cfsr_path1 = 'data/outputs/'
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
    icthin = nb.get_gridded_vars(['Concentration_thin_ice'],x,y)['Concentration_thin_ice']
    fyi = nb.get_gridded_vars(['Fyi_fraction'],x,y)['Fyi_fraction']
    myi = ic-icthin-fyi
    
    snow = nb.get_gridded_vars(['Snow'],x,y)['Snow']
    rr = nb.get_gridded_vars(['Ridge_ratio'],x,y)['Ridge_ratio']
    
    #OSI-SAF data
    tmp = year+'04301200'
    outpath = 'data/outputs/'
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
    cfsr = np.load(cfsr_path1+'cfsr_warm_freq_'+str(yr))
    cfsr_smooth = smooth_data(cfsr,lon_fm,lat_fm,lon_gm,lat_gm)
    
    plot_contour_bg(lon_g,lat_g,cfsr_smooth,data=[myi_smooth,myi_osi_smooth],levels=[.05,.1],colors=['red','purple'], lw=[3,3], \
                 labels=['neXtSIM MYI extent','OSI-SAF MYI extent'],bg_label='warm instrusions freq',outname=outpath_plots+'it_contour_bg_cfsr_'+year+'.png')



    ##maps of snow and ridge ratio
    #figname = outpath_plots+'snow_'+year+'.png'
    #nb.plot_var('Snow',figname=figname,clim=[0,1])
    ##for time series of mean snow depth in the areas not classified as MYI in OSI-SAF
    #snow = nb.get_var('Snow')
    #diff_mask = np.where(it_diff_cumul>0,1,0)
    #snow_bias.append(np.mean(snow*diff_mask))
    #dates.append(date)
    #print(np.mean(snow*diff_mask))
    
    ##also save the difference area!!!!
    #ea = nb.get_var('Element_area')
    #area = np.sum(diff_mask*ea)/1e9/1e3
    #area_bias.append(area)
    
    #figname = outpath_plots+'ridges_'+year+'.png'
    #nb.plot_var('Ridge_ratio',figname=figname,clim=[0,1])
    ##for time series of mean ridge ratio in the areas not classified as MYI in OSI-SAF
    #rr = nb.get_var('Ridge_ratio')
    #ridge_bias.append(np.mean(rr*diff_mask))
    #print(np.mean(rr*diff_mask))
    
    ##maps of ice age
    #age = nb.get_var('Age_d')
    #aoy = age/60/60/24/365
    #nb.add_var('Age_d_year', f, aoy)
    #figname = outpath_plots+'age_'+year+'.png'
    #nb.plot_var('Age_d_year',figname=figname,clim=[0,6])
    
    #age_bin = np.where(aoy>1.33,1,0)
    #age_diff = age_bin - myi_osi_cumul
    #age_diff = np.where(landmask,0,age_diff)
    #nb.add_var('Age_diff', f, age_diff)
    #figname = outpath_plots+'age_diff'+year+'.png'
    #nb.plot_var('Age_diff',cmap='bwr',figname=figname,clim=[-1,1])

##save for time series    
#np.save(outpath+'snow_bias',np.array(snow_bias))   
#np.save(outpath+'ridge_bias',np.array(ridge_bias))
#np.save(outpath+'area_bias',np.array(area_bias))
#np.save(outpath+'dates_bias',np.array(dates))

##load
#dates = np.load(outpath+'dates_bias.npy')
#snow_bias = np.load(outpath+'snow_bias.npy')
#ridge_bias = np.load(outpath+'ridge_bias.npy')
#area_bias = np.load(outpath+'area_bias.npy')


#df = pd.DataFrame({ 'mean snow depth' : snow_bias,
                    #'mean ridge ratio' : ridge_bias,
                    #'diff area size' : area_bias/10}, index=dates)

#fig2 = df.plot().get_figure()
#fig2.savefig(outpath_plots+'bias_test.png')



