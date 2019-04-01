import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample

from age_func import *

#cdo mergetime -selmon,1 -selday,1,2,3 Moorings4.nc -selmon,1 -selday,1,2,3 Moorings5.nc -selmon,1 -selday,1,2,3 Moorings6.nc ~/data/merged_Jan.nc
#scp ~/data/merged_Jan.nc johansen:/Data/sim/data/FRASIL/run01_moorings/
#scp ~/data/run04_Jan*.nc johansen:/Data/sim/data/FRASIL/



#cdo mergetime -selyear,2007 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~4~   ~/src/nextsim_age/data/run04_Jan2007.nc

#cdo mergetime -selyear,2008 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~5~   ~/src/nextsim_age/data/run04_Jan2008.nc
#cdo mergetime -selyear,2009 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~5~   ~/src/nextsim_age/data/run04_Jan2009.nc
#cdo mergetime -selyear,2010 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~5~   ~/src/nextsim_age/data/run04_Jan2010.nc

#cdo mergetime -selyear,2011 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~6~   ~/src/nextsim_age/data/run04_Jan2011.nc
#cdo mergetime -selyear,2012 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~6~   ~/src/nextsim_age/data/run04_Jan2012.nc
#cdo mergetime -selyear,2013 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~6~   ~/src/nextsim_age/data/run04_Jan2013.nc

#cdo mergetime -selyear,2014 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~7~   ~/src/nextsim_age/data/run04_Jan2014.nc
#cdo mergetime -selyear,2015 -selmon,1 ~/src/nextsim_age/data/Moorings.nc.~7~   ~/src/nextsim_age/data/run04_Jan2015.nc





#make PDF of Jan 1-2 average velocities and compare to OSI-SAF (2009-2016)
inpath = '/input_obs_data/data/FRASIL/'
inpath = 'data/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
drosi_path = '/input_obs_data/data/OSISAF_ice_drift/'
outpath_plots = 'plots/'

years = range(2007,2016)
#years = range(2007,2008)

slist = []
slist_gauss = []
slist_uniform = []
slist_full = []
slist_coarse = []

slist_osi = []

#get model data
fn = inpath+'run04_Jan2007.nc'
f = Dataset(fn)
lon_m = f.variables['longitude'][:]
lat_m = f.variables['latitude'][:]

age_mask_m = get_poly_mask(lon_m,lat_m)

#get the OSI-SAF grid
fl = sorted(glob(drosi_path+'2008/01/*01031200.nc'))
f = Dataset(fl[0])
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]

##check spacing
#yc = f.variables['yc'][:]
#xc = f.variables['xc'][:]
#print(yc)
#print(xc)
#exit()

cum_speed = np.zeros_like(lat)
cum_speed_osi = np.zeros_like(lat)

#get all the infromation for the re-gridding
orig_def = pyresample.geometry.SwathDefinition(lons=lon_m, lats=lat_m)
targ_def = pyresample.geometry.SwathDefinition(lons=lon, lats=lat)
coar_def = pyresample.geometry.SwathDefinition(lons=lon[::3,::3], lats=lat[::3,::3])  #take every 3rd grid cell to get ~140km large grid cells


age_mask = pyresample.kd_tree.resample_nearest(orig_def, age_mask_m, \
            targ_def, radius_of_influence=10000, fill_value=None)

o=0
for yr in years:
    fn = inpath+'run04_Jan'+str(yr)+'.nc'
    print(fn)
    
    fl = sorted(glob(drosi_path+str(yr)+'/01/*.nc'))
    
    for i in range(0,29):       #just for the first 30 days (then no averages possible for 2 days in Jan)
        f = Dataset(fn)
        u = f.variables['siu'][:]
        v = f.variables['siv'][:]
        sic = f.variables['sic'][:]

        #make averages over 2 days
        #this is not perfect as the image is taken not at midnight, but sometime during the day
        meanu = (u[i,:,:] + u[i+1,:,:])/2
        meanv = (v[i,:,:] + v[i+1,:,:])/2
        meanic = (sic[i,:,:] + sic[i+1,:,:])/2
        speed = np.sqrt(meanu**2 +meanv**2)
        
        speed_nearest = pyresample.kd_tree.resample_nearest(orig_def, speed, \
            targ_def, radius_of_influence=30000, fill_value=None)               #30km
        
        sic_nearest = pyresample.kd_tree.resample_nearest(orig_def, meanic, \
            targ_def, radius_of_influence=30000, fill_value=None)
        
        #gauss reasmpling
        speed_gauss = pyresample.kd_tree.resample_gauss(orig_def, speed, \
            targ_def, radius_of_influence=30000, sigmas=15000)
        
        sic_gauss = pyresample.kd_tree.resample_gauss(orig_def, meanic, \
            targ_def, radius_of_influence=30000, sigmas=15000)
        
        #equal weights resampling
        wf = lambda r: 1
        speed_uniform = pyresample.kd_tree.resample_custom(orig_def, speed, \
            targ_def, radius_of_influence=30000, weight_funcs=wf)
        
        sic_uniform = pyresample.kd_tree.resample_custom(orig_def, meanic, \
            targ_def, radius_of_influence=30000, weight_funcs=wf)
        
        #sub-image area and back
        tmp = pyresample.kd_tree.resample_custom(orig_def, speed, \
            coar_def, radius_of_influence=30000, weight_funcs=wf)
        speed_coarse = pyresample.kd_tree.resample_custom(coar_def, tmp, \
            targ_def, radius_of_influence=140000, weight_funcs=wf)
        
        #OSI-SAF data
        print(fl[i])
        f = Dataset(fl[i])
        dX = f.variables['dX'][0,:,:]*1000
        dY = f.variables['dY'][0,:,:]*1000
        sf = f.variables['status_flag'][0,:,:]
        
        speed_osi = np.sqrt(dX**2 + dY**2)/2/24/60/60
        
        
        
        
        
        #match the area and dump the data
        #no MIZ and NP hole/landmask (use sf<30 for just best OSI-SAF quality)
        mask = (sic_gauss<.15) | (sf<30) | (age_mask==0)
        slist_gauss.extend(speed_gauss[~mask].flatten())
        
        mask_uni = (sic_uniform<.15) | (sf<30) | (age_mask==0)
        slist_uniform.extend(speed_uniform[~mask_uni].flatten())
        
        mask_uni = (sic_uniform<.15) | (sf<30) | (age_mask==0)
        slist_coarse.extend(speed_coarse[~mask_uni].flatten())
        
        #mask = (sic_nearest<.15) | (sf<30)
        #slist.extend(speed_nearest[~mask].flatten())
        
        slist_osi.extend(speed_osi[~mask_uni].flatten())

        sf_m = pyresample.kd_tree.resample_nearest(targ_def, sf, \
            orig_def, radius_of_influence=30000, fill_value=None)
        mask = (meanic<.15) | (sf_m<20) | (age_mask_m==0)
        slist_full.extend(speed[~mask].flatten())
        
        
        
        
        
        #for January means
        cum_speed = cum_speed+np.ma.array(speed_uniform,mask=mask_uni)
        #cum_speed = cum_speed+np.ma.array(speed_coarse,mask=mask_uni)
        cum_speed_osi = cum_speed_osi+np.ma.array(speed_osi,mask=mask_uni)
        o = o+1

#plot a PDF
bl = np.arange(0,.41,.01)
#n, bins, patches = plt.hist(slist, bl, normed=True, histtype='step', color='m', alpha=.8, label='neXtSIM', lw = 3)
#n, bins, patches = plt.hist(slist_gauss, bl, normed=True, histtype='step', color='r', alpha=.8, label='neXtSIM', lw = 3)
n, bins, patches = plt.hist(np.clip(slist_full, bl[0], bl[-1]), bl, normed=True, histtype='step', color='darkred', alpha=.8, label='neXtSIM fine', lw = 2)
n, bins, patches = plt.hist(np.clip(slist_uniform, bl[0], bl[-1]), bl, normed=True, histtype='step', color='purple', alpha=.8, label='neXtSIM mean', lw = 3)
n, bins, patches = plt.hist(np.clip(slist_coarse, bl[0], bl[-1]), bl, normed=True, histtype='step', color='m', alpha=.8, label='neXtSIM coarse', lw = 2)
n, bins, patches = plt.hist(np.clip(slist_osi, bl[0], bl[-1]), bl, normed=True, histtype='step', alpha=.8, label='OSI-SAF', lw = 3)


plt.xlabel('Speed (m/s)')
plt.ylabel('Probability')
plt.title('Probability distribution of mean 2-day speed \nfor January 2007-2015')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.legend(loc='upper right',prop={'size':16})
plt.grid(True)
plt.savefig(outpath_plots+'drift_pdf_run04_age_mask.png')



#January means
print(o)
jan_speed = cum_speed/o
jan_speed_osi = cum_speed_osi/o

#plot maps
map_name = outpath_plots+'jan_speed_map_age_mask.png'
plot_pcolormesh(lon,lat,jan_speed,map_name,vmin=0,vmax=.1,cmap='jet',label='Mean January Speed (m/s)')

map_name = outpath_plots+'jan_speed_map_osi_age_mask.png'
plot_pcolormesh(lon,lat,jan_speed_osi,map_name,vmin=0,vmax=.1,cmap='jet',label='Mean January Speed (m/s)')
