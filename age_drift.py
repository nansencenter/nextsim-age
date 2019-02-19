import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample

#cdo mergetime -selmon,1 -selday,1,2,3 Moorings4.nc -selmon,1 -selday,1,2,3 Moorings5.nc -selmon,1 -selday,1,2,3 Moorings6.nc ~/data/merged_Jan.nc
#scp ~/data/merged_Jan.nc johansen:/Data/sim/data/FRASIL/run01_moorings/


#cdo mergetime -selyear,2007 -selmon,1 Moorings.nc.~3~   ~/data/run04_Jan2007.nc
#cdo mergetime -selyear,2008 -selmon,1 Moorings.nc.~3~   ~/data/run04_Jan2008.nc

#scp ~/data/run04_Jan*.nc johansen:/Data/sim/data/FRASIL/



#make PDF of Jan 1-2 average velocities and compare to OSI-SAF (2009-2016)
inpath = '/input_obs_data/data/FRASIL/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
drosi_path = '/input_obs_data/data/OSISAF_ice_drift/'
outpath_plots = 'plots/'

years = range(2007,2009)

slist = []
slist_osi = []

#get model data
fn = inpath+'run04_Jan2007.nc'
f = Dataset(fn)
lon_m = f.variables['longitude'][:]
lat_m = f.variables['latitude'][:]

#get the OSI-SAF grid
fl = sorted(glob(drosi_path+'2008/01/*01031200.nc'))
f = Dataset(fl[0])
lat = f.variables['lat'][:]
lon = f.variables['lon'][:]

#get all the infromation for the re-gridding
orig_def = pyresample.geometry.SwathDefinition(lons=lon_m, lats=lat_m)
targ_def = pyresample.geometry.SwathDefinition(lons=lon, lats=lat)

o = 0
for yr in years:
    fn = inpath+'run04_Jan'+str(yr)+'.nc'
    print(fn)
    
    fl = sorted(glob(drosi_path+str(yr)+'/01/*.nc'))
    
    for i in range(0,29):       #just for the first 29 days (then no averages possible for 3 days in Jan)
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
            targ_def, radius_of_influence=500000, fill_value=None)
        
        sic_nearest = pyresample.kd_tree.resample_nearest(orig_def, meanic, \
            targ_def, radius_of_influence=500000, fill_value=None)
        
        #OSI-SAF data
        print(fl[i])
        f = Dataset(fl[i])
        dX = f.variables['dX'][:]*1000
        dY = f.variables['dY'][:]*1000
        sf = f.variables['status_flag'][:]
        
        speed_osi = np.sqrt(dX**2 + dY**2)/3/24/60/60
        
        #matching the area
        mask = (sic_nearest<.15) | (sf<20)       #no low SIC data and just best quality
        speed_nearest = np.ma.array(speed_nearest,mask=mask)
        speed_osi = np.ma.array(speed_osi,mask=mask)

        #dump all values into a list
        slist.extend(speed_nearest.flatten())
        slist_osi.extend(speed_osi.flatten())
    
    #OSI-SAF counter
    o=o+1

bl = np.arange(0,.6,.01)
#plot a PDF
n, bins, patches = plt.hist(slist, bl, normed=True, histtype='step', color='purple', alpha=.8, label='neXtSIM', lw = 3)
n, bins, patches = plt.hist(slist_osi, bl, normed=True, histtype='step', alpha=.8, label='OSI-SAF', lw = 3)


plt.xlabel('Speed (m/s)')
plt.ylabel('Probability')
plt.title('Probability distribution of mean 2-day speed \nfor January 2007-2008')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.legend(loc='upper right',prop={'size':16})
plt.grid(True)
plt.savefig(outpath_plots+'drift_pdf_run04.png')
