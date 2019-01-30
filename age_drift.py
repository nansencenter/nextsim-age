import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset

#cdo mergetime -selmon,1 -selday,1,2,3 Moorings4.nc -selmon,1 -selday,1,2,3 Moorings5.nc -selmon,1 -selday,1,2,3 Moorings6.nc ~/data/merged_Jan.nc
#scp ~/data/merged_Jan.nc johansen:/Data/sim/data/FRASIL/run01_moorings/

#make PDF of Jan 1-2 average velocities and compare to OSI-SAF (2009-2016)
inpath = '/input_obs_data/FRASIL/run01_moorings/'
icosi_path = '/input_obs_data/OSISAF_ice_conc/polstere/'
drosi_path = '/input_obs_data/OSISAF_ice_drift/'
outpath_plots = 'plots/'

#get model data
f = inpath+'merged_Jan.nc'
f = Dataset(f)
u = f.variables['siu'][:]
v = f.variables['siv'][:]
sic = f.variables['sic'][:]
time = f.variables['time'][:]   #days since 1900-01-01 00:00:00

slist = []
slist_osi = []

for i in range(0,time.shape[0],3):
    diff = time[i]*24
    date = datetime(1900,1,1,0,0,0)+timedelta(hours=diff)
    print(date)

    #make averages over 3 days
    meanu = (u[i,:,:] + u[i+1,:,:] + u[i+2,:,:])/3
    meanv = (v[i,:,:] + v[i+1,:,:] + v[i+2,:,:])/3
    meansic = (sic[i,:,:] + sic[i+1,:,:]+ sic[i+2,:,:])/3
    
    speed = np.sqrt(meanu**2 +meanv**2)
    
    ##mask
    mask = (meansic<0.8)|(speed<0.01)
    speed = np.ma.array(speed,mask=mask)
    #print(speed.compressed())
    #exit()

    #dump all values into a list
    slist.extend(speed.flatten())
    
#get the OSI-SAF data
fl = sorted(glob(drosi_path+'/**/01/ice_*01031200.nc', recursive=True))

for f in fl:
    print(f)
    f = Dataset(f)
    dX = f.variables['dX'][:]*1000
    dY = f.variables['dY'][:]*1000
    
    speed = np.sqrt(dX**2 + dY**2)/3/24/60/60
    
    #this needs to be masked in a similar way as the model!

    slist_osi.extend(speed.flatten())

#plot a PDF
n, bins, patches = plt.hist(slist, 50, normed=True, histtype='step', color='purple', alpha=.8, label='neXtSIM', lw = 3)
n, bins, patches = plt.hist(slist_osi, 50, normed=True, histtype='step', alpha=.8, label='OSI-SAF', lw = 3)


plt.xlabel('Speed (m/s)')
plt.ylabel('Probability')
plt.title('Probability distribution of mean 3-day speed \nfor 1-3 January 2009-2015')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.legend(loc='upper right',prop={'size':16})
plt.grid(True)
plt.savefig(outpath_plots+'drift_pdf.png')
