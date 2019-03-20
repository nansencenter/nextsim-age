import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample
import pandas as pd

from age_func import *


#daily export of sea ice thickness, ice type and ice age (detectable and volumetric) through a gate in the Fram Strait
inpath = '/input_obs_data/data/FRASIL/'
inpath = 'data/'
outpath_plots = 'plots/'

#get all mooring files
fl = sorted(glob(inpath+'Moorings.nc'))


for fn in fl:
    print(fn)
    #read ice type, age_d and age & thickness
    f = Dataset(fn)
    fyi = f.variables['fyi_fraction'][:]
    age_d = f.variables['sia_det'][:]/60/60/24/365    #from seconds to years
    age = f.variables['sia'][:]/60/60/24/365
    sit = f.variables['sit'][:]
    sic = f.variables['sic'][:]
    u = f.variables['siu'][:]
    v = f.variables['siv'][:]
    
    #print(u.shape)
    #exit()
    
    #store date for timeseries
    time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
    base = datetime(1900,1,1)
    dates = np.array([base + timedelta(days=i) for i in time])
    print(dates)

    #define gate in the Fram Strait (consider adding Nares Strait and other Canadian Archipelago straights later to close the budget)
    #base definition on the neXtSIM grid alignment! i=115, j=[270:340]
    #only u
    
    ##plot gate!
    #lons = f.variables['longitude'][:]
    #lats = f.variables['latitude'][:]
    #print(tmp.shape)
    #print(lons.shape)
    #plot_pcolormesh(lons[270:340,90:95],lats[270:340,90:95],np.zeros_like(lons[270:340,90:95]),outpath_plots+'gate.png')
    #exit()

    #all ice
    tmp = u[:,270:340,90]*sit[:,270:340,90]*10000*sic[:,270:340,90]
    export = np.sum(tmp,axis=1)/1e9*60*60*24
    print(export)
    print(export.shape)
    
    #myi ice
    myi = sic-fyi                                   #problem!: there is no thin ice concentration in Moorings!!!    #workaround: mask out very thin ice < 0.5m
    myi = np.ma.array(myi,mask=sit<.0005)
    
    tmp = u[:,270:340,90]*sit[:,270:340,90]*10000*myi[:,270:340,90]
    export_myi = np.sum(tmp,axis=1)/1e9*60*60*24
    
   
    
    
    
    
#make pandas time series
df = pd.DataFrame({ 'export' : export,   
                    'export MYI' : export_myi,
                    }, index=dates)

dfmon = df.resample('M').mean()


#plot
fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(10,8))

#type
ax = df.iloc[:,0].plot(ax=axes[0],title='export through the Fram Strait',ylim=[-15,15])
ax = dfmon.iloc[:,0].plot(ax=axes[0])
ax.plot(dates,np.zeros_like(export))
ax.set_ylabel(r'volume (km$^3$/day)')


bx = df.iloc[:,1].plot(ax=axes[1],ylim=[-15,15])
bx = dfmon.iloc[:,1].plot(ax=axes[1])
bx.plot(dates,np.zeros_like(export))
bx.set_ylabel(r'volume (km$^3$/day)')


fig.savefig(outpath_plots+'export.png',bbox_inches='tight')


    
#compare to values in literature (simple check if the values are reasonable)
    
