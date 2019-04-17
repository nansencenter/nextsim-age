import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample
import pandas as pd

from age_func import *


#daily export of sea ice thickness, ice type and ice age (detectable and volumetric) through a gate in the Fram Strait
inpath = 'data/'
outpath = 'data/'
outpath_plots = 'plots/new/'

##get all mooring files
#fl = sorted(glob(inpath+'Moorings.nc.~*'))

#dts = []
#ex = []
#exa = []
#ex_myi = []
#exa_myi = []

#for fn in fl:
    #print(fn)
    ##read ice type, age_d and age & thickness
    #f = Dataset(fn)
    #fyi = f.variables['fyi_fraction'][:]
    #siad = f.variables['sia_det'][:]/60/60/24/365    #from seconds to years
    #sit = f.variables['sit'][:]
    #sic = f.variables['sic'][:]
    #u = f.variables['siu'][:]
    #v = f.variables['siv'][:]
    
    ##print(u.shape)
    ##exit()
    
    ##store date for timeseries
    #time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
    #base = datetime(1900,1,1)
    #dt = np.array([base + timedelta(days=i) for i in time])
    ##print(dt)

    ##define gate in the Fram Strait (consider adding Nares Strait and other Canadian Archipelago straights later to close the budget)
    ##base definition on the neXtSIM grid alignment! i=115, j=[270:340]
    ##only u
    
    ###plot gate!
    ##lons = f.variables['longitude'][:]
    ##lats = f.variables['latitude'][:]
    ##print(tmp.shape)
    ##print(lons.shape)
    ##plot_pcolormesh(lons[270:340,90:95],lats[270:340,90:95],np.zeros_like(lons[270:340,90:95]),outpath_plots+'gate.png')
    ##exit()

    ##all ice
    #tmp = u[:,270:340,90]*sit[:,270:340,90]*10000
    #export = np.sum(tmp,axis=1)/1e9*60*60*24
    
    #tmp = u[:,270:340,90]*10000*sic[:,270:340,90]
    #export_area = np.sum(tmp,axis=1)/1e6*60*60*24
    ##print(export)
    ##print(export.shape)
    
    ###myi ice
    ##myi = sic-fyi                                   #problem!: there is no thin ice concentration in Moorings!!!    #workaround: mask out very thin ice < 0.5m
    ##myi = np.ma.array(myi,mask=sit<.0005)
    
    ##MYI from 'surface age' (detectable from space)
    #diff = np.zeros_like(sit)
    #for i in range(0,sit.shape[0]):
        #if (dt[i].month>9) | ((dt[i].month==9) & (dt[i].day>15)):
            #pyr=dt[i].year
        #else:
            #pyr=dt[i].year-1
        #diff[i,:,:] = (dt[i]-datetime(pyr,9,15)).total_seconds()/60/60/24/356
    
    #myi_age = np.where(siad>diff,1,0)
    
    #tmp = u[:,270:340,90]*sit[:,270:340,90]*10000*myi_age[:,270:340,90]
    #export_myi = np.sum(tmp,axis=1)/1e9*60*60*24
    
    #tmp = u[:,270:340,90]*10000*sic[:,270:340,90]*myi_age[:,270:340,90]
    #export_myi_area = np.sum(tmp,axis=1)/1e6*60*60*24    
    
    #ex.extend(export)
    #exa.extend(export_area)
    #ex_myi.extend(export_myi)
    #exa_myi.extend(export_myi_area)
    #dts.extend(dt)
    
##save data
#outfile = outpath+'export_ts' 
#np.savez(outfile, dates = np.array(dts),\
    #export = np.array(ex), export_myi = np.array(ex_myi) ,\
    #export_area = np.array(exa), export_myi_area = np.array(exa_myi) )

#load data
container = np.load(outpath+'export_ts.npz')
dates = container['dates']
export = container['export']   
export_area = container['export_area']
export_myi = container['export_myi']  
export_myi_area = container['export_myi_area']  
    
#make pandas time series
df = pd.DataFrame({ 'export' : export,  
                    'export MYI' : export_myi,
                    'export (area)' : export_area,
                    'export MYI (area)' : export_myi_area,
                    }, index=dates)

dfmon = df.resample('M').mean()

dfyr = df.resample('Y').mean()
dfyr_std = df.resample('Y').std()

#plot
fig, axes = plt.subplots(nrows=3, ncols=1,figsize=(8,8))

#type
ax = df.iloc[:,:2].plot(ax=axes[0],title='export through the Fram Strait',ylim=[-15,15])
ax = dfmon.iloc[:,:2].plot(ax=axes[0])
ax.plot(dates,np.zeros_like(export),c='k')
ax.set_ylabel(r'volume (km$^3$/day)')

bx = df.iloc[:,2:].plot(ax=axes[1])
bx = dfmon.iloc[:,2:].plot(ax=axes[1])
bx.plot(dates,np.zeros_like(export),c='k')
bx.set_ylabel(r'area (km$^2$/day)')

cx = dfyr.iloc[:,:2].plot(ax=axes[2],yerr=dfyr_std)
cx.plot(dates,np.zeros_like(export),c='k')
cx.set_ylabel(r'volume (km$^3$/day)')


fig.savefig(outpath_plots+'export.png',bbox_inches='tight')


    
#compare to values in literature (simple check if the values are reasonable)
    
