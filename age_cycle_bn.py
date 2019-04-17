import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pandas as pd

from pynextsim.nextsim_bin import NextsimBin
from age_func import *

#seasonal cycle of ice type and ice age (detectable and volumetric)
inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
outpath = 'data/'
outpath_plots = 'plots/'

#get all daily files
fl = sorted(glob(inpath+'field*T000000Z.bin'))

dates = []
volume_myi = []
volume_age0 = []
volume_age1 = []
volume_age2 = []
volume_age3 = []
volume_age4 = []
volume_age0v = []
volume_age1v = []
volume_age2v = []
volume_age3v = []
volume_age4v = []

for fn in fl:
    print(fn)
    
    #get date
    tmp = fn.split('_')[-1].split('T')[0]
    date = datetime.strptime(tmp, "%Y%m%d")
    year = str(date.year)
    if int(year)<2014: continue

    dates.append(date)

    #get data
    nb = NextsimBin(fn)
    ea = nb.get_var('Element_area')/1e6              #from m^2 to km^2
    fyi = nb.get_var('Fyi_fraction')
    sic = nb.get_var('Concentration')
    sicthin = nb.get_var('Concentration_thin_ice')
    sit = nb.get_var('Thickness')/1000               #from m to km
    siad = nb.get_var('Age_d')/60/60/24/365          #from seconds to years
    sia = nb.get_var('Age')/60/60/24/365
    
    #make a sum of volume over the whole model domain
    #myi volume
    #myi = sic-sicthin-fyi              #this does not work, but gives similar results
    myi = 1-fyi
    vmyi = np.sum(myi*sit*ea)/1e3
    #print(vmyi)
    #total_vol = np.sum(sit*ea)/1e3
    #print(total_vol)
    volume_myi.append(vmyi)
    
    #reclassify the age_d into myi
    #MYI: age>timedelta(today,15. Sept)
    if (date.month>9) | ((date.month==9) & (date.day>15)):
        pyr=date.year
    else:
        pyr=date.year-1
    diff = (date-datetime(pyr,9,15)).total_seconds()/60/60/24/356
    #print(diff)
    
    vmyi = np.sum(np.ma.array(sit*ea,mask=siad<diff)             )/1e3                               #sit is effective sea ice thickness - thickness*concentration
    #print(vmyi)
    volume_myi.append(vmyi)
    
    #ice age (detectable from space/surface) volume
    vage0 = np.sum(np.ma.array(sit*ea,mask=siad>diff)             )/1e3                      #FYI
    vage1 = np.sum(np.ma.array(sit*ea,mask=(siad<diff)   | (siad>diff+1)))/1e3               #SYI
    vage2 = np.sum(np.ma.array(sit*ea,mask=(siad<diff+1) | (siad>diff+2)))/1e3               #3YI
    vage3 = np.sum(np.ma.array(sit*ea,mask=(siad<diff+2) | (siad>diff+3)))/1e3               #4YI
    vage4 = np.sum(np.ma.array(sit*ea,mask=siad<diff+3)             )/1e3                    #5YI
    volume_age0.append(vage0); volume_age1.append(vage1); volume_age2.append(vage2); volume_age3.append(vage3); volume_age4.append(vage4)
    #print(vage0+vage1+vage2+vage3+vage4)
    #exit()

    #ice age volume
    vage0v = np.sum(np.ma.array(sit*ea,mask=sia>diff)             )/1e3
    vage1v = np.sum(np.ma.array(sit*ea,mask=(sia<diff)   | (sia>diff+1)))/1e3
    vage2v = np.sum(np.ma.array(sit*ea,mask=(sia<diff+1) | (sia>diff+2)))/1e3
    vage3v = np.sum(np.ma.array(sit*ea,mask=(sia<diff+2) | (sia>diff+3)))/1e3
    vage4v = np.sum(np.ma.array(sit*ea,mask=sia<diff+3)             )/1e3
    volume_age0v.append(vage0v); volume_age1v.append(vage1v); volume_age2v.append(vage2v); volume_age3v.append(vage3v); volume_age4v.append(vage4v)
    
#save data
outfile = outpath+'age_volume_ts_bn' 
np.savez(outfile, dates = np.array(dates),\
    vmyi = np.array(volume_myi),\
    vage0 = np.array(volume_age0), vage1 = np.array(volume_age1),vage2 = np.array(volume_age2),vage3 = np.array(volume_age3),vage4 = np.array(volume_age4),\
    vage0v = np.array(volume_age0v),vage1v = np.array(volume_age1v),vage2v = np.array(volume_age2v),vage3v = np.array(volume_age3v),vage4v = np.array(volume_age4v) )

#load data
container = np.load(outpath+'age_volume_ts_bn.npz')
dates = container['dates']
volume_myi = container['vmyi'][::2]
volume_myi1 = container['vmyi'][1::2]
volume_age0 = container['vage0'];volume_age1 = container['vage1'];volume_age2 = container['vage2'];volume_age3 = container['vage3'];volume_age4 = container['vage4']
volume_age0v = container['vage0v'];volume_age1v = container['vage1v'];volume_age2v = container['vage2v'];volume_age3v = container['vage3v'];volume_age4v = container['vage4v']

#make pandas time series
df = pd.DataFrame({ 'MYI type' : volume_myi,
                    'MYI type from age' : volume_myi1,
                    '5YI+' : volume_age4,
                    '4YI' : volume_age3,
                    '3YI' : volume_age2,
                    'SYI' : volume_age1,
                    'FYI' : volume_age0,
                    '5YI+ (vol)' : volume_age4v,
                    '4YI (vol)' : volume_age3v,
                    '3YI (vol)' : volume_age2v,
                    'SYI (vol)' : volume_age1v,
                    'FYI (vol)' : volume_age0v,                  
                    }, index=dates)

#plot
fig, axes = plt.subplots(nrows=3, ncols=1,figsize=(10,8))

#type
ax = df.iloc[:,:2].plot(ax=axes[0],title='ice volume seasonal cycle',ylim=[0,20])
ax.set_ylabel(r'volume (10$^3$km$^3$)')

#surface age
bx = df.iloc[:,2:7].plot.area(ax=axes[1],ylim=[0,20])
bx.set_ylabel(r'volume (10$^3$km$^3$)')

#volume age
cx = df.iloc[:,7:].plot.area(ax=axes[2],ylim=[0,20])
cx.set_ylabel(r'volume (10$^3$km$^3$)')



fig.savefig(outpath_plots+'cycle_bn.png',bbox_inches='tight')



#sample reader for Moorings.nc
##read ice type, siad and age & thickness
#from netCDF4 import Dataset
#f = Dataset(fn)
#fyi = f.variables['fyi_fraction'][:]
#siad = f.variables['sia_det'][:]/60/60/24/365    #from seconds to years
#age = f.variables['sia'][:]/60/60/24/365
#sit = f.variables['sit'][:]/1000                #from m to km
#sic = f.variables['sic'][:]

##store date for timeseries
#time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
#base = datetime(1900,1,1)
#dates = np.array([base + timedelta(days=i) for i in time])
#print(dates)
##exit()
