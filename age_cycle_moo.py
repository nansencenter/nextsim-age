import numpy as np
from glob import glob
import regionmask
import cartopy.crs as ccrs
from netCDF4 import Dataset
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import pandas as pd

from pynextsim.nextsim_bin import NextsimBin
from age_func import *

#seasonal cycle of ice type and ice age (detectable and volumetric)
inpath = 'data/'
outpath = 'data/'
outpath_plots = 'plots/'

#get all daily files
fl = sorted(glob(inpath+'Moorings.nc.~7*'))

##mask out the region with thick/immobile/old ice at the end of the run
fn = fl[-1]
f = Dataset(fn)
lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
age_mask = get_poly_mask(lons,lats)

#lists to collect the time series
dates = []
area_myi = []
area_myi1 = []
volume_myi = []
volume_myi1 = []
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
    
    #get data
    f = Dataset(fn)
    fyi = f.variables['fyi_fraction'][:]
    siad = f.variables['sia_det'][:]/60/60/24/365    #from seconds to years
    sia = f.variables['sia'][:]/60/60/24/365
    sit = f.variables['sit'][:]/1000                #from m to km
    sic = f.variables['sic'][:]
    
    #store date for timeseries
    time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
    base = datetime(1900,1,1)
    dt = np.array([base + timedelta(days=i) for i in time])
    print(dt[0])
    print(dt[-1])
    dates.extend(dt)

    #extend the age_mask
    #prepare time difference array
    #summer mask
    age_mask_tm = np.zeros_like(sit)
    diff = np.zeros_like(sit)
    #summer = np.zeros_like(sit)
    for i in range(0,sit.shape[0]):
        age_mask_tm[i,:,:]=age_mask
        
        #for reclassifying the age_d into myi
        #MYI: age>timedelta(today,15. Sept)
        if (dt[i].month>9) | ((dt[i].month==9) & (dt[i].day>15)):
            pyr=dt[i].year
        else:
            pyr=dt[i].year-1
        diff[i,:,:] = (dt[i]-datetime(pyr,9,15)).total_seconds()/60/60/24/356
        
        ##summer
        #if (dt[i].month>5)  & (dt[i].month<10):
            #summer[i,:,:]=1
        
    #mask sea ice concentration/thickness   
    sic = sic*age_mask_tm
    sit = sit*age_mask_tm
    
    #MYI
    thin_ice_mask = (sit<.0005)# & (summer==0)  #mask thin ice (not accounted in FYI fraction), <.5m ---only in the growth season? (otherwise there will be a dip in summer)
    myi = (1-(sic*fyi))*age_mask_tm
    myi = np.ma.array(myi,mask=thin_ice_mask)
    #plot for checking
    #map_name = outpath_plots+'myi_test.png'
    #plot_pcolormesh(lons,lats,myi[90,:,:],map_name,vmin=0,vmax=1,cmap='jet',label='MYI')
    #print(dt[120])

    #make a sum of area/volume over the whole model domain
    amyi = np.sum(np.sum(myi*100,axis=1),axis=1)/1e6
    vmyi = np.sum(np.sum(myi*sit*100,axis=1),axis=1)/1e3
    print(vmyi[0])
    area_myi.extend(amyi)
    volume_myi.extend(vmyi)
    
    #reclassify the age_d into myi
    #myi=np.ma.array(sic,mask=siad<diff)
    #map_name = outpath_plots+'myi_test_age.png'
    #plot_pcolormesh(lons,lats,myi[90,:,:],map_name,vmin=0,vmax=1,cmap='jet',label='MYI')
    #print(dt[120])
    #exit()
    
    amyi = np.sum(np.sum(np.ma.array(sic*100,mask=siad<diff),axis=1),axis=1)/1e6
    vmyi = np.sum(np.sum(np.ma.array(sit*100,mask=siad<diff),axis=1),axis=1)/1e3                               #sit is effective sea ice thickness - thickness*concentration
    print(vmyi[0])
    area_myi1.extend(amyi)
    volume_myi1.extend(vmyi)
    
    #ice age (detectable from space/surface) volume
    vage0 = np.sum(np.sum(np.ma.array(sit*100,mask=siad>diff),axis=1),axis=1)             /1e3                      #FYI
    vage1 = np.sum(np.sum(np.ma.array(sit*100,mask=(siad<diff)   | (siad>diff+1)),axis=1),axis=1)/1e3               #SYI
    vage2 = np.sum(np.sum(np.ma.array(sit*100,mask=(siad<diff+1) | (siad>diff+2)),axis=1),axis=1)/1e3               #3YI
    vage3 = np.sum(np.sum(np.ma.array(sit*100,mask=(siad<diff+2) | (siad>diff+3)),axis=1),axis=1)/1e3               #4YI
    vage4 = np.sum(np.sum(np.ma.array(sit*100,mask=siad<diff+3),axis=1),axis=1)           /1e3                      #5YI
    volume_age0.extend(vage0); volume_age1.extend(vage1); volume_age2.extend(vage2); volume_age3.extend(vage3); volume_age4.extend(vage4)
    #print(vage0+vage1+vage2+vage3+vage4)

    #ice age volume
    vage0v = np.sum(np.sum(np.ma.array(sit*100,mask=sia>diff),axis=1),axis=1             )/1e3
    vage1v = np.sum(np.sum(np.ma.array(sit*100,mask=(sia<diff)   | (sia>diff+1)),axis=1),axis=1)/1e3
    vage2v = np.sum(np.sum(np.ma.array(sit*100,mask=(sia<diff+1) | (sia>diff+2)),axis=1),axis=1)/1e3
    vage3v = np.sum(np.sum(np.ma.array(sit*100,mask=(sia<diff+2) | (sia>diff+3)),axis=1),axis=1)/1e3
    vage4v = np.sum(np.sum(np.ma.array(sit*100,mask=sia<diff+3),axis=1),axis=1             )/1e3
    volume_age0v.extend(vage0v); volume_age1v.extend(vage1v); volume_age2v.extend(vage2v); volume_age3v.extend(vage3v); volume_age4v.extend(vage4v)
    
#save data
outfile = outpath+'age_volume_ts' 
np.savez(outfile, dates = np.array(dates),\
    amyi = np.array(area_myi), amyi1 = np.array(area_myi1),\
    vmyi = np.array(volume_myi), vmyi1 = np.array(volume_myi1),\
    vage0 = np.array(volume_age0), vage1 = np.array(volume_age1),vage2 = np.array(volume_age2),vage3 = np.array(volume_age3),vage4 = np.array(volume_age4),\
    vage0v = np.array(volume_age0v),vage1v = np.array(volume_age1v),vage2v = np.array(volume_age2v),vage3v = np.array(volume_age3v),vage4v = np.array(volume_age4v) )

#load data
container = np.load(outpath+'age_volume_ts.npz')
dates = container['dates']
area_myi = container['amyi']
area_myi1 = container['amyi1']
volume_myi = container['vmyi']
volume_myi1 = container['vmyi1']
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
                    'MYI type area' : area_myi,
                    'MYI type area from age' : area_myi1,
                    }, index=dates)

#plot
fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(10,8))

#type
ax = df.iloc[:,12:].plot(ax=axes[0],title='Ice area/volume seasonal cycle')
ax.set_ylabel(r'area (10$^6$km$^2$)')

ax = df.iloc[:,:2].plot(ax=axes[1],ylim=[0,20])
ax.set_ylabel(r'volume (10$^3$km$^3$)')

#surface age
bx = df.iloc[:,2:7].plot.area(ax=axes[2],ylim=[0,20])
bx.set_ylabel(r'volume (10$^3$km$^3$)')

#volume age
cx = df.iloc[:,7:12].plot.area(ax=axes[3],ylim=[0,20])
cx.set_ylabel(r'volume (10$^3$km$^3$)')

fig.savefig(outpath_plots+'cycle.png',bbox_inches='tight')
