import numpy as np
import pandas as pd
from glob import glob
from netCDF4 import Dataset
from datetime import datetime, timedelta
from pynextsim.nextsim_bin import NextsimBin
import matplotlib.pyplot as plt

from age_func import *

##new simulations 
#cd data/run04_sept
#ln -s /input_obs_data/einar/age_datarmor/*0915T000000Z* .


inpath_ps = 'data/'
#inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
inpath = 'data/'
inpath2 ='/input_obs_data/einar/datarmor/age_95_noice/'
outpath = 'data/'
outpath_plots = 'plots/new/'

#read the PIOMAS sea ice volume
fn = inpath_ps+'PIOMAS.vol.daily.1979.2018.Current.v2.1.dat'

data = pd.read_fwf(fn, colspecs='infer', widths=[4,5,7])
year = data['Year'].values
day = data['#day'].values
vol = data['Vol'].values

dates = []
volume = []
volume_model = []
dates_model = []

volume_model2 = []
dates_model2 = []


for i in range(0,len(year)):
    date = datetime(year[i],1,1)+timedelta(days=int(day[i]))
    print(date)
    
    dates.append(date)
    volume.append(vol[i])

#model data
#mask out the region with thick/immobile/old ice at the end of the run
fn = inpath+'Moorings.nc.~7~'
f = Dataset(fn)
lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
age_mask = get_poly_mask(lons,lats)

#fl = sorted(glob(inpath+'Moorings.nc.~*'))
#print(fl)
#for fn in fl:
    #print(fn)
    #f = Dataset(fn)
    #sit = f.variables['sit'][:]
    
    ##mask with the age_mask
    #age_mask_tm = np.zeros_like(sit)
    #for i in range(0,sit.shape[0]):
        #age_mask_tm[i,:,:]=age_mask
    #sit = sit*age_mask_tm/1000               #convert from m to km
    #vol = np.sum(np.sum(sit*100,axis=1),axis=1)/1e3         #sit is effective sea ice thickness - thickness*concentration

    ##store date for timeseries
    #time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
    #base = datetime(1900,1,1)
    #dt = np.array([base + timedelta(days=i) for i in time])
    #print(dt[0])
    #print(dt[-1])

    #volume_model.extend(vol)
    #dates_model.extend(dt)
    ##print(volume_model)
    ##exit()

#new run
fl = sorted(glob(inpath2+'Moorings*.nc'))
print(fl)
for fn in fl:
    print(fn)
    f = Dataset(fn)
    sit = f.variables['sit'][:]
    
    #mask with the age_mask
    age_mask_tm = np.zeros_like(sit)
    for i in range(0,sit.shape[0]):
        age_mask_tm[i,:,:]=age_mask
    sit = sit*age_mask_tm/1000               #convert from m to km
    vol = np.sum(np.sum(sit*100,axis=1),axis=1)/1e3         #sit is effective sea ice thickness - thickness*concentration

    #store date for timeseries
    time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
    base = datetime(1900,1,1)
    dt = np.array([base + timedelta(days=i) for i in time])
    print(dt[0])
    print(dt[-1])

    volume_model2.extend(vol)
    dates_model2.extend(dt)


##save data
#outfile = outpath+'volume_ts' 
#np.savez(outfile, dt = np.array(dates_model), vm = np.array(volume_model) )

outfile = outpath+'volume_ts2' 
np.savez(outfile, dt = np.array(dates_model2), vm = np.array(volume_model2) )


#load data
container = np.load(outpath+'volume_ts.npz')
dates_model = container['dt']
volume_model = container['vm']

container = np.load(outpath+'volume_ts2.npz')
dates_model2 = container['dt']
volume_model2 = container['vm']


#import to pandas and plot
df1 = pd.DataFrame({ 'PIOMAS' : volume}, index=dates)
dfmon1 = df1.resample('M').mean()
dfmon1_std = df1.resample('M').std()

df2 = pd.DataFrame({ 'neXtSIM' : volume_model}, index=dates_model)
dfmon2 = df2.resample('M').mean()
dfmon2_std = df2.resample('M').std()

df3 = pd.DataFrame({ 'neXtSIM 95_noice' : volume_model2}, index=dates_model2)
dfmon3 = df3.resample('M').mean()
dfmon3_std = df3.resample('M').std()

#merge data
tmp = dfmon1.join(dfmon2, how='outer')
dfmon = tmp.join(dfmon3, how='outer')

plt.figure()
dfmon.loc[dfmon.index.month==9].iloc[15:,:].plot(title='September sea ice volume',lw=3,yerr=dfmon1_std)
plt.ylabel(r'Volume (10$^3$ km$^3$)')
plt.savefig(outpath_plots+'piomas_nextsim.png')


