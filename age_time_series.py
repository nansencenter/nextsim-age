import scipy.io
from glob import glob
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from age_func import *

outpath = 'data/outputs/'
outpath_plots = 'plots/new/'

#also check: http://www.scp.byu.edu/data/Quikscat/iceage_v2/Quikscat_MYFY.html
#and: http://scp.byu.edu/data/OSCAT/iceage_v2/Oscat_MYFY.html

#also use Quickscat MYI fraction (1999-2009)
#data/sim/data/QUICKSCAT_ice_type/MYIF/
inpath = '/input_obs_data/data/QUICKSCAT_ice_type/MYIF/'
#inpath_sir = 'data/QSCAT_BYU/'

#NSIDC ice age data
inpath_nsidc = 'data/nsidc/'

#get model grid from moorings for age mask
fn = 'data/Moorings.nc.~7~'
f = Dataset(fn)
lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
age_mask = get_poly_mask(lons,lats)

#get QSCAT grid
mat = scipy.io.loadmat('/input_obs_data/data/QUICKSCAT_ice_type/QSgrid.mat')
lat = mat['lat']
lon = mat['lon']
print(lon.shape)

#get EASE grid
fn = inpath_nsidc+'iceage_nh_12.5km_1995.nc'
f = Dataset(fn)
lone = f.variables['longitude'][:]
late = f.variables['latitude'][:]

#landmask
mat = scipy.io.loadmat('/input_obs_data/data/QUICKSCAT_ice_type/sirareamask.mat')
am = mat['areamask']
mat = scipy.io.loadmat('/input_obs_data/data/QUICKSCAT_ice_type/sirlandmask.mat')
lm = mat['sirlandmask']

lm = am+lm
print(lm.shape)

fl = sorted(glob(inpath+'**/MYfraction*106.mat'))

myi_qs = []
dts_qs = []

#for fn in fl:
    #print(fn)
    #mat = scipy.io.loadmat(fn)
    #myi = mat['MYIfraction']
    #myi = myi*np.where(lm>0,1,0)
    
    ###Plot data
    ##outpath_plots = 'plots/new/'
    ##outname = outpath_plots+'qscat_test.png'
    ##plot_pcolormesh(lon,lat,myi,outname,cmap='viridis',label='MYI fraction')
    ##exit()
    
    ##regrid and mask
    #myi_10 = regrid_data(myi,lon,lat,lons,lats)
    #myi_10 = np.where(myi_10>.3,1,0)
    #nph = np.where(lats > 88.,0,1)
    #myi_10 = myi_10*age_mask*nph
        
    #myia = np.sum(myi_10*100)/1e3 #10^3 km^2
    #print(myia)
    
    #myi_qs.append(myia)
    
    #tmp = fn.split('fraction')[-1].split('.')[0]
    #date = datetime.strptime(tmp, "%Y%j")
    #print(date)
    #dts_qs.append(date)
    
    ##Plot data
    #outpath_plots = 'plots/new/'
    #outname = outpath_plots+'qscat_'+tmp+'.png'
    #plot_pcolormesh(lons,lats,myi_10,outname,cmap='viridis',label='MYI QSCAT')
    ##exit()    

##get BYU QSCAT and OSCAT data
##netcdf files: ftp://podaac.jpl.nasa.gov/allData/quikscat/preview/L3/byu_scp/sea_ice_age/arctic/v1
#fl = sorted(glob(inpath_sir+'*.sir'))
#print(fl)
#for fn in fl:
    #print(fn)
    #sir = read_sir(fn)
#exit()

#get ice age from NSIDC
fl = sorted(glob(inpath_nsidc+'iceage*2005.nc'))

myi_nsidc = []
dts_nsidc = []


for fn in fl:
    print(fn)
    f = Dataset(fn)
    aosi = f.variables['age_of_sea_ice'][:]
    #print(aosi.shape)
    
    #pick 2nd week of April and convert to MYI
    wn = int(106/7)
    mask = (aosi[wn,:,:]>1) & (aosi[wn,:,:]<20)   #value 20 is for the landmask
    myie = np.where(mask,1,0)

    #Plot data
    outpath_plots = 'plots/new/'
    outname = outpath_plots+'nsidc_test.png'
    field = np.where(aosi[wn,:,:]>10,0,aosi[wn,:,:])
    plot_pcolormesh(lone,late,field,outname,cmap='viridis',label='ice age',vmin=0,vmax=10)
    exit()

    #regrid and mask
    myi_10 = regrid_data(myie,lone,late,lons,lats)
    #myi_10 = np.where(myi_10>.3,1,0)
    nph = np.where(lats > 88.,0,1)
    myi_10 = myi_10*age_mask*nph
        
    myia = np.sum(myi_10*100)/1e3 #10^3 km^2
    print(myia)
    
    myi_nsidc.append(myia)
    
    yyyy = fn.split('_')[-1].split('.')[0]
    date = datetime(int(yyyy),4,15)
    print(date)
    dts_nsidc.append(date)
    
    #Plot data
    outpath_plots = 'plots/new/'
    outname = outpath_plots+'nsidc_'+yyyy+'.png'
    plot_pcolormesh(lons,lats,myi_10,outname,cmap='viridis',label='MYI NSIDC')
    #exit()    


#load ice type data produced by age_maps_cont.py
container = np.load(outpath+'age_ts.npz')
dates = container['dates']
myi_area = container['myit']/1e3 #10^3 km^2
myi_age = container['myia']/1e3
myi_osi = container['myio']/1e3
myi_it = container['myi_sit'];fyi_it = container['fyi_sit']
myi_rr = container['myi_rr'];fyi_rr = container['fyi_rr']
myi_sd = container['myi_sd'];fyi_sd = container['fyi_sd']
snow_bias = container['sb']
ridge_bias = container['rb']
warm_bias = container['wb']

myi_osi = np.ma.array(myi_osi, mask=myi_osi==0.)    #mask out dummy values before 2005

#import to pandas
df = pd.DataFrame({ 'MYI area OSI-SAF' : myi_osi,
                    'MYI area' : myi_area,
                    'MYI area (age)' : myi_age,
                    'MYI ice thickness' : myi_it,
                    'FYI ice thickness' : fyi_it,
                    'MYI ridge fraction' : myi_rr,
                    'FYI ridge fraction' : fyi_rr,
                    'MYI snow depth' : myi_sd,
                    'FYI snow depth' : fyi_sd}, index=dates)

#make winter (April) averages
dfmon = df.resample('M').mean()
dfmon_std = df.resample('M').std()
dfapr = dfmon.loc[dfmon.index.month==4]
dfapr_std = dfmon_std.loc[dfmon.index.month==4]

print(dfapr)

#QuickSCAT
df_qs = pd.DataFrame({'MYI area QSCAT' : myi_qs}, index=dts_qs)
dfmon_qs = df_qs.resample('M').mean()
dfmon_qs_std = df_qs.resample('M').std()
dfapr_qs = dfmon_qs.loc[dfmon_qs.index.month==4]
dfapr_qs_std = dfmon_qs_std.loc[dfmon_qs.index.month==4]
print(dfapr_qs)

#NSIDC
df_nsidc = pd.DataFrame({'MYI area NSIDC' : myi_nsidc}, index=dts_nsidc)
dfmon_nsidc = df_nsidc.resample('M').mean()
dfmon_nsidc_std = df_nsidc.resample('M').std()
dfapr_nsidc = dfmon_nsidc.loc[dfmon_nsidc.index.month==4]
dfapr_nsidc_std = dfmon_nsidc_std.loc[dfmon_nsidc.index.month==4]
print(dfapr_nsidc)

dfapr_sat = dfapr_qs.join(dfapr_nsidc, how='outer')
dfapr_sat_std = dfapr_qs_std.join(dfapr_nsidc_std, how='outer')


result = dfapr_sat.join(dfapr, how='outer')
result_std = dfapr_sat_std.join(dfapr_std, how='outer')

print(result)

start = dfapr.iloc[:,3:].index[0]
end = dfapr.iloc[:,3:].index[-1]

fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(8,8))

#area
ax = result.iloc[:,:5].plot(ax=axes[0],xlim=(start,end),yerr=result_std,title='April ice type/age',lw=2,color = ['darkblue','royalblue','k','purple', 'salmon'])
ax.set_ylabel(r'Area (10$^3$ km$^2$)')

#thickness
bx = result.iloc[:,5:7].plot(ax=axes[1],xlim=(start,end),yerr=result_std,lw=2,color = ['purple', 'salmon'])
bx.set_ylabel(r'Thickness (m)')

#ridge fraction
cx = result.iloc[:,7:9].plot(ax=axes[2],xlim=(start,end),yerr=result_std,lw=2,color = ['purple', 'salmon'])
cx.set_ylabel(r'Fraction')

#snow depth
dx = result.iloc[:,9:].plot(ax=axes[3],xlim=(start,end),yerr=result_std,lw=2,color = ['purple', 'salmon'])
dx.set_ylabel(r'Depth (m)')


fig.savefig(outpath_plots+'osi.png',bbox_inches='tight')

##make scatter plots of area difference vs. ridges, snow and warm intrusions
diff = myi_age-myi_osi

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
ax[0].scatter(diff, ridge_bias)
ax[0].set_title('Ridges')
ax[1].scatter(diff, snow_bias)
ax[1].set_title('Snow')
ax[2].scatter(diff, warm_bias)
ax[2].set_title('Crust')

#aa.set_xlabel(r'Area difference (10$^3$ km$^2$)')

fig.savefig(outpath_plots+'scatter_rest.png',bbox_inches='tight')



