import numpy as np
import pandas as pd
from glob import glob
from netCDF4 import Dataset
import pyresample
from datetime import datetime, timedelta
from pynextsim.nextsim_bin import NextsimBin
import matplotlib.pyplot as plt

from age_func import *

##new simulations 
#cd data/run04_sept
#ln -s /input_obs_data/einar/age_datarmor/*0915T000000Z* .


inpath_is = 'data/icesat/'
inpath_cs = '/input_obs_data/data/CS2_SMOS_v2.0/'
#inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
inpath = 'data/'
outpath = 'data/'
outpath_plots = 'plots/'

##read all the girds and construct a mask
fl = sorted(glob(inpath_is+'icesat_icethk*.nc'))
f = Dataset(fl[0])
loni = f.variables['Lon'][:]
lati = f.variables['Lat'][:]
sit = f.variables['Th'][:]
loni = np.where(loni<-180,180+(loni+180),loni)
lati = np.where(lati>0,lati,0)
mask_is = np.where(sit<0,0,1)

#outpath_plots = 'plots/run04/'
#outname = outpath_plots+'age_mask_comb.png'
#plot_pcolormesh(loni,lati,sit,outname,vmin=0,vmax=400,cmap='viridis',label='Central Arctic Mask=1')
#exit()

#mask out the region with thick/immobile/old ice at the end of the run
fn = inpath+'Moorings.nc.~7~'
f = Dataset(fn)
lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
mask_m = get_poly_mask(lons,lats)

orig_def = pyresample.geometry.SwathDefinition(lons=loni, lats=lati)
targ_def = pyresample.geometry.SwathDefinition(lons=lons, lats=lats)
mask_is = pyresample.kd_tree.resample_nearest(orig_def, mask_is, \
            targ_def, radius_of_influence=30000)

mask_comb = mask_is * mask_m

##Plot mask
#outpath_plots = 'plots/run04/'
#outname = outpath_plots+'age_mask_comb.png'
#plot_pcolormesh(lons,lats,mask_comb,outname,cmap='viridis',label='Central Arctic Mask=1')
#exit()

mask_comb = mask_comb==0

##read the IS sea ice thickness
dates_is = []
thick_is =[]

for fn in fl:
    print(fn)
    
    tmp = fn.split('_')[2]
    yr = 2000+int(tmp.split('0')[1])
    tmp1 = tmp.split('0')[0]
    if tmp1=='FM':
        date = datetime(yr,3,1)
    elif tmp1=='ON':
        date = datetime(yr,11,1)
    print(date)
    
    f = Dataset(fn)
    sit = f.variables['Th'][:]/100 #from cm to m
        
    sit = pyresample.kd_tree.resample_nearest(orig_def, sit, \
                targ_def, radius_of_influence=30000)    
    
    sit = np.ma.array(sit,mask=mask_comb)
    
    ##Plot data
    #outpath_plots = 'plots/run04/'
    #outname = outpath_plots+'icesat_test1.png'
    #plot_pcolormesh(lons,lats,sit,outname,vmin=0,vmax=4,cmap='viridis',label='Central Arctic Mask=1')
    #exit()
    
    msit = np.mean(sit)
    
    dates_is.append(date)
    thick_is.append(msit)

print(dates_is)
print(thick_is)
#exit()

###read the CS-2 sea ice thickness
#fl = sorted(glob(inpath_cs+'awi-cs2smos-l4-sithick-cryosat2_smos_merged-rep-nh25km_ease2-*.nc'))

#f = Dataset(fl[0])
#lonc = f.variables['lon'][:]
#latc = f.variables['lat'][:]
##sit = f.variables['cs2_ice_thickness'][:]
##mask_cs = np.where(sit>0,1,0)

#orig_def_cs = pyresample.geometry.SwathDefinition(lons=lonc, lats=latc)
##mask_cs = pyresample.kd_tree.resample_nearest(orig_def_cs, mask_cs, \
            ##targ_def, radius_of_influence=30000)

##mask_comb_cs = mask_cs * mask_m

####Plot mask
###outpath_plots = 'plots/run04/'
###outname = outpath_plots+'age_mask_comb_cs.png'
###plot_pcolormesh(lons,lats,mask_comb_cs,outname,cmap='viridis',label='Central Arctic Mask=1')
###exit()

#dates_cs = []
#thick_cs =[]
#unc_cs =[]

#for fn in fl:
    #print(fn)

    #tmp = fn.split('-')[-2].split('_')[0]
    #date = datetime.strptime(tmp, "%Y%m%d")                              
    #print(date)
    ##exit()

    #f = Dataset(fn)
    #sit = f.variables['analysis_ice_thickness'][:]
    #unc = f.variables['analysis_thickness_unc'][:]
    
    ##mask_cs = np.where((sit>0)&(unc<.4),1,0)
    #mask_cs = np.where((sit>0),1,0)
    
    #sit = pyresample.kd_tree.resample_nearest(orig_def_cs, sit, \
            #targ_def, radius_of_influence=30000)
    #mask_cs = pyresample.kd_tree.resample_nearest(orig_def_cs, mask_cs, \
                #targ_def, radius_of_influence=30000)
    #unc = pyresample.kd_tree.resample_nearest(orig_def_cs, unc, \
            #targ_def, radius_of_influence=30000)
    
    
    #mask_comb_cs = mask_cs * mask_m * mask_is
    #sit = np.ma.array(sit,mask=mask_comb_cs)
    
    #msit = np.mean(sit)
    #munc = np.mean(unc)
    #print(msit)
    
    #dates_cs.append(date)
    #thick_cs.append(msit)
    #unc_cs.append(munc)


##model data
#dates_m = []
#thick_m = []

#fl = sorted(glob(inpath+'Moorings.nc.~*'))
#print(fl)
#for fn in fl:
    #print(fn)
    #f = Dataset(fn)
    #sit = f.variables['sit'][:]
    
    ##mask with the age_mask
    #mask_comb_tm = np.zeros_like(sit)
    #for i in range(0,sit.shape[0]):
        #mask_comb_tm[i,:,:]=mask_comb
    #sit = np.ma.array(sit,mask=mask_comb_tm)
    #msit = np.mean(np.mean(sit,axis=1),axis=1)

    ##store date for timeseries
    #time = f.variables['time'][:]   #days since 1900-01-01 00:00:00
    #base = datetime(1900,1,1)
    #dt = np.array([base + timedelta(days=i) for i in time])
    #print(dt[0])
    #print(dt[-1])

    #thick_m.extend(msit)
    #dates_m.extend(dt)

##save data
#outfile = outpath+'thickness_ts' 
#np.savez(outfile, dt = np.array(dates_m), vm = np.array(thick_m) )

#outfile = outpath+'thickness_ts_cs' 
#np.savez(outfile, dt = np.array(dates_cs), vm = np.array(thick_cs), uc = np.array(unc_cs) )

#load data
container = np.load(outpath+'thickness_ts.npz')
dates_m = container['dt']
thick_m = container['vm']

container = np.load(outpath+'thickness_ts_cs.npz')
dates_cs = container['dt']
thick_cs = container['vm']
unc_cs = container['uc']

#import to pandas and plot
df1 = pd.DataFrame({ 'ICESat' : thick_is}, index=dates_is)
dfmon1 = df1.resample('M').mean()
dfmon1_std = df1.resample('M').std()

df2 = pd.DataFrame({ 'neXtSIM' : thick_m}, index=dates_m)
dfmon2 = df2.resample('M').mean()
dfmon2_std = df2.resample('M').std()

df3 = pd.DataFrame({ 'CS-2' : thick_cs,
                     'CS-2 unc' : unc_cs}, index=dates_cs)
dfmon3 = df3.resample('M').mean()
dfmon3_std = df3.resample('M').std()


#print(dfmon3)
dfmon3.plot()
plt.figure()
plt.plot(dfmon3.index, dfmon3['CS-2'], 'k')
plt.fill_between(dfmon3.index, dfmon3['CS-2'] - dfmon3['CS-2 unc'], dfmon3['CS-2'] + dfmon3['CS-2 unc'],
                  color='b', alpha=0.2)
#dfmon3.loc[dfmon3.index.month==11].plot(title='Nov',ylim=[0,2.5])
#dfmon3.loc[dfmon3.index.month==12].plot(title='Dec',ylim=[0,2.5])
#dfmon3.loc[dfmon3.index.month==1].plot(title='Jan',ylim=[0,2.5])
#dfmon3.loc[dfmon3.index.month==2].plot(title='Feb',ylim=[0,2.5])
#dfmon3.loc[dfmon3.index.month==3].plot(title='Mar',ylim=[0,2.5])
#dfmon3.loc[dfmon3.index.month==4].plot(title='Apr',ylim=[0,2.5])
plt.show()
#exit()



dfmon = dfmon1.join(dfmon2, how='outer')
dfmon = dfmon.join(dfmon3, how='outer')

print(dfmon.loc[dfmon.index.month==11])
print(dfmon.loc[dfmon.index.month==3])

#plot
fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(10,8))

#dfmon.plot(title='Ice Thickness',lw=3,yerr=dfmon2_std)
ax = dfmon.iloc[:,:-1].loc[dfmon.index.month==11].plot(ax=axes[0],lw=3,yerr=dfmon2_std,ylim=[0,2.5])
bx = dfmon.iloc[:,:-1].loc[dfmon.index.month==3].plot(ax=axes[1],lw=3,yerr=dfmon2_std,ylim=[0,2.5])


#November uncertany
nov = dfmon3.loc[dfmon3.index.month==11]
ax.fill_between(nov.index, nov['CS-2'] - nov['CS-2 unc'], nov['CS-2'] + nov['CS-2 unc'],
                  color='b', alpha=0.2)

#March uncertany
mar = dfmon3.loc[dfmon3.index.month==3]
bx.fill_between(mar.index, mar['CS-2'] - mar['CS-2 unc'], mar['CS-2'] + mar['CS-2 unc'],
                  color='b', alpha=0.2)

ax.set_ylabel(r'November mean Thickness (m)')
bx.set_ylabel(r'March mean Thickness (m)')
plt.savefig(outpath_plots+'thickness_ts.png')


#plot with uncertanty


#In [188]: price = pd.Series(np.random.randn(150).cumsum(),
   #.....:                   index=pd.date_range('2000-1-1', periods=150, freq='B'))
   #.....: 

#In [189]: ma = price.rolling(20).mean()

#In [190]: mstd = price.rolling(20).std()

#In [191]: plt.figure()
#Out[191]: <Figure size 640x480 with 0 Axes>

#In [192]: plt.plot(price.index, price, 'k')
#Out[192]: [<matplotlib.lines.Line2D at 0x7f2b22659f60>]

#In [193]: plt.plot(ma.index, ma, 'b')
#Out[193]: [<matplotlib.lines.Line2D at 0x7f2b49442ac8>]

#In [194]: plt.fill_between(mstd.index, ma - 2 * mstd, ma + 2 * mstd,
   #.....:                  color='b', alpha=0.2)
   #.....: 
#Out[194]: <matplotlib.collections.PolyCollection at 0x7f2b2245bef0>
