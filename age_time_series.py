import numpy as np
from glob import glob
from netCDF4 import Dataset
from datetime import datetime,timedelta
import pandas as pd
import matplotlib.pyplot as plt

from pynextsim.nextsim_bin import NextsimBin


from age_func import *

#-------------------------------------------------------------------
#inpath = '/input_obs_data/FRASIL/run01_part3/nextsim_outputs/'

##to get all the files in one folder creat symbolic links as sim data user (symlinks dont realyl work accross the users!)
#ulimit -s 65536
#mkdir data/run01
#cd data/run01
#ln -s /input_obs_data/data/FRASIL/run01_part1/nextsim_outputs/*0101T000000Z* .   #and same for all the other parts

##new simulations (with new dynamical core) is called run04
##best comparision point might be April and not January (end of winter/growth, most ice, least motion)
#cd data/run04
#ln -s /input_obs_data/einar/age_datarmor/*0401T000000Z* .


###scp just 2015 files!!!
#output older than 2015 can be removed from hexagon!!!
#then copy next output into part7!!!!

inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
outpath = 'data/outputs/'
outpath_plots = 'plots/run04/'

##get all daily files
#fl = sorted(glob(inpath+'field*201*04*[0-32]*T000000Z.bin'))

###mask out the region with thick/immobile/old ice at the end of the run
#fn = inpath+'Moorings.nc.~7~'
#f = Dataset(fn)
#lons = f.variables['longitude'][:]
#lats = f.variables['latitude'][:]
#age_mask = get_poly_mask(lons,lats)

##prepare OSI-SAF land and NP hole mask
##get the land and NP-hole mask (north pole hole got smaller over time!!!)
#f = Dataset('/input_obs_data/data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200604011200.nc')
#sf = f.variables['status_flag'][0,:,:]
#osilats = f.variables['lat'][:]
#osilons = f.variables['lon'][:]
#landmask = sf==0

##lists to collect the time series
#dates = []
#myi_area_l = []; fyi_area_l = []
#myi_it_l = []; fyi_it_l = []
#myi_rr_l = []; fyi_rr_l = []
#myi_sd_l = []; fyi_sd_l = []

#myi_area_age_l = []

#for fn in fl:
    #print(fn)
    
    ##get date
    #tmp = fn.split('_')[-1].split('T')[0]
    #dt = datetime.strptime(tmp, "%Y%m%d")
    #dates.append(dt)

    ##mesh coordinates
    #nb = NextsimBin(fn)
    #x,y = nb.get_xy_grids()
    #lon,lat=nb.mapping(x,y,inverse=True)

    ##get data
    #sic = nb.get_gridded_vars(['Concentration'],x,y)['Concentration']
    #sicthin = nb.get_gridded_vars(['Concentration_thin_ice'],x,y)['Concentration_thin_ice']
    #fyi = nb.get_gridded_vars(['Fyi_fraction'],x,y)['Fyi_fraction']
    #myi = sic-sicthin-fyi
    #sd = nb.get_gridded_vars(['Snow'],x,y)['Snow']
    #rr = nb.get_gridded_vars(['Ridge_ratio'],x,y)['Ridge_ratio']
    #sit = nb.get_gridded_vars(['Thickness'],x,y)['Thickness']
    #siad = nb.get_gridded_vars(['Age_d'],x,y)['Age_d']/60/60/24/365          #from seconds to years

    ##regrid everything to the moorings grid
    ##should we regrid everything to OSI-SAF grid instead?
    #sic = regrid_data(sic,lon,lat,lons,lats)
    #fyi = regrid_data(fyi,lon,lat,lons,lats)
    #myi = regrid_data(myi,lon,lat,lons,lats)
    #sd = regrid_data(sd,lon,lat,lons,lats)
    #rr = regrid_data(rr,lon,lat,lons,lats)
    #sit = regrid_data(sit,lon,lat,lons,lats)
    #siad = regrid_data(siad,lon,lat,lons,lats)
    
    ##mask with the age_mask and OSI-SAF mask
    #osi_mask = regrid_data(landmask,osilons,osilats,lons,lats)
    #mask = age_mask | osi_mask
    
    #sic = np.ma.array(sic,mask=mask)
    #fyi = np.ma.array(fyi,mask=mask)
    #myi = np.ma.array(myi,mask=mask)
    #sd = np.ma.array(sd,mask=mask)
    #rr = np.ma.array(rr,mask=mask)
    #sit = np.ma.array(sit,mask=mask)
    #siad = np.ma.array(siad,mask=mask)
    
    ##'binaries'
    #myi_bin = np.where(myi>.5,1,0)
    #fyi_bin = np.where(fyi>.5,1,0)
    
    ##area in km^2!!!
    #myi_area = np.sum(myi_bin*100*sic)
    #fyi_area = np.sum(fyi_bin*100*sic)
    #myi_area_l.append(myi_area); fyi_area_l.append(fyi_area)
    
    ##mean sea ice thickness
    #myi_it = np.mean(myi_bin*sit)
    #fyi_it = np.mean(fyi_bin*sit)
    #myi_it_l.append(myi_it); fyi_it_l.append(fyi_it)
    
    ##mean ridge ratio
    #myi_rr = np.mean(myi_bin*rr*sic)
    #fyi_rr = np.mean(fyi_bin*rr*sic)
    #myi_rr_l.append(myi_rr); fyi_rr_l.append(fyi_rr)
    
    ##mean snow depth
    #myi_sd = np.mean(myi_bin*sd*sic)
    #fyi_sd = np.mean(fyi_bin*sd*sic)
    #myi_sd_l.append(myi_sd); fyi_sd_l.append(fyi_sd)
    
    ##classify ice age into year classes and calculate volume
    ##for reclassifying the age_d into myi
    ##MYI: age>timedelta(today,15. Sept)
    #if (dt.month>9) | ((dt.month==9) & (dt.day>15)):
        #pyr=dt.year
    #else:
        #pyr=dt.year-1
    #diff = (dt-datetime(pyr,9,15)).total_seconds()/60/60/24/356
    
    #diff_mask = np.where(siad<diff,0,1)
    #myi_area_age = np.sum(sic*100*diff_mask)
    #myi_area_age_l.append(myi_area_age)
    
##save ice type data 
#outfile = outpath+'april_myi' 
#np.savez(outfile,dates = np.array(dates),myi_area = np.array(myi_area_l), fyi_area = np.array(fyi_area_l), myi_it = np.array(myi_it_l),fyi_it = np.array(fyi_it_l), \
    #myi_rr = np.array(myi_rr_l), fyi_rr = np.array(fyi_rr_l), myi_sd = np.array(myi_sd_l), fyi_sd = np.array(fyi_sd_l), myi_age =  np.array(myi_area_age_l))    


#load ice type data
container = np.load(outpath+'april_myi.npz')
dates = container['dates']
myi_area = container['myi_area']/1e9 #10^3 km^2
fyi_area = container['fyi_area']/1e9 #10^3 km^2
myi_it = container['myi_it'];fyi_it = container['fyi_it']
myi_rr = container['myi_rr'];fyi_rr = container['fyi_rr']
myi_sd = container['myi_sd'];fyi_sd = container['fyi_sd']
myi_age = container['myi_age']/1e9

#load OSI-SAF cumulated ice type
dates_osi_cumul = np.load(outpath+'dates_osi_jan_cumul.npy')-timedelta(hours=12)
myi_osi_cumul = np.load(outpath+'myi_osi_jan_cumul.npy')/1e3  #10^3 km^2

#import to pandas
df_osi = pd.DataFrame({ 'MYI area OSI-SAF' : myi_osi_cumul}, index=dates_osi_cumul)

#make winter (January) averages
dfmon_osi = df_osi.resample('M').mean()
dfmon_osi_std = df_osi.resample('M').std()
dfapr_osi = dfmon_osi.loc[dfmon_osi.index.month==4]
dfapr_osi_std = dfmon_osi_std.loc[dfmon_osi_std.index.month==4]


#import to pandas
df = pd.DataFrame({ 'MYI area' : myi_area,
                    'MYI area (age)' : myi_age,
                    'MYI ice thickness' : myi_it,
                    'FYI ice thickness' : fyi_it,
                    'MYI ridge fraction' : myi_rr,
                    'FYI ridge fraction' : fyi_rr,
                    'MYI snow depth' : myi_sd,
                    'FYI snow depth' : fyi_sd}, index=dates)

#make winter (April) averages
dfmon = df.resample('M').mean()
dfapr = dfmon.loc[dfmon.index.month==4]


#join obs and simulations into one table
result = dfapr_osi.join(dfapr, how='outer')
print(result)



start = result.iloc[:,3:].index[2]
end = result.iloc[:,3:].index[-1]

fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(10,8))

#area
ax = result.iloc[:,:3].plot(ax=axes[0],xlim=(start,end),yerr=dfapr_osi_std,title='1. April ice type',lw=2,color = ['royalblue','purple', 'salmon'])
ax.set_ylabel(r'Area (10$^3$ km$^2$)')

#thickness
bx = result.iloc[:,3:5].plot(ax=axes[1],xlim=(start,end),yerr=dfapr_osi_std,lw=2,color = ['purple', 'salmon'])
bx.set_ylabel(r'Thickness (m)')

#ridge fraction
cx = result.iloc[:,5:7].plot(ax=axes[2],xlim=(start,end),yerr=dfapr_osi_std,lw=2,color = ['purple', 'salmon'])
cx.set_ylabel(r'Fraction')

#snow depth
dx = result.iloc[:,7:].plot(ax=axes[3],xlim=(start,end),yerr=dfapr_osi_std,lw=2,color = ['purple', 'salmon'])
dx.set_ylabel(r'Depth (m)')


fig.savefig(outpath_plots+'osi.png',bbox_inches='tight')

##make scatter plots of area difference vs. ridges, snow and warm intrusions
##load
#dates = np.load(outpath+'dates_bias.npy')
#snow_bias = np.load(outpath+'snow_bias.npy')
#ridge_bias = np.load(outpath+'ridge_bias.npy')
#warm_bias = np.load(outpath+'warm_bias.npy')

#diff = myi_area[4:]-myi_osi_cumul[:10]


#fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
#ax[0].scatter(diff, ridge_bias)
#ax[0].set_title('Ridges')
#ax[1].scatter(diff, snow_bias)
#ax[1].set_title('Snow')
#ax[2].scatter(diff, warm_bias)
#ax[2].set_title('Crust')

##aa.set_xlabel(r'Area difference (10$^3$ km$^2$)')

#fig.savefig(outpath_plots+'scatter.png',bbox_inches='tight')



