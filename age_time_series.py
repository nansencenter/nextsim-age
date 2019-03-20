import numpy as np
from glob import glob
from datetime import datetime,timedelta
import pandas as pd
import matplotlib.pyplot as plt

from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

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

inpath = 'data/run04/'
outpath = 'data/outputs/'
outpath_plots = 'plots/run04/'

#date_l = []
#myi_area_l = []; fyi_area_l = []
#myi_it_l = []; fyi_it_l = []
#myi_rr_l = []; fyi_rr_l = []
#myi_sd_l = []; fyi_sd_l = []

#aoy0_l = []
#aoy1_l = []
#aoy2_l = []


#f = FileList(inpath)
#for nb, fn, date in zip(f.objects[:], f.filelist[:], f.datetimes[:]):
    #print(date)
    #date_l.append(date)
    
    #ea = nb.get_var('Element_area')
    #fyi = nb.get_var('Fyi_fraction')
    #ic = nb.get_var('Concentration')
    #icthin = nb.get_var('Concentration_thin_ice')
    #it = nb.get_var('Thickness')
    #rr = nb.get_var('Ridge_ratio')
    #sd = nb.get_var('Snow')
    #ao = nb.get_var('Age_d')
    
    #myi = ic-icthin-fyi
    #myi_bin = np.where(myi>.1,1,0)
    #fyi_bin = np.where(fyi>.9,1,0)
    
    ##masking with OSI-SAF mask
    #sf_osi = nb.get_external_data('/input_obs_data/data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200601011200.nc',
                        #'status_flag')
    #mask=sf_osi>20                      #this will be land mask and north pole hole
    #ea = np.where(mask,0,ea)
        
    ##area
    #myi_area = np.sum(myi_bin*ea*ic)
    #fyi_area = np.sum(fyi_bin*ea*ic)
    #myi_area_l.append(myi_area); fyi_area_l.append(fyi_area)
    
    ##mean sea ice thickness
    #myi_it = np.mean(myi_bin*it*ic)
    #fyi_it = np.mean(fyi_bin*it*ic)
    #myi_it_l.append(myi_it); fyi_it_l.append(fyi_it)
    
    ##mean ridge ratio
    #myi_rr = np.mean(myi_bin*rr*ic)
    #fyi_rr = np.mean(fyi_bin*rr*ic)
    #myi_rr_l.append(myi_rr); fyi_rr_l.append(fyi_rr)
    
    ##mean snow depth
    #myi_sd = np.mean(myi_bin*sd*ic)
    #fyi_sd = np.mean(fyi_bin*sd*ic)
    #myi_sd_l.append(myi_sd); fyi_sd_l.append(fyi_sd)
    
    
    ##ice age
    ##classify ice age into year classes and calculate volume
    #aoy = ao/60/60/24/365
    
    ##in January, FYI is up to 3.5 months old (1/3 of year)
    ##1st year ice area (area composed purely of FYI)
    #mask = (aoy>.33)
    #aoy0 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    ##mixed ice area (area composed of ice that is a mixture of FYI and older ice - e.g. lead ice in the central Arctic)
    #mask = (aoy<=.33) | (aoy>1.33)
    #aoy1 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    ##MYI area
    #mask = (aoy<=1.33)
    #aoy2 = np.sum(np.ma.array(ic,mask=mask)*ea)
            
    #aoy0_l.append(aoy0)
    #aoy1_l.append(aoy1)
    #aoy2_l.append(aoy2)

##save ice type data 
#outfile = outpath+'april_myi' 
#np.savez(outfile,dates = np.array(date_l),myi_area = np.array(myi_area_l), fyi_area = np.array(fyi_area_l), myi_it = np.array(myi_it_l),fyi_it = np.array(fyi_it_l), \
    #myi_rr = np.array(myi_rr_l), fyi_rr = np.array(fyi_rr_l), myi_sd = np.array(myi_sd_l), fyi_sd = np.array(fyi_sd_l) )    

##save ice age data
#outfile = outpath+'april_aoy' 
#np.savez(outfile,dates = np.array(date_l),aoy0 = np.array(aoy0_l), aoy1 = np.array(aoy1_l), aoy2 = np.array(aoy2_l) )

#load ice type data
container = np.load(outpath+'april_myi.npz')
dates = container['dates']
myi_area = container['myi_area']/1e9 #10^3 km^2
fyi_area = container['fyi_area']/1e9 #10^3 km^2
myi_it = container['myi_it'];fyi_it = container['fyi_it']
myi_rr = container['myi_rr'];fyi_rr = container['fyi_rr']
myi_sd = container['myi_sd'];fyi_sd = container['fyi_sd']

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
                    'FYI area' : fyi_area,
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

#make scatter plots of area difference vs. ridges, snow and warm intrusions
#load
dates = np.load(outpath+'dates_bias.npy')
snow_bias = np.load(outpath+'snow_bias.npy')
ridge_bias = np.load(outpath+'ridge_bias.npy')
warm_bias = np.load(outpath+'warm_bias.npy')

diff = myi_area[4:]-myi_osi_cumul[:10]


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
ax[0].scatter(diff, ridge_bias)
ax[0].set_title('Ridges')
ax[1].scatter(diff, snow_bias)
ax[1].set_title('Snow')
ax[2].scatter(diff, warm_bias)
ax[2].set_title('Crust')

#aa.set_xlabel(r'Area difference (10$^3$ km$^2$)')

fig.savefig(outpath_plots+'scatter.png',bbox_inches='tight')



