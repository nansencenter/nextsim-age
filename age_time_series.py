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
#ln -s /input_obs_data/FRASIL/run01_part1/nextsim_outputs/*0101T000000Z* .   #and same for all the other parts

###scp just 2015 files!!!
#output older than 2015 can be removed from hexagon!!!
#then copy next output into part7!!!!

inpath = 'data/run01/'
outpath = 'data/outputs/'
outpath_plots = 'plots/'

#date_list = []
#myi_area_list = []
#myi_area_bin_list = []
#area_list = []
#myi_vol_list = []
#vol_list = []

#aoy0_list = []
#aoy1_list = []
#aoy2_list = []


#f = FileList(inpath)
#for nb, fn, date in zip(f.objects[:], f.filelist[:], f.datetimes[:]):
    #print(date)
    #date_list.append(date)
    
    #ea = nb.get_var('Element_area')
    #fyi = nb.get_var('Fyi_fraction')
    #ic = nb.get_var('Concentration')
    #icthin = nb.get_var('Concentration_thin_ice')
    ##it = nb.get_var('Thickness')
    #ao = nb.get_var('Age_o')
    
    #ic = ic - icthin
    
    ##masking with OSI-SAF mask
    #sf_osi = nb.get_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200601011200.nc',
                        #'status_flag')
    #mask=sf_osi>20                      #this will be land mask and north pole hole
    ##ea = np.ma.array(ea,mask=mask)
    #ea = np.where(mask,0,ea)
    
    ##MIZ
    #miz = ic<.8
    
    ##MYI area
    #fyi_area = np.sum(fyi*ea)
    #area = np.sum(ea*ic)
    #myi_area = area - fyi_area
    
    ##MYI binary values
    #myi = np.where((fyi<.5)&~miz,1,0)
    #myi_area_bin = np.sum(myi*ic*ea)

    ###MYI volume
    ##fyi_vol = np.sum(fyi*ea*it)
    ##vol = np.sum(ea*ic*it)
    ##myi_vol = vol - fyi_vol
    
    #area_list.append(area)
    ##vol_list.append(vol)
    #myi_area_list.append(myi_area)
    #myi_area_bin_list.append(myi_area_bin)
    ##myi_vol_list.append(myi_vol)
    
    ##ice age
    ##classify ice age into year classes and calculate volume
    #aoy = ao/60/60/24/365
    
    ##in January, FYI is up to 3.5 months old (1/3 of year)
    ##1st year ice area (area composed purely of FYI)
    #mask = (aoy>.33) & ~miz
    #aoy0 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    ##mixed ice area (area composed of ice that is a mixture of FYI and older ice - e.g. lead ice in the central Arctic)
    #mask = (aoy<=.33) | (aoy>1.33) & ~miz
    #aoy1 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    ##MYI area
    #mask = (aoy<=1.33) & ~miz
    #aoy2 = np.sum(np.ma.array(ic,mask=mask)*ea)
            
    #aoy0_list.append(aoy0)
    #aoy1_list.append(aoy1)
    #aoy2_list.append(aoy2)
    
    

#np.save(outpath+'dates',np.array(date_list))
#np.save(outpath+'myi_area',np.array(myi_area_list))
#np.save(outpath+'myi_area_bin',np.array(myi_area_bin_list))
##np.save('myi_vol',np.array(myi_vol_list))
##np.save('vol',np.array(vol_list))
#np.save(outpath+'area',np.array(area_list))
#np.save(outpath+'aoy0',np.array(aoy0_list))
#np.save(outpath+'aoy1',np.array(aoy1_list))
#np.save(outpath+'aoy2',np.array(aoy2_list))

dates = np.load(outpath+'dates.npy')
myi_area = np.load(outpath+'myi_area.npy')/1e9 #10^3 km^2
myi_area_bin = np.load(outpath+'myi_area_bin.npy')/1e9 #10^3 km^2
#myi_vol = np.load('myi_vol.npy')/1e9 # km^3
#vol = np.load('vol.npy')/1e9 # km^3
area = np.load(outpath+'area.npy')/1e9 #10^3 km^2

aoy0 = np.load(outpath+'aoy0.npy')/1e9 #10^3 km^2
aoy1 = np.load(outpath+'aoy1.npy')/1e9
aoy2 = np.load(outpath+'aoy2.npy')/1e9

#load OSI-SAF sea ice type data
dates_osi = np.load(outpath+'dates_osi.npy')
fyi_osi = np.load(outpath+'fyi_osi.npy')/1e3  #10^3 km^2
myi_osi = np.load(outpath+'myi_osi.npy')/1e3
oi_osi = np.load(outpath+'oi_osi.npy')/1e3

#load OSI-SAF cumulated ice type
dates_osi_cumul = np.load(outpath+'dates_osi_jan_cumul.npy')-timedelta(hours=12)
myi_osi_cumul = np.load(outpath+'myi_osi_jan_cumul.npy')/1e3  #10^3 km^2

#import to pandas
df_osi = pd.DataFrame({ 'MYI OSI-SAF' : myi_osi+oi_osi}, index=dates_osi)

tmp = pd.DataFrame({ 'MYI OSI-SAF cumul' : myi_osi_cumul}, index=dates_osi_cumul)

#make winter (January) averages
dfmon_osi = df_osi.resample('M').mean()
dfmon_osi_std = df_osi.resample('M').std()
dfjan_osi = dfmon_osi.loc[dfmon_osi.index.month==1]
dfjan_osi_std = dfmon_osi_std.loc[dfmon_osi_std.index.month==1]
result_osi = dfjan_osi.join(tmp, how='outer')

#import to pandas
df = pd.DataFrame({ 'total' : area,
                    'FYI' : area - myi_area,
                    'FYI age' : aoy0,
                    'MYI' : myi_area,
                    'MYI bin' : myi_area_bin,
                    'MYI age' : aoy2}, index=dates)

#make winter (January) averages
dfmon = df.resample('M').mean()
dfjan = dfmon.loc[dfmon.index.month==1]


#join obs and simulations into one table
result = dfjan.join(result_osi, how='outer')
print(result)


#load
dates_bias = np.load(outpath+'dates_bias.npy')
snow_bias = np.load(outpath+'snow_bias.npy')
ridge_bias = np.load(outpath+'ridge_bias.npy')


bias = pd.DataFrame({ 'mean snow depth' : snow_bias,
                    'mean ridge ratio' : ridge_bias}, index=dates_osi_cumul[:-2])


result_bias = result.join(bias, how='outer')
print(result_bias)



start = result.iloc[:,3:].index[5]
end = result.iloc[:,3:].index[-1]

fig, axes = plt.subplots(nrows=2, ncols=1,figsize=(10,5))

#plt.figure()
result.iloc[:,3:].plot(ax=axes[0],xlim=(start,end),yerr=dfjan_osi_std,title='1. January ice type/age area',lw=2)
plt.ylabel(r'Area (10$^3$ km$^2$)')


result_bias.iloc[:,8:].plot(ax=axes[1],xlim=(start,end),lw=2)

fig.savefig('osi.png')





