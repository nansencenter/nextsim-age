import numpy as np
from glob import glob
from datetime import datetime
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

inpath = 'data/run01/'
outpath = 'data/outputs/'
outpath_plots = 'plots/'

#fl = sorted(glob(inpath+'*.bin'))
#print(fl)



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
    #ea = np.ma.array(ea,mask=mask)
    
    ##MYI area
    #fyi_area = np.sum(fyi*ea)
    #area = np.sum(ea*ic)
    #myi_area = area - fyi_area
    
    ##MYI binary values
    #myi = np.where((fyi<.3)&(ic>.5),1,0)
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
    ##aoy = aoy.astype(int)
    
    ##print(aoy)
    ##exit()
    
    ##in January, FYI is up to 4 months old (1/3 of year)
    ##1st year ice area
    #mask = aoy>.3
    #aoy0 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    ##2nd year ice area
    #mask = (aoy<=.3) | (aoy>1.3)
    #aoy1 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    ##MYI area
    #mask = aoy<=1.3
    #aoy2 = np.sum(np.ma.array(ic,mask=mask)*ea)
            
    #aoy0_list.append(aoy0)
    #aoy1_list.append(aoy1)
    #aoy2_list.append(aoy2)
    
    

#np.save('dates',np.array(date_list))
#np.save('myi_area',np.array(myi_area_list))
#np.save('myi_area_bin',np.array(myi_area_bin_list))
##np.save('myi_vol',np.array(myi_vol_list))
##np.save('vol',np.array(vol_list))
#np.save('area',np.array(area_list))
#np.save('aoy0',np.array(aoy0_list))
#np.save('aoy1',np.array(aoy1_list))
#np.save('aoy2',np.array(aoy2_list))

dates = np.load('dates.npy')
myi_area = np.load('myi_area.npy')/1e9 #10^3 km^2
myi_area_bin = np.load('myi_area_bin.npy')/1e9 #10^3 km^2
#myi_vol = np.load('myi_vol.npy')/1e9 # km^3
#vol = np.load('vol.npy')/1e9 # km^3
area = np.load('area.npy')/1e9 #10^3 km^2

aoy0 = np.load('aoy0.npy')/1e9 #10^3 km^2
aoy1 = np.load('aoy1.npy')/1e9
aoy2 = np.load('aoy2.npy')/1e9

#load OSI-SAF sea ice type data
dates_osi = np.load('data/outputs/dates_osi.npy')
fyi_osi = np.load('data/outputs/fyi_osi.npy')/1e3  #10^3 km^2
myi_osi = np.load('data/outputs/myi_osi.npy')/1e3
oi_osi = np.load('data/outputs/oi_osi.npy')/1e3


dates_osi = np.load('dates_osi_jan.npy')
fyi_osi = np.load('fyi_osi_jan.npy')/1e3  #10^3 km^2
myi_osi = np.load('myi_osi_jan.npy')/1e3
oi_osi = np.load('oi_osi_jan.npy')/1e3




#import to pandas and plot
df_osi = pd.DataFrame({ 'MYI OSI-SAF' : myi_osi+oi_osi,
                        'FYI OSI-SAF' : fyi_osi,
                        'total OSI-SAF' : fyi_osi+myi_osi+oi_osi}, index=dates_osi)

#make winter (January) averages
dfmon_osi = df_osi.resample('M').mean()
#print(dfmon.index.month)
#print(dfmon['2008'])
dfjan_osi = dfmon_osi.loc[dfmon_osi.index.month==1]
print(dfjan_osi)



#import to pandas and plot
df = pd.DataFrame({ 'total' : area,
                    'FYI' : area - myi_area,
                    'FYI age' : aoy0,
                    'MYI' : myi_area,
                    'MYI bin' : myi_area_bin,
                    'MYI age' : aoy2}, index=dates)

#make winter (January) averages
dfmon = df.resample('M').mean()
#print(dfmon.index.month)
#print(dfmon['2008'])
dfjan = dfmon.loc[dfmon.index.month==1]
print(dfjan)

#dfmon.loc[(dfmon[index=])]


result = dfjan.join(dfjan_osi, how='outer')
print(result)


#fig1 = df.iloc[:, 4:].plot.area().get_figure()
#fig1.savefig('test1.png')


##fig2 = df.iloc[:, 0].plot().get_figure()
##fig2.savefig('test.png')


#fig3 = dfjan.iloc[:, 0].plot().get_figure()
#fig3.savefig('test2.png')


#fig4 = result.iloc[4:,3:7].plot(title='1. January ice type/age area').get_figure()
#fig4.savefig('osi.png')


fig5 = result.iloc[:, :].plot(title='1. January ice type/age area').get_figure()
fig5.savefig('total_area.png')



