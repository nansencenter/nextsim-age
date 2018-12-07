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
#ln -s /input_obs_data/FRASIL/run01_part1/nextsim_outputs/* .   #and same for all the other parts

inpath = '/home/sim/nextsim_age/data/run01/'

#fl = sorted(glob(inpath+'*.bin'))
#print(fl)

date_list = []
myi_area_list = []
area_list = []
myi_vol_list = []
vol_list = []

aoy0_list = []
aoy1_list = []
aoy2_list = []
aoy3_list = []

f = FileList(inpath)
for nb, fn, date in zip(f.objects, f.filelist, f.datetimes):
    print(date)
    date_list.append(date)
    
    ea = nb.get_var('Element_area')
    fyi = nb.get_var('Fyi_fraction')
    ic = nb.get_var('Concentration')
    it = nb.get_var('Thickness')
    ao = nb.get_var('Age_o')
    
    #MYI area
    fyi_area = np.sum(fyi*ea)
    area = np.sum(ea*ic)
    myi_area = area - fyi_area

    #MYI volume
    fyi_vol = np.sum(fyi*ea*it)
    vol = np.sum(ea*ic*it)
    myi_vol = vol - fyi_vol
    
    area_list.append(area)
    vol_list.append(vol)
    myi_area_list.append(myi_area)
    myi_vol_list.append(myi_vol)
    
    #ice age
    #classify ice age into year classes and calculate volume
    aoy = ao/60/60/24/365
    #aoy = aoy.astype(int)
    
    #print(aoy)
    #exit()
    
    
    #1st year ice area
    mask = aoy>1.0
    aoy0 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    #2nd year ice area
    mask = (aoy<=1.0) | (aoy>2.0)
    aoy1 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    #MYI area
    mask = (aoy<=2.0) | (aoy>3.0)
    aoy2 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    mask = aoy<=3.0
    aoy3 = np.sum(np.ma.array(ic,mask=mask)*ea)
    
    #print(area)
    #print(aoy0)
    #print(aoy1)
    #print(aoy2)
    #print(aoy3)
    #exit()
    
    
    aoy0_list.append(aoy0)
    aoy1_list.append(aoy1)
    aoy2_list.append(aoy2)
    aoy3_list.append(aoy3)
    
    
    
#print(date_list)
#print(myi_area_list)
#print(myi_vol_list)

np.save('dates',np.array(date_list))
np.save('myi_area',np.array(myi_area_list))
np.save('myi_vol',np.array(myi_vol_list))
np.save('vol',np.array(vol_list))
np.save('area',np.array(area_list))
np.save('aoy0',np.array(aoy0_list))
np.save('aoy1',np.array(aoy1_list))
np.save('aoy2',np.array(aoy2_list))
np.save('aoy3',np.array(aoy3_list))

dates = np.load('dates.npy')
myi_area = np.load('myi_area.npy')/1e9 #10^3 km^2
myi_vol = np.load('myi_vol.npy')/1e9 # km^3
vol = np.load('vol.npy')/1e9 # km^3
area = np.load('area.npy')/1e9 #10^3 km^2

aoy0 = np.load('aoy0.npy')/1e9 #10^3 km^2
aoy1 = np.load('aoy1.npy')/1e9
aoy2 = np.load('aoy2.npy')/1e9
aoy3 = np.load('aoy3.npy')/1e9

#load OSI-SAF sea ice type data



#import to pandas and plot
df = pd.DataFrame({ 'MYI area' : myi_area,
                     'MYI volume' : myi_vol,
                     'total volume' : vol,
                     'total area' : area,
                     'YI area' : aoy0,
                     'FYI area' : aoy1,
                     'SYI area' : aoy2,
                     'MYI+ area' : aoy3}, index=dates)

#make winter (January) averages
dfmon = df.resample('M').mean()
#print(dfmon.index.month)
#print(dfmon['2008'])
dfjan = dfmon.loc[dfmon.index.month==1]
print(dfjan)

#dfmon.loc[(dfmon[index=])]


fig1 = df.iloc[:, 4:].plot.area().get_figure()
fig1.savefig('test1.png')


##fig2 = df.iloc[:, 0].plot().get_figure()
##fig2.savefig('test.png')


#fig3 = dfjan.iloc[:, 0].plot().get_figure()
#fig3.savefig('test2.png')


