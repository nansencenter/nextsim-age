import numpy as np
from netCDF4 import Dataset
from glob import glob
from datetime import datetime
import matplotlib.pyplot as plt
import pandas as pd

inpath='/input_obs_data/'
outpath = 'data/outputs/'
outpath_plots = 'plots/'

#outpath='../plots/'

#OSI-SAF sea ice type
#"  1 -> no ice or very open ice \n",
#"  2 -> relatively young ice\n",
#"  3 -> ice that survived a summer melt\n",
#"  4 -> ambiguous ice type" ;


#make the list of all files
fl = sorted(glob(inpath+'OSISAF_ice_type/**/01/ice_type_nh_polstere-100_multi_*1200.nc', recursive=True))
print(fl)

date_list = []
fyi_list = []
myi_list = []
oi_list = []

#get the land and NP-hole mask (north pole hole got smaller over time!!!)
f = Dataset(fl[0])
sf = f.variables['status_flag'][:]
landmask = sf==0

for i in range(0,len(fl[:])):
    f = fl[i]
    print(f)
    
    #date
    tmp = f.split('_')[-1].split('.')[0]
    date = datetime.strptime(tmp, "%Y%m%d%H%M")
    #print(date)
    
    f = Dataset(f)
    itype = f.variables['ice_type'][:]
    cl = f.variables['confidence_level'][:]
    #sf = f.variables['status_flag'][:]
    
    #sea ice concentration
    #Attention!!! First 1. and 3. day of 2007 are missing and replaced by 2. Jan 2007!!! (files are copies)
    #similar ice_conc_nh_polstere-100_multi_200601291200.nc is a copy of 30
    #ice_conc_nh_polstere-100_multi_200701071200.nc is a copy of 08
    year = str(date.year)
    f = inpath+'OSISAF_ice_conc/polstere/'+year+'_nh_polstere/ice_conc_nh_polstere-100_multi_'+tmp+'.nc'
    print(f)
    f = Dataset(f)
    iconc = f.variables['ice_conc'][:]/100 #scale from % to fraction
    
    miz = iconc<.8
        
    ones = np.ones_like(itype)
    a = 100 #each grid box is 10x10km2

    #FYI area (including lakes and sub-arctic seas)
    mask = (itype==2)&landmask&~miz
    fyi = np.sum(np.ma.array(ones,mask=~mask)*iconc)*a

    #MYI area
    mask = (itype==3)&landmask&~miz
    myi = np.sum(np.ma.array(ones,mask=~mask)*iconc)*a

    #other ice area
    mask = (itype==4)&landmask&~miz
    oi = np.sum(np.ma.array(ones,mask=~mask)*iconc)*a

    date_list.append(date)
    fyi_list.append(fyi)
    myi_list.append(myi)
    oi_list.append(oi)

#save all the data   
np.save(outpath+'dates_osi_jan',np.array(date_list))
np.save(outpath+'fyi_osi_jan',np.array(fyi_list))
np.save(outpath+'myi_osi_jan',np.array(myi_list))
np.save(outpath+'oi_osi_jan',np.array(oi_list))

#load OSI-SAF sea ice type data
dates = np.load(outpath+'dates_osi_jan.npy')
fyi = np.load(outpath+'fyi_osi_jan.npy')/1e3  #10^3 km^2
myi = np.load(outpath+'myi_osi_jan.npy')/1e3
oi = np.load(outpath+'oi_osi_jan.npy')/1e3

#import to pandas and plot
df = pd.DataFrame({ 'FYI area' : fyi,
                     'MYI area' : myi,
                     'other ice area' : oi}, index=dates)

#fig1 = df.iloc[:, :].plot.area().get_figure()
#fig1.savefig('test3.png')


#make winter (January) averages
dfmon = df.resample('M').mean()
dfmon_std = df.resample('M').std()
dfjan = dfmon.loc[dfmon.index.month==1]
print(dfjan)

fig2 = dfmon.loc[dfmon.index.month==1].plot(yerr=dfmon_std).get_figure()
fig2.savefig('osi_test.png')




