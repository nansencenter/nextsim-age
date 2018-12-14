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
fl = sorted(glob(inpath+'OSISAF_ice_type/**/01/ice_type_nh_polstere-100_multi_*01011200.nc', recursive=True))
print(fl)

#sea ice concentration files
fl_conc = sorted(glob(inpath+'OSISAF_ice_conc/**/ice_conc_nh_polstere-100_multi_*01011200.nc', recursive=True))
print(fl_conc)


date_list = []
fyi_list = []
myi_list = []
oi_list = []



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
    sf = f.variables['status_flag'][:]
    
    #sea ice concentration
    f = fl_conc[i]
    print(f)
    f = Dataset(f)
    iconc = f.variables['ice_conc'][:]/100 #scale from % to fraction
    
    #mask
    mask=sf>20                      #this will be land mask and north pole hole
    iconc = np.ma.array(iconc,mask=mask)
    
    ones = np.ones_like(itype)
    a = 100 #each grid box is 10x10km2

    #FYI area
    mask = itype==2
    fyi = np.sum(np.ma.array(ones,mask=~mask)*iconc)*a

    #MYI area
    mask = itype==3
    myi = np.sum(np.ma.array(ones,mask=~mask)*iconc)*a

    #other ice area
    mask = itype==4
    oi = np.sum(np.ma.array(ones,mask=~mask)*iconc)*a

    date_list.append(date)
    fyi_list.append(fyi)
    myi_list.append(myi)
    oi_list.append(oi)
   
    #print(fyi)
    #print(myi)
    #print(oi)
    #exit()
    

np.save('dates_osi_jan',np.array(date_list))
np.save('fyi_osi_jan',np.array(fyi_list))
np.save('myi_osi_jan',np.array(myi_list))
np.save('oi_osi_jan',np.array(oi_list))

#load OSI-SAF sea ice type data
dates = np.load('dates_osi_jan.npy')
fyi = np.load('fyi_osi_jan.npy')/1e3  #10^3 km^2
myi = np.load('myi_osi_jan.npy')/1e3
oi = np.load('oi_osi_jan.npy')/1e3





#import to pandas and plot
df = pd.DataFrame({ 'FYI area' : fyi,
                     'MYI area' : myi,
                     'other ice area' : oi}, index=dates)

#fig1 = df.iloc[:, :].plot.area().get_figure()
#fig1.savefig('test3.png')



#make winter (January) averages
dfmon = df.resample('M').mean()
#print(dfmon.index.month)
#print(dfmon['2008'])
dfjan = dfmon.loc[dfmon.index.month==1]
print(dfjan)



fig2 = dfjan.iloc[:, :].plot().get_figure()
fig2.savefig('osi_test.png')




