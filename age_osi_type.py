import numpy as np
from netCDF4 import Dataset
from glob import glob
from datetime import datetime
import matplotlib.pyplot as plt

inpath='/input_obs_data/OSISAF_ice_type/'
outpath='../plots/'

#OSI-SAF sea ice type
#"  1 -> no ice or very open ice \n",
#"  2 -> relatively young ice\n",
#"  3 -> ice that survived a summer melt\n",
#"  4 -> ambiguous ice type" ;


#make the list of all files
fl = sorted(glob(inpath+'**/ice_type_nh_polstere-100_multi_*.nc', recursive=True))
print(fl)


date_list = []
fyi_list = []
myi_list = []
oi_list = []



for f in fl[:10]:
    print(f)
    
    #date
    tmp = f.split('_')[-1].split('.')[0]
    date = datetime.strptime(tmp, "%Y%m%d%H%M")
    #print(date)
    
    f = Dataset(f)
    itype = f.variables['ice_type'][:]
    cl = f.variables['confidence_level'][:]
    sf = f.variables['status_flag'][:]

    ones = np.ones_like(itype)
    a = 100 #each grid box is 10x10km2

    #FYI area
    mask = itype==2
    fyi = np.sum(np.ma.array(ones,mask=~mask))*a

    #MYI area
    mask = itype==3
    myi = np.sum(np.ma.array(ones,mask=~mask))*a

    #other ice area
    mask = (itype==4) | ((itype==-1) & (sf==101) )
    oi = np.sum(np.ma.array(ones,mask=~mask))*a

    date_list.append(date)
    fyi_list.append(fyi)
    myi_list.append(myi)
    oi_list.append(oi)
   
    #print(fyi)
    #print(myi)
    #print(oi)
    #exit()
    

np.save('dates_osi',np.array(date_list))
np.save('fyi_osi',np.array(fyi_list))
np.save('myi_osi',np.array(myi_list))
np.save('oi_osi',np.array(oi_list))

#load OSI-SAF sea ice type data
dates = np.load('dates_osi.npy')
fyi = np.load('fyi_osi.npy')/1e3  #10^3 km^2
myi = np.load('myi_osi.npy')/1e3
oi = np.load('oi_osi.npy')/1e3





#import to pandas and plot
df = pd.DataFrame({ 'FYI area' : fyi,
                     'MYI area' : myi,
                     'other ice area' : oi}, index=dates)

##make winter (January) averages
#dfmon = df.resample('M').mean()
##print(dfmon.index.month)
##print(dfmon['2008'])
#dfjan = dfmon.loc[dfmon.index.month==1]
#print(dfjan)

#dfmon.loc[(dfmon[index=])]


fig1 = df.iloc[:, :].plot.area().get_figure()
fig1.savefig('test3.png')

