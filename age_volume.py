import numpy as np
import pandas as pd
from glob import glob
from datetime import datetime, timedelta
from pynextsim.nextsim_bin import NextsimBin
import matplotlib.pyplot as plt

##new simulations 
#cd data/run04_sept
#ln -s /input_obs_data/einar/age_datarmor/*0915T000000Z* .


inpath_ps = 'data/'
inpath = '/input_obs_data/polona/FRASIL/age_datamor_long/'
outpath = 'data/'
outpath_plots = 'plots/'

#read the PIOMAS sea ice volume
fn = inpath_ps+'PIOMAS.vol.daily.1979.2018.Current.v2.1.dat'

data = pd.read_fwf(fn, colspecs='infer', widths=[4,5,7])
year = data['Year'].values
day = data['#day'].values
vol = data['Vol'].values

dates = []
volume = []
volume_model = []
dates_model = []

for i in range(0,len(year)):
    date = datetime(year[i],1,1)+timedelta(days=int(day[i]))
    print(date)
    
    dates.append(date)
    volume.append(vol[i])

#model data
fl = sorted(glob(inpath+'field*0915*T000000Z.bin'))
print(fl)
for f in fl:
    
    tmp = f.split('_')[-1].split('T')[0]
    date = datetime.strptime(tmp, "%Y%m%d")
    print(date)

    #neXtSIM data
    nb = NextsimBin(f)
    #ic = nb.get_var('Concentration')
    it = nb.get_var('Thickness')
    ea = nb.get_var('Element_area')
    sia = nb.get_var('Age_d')/60/60/24/365          #from seconds to years
    
    #mask out the ice that older than run-1 year
    mask=sia>date.year-1996
    mask=sia>7
    it = np.ma.array(it,mask=mask)
    vol = np.sum(it*ea)/1e9/1e3                               #sit is effective sea ice thickness - thickness*concentration

    volume_model.append(vol)
    dates_model.append(date)


#save data
outfile = outpath+'volume_ts' 
np.savez(outfile, dates = np.array(dates_model), vm = np.array(volume_model) )

#load data
container = np.load(outpath+'volume_ts.npz')
dates_model = container['dates']
volume_model = container['vm']

#import to pandas and plot
df1 = pd.DataFrame({ 'PIOMAS' : volume}, index=dates)
dfmon1 = df1.resample('M').mean()
dfmon1_std = df1.resample('M').std()

df2 = pd.DataFrame({ 'neXtSIM' : volume_model}, index=dates_model)
dfmon2 = df2.resample('M').mean()
dfmon2_std = df2.resample('M').std()

dfmon = dfmon1.join(dfmon2, how='outer')

plt.figure()
dfmon.loc[dfmon.index.month==9].iloc[15:,:].plot(title='September sea ice volume',lw=3,yerr=dfmon1_std)
plt.ylabel(r'Volume (10$^3$ km$^3$)')
plt.savefig(outpath_plots+'piomas_run04.png')


