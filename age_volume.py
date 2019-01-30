import numpy as np
import pandas as pd
from glob import glob
from datetime import datetime, timedelta
from pynextsim.nextsim_bin import NextsimBin
import matplotlib.pyplot as plt


inpath_ps = 'data/'
inpath = 'data/run01_sept/'
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
fl = sorted(glob(inpath+'field*0915T000000Z.bin'))
print(fl)

for f in fl:
    #neXtSIM data
    nb = NextsimBin(f)
    ic = nb.get_var('Concentration')
    it = nb.get_var('Thickness')
    ea = nb.get_var('Element_area')
    
    vol = np.sum(it*ea*ic)/1e9/1e3

    tmp = f.split('_')[-1].split('T')[0]
    date = datetime.strptime(tmp, "%Y%m%d")
    print(date)

    volume_model.append(vol)
    dates_model.append(date)

    
#import to pandas and plot
df1 = pd.DataFrame({ 'PIOMAS' : volume}, index=dates)
dfmon1 = df1.resample('M').mean()

df2 = pd.DataFrame({ 'neXtSIM' : volume_model}, index=dates_model)
dfmon2 = df2.resample('M').mean()

dfmon = dfmon1.join(dfmon2, how='outer')

plt.figure()
dfmon.loc[dfmon.index.month==9].iloc[27:,:].plot(title='September sea ice volume',lw=3)
plt.ylabel(r'Volume (10$^3$ km$^3$)')
plt.savefig(outpath_plots+'piomas.png')


