import scipy.io
from glob import glob
from netCDF4 import Dataset
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from age_func import *

outpath = 'data/outputs/'
outpath_plots = 'plots/new/'



#load ice type data produced by age_maps_cont.py
container = np.load(outpath+'age_ts_winter.npz')
dates = container['dates']
myi_area = container['myit']/1e3 #10^3 km^2
myi_age = container['myia']/1e3
myi_osi = container['myio']/1e3
myi_it = container['myi_sit'];fyi_it = container['fyi_sit']
myi_rr = container['myi_rr'];fyi_rr = container['fyi_rr']
myi_sd = container['myi_sd'];fyi_sd = container['fyi_sd']
snow_bias = container['sb']
ridge_bias = container['rb']
warm_bias = container['wb']

myi_osi = np.ma.array(myi_osi, mask=myi_osi==0.)    #mask out dummy values before 2005

#import to pandas
df = pd.DataFrame({ 'MYI area OSI-SAF' : myi_osi,
                    'MYI area' : myi_area,
                    'MYI area (age)' : myi_age,
                    'MYI ice thickness' : myi_it,
                    'FYI ice thickness' : fyi_it,
                    'MYI ridge fraction' : myi_rr,
                    'FYI ridge fraction' : fyi_rr,
                    'MYI snow depth' : myi_sd,
                    'FYI snow depth' : fyi_sd}, index=dates)

#make winter (April) averages
dfmon = df.resample('M').mean()
dfmon_std = df.resample('M').std()
dfapr = dfmon.loc[dfmon.index.month==4]
dfapr_std = dfmon_std.loc[dfmon.index.month==4]

print(dfapr)


#start = dfapr.iloc[:,3:].index[0]
#end = dfapr.iloc[:,3:].index[-1]

#fig, axes = plt.subplots(nrows=4, ncols=1,figsize=(8,8))

##area
#ax = dfapr.iloc[:,:5].plot(ax=axes[0],xlim=(start,end),yerr=dfapr_std,title='April ice type/age',lw=2,color = ['darkblue','royalblue','k','purple', 'salmon'])
#ax.set_ylabel(r'Area (10$^3$ km$^2$)')

##thickness
#bx = dfapr.iloc[:,5:7].plot(ax=axes[1],xlim=(start,end),yerr=dfapr_std,lw=2,color = ['purple', 'salmon'])
#bx.set_ylabel(r'Thickness (m)')

##ridge fraction
#cx = dfapr.iloc[:,7:9].plot(ax=axes[2],xlim=(start,end),yerr=dfapr_std,lw=2,color = ['purple', 'salmon'])
#cx.set_ylabel(r'Fraction')

##snow depth
#dx = dfapr.iloc[:,9:].plot(ax=axes[3],xlim=(start,end),yerr=dfapr_std,lw=2,color = ['purple', 'salmon'])
#dx.set_ylabel(r'Depth (m)')


#fig.savefig(outpath_plots+'osi_winter.png',bbox_inches='tight')

##make scatter plots of area difference vs. ridges, snow and warm intrusions
diff = myi_age-myi_osi

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,4))
ax[0].scatter(diff, ridge_bias)
ax[0].set_title('Ridges')
ax[1].scatter(diff, snow_bias)
ax[1].set_title('Snow')
ax[2].scatter(diff, warm_bias)
ax[2].set_title('Crust')

ax[0].set_xlim(-1000,1000)
ax[1].set_xlim(-1000,1000)
ax[2].set_xlim(-1000,1000)

ax[0].set_ylim(.2,.5)
ax[2].set_ylim(0,.08)

#aa.set_xlabel(r'Area difference (10$^3$ km$^2$)')

fig.savefig(outpath_plots+'scatter_rest_winter.png',bbox_inches='tight')



