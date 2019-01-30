import numpy as np
from netCDF4 import Dataset
from glob import glob
from datetime import datetime,timedelta
import matplotlib.pyplot as plt
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy
import cartopy.crs as ccrs
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
#print(fl)

date_list = []
myi_list = []

#get the land and NP-hole mask (north pole hole got smaller over time!!!)
f = Dataset(fl[0])
sf = f.variables['status_flag'][0,:,:]
lats = f.variables['lat'][:]
lons = f.variables['lon'][:]
xc_osi = f.dimensions['xc']
yc_osi = f.dimensions['yc']

landmask = sf==0

#loop over all the years
#accummulate ice type on maps for each January
years = range(2006,2019)
for yr in years:
    fl = sorted(glob(inpath+'OSISAF_ice_type/'+str(yr)+'/01/ice_type_nh_polstere-100_multi_*1200.nc'))
    print(yr)
    
    freq = np.zeros_like(sf)
    
    for i in range(0,len(fl[:])):
        f = fl[i]
        #print(f)
        
        #date
        tmp = f.split('_')[-1].split('.')[0]
        date = datetime.strptime(tmp, "%Y%m%d%H%M")
        #print(date)
        
        f = Dataset(f)
        itype = f.variables['ice_type'][0,:,:]
        cl = f.variables['confidence_level'][0,:,:]
        #sf = f.variables['status_flag'][:]
        
        #sea ice concentration
        #Attention!!! First 1. and 3. day of 2007 are missing and replaced by 2. Jan 2007!!! (files are copies)
        #similar ice_conc_nh_polstere-100_multi_200601291200.nc is a copy of 30
        #ice_conc_nh_polstere-100_multi_200701071200.nc is a copy of 08
        year = str(date.year)
        f = inpath+'OSISAF_ice_conc/polstere/'+year+'_nh_polstere/ice_conc_nh_polstere-100_multi_'+tmp+'.nc'
        #print(f)
        f = Dataset(f)
        iconc = f.variables['ice_conc'][0,:,:]/100 #scale from % to fraction
        
        miz = iconc<.8
        
        #MYI area
        mask = (itype>2)&landmask&~miz
        myi = np.where(mask,1,0)
        freq = freq + myi
        
    freq = freq/len(fl)
    
    #export data to netCDF4
    netcdf_name = 'OSI-SAF_ice_type_cumul_'+tmp+'.nc'
    print(netcdf_name)
    dataset = Dataset(outpath+netcdf_name, 'w', format='NETCDF4')
    dataset.description = 'Example data'

    # dimensions
    time = dataset.createDimension('time', None)
    xc = dataset.createDimension('xc', len(xc_osi))
    yc = dataset.createDimension('yc', len(yc_osi))

    # variables
    time = dataset.createVariable('time', 'f8', ('time',))
    lat = dataset.createVariable('lat', 'f4', ('yc','xc'))
    lon = dataset.createVariable('lon', 'f4', ('yc','xc'))
    ice_age_cumul = dataset.createVariable('myi_freq', 'f8', ('time', 'yc', 'xc',))

    # data
    lat[:] =  lats
    lon[:] =  lons
    ice_age_cumul[0,:,:] = freq
    
    # Variable Attributes  
    lat.units='degree_north'  
    lon.units='degree_east'  
    ice_age_cumul.units= '%' 
    time.units='seconds since 1978-01-01 00:00:00'  
    time.calendar='standard' 
    
    #time
    from netCDF4 import num2date, date2num
    time[:] = date2num(date, units = time.units,   
                          calendar = time.calendar) 
    print('time values (in units %s): ' % time.units + 
                          '\n', time[:])
    
    #close/save file
    dataset.close()
    

    
    #plot a map
    #%% Make the figure

    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 60, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    pp = plt.pcolormesh(lons,lats,freq,vmin=0,vmax=1, cmap='jet', transform=ccrs.PlateCarree())


    # add the colourbar to the bottom of the plot.
    # The first moves the bottom of the map up to 15% of total figure height, 
    # the second makes the new axes for the colourbar, 
    # the third makes the colourbar, and the final adds the label
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(0,1.1,0.1))
    cbar.set_label(label='MYI frequency',size=14, family='serif')

    output = 'osi_myi_freq_'+year+'.png'
    plt.savefig(outpath_plots+output,bbox_inches='tight')
    
    #there are any large differences, then the artifacts are changing in a very short time-scale: e.g. atmospheric affects, processing errors
    #in areas where the boundary between MYI-FYI is stable: do they correlate with snow-covered regions? areas where there could be ice in snow (areas where temperatues were episodically high in fall/early winter - new tracer for ice in snow needs to be implemented if this lookes interesting)?
    
    a = 100 #each grid box is 10x10km2
    area = np.sum(np.where(freq>.1,1,0)*iconc)*a
    date_list.append(date)
    myi_list.append(area)
        

#save all the data   
np.save(outpath+'dates_osi_jan_cumul',np.array(date_list))
np.save(outpath+'myi_osi_jan_cumul',np.array(myi_list))

#load OSI-SAF sea ice type data
dates_cumul = np.load(outpath+'dates_osi_jan_cumul.npy')-timedelta(hours=12)
myi_cumul = np.load(outpath+'myi_osi_jan_cumul.npy')/1e3  #10^3 km^2


#load OSI-SAF sea ice type data
dates = np.load(outpath+'dates_osi_jan.npy')
myi = np.load(outpath+'myi_osi_jan.npy')/1e3
#import to pandas and plot
df = pd.DataFrame({ 'MYI area' : myi}, index=dates)
#make winter (January) averages
dfmon = df.resample('M').mean()
dfmon_std = df.resample('M').std()
dfjan = dfmon.loc[dfmon.index.month==1]
dfjan_std = dfmon_std.loc[dfmon_std.index.month==1]

#import to pandas and plot
tmp = pd.DataFrame({ 'MYI area cumul' : myi_cumul}, index=dates_cumul)
result = dfjan.join(tmp, how='outer')


print(dfjan)
print(tmp)
print(result)

fig1 = result.iloc[:, :].plot(yerr=dfjan_std).get_figure()
fig1.savefig('osi_test_cumul.png')



