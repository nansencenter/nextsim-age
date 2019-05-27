import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample
import pyproj

from age_func import *

#evaluate neXtSIM drift based on the output in the drifters. Compare to OSI-SAF

#compare April daily means

inpath='/input_obs_data/polona/FRASIL/age_datamor_long/'
inpath='data/drifters/'
outpath = 'data/outputs/'
outpath_plots = 'plots/new/'

icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
drosi_path = '/input_obs_data/data/OSISAF_ice_drift/'

#get OSI-SAF grid
fn = drosi_path+'2007/01/ice-drift_ice_drift_nh_polstere-625_multi-oi_200701011200-200701031200.nc'
f = Dataset(fn)
lat_osi = f.variables['lat'][:]
lon_osi = f.variables['lon'][:]
xc = f.variables['xc'][:]*1000      #change from km to m
yc = f.variables['yc'][:]*1000

#for every year (2007-) collect all winter data (November-April)


d=0                                     #day counter
years = range(2011,2016)
for yr in years:
    print(yr)
    
    #make empty arrays
    wdays = 30+31+31+28+31+30           #number of winter days
    ws = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wa = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wso = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wao = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    
    #get all model files
    fl = sorted(
            glob(inpath+'OSISAF_'+str(yr-1)+'11*.nc')+ \
            glob(inpath+'OSISAF_'+str(yr-1)+'12*.nc')+ \
        
            glob(inpath+'OSISAF_'+str(yr)+'01*.nc')+ \
            glob(inpath+'OSISAF_'+str(yr)+'02*.nc')+ \
            glob(inpath+'OSISAF_'+str(yr)+'03*.nc')+ \
            glob(inpath+'OSISAF_'+str(yr)+'04*.nc')  )
    #print(fl)
    
    for fn in fl:
        print(fn)
        f = Dataset(fn)
        time = f.variables['time'][:]
        base = datetime(1900,1,1)
        dt = base + timedelta(days=int(time[0]))
    
        year = str(dt.year)
        mon = dt.strftime("%m")
        day = dt.strftime("%d")
        
        lats0 = f.variables['latitude'][0,:,0]
        lons0 = f.variables['longitude'][0,:,0]
        lats1 = f.variables['latitude'][1,:,0]
        lons1 = f.variables['longitude'][1,:,0]
        #index = f.variables['index'][0,:,0]
        #sic = f.variables['sic'][0,:,0]
        
        #project lat,lon coordinates and calculate displacements
        #use OSI-SAF projection: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
        wgs84=pyproj.Proj("+init=EPSG:4326") 
        nh_stere=pyproj.Proj("+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45")
        x0,y0 = pyproj.transform(wgs84, nh_stere,lons0,lats0)
        x1,y1 = pyproj.transform(wgs84, nh_stere,lons1,lats1)
        dx = x0-x1
        dy = y0-y1
        
        #put displacements on a regular grid (pyresample) - they should be very close in space and no information lost by interpolation
        swath_def = pyresample.geometry.SwathDefinition(lons=lons0, lats=lats0)
        targ_def = pyresample.geometry.SwathDefinition(lons=lon_osi, lats=lat_osi)
        
        dx_g = pyresample.kd_tree.resample_nearest(swath_def, dx, targ_def, radius_of_influence=65000, fill_value=None)     #undefined pixles are masked
        dy_g = pyresample.kd_tree.resample_nearest(swath_def, dy, targ_def, radius_of_influence=65000, fill_value=None)
        #sic_g = pyresample.kd_tree.resample_nearest(swath_def, sic, targ_def, radius_of_influence=62500)

        #get velocities
        tm = (time[0]-time[1])*24*60*60 #this is exactly 2 days anyway
        u = dx_g/tm
        v = dy_g/tm
                
        ##quiver plot
        #outname = outpath_plots+'nextsim_drifters_test_arrows.png'
        #plot_quiver(xc,yc,u,v,outname,cmap='viridis',label='speed',vmin=0,vmax=.2)
        ##exit()

        #make corresponding OSI-SAF maps
        #OSI-SAF data
        try: netcdf_name = glob(drosi_path+year+'/'+mon+'/ice_drift_nh_polstere-625_multi-oi_'+year+mon+day+'*.nc')[0]
        except:
            #the last days of the month are already in the next month folder
            mon1 = (dt + timedelta(weeks=4)).strftime("%m")
            year1 = year
            if int(mon)==12: year1=str(int(year)+1); mon1='01'
            netcdf_name = glob(drosi_path+year1+'/'+mon1+'/ice_drift_nh_polstere-625_multi-oi_'+year+mon+day+'*.nc')[0]
        print(netcdf_name)
        
        f = Dataset(netcdf_name) 
        dX = f.variables['dX'][0,:,:]*1000
        dY = f.variables['dY'][0,:,:]*1000
        
        uo = dX/tm
        vo = dY/tm
                
        ##quiver plot
        #outname = outpath_plots+'nextsim_drifters_test_arrows_osi.png'
        #plot_quiver(xc,yc,uo,vo,outname,cmap='viridis',label='speed',vmin=0,vmax=.2)
        #exit()
        
        #collect all velocites for the month/winter and correlate with OSI-SAF
        #first combine the mask for both datasets
        cmsk = (u.mask) + (uo.mask)
        u = np.ma.array(u.data,mask=cmsk)
        uo = np.ma.array(uo.data,mask=cmsk)
        
        ##quiver plot
        #outname = outpath_plots+'nextsim_drifters_cmsk.png'
        #plot_quiver(xc,yc,u,v,outname,cmap='viridis',label='speed',vmin=0,vmax=.2)
        ##exit()
       
        #make speed,dir arrays for each dataset and add one layer for each day, keep -999 as missing values
        speed = np.sqrt(u**2+v**2)
        ang = np.arccos(u/speed)
        ws[d,:,:] = speed
        wa[d,:,:] = ang             #angles in radians
        
        speedo = np.sqrt(uo**2+vo**2)
        ango = np.arccos(uo/speedo)
        #diro = 
        wso[d,:,:] = speedo
        wao[d,:,:] = ango
        
        #print(wso[d,:,:])
        #print(wao[d,:,:])
        #exit()
        
        d = d+1  
        #exit()
        
        
        
    #calculate correlation for each grid cell for each winter
    ws = np.ma.array(ws,mask=ws==0)
    wso = np.ma.array(wso,mask=wso==0)
    wa = np.ma.array(wa,mask=wa==0)          #masked values became 0 in division
    wao = np.ma.array(wao,mask=wao==0)
    
    scorr = corr_pearson(wso,ws)
    
    acorr = corr_pearson_circ(wao,wa)
    
    print(scorr)    
    print(acorr)
    
    #make correlation maps (for speed and direction) for each winter
    outname = outpath_plots+'drifters_corr_speed.png'
    plot_pcolormesh(lon_osi,lat_osi,scorr,outname,vmin=-1,vmax=1,cmap='bwr',label='correlation')
    
    outname = outpath_plots+'drifters_corr_angle.png'
    plot_pcolormesh(lon_osi,lat_osi,acorr,outname,vmin=-1,vmax=1,cmap='bwr',label='correlation')

        
    #make also PDF for speeds
    #for every year
    
    slist = ws[~mask].flatten()
    slisto = wso[~mask].flatten()
    
    outname = outpath_plots+'drifters_pdf.png'
    plot_pdf(slist,slisto,outname)

    exit()    
    
#PDF for the whole period




