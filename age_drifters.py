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
sf = f.variables['status_flag'][:]

#lists for PDFs
sl_all = []
slo_all = []

#for every year (2011-) collect all winter data (November-April)
years = range(2011,2016)
for yr in years:
    print(yr)
    d = 0                                                   #reset day-counter
    
    #make empty arrays
    wdays = 30+31+31+28+31+30                               #number of winter days
    if yr%4 ==0 : wdays = wdays+1; print('leap year')       #leap year
    ws = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wa = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wso = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wao = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    
    wex = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wey = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    wic = np.ones([wdays,yc.shape[0],xc.shape[0]], dtype=float)*-999
    
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
        sic = f.variables['sic'][0,:,0]
        
        #project lat,lon coordinates and calculate displacements
        #use OSI-SAF projection: proj4_string = "+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
        wgs84=pyproj.Proj("+init=EPSG:4326") 
        nh_stere=pyproj.Proj("+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45")
        x0,y0 = pyproj.transform(wgs84, nh_stere,lons0,lats0)
        x1,y1 = pyproj.transform(wgs84, nh_stere,lons1,lats1)
        dx = x1-x0
        dy = y1-y0
        
        #put displacements on a regular grid (pyresample) - they should be very close in space and no information lost by interpolation
        swath_def = pyresample.geometry.SwathDefinition(lons=lons0, lats=lats0)
        targ_def = pyresample.geometry.SwathDefinition(lons=lon_osi, lats=lat_osi)
        
        dx_g = pyresample.kd_tree.resample_nearest(swath_def, dx, targ_def, radius_of_influence=65000, fill_value=None)     #undefined pixles are masked
        dy_g = pyresample.kd_tree.resample_nearest(swath_def, dy, targ_def, radius_of_influence=65000, fill_value=None)
        sic_g = pyresample.kd_tree.resample_nearest(swath_def, sic, targ_def, radius_of_influence=62500)

        #get velocities
        tm = (time[0]-time[1])*24*60*60 #this is exactly 2 days anyway
        u = dx_g/tm
        v = dy_g/tm
        
        #when compared to OSI-SAF dX,dY, we notice a difference in grid orientation, therefore we invert the sign here (see bellow):
        dy_g = -dy_g
        v = -v
                
        ##quiver plot
        #outname = outpath_plots+'drifters_test_arrows.png'
        #plot_quiver(xc,yc,u,v,outname,cmap='viridis',label='speed',vmin=0,vmax=.2)
        #exit()

        #make corresponding OSI-SAF maps
        #OSI-SAF data
        try: netcdf_name = glob(drosi_path+year+'/'+mon+'/ice_drift_nh_polstere-625_multi-oi_'+year+mon+day+'*.nc')[0]
        except:
            #the last days of the month are already in the next month folder
            mon1 = (dt + timedelta(weeks=4)).strftime("%m")
            year1 = year
            if int(mon)==12: year1=str(int(year)+1); mon1='01'
            try: netcdf_name = glob(drosi_path+year1+'/'+mon1+'/ice_drift_nh_polstere-625_multi-oi_'+year+mon+day+'*.nc')[0]
            except:
                #some files are simply missing, those should days should not be analysed
                d = d+1
                continue
        print(netcdf_name)
        
        f = Dataset(netcdf_name) 
        dX = f.variables['dX'][0,:,:]*1000
        dY = f.variables['dY'][0,:,:]*1000
        
        ##re-calculate from coordinates and compare
        #lat1_osi = f.variables['lat1'][0,:,:]
        #lon1_osi = f.variables['lon1'][0,:,:]
        #x0,y0 = pyproj.transform(wgs84, nh_stere,lon_osi,lat_osi)
        #x1,y1 = pyproj.transform(wgs84, nh_stere,lon1_osi,lat1_osi)
        #dx = x1-x0
        #dy = y1-y0       
        
        #print(dX[85,55])
        #print(dY[85,55])
        
        #print(dx[85,55])
        #print(dy[85,55])
        
        #the coordintes are same in magnitude, but the dY has swapped sign. To compensate for this we turn around v in nextsim (see above)
        
        uo = dX/tm
        vo = dY/tm
                
        ##quiver plot
        #outname = outpath_plots+'drifters_test_arrows_osi.png'
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
        #arctan2 gives the angle from origin at [1,0] - that is 90 deg. clockwise of the N direction
        ang = 90 - np.degrees(np.arctan2(v,u))
        ws[d,:,:] = speed
        wa[d,:,:] = ang
        
        speedo = np.sqrt(uo**2+vo**2)
        ango = 90 - np.degrees(np.arctan2(vo,uo))
        wso[d,:,:] = speedo
        wao[d,:,:] = ango
        
        
        ######################3
        #instead of this i could simply make mean error vectors
        x_error = dX - dx_g
        y_error = dY - dy_g
        
        
        #mask the data where one of the datasets is missing
        mask = dX.mask | dx_g.mask
        x_error = np.ma.array(x_error,mask=mask,fill_value=-999)
        y_error = np.ma.array(y_error,mask=mask,fill_value=-999)
        
        #outname = outpath_plots+'drifters_test_'+str(yr)+'.png'
        #plot_pcolormesh(lon_osi,lat_osi,x_error,outname,vmin=-1000,vmax=1000,cmap='bwr',label='displacement error')        
        #exit()
                
        wex[d,:,:] = x_error/1000           #from m to km
        wey[d,:,:] = y_error/1000
        wic[d,:,:] = sic_g
        
        #outname = outpath_plots+'drifters_test_'+str(yr)+'.png'
        #plot_pcolormesh(lon_osi,lat_osi,mask,outname,vmin=-10,vmax=10,cmap='bwr',label='displacement error')        
        #exit()
        
        ####END OF INNER CYCLE
        #increase day-counter
        d = d+1  
        
    #calculate correlation for each grid cell for each winter
    mask = (ws==0)|(ws==-999)
    ws = np.ma.array(ws,mask=mask)
    wso = np.ma.array(wso,mask=mask)
    wa = np.ma.array(wa,mask=mask)          #masked values became 0 in division
    wao = np.ma.array(wao,mask=mask)
    
    scorr = corr_pearson(wso,ws)
    
    ma,mb,diff,acorr = corr_pearson_circ(wao,wa)
    
    #make correlation maps (for speed and direction) for each winter
    outname = outpath_plots+'drifters_corr_speed_'+str(yr)+'.png'
    plot_pcolormesh(lon_osi,lat_osi,scorr,outname,vmin=0,vmax=1,cmap='Reds',label='correlation')

    outname = outpath_plots+'drifters_corr_angle_'+str(yr)+'.png'
    plot_pcolormesh(lon_osi,lat_osi,acorr,outname,vmin=-1,vmax=1,cmap='bwr',label='correlation')

    #and angle differences
    #import matplotlib.colors
    #circ_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["purple","red","white","blue","purple"])
    #outname = outpath_plots+'drifters_diff_angle1_'+str(yr)+'.png'
    #plot_pcolormesh(lon_osi,lat_osi,ma,outname,vmin=0,vmax=360,cmap='bwr',label='mean angle')
    #outname = outpath_plots+'drifters_diff_angle2_'+str(yr)+'.png'
    #plot_pcolormesh(lon_osi,lat_osi,mb,outname,vmin=0,vmax=360,cmap='bwr',label='mean angle')    
    outname = outpath_plots+'drifters_diff_angle_'+str(yr)+'.png'
    plot_pcolormesh(lon_osi,lat_osi,diff,outname,vmin=-180,vmax=180,cmap='bwr',label='angle diff')

    #and error/residual maps
    mask = (wic<.15)|(sf<10)|(wex==wey)
    wex = np.ma.array(wex,mask=mask)
    wey = np.ma.array(wey,mask=mask)
    
    #mean daily 2-daily displacements
    mwex = np.mean(wex,axis=0)/2        #get km/day
    mwey = np.mean(wey,axis=0)/2
    
    #quiver plot
    outname = outpath_plots+'drifters_residual_'+str(yr)+'.png'
    plot_quiver(xc,yc,mwex,mwey,outname,cmap='viridis',label='mean error motion (km/day)',vmin=0,vmax=3,scale=150)
        
    #make also PDF for speeds for every year
    mask=ws==0
    slist = ws[ws.mask == False].flatten()
    slisto = wso[ws.mask == False].flatten()
    
    outname = outpath_plots+'drifters_pdf_'+str(yr)+'.png'
    plot_pdf(slist,slisto,outname)
    
    #print(slist)
    #exit()
   
    sl_all.extend(slist)
    slo_all.extend(slist)
    
#PDF for the whole period
outname = outpath_plots+'drifters_pdf_all.png'
plot_pdf(sl_all,slo_all,outname)



