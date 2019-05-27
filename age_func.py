import numpy as np
import matplotlib.pyplot as plt
import regionmask
import cartopy
import cartopy.crs as ccrs
import pyresample as pr
import scipy.ndimage as ndimage
from pyproj import Proj, transform



def plot_pcolormesh(lons,lats,var,outname,vmin=0,vmax=1,cmap='jet',label='Variable'):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 66, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    pp = plt.pcolormesh(lons,lats,var,vmin=vmin,vmax=vmax, cmap=cmap, transform=ccrs.PlateCarree())


    # add the colourbar to the bottom of the plot.
    # The first moves the bottom of the map up to 15% of total figure height, 
    # the second makes the new axes for the colourbar, 
    # the third makes the colourbar, and the final adds the label
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    stp = (vmax-vmin)/10.
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(vmin,vmax+stp,stp))
    cbar.set_label(label=label,size=14, family='serif')
    
    plt.savefig(outname,bbox_inches='tight')
    
def plot_contour(lons,lats,data,levels=[.15],colors=['purple'],lw=[1.],labels=['Variable'],outname='test.png'):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 66, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    for i in range(len(levels)):   
        cs = plt.contour(lons,lats,data[i],levels=[levels[i]], colors=colors[i], linewidths=lw[i], transform=ccrs.PlateCarree())
        cs.collections[0].set_label(labels[i])
        
    ax1.legend(loc='upper left')
    
    plt.savefig(outname,bbox_inches='tight')
    
def plot_contour_bg(lons,lats,bg,data,levels=[.15],colors=['purple'],lw=[1.],labels=['Variable'],bg_label='Snow_depth',outname='test.png',vmin=0,vmax=1):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy Orthographic projection, selecting the South Pole
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo())
    ax1.set_extent([-180, 180, 66, 90], ccrs.PlateCarree())

    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND, zorder=1,facecolor=cartopy.feature.COLORS['land_alt1'])

    #plot 'background'
    pp = plt.pcolormesh(lons,lats,bg,vmin=vmin,vmax=vmax, cmap='jet', transform=ccrs.PlateCarree())
    # add the colourbar to the bottom of the plot.
    
    # plot sea ice field
    for i in range(len(levels)):   
        cs = plt.contour(lons,lats,data[i],levels=[levels[i]], colors=colors[i], linewidths=lw[i], transform=ccrs.PlateCarree())
        cs.collections[0].set_label(labels[i])
        
    ax1.legend(loc='upper left')
    
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(0,1.1,0.1))
    cbar.set_label(label=bg_label,size=14, family='serif')

    
    plt.savefig(outname,bbox_inches='tight')

def plot_quiver(x,y,u,v,outname,vmin=0,vmax=1,cmap='jet',label='Variable'):
    # create the figure panel 
    fig = plt.figure(figsize=(10,10), facecolor='w')

    # create the map using the cartopy NorthPoleStereo
    # +proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45"
    globe = cartopy.crs.Globe(semimajor_axis=6378273, semiminor_axis=6356889.44891)
    ax1 = plt.subplot(1,1,1, projection=ccrs.NorthPolarStereo(central_longitude=-45, true_scale_latitude=70, globe=globe))
    ax1.set_extent([15, -180, 72, 62], crs=ccrs.PlateCarree())
    
    # add coastlines, gridlines, make sure the projection is maximised inside the plot, and fill in the land with colour
    ax1.coastlines(resolution='110m', zorder=3) # zorder=3 makes sure that no other plots overlay the coastlines
    #ax1.gridlines(crs=ccrs.PlateCarree(),xlocs=range(0,370,10),ylocs=range(60,90,5))
    ax1.gridlines()
    ax1.add_feature(cartopy.feature.LAND,facecolor=cartopy.feature.COLORS['land_alt1'])

    # plot sea ice field
    speed = np.sqrt(u**2+v**2)
    pp = plt.pcolormesh(x,y,speed,vmin=vmin,vmax=vmax, cmap=cmap)
            
    #but our northings and eastings are in the projcted grid and not in lat, lon!!!!
    ax1.quiver(x, y, u, v, scale=5)
        
    # add the colourbar to the bottom of the plot.
    # The first moves the bottom of the map up to 15% of total figure height, 
    # the second makes the new axes for the colourbar, 
    # the third makes the colourbar, and the final adds the label
    fig.subplots_adjust(bottom=0.15)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.625, 0.033])
    stp = (vmax-vmin)/10.
    cbar = plt.colorbar(pp, cax=cbar_ax, orientation='horizontal', ticks=np.arange(vmin,vmax+stp,stp))
    cbar.set_label(label=label,size=14, family='serif')
    
    plt.savefig(outname,bbox_inches='tight')

def smooth_data(data,lon,lat,coarse_lon,coarse_lat):    
    #smoothen the data for nicer contours with a lowpass filter
    data = ndimage.gaussian_filter(data, 3)    #sigma
    
    
    #regrid to equally spaced grid in latlon - otherwise there will be problems with cyclic point in contour plots
    orig_def = pr.geometry.SwathDefinition(lons=lon, lats=lat)
    targ_def = pr.geometry.SwathDefinition(lons=coarse_lon, lats=coarse_lat)
    #coarse_def = pr.geometry.SwathDefinition(lons=coarse_lon[::5,::5], lats=coarse_lat[::5,::5])
    #data_coarse = pr.kd_tree.resample_nearest(orig_def, data, coarse_def, radius_of_influence=50000, fill_value=0)
    ##fill all nans with 0 >> closed contours
    #data_coarse = np.nan_to_num(data_coarse)
    data_smooth = pr.kd_tree.resample_gauss(orig_def, data, targ_def, radius_of_influence=500000, neighbours=10, sigmas=250000, fill_value=0)
    #data_smooth = pr.kd_tree.resample_nearest(coarse_def, data_coarse, targ_def, radius_of_influence=500000, fill_value=0)
    #wf = lambda r: 1
    #data_smooth = pr.kd_tree.resample_custom(coarse_def, data_coarse, targ_def, radius_of_influence=100000, weight_funcs=wf)
    
    
    data_smooth = np.nan_to_num(data_smooth)
    
    
    #plot_pcolormesh(lon,lat,myi,'test.png',vmin=0,vmax=1,label='MYI fraction') 
    #plot_pcolormesh(lon_g[::5],lat_g[::5],myi_coarse,'test1.png',vmin=0,vmax=1,label='MYI fraction')    
    #plot_pcolormesh(lon_g,lat_g,myi_smooth,'test2.png',vmin=0,vmax=1,label='MYI fraction')
    #plot_contour(lon_g,lat_g,myi_smooth,'test3.png',levels=[.1], lw=[10], label='MYI extent')

    return(data_smooth)

def regrid_data(data,inlon,inlat,outlon,outlat):    
    #regrid to equally spaced grid in latlon - otherwise there will be problems with cyclic point in contour plots
    orig_def = pr.geometry.SwathDefinition(lons=inlon, lats=inlat)
    targ_def = pr.geometry.SwathDefinition(lons=outlon, lats=outlat)
    
    data = pr.kd_tree.resample_nearest(orig_def, data, targ_def, radius_of_influence=50000, fill_value=0)
    #fill all nans with 0 >> closed contours
    data = np.nan_to_num(data)
    
    return(data)

def get_poly_mask(lons,lats):
    #create a geographical polygon for the Central Arctic (without the narrow band off the CAA)
    #https://regionmask.readthedocs.io/en/stable/_static/notebooks/create_own_regions.html
    #make two masks - one for W and one for E Arctic
    
    #regionmask does not handle well the circular polygons around the NP
    lon360 = np.where(lons<0,360+lons,lons)
    #print(lon360)

    #i,j coordinates of corner points can be found by exploring display in ncview
    #W Arctic
    poly1 = []
    pt = [360,90];poly1.append(pt)
    pt = [360,lats[273,115]];poly1.append(pt)
    pt = [lon360[273,115],lats[273,115]];poly1.append(pt)
    
    pt = [lon360[260,128],lats[260,128]];poly1.append(pt)
    pt = [lon360[239,136],lats[239,136]];poly1.append(pt)
    pt = [lon360[228,145],lats[228,145]];poly1.append(pt)
    pt = [lon360[210,148],lats[210,148]];poly1.append(pt)
    
    pt = [lon360[194,147],lats[194,147]];poly1.append(pt)
    pt = [lon360[157,156],lats[157,156]];poly1.append(pt)
    pt = [lon360[113,174],lats[113,174]];poly1.append(pt)        
    pt = [lon360[89,157],lats[89,157]];poly1.append(pt)
    pt = [lon360[29,123],lats[29,123]];poly1.append(pt)

    ##more radical (even further away from the CAA coast)
    #pt = [lon360[260,132],lats[260,132]];poly1.append(pt)
    #pt = [lon360[239,140],lats[239,140]];poly1.append(pt)
    #pt = [lon360[228,149],lats[228,149]];poly1.append(pt)
    #pt = [lon360[210,152],lats[210,152]];poly1.append(pt)
    
    #pt = [lon360[194,151],lats[194,151]];poly1.append(pt)
    #pt = [lon360[157,160],lats[157,160]];poly1.append(pt)
    #pt = [lon360[113,178],lats[113,178]];poly1.append(pt)
    #pt = [lon360[65,160],lats[65,160]];poly1.append(pt)
    #pt = [lon360[24,162],lats[29,162]];poly1.append(pt)

    pt = [lon360[3,194],lats[3,194]];poly1.append(pt)
    pt = [lon360[3,344],lats[3,344]];poly1.append(pt)
    pt = [180,65];poly1.append(pt)
    pt = [180,90];poly1.append(pt)
    pt = [270,90];poly1.append(pt)
    pt = [360,90];poly1.append(pt)
    #print(poly1)

    #E Arctic
    poly2 = []
    pt = [0,90];poly2.append(pt)
    pt = [90,90];poly2.append(pt)
    pt = [180,90];poly2.append(pt)
    pt = [180,65];poly2.append(pt)
    pt = [lon360[135,386],lats[135,386]];poly2.append(pt)
    pt = [lon360[238,390],lats[238,390]];poly2.append(pt)
    pt = [lon360[310,344],lats[310,344]];poly2.append(pt)
    pt = [lon360[449,301],lats[449,301]];poly2.append(pt)
    pt = [lon360[350,122],lats[350,122]];poly2.append(pt)
    pt = [0,lats[273,115]];poly2.append(pt)
    pt = [0,90];poly2.append(pt)
    #print(poly2)

    numbers = [0, 1]
    names = ['Arctic_west', 'Arctic_east']
    abbrevs = ['Aw', 'Ae']
    Arctic_mask = regionmask.Regions_cls('Arctic_mask', numbers, names, abbrevs, [poly1, poly2])

    ##Plot polygons in Mercator projection
    #ax=Arctic_mask.plot()
    #ax.set_extent([-180, 180, 45, 90], ccrs.PlateCarree())
    #plt.show()

    #Make raster
    mask = Arctic_mask.mask(lons, lats, wrap_lon=True)
    #Merge mask
    mask = np.where(mask>=0,1,0)
    # pcolormesh does not handle NaNs, requires masked array
    age_mask = np.ma.masked_invalid(mask)
    
    ##Plot mask
    #outpath_plots = 'plots/run04/'
    #outname = outpath_plots+'age_mask_rest.png'
    #plot_pcolormesh(lons,lats,age_mask,outname,cmap='viridis',label='Central Arctic Mask=1')
    #exit()
    
    return(age_mask)

def get_dra_mask(lons,lats):
    #get 'Data Release Area' mask for the Central Arctic, published by Rothrock et al, 2008
    #this is the area for which submarine draft data is available (1979-2000)
    #and all Kwok papers use this area to show trend extended by IS and CS-2 data
    
    #regionmask does not handle well the circular polygons around the NP
    lon360 = np.where(lons<0,360+lons,lons)
    
    poly1 =[[360.,90.],
            [360.,87.],
            [345.,87.],
            [300.,86.58],
            [230.,80.],
            [219.,80.],
            [219.,70.],
            [205.,72.],
            [180.,74.],
            [180.,90.],
            [360.,90.]]
            
    poly2 =[[  0.,86.],
            [  0.,90.],
            [180.,90.],
            [180.,74.],
            [175.,75.50],
            [172.,78.50],
            [163.,80.50],
            [126.,78.50],
            [110.,84.33],
            [ 80.,84.42],
            [ 57.,85.17],
            [ 33.,83.83],
            [  8.,84.08],
            [  0.,86.]]

    numbers = [0, 1]
    names = ['Arctic_west', 'Arctic_east']
    abbrevs = ['Aw', 'Ae']
    Arctic_mask = regionmask.Regions_cls('DRA_mask', numbers, names, abbrevs, [poly1,poly2])

    #Make raster
    mask = Arctic_mask.mask(lons, lats, wrap_lon=True)
    #Merge mask
    mask = np.where(mask>=0,1,0)
    # pcolormesh does not handle NaNs, requires masked array
    age_mask = np.ma.masked_invalid(mask)
    
    ##Plot mask
    #outpath_plots = 'plots/new/'
    #outname = outpath_plots+'age_mask_DRA.png'
    #plot_pcolormesh(lons,lats,age_mask,outname,cmap='viridis',label='Central Arctic Mask=1')
    #exit()
    
    return(age_mask)

def read_sir(sirfile):
    #Matlab code
    #fid=fopen(filename,'r','ieee-be');     %ieee-be is big endian integer in python: <i2
    #head=fread(fid,[256],'short');  % read header %usage:A = fread(fileID,sizeA,precision) %short: signed integers, 16bit, 2 byte (same as int16)

    data = np.fromfile(sirfile, dtype='<f')
    print(data)
    print(data.shape)
    print(data[0])
    head = data[:256]
    #print(head)
    
    with open(sirfile,'rb') as fin:
        header = fin.read(256)
        print(header)
        hh = np.fromstring(header, dtype=np.int16)
        nhead = hh[40]  #number of data blocks
        ipol = hh[44]  #polarisation (valid options: 0=n/a,1=H,2=V)
        idatatype = hh[47]  #head(48) = idatatype            ! data type code 0,2=i*2,1=i*1,4=f
        print(nhead)
        print(ipol)
        print(idatatype)
        
    exit()
    return()

def corr_pearson(x, y):

    """
    Compute Pearson correlation.
    """

    x_mean = np.mean(x, axis=0)
    x_stddev = np.std(x, axis=0)

    y_mean = np.mean(y, axis=0)
    y_stddev = np.std(y, axis=0)

    x1 = (x - x_mean)/x_stddev
    y1 = (y - y_mean)/y_stddev

    x1y1mult = x1 * y1

    x1y1sum = np.sum(x1y1mult, axis=0)

    corr = x1y1sum/x.shape[0]

    return corr


def corr_pearson_circ(x, y):

    """
    Compute Pearson correlation for circular data (angles).
    """

    #calculate means
    pirad = np.radians(np.pi)
    ssx = np.sum(np.sin(x),axis=0)
    scx = np.sum(np.cos(x),axis=0)
    x_mean = np.where((ssx>0)&(scx>0),  np.arctan(ssx/scx)          ,0)
    x_mean = np.where((scx<0),          np.arctan(ssx/scx)+pirad    ,x_mean)
    x_mean = np.where((ssx<0)&(scx>0),  np.arctan(ssx/scx)+2*pirad  ,x_mean)
    
    print(x.shape)
    print(x_mean.shape)
    
    
    ssy = np.sum(np.sin(y),axis=0)
    scy = np.sum(np.cos(y),axis=0)
    y_mean = np.where((ssy>0)&(scy>0),  np.arctan(ssy/scy)           ,0)
    y_mean = np.where((scy<0),          np.arctan(ssy/scy)+pirad     ,y_mean)
    y_mean = np.where((ssy<0)&(scy>0),  np.arctan(ssy/scy)+2*pirad   ,y_mean)
    
    #calculate residuals
    resx = np.sin(x-x_mean)
    resy = np.sin(y-y_mean)

    #calculate Pearson correlation coefficient
    corr = np.sum(resx*resy,axis=0)/np.sqrt(np.sum(resx**2*resy**2,axis=0))
    
    print(corr.shape)
    
    return corr

def plot_pdf(l1,l2,outname):
    
    fig = plt.figure(figsize=(8,8), facecolor='w')
    ax = plt.subplot(1,1,1)
    
    
    #plot a PDF
    bl = np.arange(0,.41,.01)
    #n, bins, patches = plt.hist(slist, bl, normed=True, histtype='step', color='m', alpha=.8, label='neXtSIM', lw = 3)
    #n, bins, patches = plt.hist(slist_gauss, bl, normed=True, histtype='step', color='r', alpha=.8, label='neXtSIM', lw = 3)
    n, bins, patches = plt.hist(np.clip(l1, bl[0], bl[-1]), bl, normed=True, histtype='step', color='darkred', alpha=.8, label='neXtSIM', lw = 3)
    n, bins, patches = plt.hist(np.clip(l2, bl[0], bl[-1]), bl, normed=True, histtype='step', alpha=.8, label='OSI-SAF', lw = 3)


    plt.xlabel('Speed (m/s)')
    plt.ylabel('Probability')
    plt.title('Probability distribution of mean 2-day speed \nfor January 2007-2015')
    #plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
    #plt.axis([40, 160, 0, 0.03])
    plt.legend(loc='upper right',prop={'size':16})
    plt.grid(True)
    plt.savefig(outname)

    


#from pynextsim.projection_info import ProjectionInfo
#mm = ProjectionInfo(f)
###def __init__(self,
        ###ecc    = 0.081816153,
        ###a      = 6378.273e3,
        ###lat_0  = 90.,
        ###lon_0  = -45.,
        ###lat_ts = 60.,
        ###proj='stere',
##print(mm.proj,mm.lat_ts,mm.lat_0,mm.lon_0,mm.a,mm.ecc)
#nextsim_proj = '+proj=%s +lat_ts=%f +lat_0=%f +lon_0=%f +a=%f +e=%f +units=m' %(mm.proj,mm.lat_ts,mm.lat_0,mm.lon_0,mm.a,0.081816153)
#inProj  = Proj(nextsim_proj,preserve_units=True)
#outProj = Proj("+init=EPSG:4326") # WGS84 in degrees
#lonc,latc=transform(inProj,outProj,xc,yc)

#Hi @loniitkina, this class can't be initialised in this way.
#To get the default nextsim projection:
#proj=ProjectionInfo()
#You can also get it from
#nbi = NextsimBin(f)
#proj = nbi.mesh_info.projection
#and if you have an mppfile
#proj==ProjectionInfo.init_from_mppfile(mppfile=...)
