import numpy as np
from glob import glob
from datetime import datetime
import matplotlib.pyplot as plt
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

#make Jan 1 maps (from 2006 on)
inpath = 'data/run01/'
icosi_path = '/input_obs_data/OSISAF_ice_conc/polstere/'
tyosi_path = '/input_obs_data/OSISAF_ice_type/'
outpath_plots = 'plots/'

fl = sorted(glob(inpath+'field*0101T000000Z.bin'))
print(fl)

for f in fl:
    tmp = f.split('_')[-1].split('T')[0]
    date = datetime.strptime(tmp, "%Y%m%d")
    print(date)
    year = str(date.year)
    if int(year)<2006: continue
    
    #neXtSIM data
    nb = NextsimBin(f)
    ic = nb.get_var('Concentration')
    icthin = nb.get_var('Concentration_thin_ice')
    fyi = nb.get_var('Fyi_fraction')
    
    icthick = ic-icthin

    #OSI-SAF data
    stamp = datetime.strftime(date, "%Y%m%d")+'1200'
    print(stamp)
    
    fosi = icosi_path+year+'_nh_polstere/ice_conc_nh_polstere-100_multi_'+stamp+'.nc'
    ic_osi = nb.get_external_data(fosi,'ice_conc')
    sf_osi = nb.get_external_data(fosi,'status_flag')
    
    fosi = tyosi_path+year+'/01/ice_type_nh_polstere-100_multi_'+stamp+'.nc'
    it_osi = nb.get_external_data(fosi,'ice_type')
  
    #masking
    ic_osi = ic_osi/100  
    landmask=sf_osi>20                      #this will be land mask and north pole hole
    #ic = np.where(landmask,0,ic)
    fyi = np.where(landmask,0,fyi)
    ic_osi = np.where(landmask,0,ic_osi)
    it_osi = np.where(landmask,0,it_osi)
    
    icthick = ic-icthin
    icthick = np.where(landmask,0,icthick)
    
    #ice concentration difference
    ic_diff = icthick - ic_osi
    
    #plotting
    nb.add_var('IC difference', f, ic_diff)
    figname = outpath_plots+'icecon_diff_'+year+'.png'
    nb.plot_var('IC difference',cmap='bwr',figname=figname,clim=[-1,1])

    #ice type difference
    myi = icthick-fyi

    miz = ic_osi<.5
    myi_osi = np.where((it_osi>2)&~miz,1,0)*ic_osi

    it_diff = myi-myi_osi
    nb.add_var('IType difference', f, it_diff)
    figname = outpath_plots+'it_diff_'+year+'.png'
    nb.plot_var('IType difference',cmap='bwr',figname=figname,clim=[-1,1])

    myi_bin = np.where((fyi<.3)&(icthick>.5),1,0)
    myi_osi_bin = np.where((it_osi>2)&~miz,1,0)
    it_diff_bin = myi_bin - myi_osi_bin
    nb.add_var('IType difference bin', f, it_diff_bin)
    figname = outpath_plots+'itbin_diff_'+year+'.png'
    nb.plot_var('IType difference bin',cmap='bwr',figname=figname,clim=[-1,1])
    
    #maps of snow and ridge ratio
    figname = outpath_plots+'snow_'+year+'.png'
    nb.plot_var('Snow',figname=figname,clim=[0,1])
    
    figname = outpath_plots+'ridges_'+year+'.png'
    nb.plot_var('Ridge_ratio',figname=figname,clim=[0,1])
    
    #maps of ice age
    age = nb.get_var('Age_o')
    aoy = age/60/60/24/365
    nb.add_var('Age_o_year', f, aoy)
    figname = outpath_plots+'age_'+year+'.png'
    nb.plot_var('Age_o_year',figname=figname,clim=[0,6])
    
    age_bin = np.where(aoy>1.33,1,0)
    age_diff = age_bin - myi_osi_bin
    age_diff = np.where(landmask,0,age_diff)
    nb.add_var('Age_diff', f, age_diff)
    figname = outpath_plots+'age_diff'+year+'.png'
    nb.plot_var('Age_diff',cmap='bwr',figname=figname,clim=[-1,1])
    
    
    
exit()

#plot just the mask
osi_mask = np.where(mask,0,1)
nb.add_var('OSI-SAF mask', f, osi_mask)
nb.plot_var('OSI-SAF mask',figname='osi_mask.png')


###for Glen and state check-up
##f = 'data/run01/field_20060401T000000Z.bin'
#f = '/input_obs_data/FRASIL/run01_part2/nextsim_outputs/field_20060401T000000Z.bin'
#nb = NextsimBin(f)
#nb.plot_var('Thickness',figname='it_for_glen.png',clim=[0,6])
#nb.plot_var('Snow',figname='sn_for_glen.png',clim=[0,1])
#nb.plot_var('Concentration',figname='ic_for_glen.png',clim=[0,1])
#nb.plot_var('Fyi_fraction',figname='fyi_for_glen.png',clim=[0,1])
#age_os = nb.get_var('Age_o')
#age_oy = age_os/60/60/24/365
#print(np.max(age_oy))
#nb.add_var('Age_oy', f, age_oy)
#nb.plot_var('Age_oy',figname='age_for_glen.png')





