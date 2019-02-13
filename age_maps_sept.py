import numpy as np
from glob import glob
from datetime import datetime
import matplotlib.pyplot as plt
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList

##to get all the files in one folder creat symbolic links as sim data user (symlinks dont realyl work accross the users!)
#ulimit -s 65536
#mkdir data/run01_sept
#cd data/run01_sept
#ln -s /input_obs_data/FRASIL/run01_part6/nextsim_outputs/*0915T000000Z* .   #and same for all the other parts

#September (aka 'switch') sea ice concentration differences
inpath = 'data/run04_sept/'
icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
outpath_plots = 'plots/run04/'


fl = sorted(glob(inpath+'field*0915T000000Z.bin'))
print(fl)

for f in fl:
    #neXtSIM data
    nb = NextsimBin(f)
    ic = nb.get_var('Concentration')
    icthin = nb.get_var('Concentration_thin_ice')

    icthick = ic-icthin
    icthick=ic

    tmp = f.split('_')[-1].split('T')[0]
    date = datetime.strptime(tmp, "%Y%m%d")
    print(date)

    #OSI-SAF data
    year = str(date.year)
    if int(year)<2007: continue
    stamp = datetime.strftime(date, "%Y%m%d")+'1200'
    print(stamp)
    
    fosi = icosi_path+year+'_nh_polstere/ice_conc_nh_polstere-100_multi_'+stamp+'.nc'
    ic_osi = nb.get_external_data(fosi,'ice_conc')
    sf_osi = nb.get_external_data(fosi,'status_flag')
      
    #masking
    ic_osi = ic_osi/100  
    landmask=sf_osi>20                      #this will be land mask and north pole hole
    ic = np.where(landmask,0,ic)
    ic_osi = np.where(landmask,0,ic_osi)

    #icthick = ic-icthin
    icthick = np.where(landmask,0,icthick)
    
    #ice concentration difference
    ic_diff = icthick - ic_osi
    
    #plotting
    nb.add_var('IC difference', f, ic_diff)
    figname = outpath_plots+'icecon_diff_sept'+year+'.png'
    nb.plot_var('IC difference',cmap='bwr',figname=figname,clim=[-1,1])
