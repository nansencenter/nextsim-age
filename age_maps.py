import matplotlib.pyplot as plt
from pynextsim.nextsim_bin import NextsimBin
from pynextsim.file_list import FileList
import numpy as np


inpath = 'data/run01/*0101T000000Z*'
outpath_plots = 'plots/'


f = 'data/run01/field_20060101T000000Z.bin'
nb = NextsimBin(f)

ic = nb.get_var('Concentration')
icthin = nb.get_var('Concentration_thin_ice')
ic = ic-icthin

ic_osi = nb.get_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200601011200.nc',
                    'ice_conc')

#apply masks to both datasets
sf_osi = nb.get_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200601011200.nc',
                    'status_flag')

mask=sf_osi>20                      #this will be land mask and north pole hole
ic = np.where(mask,0,ic)

ic_diff = ic - ic_osi/100

nb.add_var('IC difference', f, ic_diff)
#nb.set_var('IC difference', ic_diff)
nb.plot_var('IC difference',cmap='bwr',figname='test_map1.png',clim=[-1,1])

ic = nb.get_var('Concentration')    #Concentration of thick ice
icthin = nb.get_var('Concentration_thin_ice')
fyi = nb.get_var('Fyi_fraction')

icthick = ic-icthin
myi = icthick-fyi

myi = np.where(mask,0,myi)
icthick = np.where(mask,0,icthick)

it_osi = nb.get_external_data('/input_obs_data/OSISAF_ice_type/2006/01/ice_type_nh_polstere-100_multi_200601011200.nc',
                'ice_type')

#myi_osi = np.where(it_osi==3,1,0)*ic_osi/100
#if other ice type is also MYI 
myi_osi = np.where(it_osi>2,1,0)*ic_osi/100

it_diff = myi-myi_osi
nb.add_var('IType difference', f, it_diff)
nb.plot_var('IType difference',cmap='bwr',figname='test_map2.png',clim=[-1,1])

myi_bin = np.where((fyi<.5)&(icthick>.5),1,0)
#myi_bin = np.where((ic==1),1,0)

#myi_osi_bin = np.where(it_osi==3,1,0)
#if other ice type is also MYI 
myi_osi_bin = np.where(it_osi>2,1,0)

it_diff_bin = myi_bin - myi_osi_bin
nb.add_var('IType difference bin', f, it_diff_bin)
nb.plot_var('IType difference bin',cmap='bwr',figname='test_map3.png',clim=[-1,1])

#switch
f = '/input_obs_data/FRASIL/run01_part1/nextsim_outputs/field_20050915T000000Z.bin'
nb = NextsimBin(f)

ic = nb.get_var('Concentration')
icthin = nb.get_var('Concentration_thin_ice')
ic = ic-icthin

ic_osi = nb.get_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2005_nh_polstere/ice_conc_nh_polstere-100_multi_200509151200.nc',
                    'ice_conc')

#apply masks to both datasets
sf_osi = nb.get_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2005_nh_polstere/ice_conc_nh_polstere-100_multi_200509151200.nc',
                    'status_flag')

mask=sf_osi>20                      #this will be land mask and north pole hole
#ic = np.ma.array(ic,mask=mask)
ic = np.where(mask,0,ic)


ic_diff = ic - ic_osi/100
nb.add_var('IC difference', f, ic_diff)
nb.plot_var('IC difference',cmap='bwr',figname='test_map4.png',clim=[-1,1])


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





