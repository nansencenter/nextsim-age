#This is a script to analyze basic model output and follows the tutorial:  nextsim-tools/python/pynextsim/tutorials/basics.ipynb
#https://github.com/nansencenter/nextsim-tools/blob/tutorials/python/pynextsim/tutorials/basics.ipynb


import matplotlib.pyplot as plt
from pynextsim.nextsim_bin import NextsimBin
import numpy as np


#read last output
nb = NextsimBin('/output/test01/mesh/field_final.bin')
f = 'data/run01/field_20060101T000000Z.bin'
nb = NextsimBin(f)

# print date of the file
nb.datetime

#list all available variables
nb.variables

#preview a variable
nb.plot_var('Fyi_fraction')

#scale observable age to years
age_os = nb.get_var('Age_o')
age_oy = age_os/60/60/24/365
print(np.max(age_oy))
nb.add_var('Age_oy', f, age_oy)
nb.plot_var('Age_oy')

#scale true (volume) age to years
age_s = nb.get_var('Age')
age_y = age_s/60/60/24/365
print(np.max(age_y))
nb.add_var('Age_y', f, age_y)
nb.plot_var('Age_y')



#compare ridge ratio of 2 runs
run nextsim-tools/python/pynextsim/scripts/plot_diff_binary_files.py /output/test01/mesh/field_20060904T000000Z.bin /output/test02/mesh/field_20060904T000000Z.bin -d /output/

run nextsim-tools/python/pynextsim/scripts/plot_diff_binary_files.py /output/test04/field_20070903T000000Z.bin /output/out_cpp/mesh/field_20070903T000000Z.bin -d /output/test04




#testing other stuff
# plot SIC from OSISAF
nb.plot_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200611011200.nc',
                     'ice_conc')
# get SIC from osisaf on neXtSIM mesh elements
c_osisaf = nb.get_external_data('/input_obs_data/OSISAF_ice_conc/polstere/2006_nh_polstere/ice_conc_nh_polstere-100_multi_200611011200.nc',
                     'ice_conc')

#OSI-SAF sea ice type
#"  1 -> no ice or very open ice \n",
#"  2 -> relatively young ice\n",
#"  3 -> ice that survived a summer melt\n",
#"  4 -> ambiguous ice type" ;
nb.plot_external_data('/input_obs_data/OSISAF_ice_type/2006/09/ice_type_nh_polstere-100_multi_200609041200.nc',
                      'ice_type')
t_osisaf = nb.get_external_data('/input_obs_data/OSISAF_ice_type/2006/09/ice_type_nh_polstere-100_multi_200609041200.nc',
                      'ice_type')


# get SIC from osisaf on neXtSIM mesh elements
d_osisaf = nb.get_external_data('/input_obs_data/OSISAF_ice_drift/2010/01/ice_drift_nh_polstere-625_multi-oi_201001011200-201001031200.nc',
                     'dX')
