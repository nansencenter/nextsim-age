import numpy as np
from glob import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import pyresample

from age_func import *

#evaluate neXtSIM drift based on the output in the drifters. Compare to OSI-SAF

inpath='/input_obs_data/polona/FRASIL/age_datamor_long/'
inpath='/data/drifters/'
outpath = 'data/outputs/'
outpath_plots = 'plots/new/'

icosi_path = '/input_obs_data/data/OSISAF_ice_conc/polstere/'
drosi_path = '/input_obs_data/data/OSISAF_ice_drift/'

#get OSI-SAF grid
fn = drosi_path+'2007/01/ice-drift_ice_drift_nh_polstere-625_multi-oi_200701011200-200701031200.nc'
lat_osi = f.variables['lat'][:]
lon_osi = f.variables['lon'][:]

#get all model drifters
fl = sorted(glob(inpath+'field*0415T000000Z.bin'))


#make PDFs






#make scatter plots

#make correlation maps
