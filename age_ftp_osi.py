##send all files to ftp
from ftplib import FTP
from glob import glob
import os

ftp = FTP('osisaf.met.no')
ftp.login('anonymous','polona.itkin@nersc.no')
ftp.cwd('archive/ice/conc/2008/01/')
#os.chdir('../data/OSI-SAF/icecon')
os.chdir('/input_obs_data/OSISAF_ice_conc/polstere/2008_nh_polstere/')
filenames = ftp.nlst() # get filenames within the directory
#print filenames

for filename in filenames:
    if filename.split('.')[-1]!='nc':continue           #download just nc files
    if filename.split('_')[2]!='nh':continue            #download just Northern H. files
    if filename.split('_')[3]!='polstere-100':continue  #download just the polstere grid files
    local_filename = filename
    file = open(local_filename, 'wb')
    ftp.retrbinary('RETR '+ filename, file.write)

    file.close()

ftp.quit()
