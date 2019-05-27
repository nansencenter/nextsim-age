##send all files to ftp
from ftplib import FTP
from glob import glob
import os

ftp = FTP('osisaf.met.no')
ftp.login('anonymous','polona.itkin@nersc.no')
##sea ice concentration
#ftp.cwd('archive/ice/conc/2015/02/')
#os.chdir('/input_obs_data/data/OSISAF_ice_conc/polstere/2015_nh_polstere/')
#filenames = ftp.nlst() # get filenames within the directory
##print filenames

#for filename in filenames:
    #if filename.split('.')[-1]!='nc':continue           #download just nc files
    #if filename.split('_')[2]!='nh':continue            #download just Northern H. files
    #if filename.split('_')[3]!='polstere-100':continue  #download just the polstere grid files
    #local_filename = filename
    ##print(filename); exit()
    #file = open(local_filename, 'wb')
    #ftp.retrbinary('RETR '+ filename, file.write)

    #file.close()

#ftp.quit()

#sea ice drift
ftp.cwd('archive/ice/drift_lr/merged/2010/02/')
os.chdir('/input_obs_data/data/OSISAF_ice_drift/2007/02/')
filenames = ftp.nlst() # get filenames within the directory
#print filenames

for filename in filenames:
    if filename.split('_')[2]!='nh':continue            #download just Northern H. files
    local_filename = filename
    file = open(local_filename, 'wb')
    ftp.retrbinary('RETR '+ filename, file.write)

    file.close()

ftp.quit()
