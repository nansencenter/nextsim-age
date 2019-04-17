##send all files to ftp
from ftplib import FTP
from glob import glob
import os

ftp = FTP('ftp.scp.byu.edu')
ftp.login('anonymous','polona.itkin@nersc.no')
#sea ice concentration
ftp.cwd('data/qscat/iceage_v2/1999/')
os.chdir('/input_obs_data/data/QSCAT_BYU/1999/')
filenames = ftp.nlst() # get filenames within the directory
print(filenames)

for filename in filenames:
    if filename.split('.')[-1]!='gz':continue           #download just mat files
    #local_filename = filename
    file = open(filename, 'wb')
    ftp.retrbinary('RETR '+ filename, file.write)

    file.close()

ftp.quit()

##sea ice drift
#ftp.cwd('archive/ice/drift_lr/merged/2016/01/')
#os.chdir('/input_obs_data/OSISAF_ice_drift/2016/01/')
#filenames = ftp.nlst() # get filenames within the directory
##print filenames

#for filename in filenames:
    #if filename.split('_')[2]!='nh':continue            #download just Northern H. files
    #local_filename = filename
    #file = open(local_filename, 'wb')
    #ftp.retrbinary('RETR '+ filename, file.write)

    #file.close()

#ftp.quit()
