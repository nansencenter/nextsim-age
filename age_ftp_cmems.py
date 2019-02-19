##send all files to ftp
from ftplib import FTP
from glob import glob
import os

ftp = FTP('my.cmems-du.eu')
ftp.login('pitkin','J@nk0vna')
#sea ice concentration
ftp.cwd('Core/ARCTIC_REANALYSIS_PHYS_002_003/dataset-ran-arc-day-myoceanv2-be/2017/12/')
os.chdir('/input_obs_data/data/TOPAZ_new/')
filenames = ftp.nlst() # get filenames within the directory
#print filenames

for filename in filenames:
    local_filename = filename
    file = open(local_filename, 'wb')
    ftp.retrbinary('RETR '+ filename, file.write)

    file.close()

ftp.quit()

