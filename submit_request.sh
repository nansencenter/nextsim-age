#! /bin/bash

# User modifiable: southern boundary of the region
slat=37

# No modifications needed below

cookie_file=.RDA_cookie

if [ "$1" == "CFSR" ]
then
	dsn="093.0"
elif [ "$1" == "CFSv2" ]
then
	dsn="094.0"
else
	dsn=""
fi

startdate="$2 $3"
enddate="$4 $5"

if [ -z "$dsn" -o -z "$5" -o "$1" == "-h" -o "$1" == "--help" ]
then
	echo "Usage: $0 dataset start_date start_time end_date end_time"
	echo "Dataset is either CFSR (for 1979 to 2010) or CFSv2 (for 2011 onwards)"
	echo "Date format is YYYY-MM-DD and time format is HH:MM"
	exit 1
fi

if [ ! -f $cookie_file ]
then
	echo Enter your RDA registered e-mail address
	read email

	echo Enter your RDA password
	read -s password

	wget -nv -O - --save-cookies $cookie_file --post-data "email=$email&passwd=$password&action=login" https://rda.ucar.edu/cgi-bin/login
fi

wget -nv -O - --load-cookies $cookie_file --post-data \
"dsid=ds$dsn&rtype=S&afmt=BZ2&rinfo=dsnum=$dsn;\
startdate=$startdate;\
enddate=$enddate;\
parameters=3!7-0.2-1:0.1.195,3!7-0.2-1:0.5.192,3!7-0.2-1:0.4.192,3!7-0.2-1:0.3.0,3!7-0.2-1:0.1.0,3!7-0.2-1:0.0.0,3!7-0.2-1:0.2.2,3!7-0.2-1:0.2.3,3!7-0.2-1:0.1.8;\
product=3,609,652;\
grid_definition=57;\
level=107,221,223;\
nlat=90;slat=$slat;wlon=-180;elon=180;\
ofmt=netCDF" https://rda.ucar.edu/php/dsrqst.php
