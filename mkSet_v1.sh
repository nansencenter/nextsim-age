#! /bin/bash

# The variables we use:
# SPF_H_L103 - 2 m speciffic humidity [kg/kg]
# TMP_L103 - 2 m temperature [K]
# DLWRF_L1_Avg_1 - Downward longwave radiation [W/m^2]
# DSWRF_L1_Avg_1 - Downward shortwave radiation [W/m^2]
# U_GRD_L103 - U component of 10 m wind
# V_GRD_L103 - V component of 10 m wind
# A_PCP_L1_Accum_1 - Accumulated total precipitation [m]
# CSNOW_L1_Avg_1 - Snow fall fraction [0-1]
# PRES_L1 - Pressure at surface [Pa]
#
# These are expected to be 6-hourly fields

year0=$1
year1=$2

fcstlen=6

for year in $(seq $year0 $year1)
do
        (( yearm1 = year - 1 ))
        if [ "$(ls pgbh0$fcstlen.gdas.*${yearm1}1231.grb2.nc.bz2 pgbh0$fcstlen.gdas.*$year*.nc.bz2 2>/dev/null)" ] 
        then
                for file in  pgbh0$fcstlen.gdas.*${yearm1}1231.grb2.nc.bz2 pgbh0$fcstlen.gdas.*$year*.nc.bz2
                do
                        [[ ! -f $file ]] && continue
                        bunzip2 $file &
                done
                wait
        fi

        for file in  pgbh0$fcstlen.gdas.*${yearm1}1231.grb2.nc pgbh0$fcstlen.gdas.*$year*.nc
        do
                [[ ! -f $file ]] && continue
                no=$(echo $file | cut -f3 -d.)
                cdo selvar,SPF_H_L103,PRES_L1,U_GRD_L103,V_GRD_L103,TMP_L103 $file forecast.$no.nc || exit 5
                cdo selvar,DLWRF_L1_Avg_1,DSWRF_L1_Avg_1,A_PCP_L1_Accum_1,CSNOW_L1_Avg_1 $file accumulated.$no.nc || exit 5
        done

        cdo cat accumulated.*.nc buffy.nc || exit 5
        cdo settaxis,$year-01-01,0${fcstlen}:00:00,${fcstlen}hour buffy.nc accumulated.nc || exit 5

        cdo cat forecast.*.nc forecast.nc || exit 5

        cdo -O merge forecast.nc accumulated.nc buffy.nc || exit 5
        cdo setreftime,1901-01-01,00:00:00,hours -selyear,$year buffy.nc buffy2.nc || exit 5
        cdo -f nc4 -r splitmon buffy2.nc cfsr.${fcstlen}h.$year || exit 5

        rm accumulated*.nc forecast*.nc buffy*.nc
done
