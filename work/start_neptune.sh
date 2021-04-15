#!/bin/bash
#
ln -sf ../bin/neptune-sa
#
# Download solar and geomagnetic activity data, when the files are older than 1 da#
if [[ `find "data/fap_day.dat" -mtime +24` ]]; then
  wget -N https://static.sdo.esoc.esa.int/SOLMAG/fap_day.dat --directory-prefix=./data
  wget -N https://static.sdo.esoc.esa.int/SOLMAG/fap_mon.dat --directory-prefix=./data
fi
if [[ `find "data/sw19571001.txt" -mtime +24` ]]; then
  wget -N https://celestrak.com/SpaceData/SW-All.txt --directory-prefix=./data
  mv ./data/SW-All.txt ./data/sw19571001.txt
fi
if [[ `find "data/eop19620101.txt" -mtime +24` ]]; then
  wget -N https://celestrak.com/SpaceData/EOP-All.txt --directory-prefix=./data
  mv ./data/EOP-All.txt ./data/eop19620101.txt
fi
#if [[ `find "data/SOLFSMY.TXT" -mtime +24` ]]; then
#  wget -N http://sol.spacenvironment.net/~JB2008/indices/SOLFSMY.TXT --directory-prefix=./data
#fi
#if [[ `find "data/DTCFILE.TXT" -mtime +24` ]]; then
#  wget -N http://sol.spacenvironment.net/~JB2008/indices/DTCFILE.TXT --directory-prefix=./data
#fi
#
# Create the output directory if it does not exist
#
if [[ ! -d "output" ]]; then
  mkdir output
fi
LD_LIBRARY_PATH=../lib ./neptune-sa
