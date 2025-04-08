#!/bin/bash
#
ln -sf ../bin/neptune-sa .
#
# Download solar and geomagnetic activity data, when the files are older than 1 da#
if [[ `find "data/fap_day.dat" -mtime +0` ]]; then
  wget -N https://static.sdo.esoc.esa.int/SOLMAG/fap_day.dat --directory-prefix=./data
  wget -N https://static.sdo.esoc.esa.int/SOLMAG/fap_mon.dat --directory-prefix=./data
fi
if [[ `find "data/sw19571001.txt" -mtime +0` ]]; then
  wget -N https://celestrak.com/SpaceData/SW-All.txt --directory-prefix=./data
  mv ./data/SW-All.txt ./data/sw19571001.txt
fi
if [[ `find "data/eop19620101.txt" -mtime +0` ]]; then
  wget -N https://celestrak.com/SpaceData/EOP-All.txt --directory-prefix=./data
  sed -i '11 i BEGIN NGA_COEFFICIENTS\n  59486.00   .145499   .000000  -.014377  -.049803   .124462  -.051964365.25\n435.00   .352522   .000000  -.114858   .050089  -.016756  -.048254365.25435.00\n  59579.00   .070781   .000019   .000000   .000000  -.022000   .006000\n   .000000   .000000   .012000  -.007000 500.0000 500.0000 365.2500 182.6250\n  37 2326 59905  59904 00000     .019212\nEND NGA_COEFFICIENTS\n#' ./data/EOP-All.txt
  mv ./data/EOP-All.txt ./data/eop19620101.txt
fi
if [[ `find "data/SOLFSMY.TXT" -mtime +0` ]]; then
  wget -N http://sol.spacenvironment.net/JB2008/indices/SOLFSMY.TXT --directory-prefix=./data
fi
if [[ `find "data/DTCFILE.TXT" -mtime +0` ]]; then
  wget -N http://sol.spacenvironment.net/JB2008/indices/DTCFILE.TXT --directory-prefix=./data
fi
#
# Create the output directory if it does not exist
#
if [[ ! -d "output" ]]; then
  mkdir output
fi
LD_LIBRARY_PATH=../lib ./neptune-sa
