#!/bin/bash
#
# Sets up the working directory properly with data files needed and retrieves the validation data
#
# Create the poeData directory if it does not exist and retrieve the sentinel poe data
if [[ ! -d "data" ]]; then
  ln -sf ../../../work/data
fi
#
ln -sf ../../../bin/neptune-valsent1b
#
# Download solar and geomagnetic activity data, when the files are older than 1 day
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
if [[ `find "data/SOLFSMY.TXT" -mtime +24` ]]; then
  wget -N http://sol.spacenvironment.net/~JB2008/indices/SOLFSMY.TXT --directory-prefix=./data
fi
if [[ `find "data/DTCFILE.TXT" -mtime +24` ]]; then
  wget -N http://sol.spacenvironment.net/~JB2008/indices/DTCFILE.TXT --directory-prefix=./data
fi
#
# Create the poeData directory if it does not exist and retrieve the sentinel poe data
if [[ ! -d "poeData" ]]; then
  mkdir poeData
  python3 getPoeData.py
fi
# Create the output directory if it does not exist
if [[ ! -d "output" ]]; then
  mkdir output
fi
# Remove old input file
if [[ -f "input/valsent.inp" ]]; then
  rm "input/valsent.inp"
fi
cp -v "input/valsent_forward.inp" "input/valsent.inp"
cat "input/valsent.inp" && echo "\n"
LD_LIBRARY_PATH=../../../lib ./neptune-valsent1b-set
unlink "input/valsent.inp"
