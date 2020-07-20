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
  wget -q https://static.sdo.esoc.esa.int/SOLMAG/fap_day.dat --directory-prefix=./data
  wget -q https://static.sdo.esoc.esa.int/SOLMAG/fap_mon.dat --directory-prefix=./data
fi
if [[ `find "data/sw19571001.txt" -mtime +24` ]]; then
  wget -q https://celestrak.com/SpaceData/SW-All.txt --directory-prefix=./data
  mv ./data//SW-All.txt ./data/sw19571001.txt
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
./neptune-valsent1b
