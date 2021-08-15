#!bin/bash
LIBSUFFIX="so"
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  LIBSUFFIX="so"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  LIBSUFFIX="dylib"
elif [[ "$OSTYPE" == "CYGWIN"* ]]; then
  LIBSUFFIX="dll"
elif [[ "$OSTYPE" == "MINGW"* ]]; then
  LIBSUFFIX="dll"
fi

ln -sf "$(pwd)"/libslam/build/include libslam/include
ln -sf "$(pwd)"/libslam/build/lib libslam/lib

mkdir lib
ln -sf "$(pwd)"/libslam/lib/libslam-Fortran.$LIBSUFFIX lib/libslam-Fortran.$LIBSUFFIX
ln -sf "$(pwd)"/OPI/lib/libOPI-Fortran.$LIBSUFFIX lib/libOPI-Fortran.$LIBSUFFIX
ln -sf "$(pwd)"/OPI/lib/libOPI.$LIBSUFFIX lib/libOPI.$LIBSUFFIX

mkdir include
ln -sf "$(pwd)"/libslam/include/SLAM include/SLAM
ln -sf "$(pwd)"/OPI/include/OPI include/OPI