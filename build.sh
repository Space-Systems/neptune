#!/bin/bash
#
# Define how to build the libraries and executables:
BUILD_TYPE=Debug
Fortran_COMPILER=gfortran
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
################################################################################
#                                                                              #
#                                  Build OPI                                   #
#                                                                              #
################################################################################
echo "Checking for OPI"
if [[ ! -d "OPI" ]]; then
  echo "Not found - cloning from  https://github.com/Space-Systems/OPI.git --branch OPI-2015"
  git clone https://github.com/Space-Systems/OPI.git --branch OPI-2015
else
  echo "Found - updating branch."
  cd OPI
  git pull
  cd ..
fi
# No need to continue when cloning did not work
if [[ $? -ne 0 ]]; then
  echo "OPI could not be found. Exiting."
  exit $?
fi
# Build OPI
cd OPI
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build
echo "Updating cmake"
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_INSTALL_PREFIX=../ -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_CXX11=ON -DENABLE_CL_SUPPORT=OFF -DENABLE_CUDA_SUPPORT=OFF -DENABLE_PYTHON=OFF ../
echo "Building OPI"
make install
if [[ $? -ne 0 ]]; then
  echo "Could not build OPI. Exiting."
  exit $?
fi
echo "Leaving OPI"
cd ../../
################################################################################
#                                                                              #
#                                Build libslam                                 #
#                                                                              #
################################################################################
echo "Checking for libslam"
if [[ ! -d "libslam" ]]; then
  echo "Not found - cloning from https://github.com/Space-Systems/libslam.git"
  git clone https://github.com/Space-Systems/libslam.git --branch v2020-12
else
  echo "Found - updating branch."
  cd libslam
  git pull
  cd ..
fi
# No need to continue when cloning did not work
if [[ $? -ne 0 ]]; then
  echo "libslam could not be found. Exiting."
  exit $?
fi
# Build libslam
cd libslam
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build
echo "Updating cmake"
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_OpenMP_SUPPORT=OFF -DENABLE_POSTGRESQL_SUPPORT=OFF ../
echo "Building libslam"
make install
if [[ $? -ne 0 ]]; then
    echo "Could not build libslam. Exiting."
    exit $?
fi
echo "Manually preparing 'lib' and 'include' directories"
cd ../
ln -sf build/include include
ln -sf build/lib lib
echo "Leaving libslam"
cd ../
################################################################################
#                                                                              #
#                                 Build NEPTUNE                                #
#                                                                              #
################################################################################
# Create the lib directory if it does not exist
if [[ ! -d "lib" ]]; then
  mkdir lib
fi
# Create the links to libraries needed by CAMP
cd lib
ln -sf ../libslam/lib/libslam-Fortran.$LIBSUFFIX .
ln -sf ../OPI/lib/libOPI-Fortran.$LIBSUFFIX .
ln -sf ../OPI/lib/libOPI.$LIBSUFFIX .
cd ..
# Create the include directory if it does not exist
if [[ ! -d "include" ]]; then
  mkdir include
fi
# Create the links to includes needed by NEPTUNE
cd include
ln -sf ../libslam/include/SLAM .
ln -sf ../OPI/include/OPI .
cd ..
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build
echo "Updating cmake"
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_OpenMP_SUPPORT=OFF -DENABLE_OPI_SUPPORT=ON ../
echo "Building NEPTUNE"
make install
if [[ $? -ne 0 ]]; then
    echo "Could not build NEPTUNE. Exiting."
    exit $?
fi
echo "Leaving NEPTUNE"
cd ../work
ln -sf ../bin/neptune-sa
cd ..
echo "Done"
