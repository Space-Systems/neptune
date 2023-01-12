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
git submodule update --init --recursive
################################################################################
#                                                                              #
#                                  Build OPI                                   #
#                                                                              #
################################################################################
# Build OPI
cd OPI || exit
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build || exit
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
# Build libslam
cd libslam || exit
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build || exit
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
#                                Build pFUnit                                  #
#                                                                              #
################################################################################
cd pFUnit || exit
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build || exit
echo "Updating cmake"
export FC=$Fortran_COMPILER
cmake -DSKIP_MPI=yes ../
echo "Building pFUnit"
make
make install
cd ../
echo "Leaving pFUnit"
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
cd lib || exit
ln -sf ../libslam/lib/libslam-Fortran.$LIBSUFFIX .
ln -sf ../OPI/lib/libOPI-Fortran.$LIBSUFFIX .
ln -sf ../OPI/lib/libOPI.$LIBSUFFIX .
cd ..
# Create the include directory if it does not exist
if [[ ! -d "include" ]]; then
  mkdir include
fi
# Create the links to includes needed by NEPTUNE
cd include || exit
ln -sf ../libslam/include/SLAM .
ln -sf ../OPI/include/OPI .
cd ..
# Create the build directory if it does not exist
if [[ ! -d "build" ]]; then
  mkdir build
else
  rm -rf build/*
fi
cd build || exit
echo "Updating cmake"
export PFUNIT_DIR=..//pFUnit/build/installed
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_OPI_SUPPORT=ON ../
echo "Building NEPTUNE"
make install
if [[ $? -ne 0 ]]; then
    echo "Could not build NEPTUNE. Exiting."
    exit $?
fi
echo "Leaving NEPTUNE"
cd ../work || exit
ln -sf ../bin/neptune-sa .
cd ..
echo "Done"
