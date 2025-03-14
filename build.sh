#!/bin/bash
#
# Define how to build the libraries and executables:
BUILD_TYPE=Debug
Fortran_COMPILER=gfortran
GENERATOR="Unix Makefiles"
ENABLE_OPI_SUPPORT=OFF
LIBSUFFIX="so"
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  LIBSUFFIX="so"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  LIBSUFFIX="dylib"
elif [[ "$OSTYPE" == "CYGWIN"* ]]; then
  LIBSUFFIX="dll"
elif [[ "$OSTYPE" == "MINGW"* ]]; then
  LIBSUFFIX="dll"
elif [[ "$OSTYPE" == "msys"* ]]; then
  LIBSUFFIX="dll"
  GENERATOR="MSYS Makefiles"
fi
git submodule update --init --recursive
################################################################################
#                                                                              #
#                                  Build OPI                                   #
#                                                                              #
################################################################################
if [[ "$ENABLE_OPI_SUPPORT" == "ON" ]]; then
  # Build OPI
  cd OPI || exit || exit
  # Create the build directory if it does not exist
  if [[ ! -d "build" ]]; then
    mkdir build
  else
    rm -rf build/*
  fi
  cd build || exit || exit
  echo "Updating cmake"
  cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_INSTALL_PREFIX=../ -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_CXX11=ON -DENABLE_CL_SUPPORT=OFF -DENABLE_CUDA_SUPPORT=OFF -DENABLE_PYTHON=OFF -G "$GENERATOR" ../
  echo "Building OPI"
  cmake --build .
  cmake --install .
  if [[ $? -ne 0 ]]; then
    echo "Could not build OPI. Exiting."
    exit $?
  fi
  echo "Leaving OPI"
  cd ../../
fi
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
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_OpenMP_SUPPORT=OFF -DENABLE_POSTGRESQL_SUPPORT=OFF -G "$GENERATOR" ../
echo "Building libslam"
cmake --build .
cmake --install .
if [[ $? -ne 0 ]]; then
    echo "Could not build libslam. Exiting."
    exit $?
fi
echo "Manually preparing 'lib' and 'include' directories"
cd ../
# Delete the include and lib folders when they exist, 
# as linking on windows will create directories which 
# leads to errors on the second run of the script
if [[ ! -d "include" ]]; then
  rm -rf include
fi
if [[ ! -d "lib" ]]; then
  rm -rf lib
fi
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
cmake -DSKIP_MPI=yes -G "$GENERATOR" ../
echo "Building pFUnit"
cmake --build .
cmake --install .
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
# Create the links to libraries needed by NEPTUNE
cd lib || exit || exit
ln -sf ../libslam/lib/libslam-Fortran.$LIBSUFFIX .
if [[ "$ENABLE_OPI_SUPPORT" == "ON" ]]; then
  ln -sf ../OPI/lib/libOPI-Fortran.$LIBSUFFIX .
  ln -sf ../OPI/lib/libOPI.$LIBSUFFIX .
fi
cd ..
# Create the include directory if it does not exist
if [[ ! -d "include" ]]; then
  mkdir include
fi
# Create the links to includes needed by NEPTUNE
cd include || exit
ln -sf ../libslam/include/SLAM .
if [[ "$ENABLE_OPI_SUPPORT" == "ON" ]]; then
  ln -sf ../OPI/include/OPI .
fi
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
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE -DCMAKE_Fortran_COMPILER=$Fortran_COMPILER -DENABLE_OPI_SUPPORT=$ENABLE_OPI_SUPPORT -DSKIP_MSIS_2=ON -G "$GENERATOR" ../
echo "Building NEPTUNE"
cmake --build .
cmake --install .
if [[ $? -ne 0 ]]; then
    echo "Could not build NEPTUNE. Exiting."
    exit $?
fi
echo "Leaving NEPTUNE"
cd ../work || exit || exit
if [[ "$OSTYPE" == "msys"* || "$OSTYPE" == "MINGW"* ]]; then
  ln -sf ../bin/neptune-sa.exe . || exit || exit
  ln -sf ../bin/libneptune.dll . || exit || exit
  ln -sf ../lib/libneptune.dll.a . || exit || exit
  ln -sf ../libslam/lib/libslam-Fortran.dll . || exit || exit
  ln -sf ../libslam/lib/libslam-Fortran.dll.a . || exit || exit
else 
  ln -sf ../bin/neptune-sa . || exit || exit
fi
cd ..
echo "Done"
