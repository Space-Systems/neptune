set(SRC_FILES
    neptuneErrorHandling.f90
    solarsystem.f90
    astroConvert.f90
    satellite.f90
    version.f90
    hwm07.f90
    maneuvers.f90
    neptuneParameters.f90
    gravity.f90
    nrlmsise00Data.f90
    nrlmsise00Class.f90
    JB2008.for
    atmosphere.f90
    thirdbody.f90
    radiation.f90
    tides.f90
    averaging.f90
    correlation.f90
    derivatives.f90
    numint.f90
    neptuneClass.f90
    neptuneOutput.f90
    neptuneClock.f90
    reentry.f90
    rdinp.f90
    neptune.f90
    JB2008module.for
)

# Build the neptune library based on the given source files
add_library(neptune    SHARED ${SRC_FILES})

# Build the neptune executable based on the given source files
add_executable(neptune-sa neptune-sa.f90)

# Build the neptune executable based on the given source files
add_executable(openmp-test-sa openmp-test-sa.f90)

# Dynamic linking against our libraries
target_link_libraries(neptune
#                        ${LIBSPICE}
                        ${LIBSLAM})

# Dynamic linking against our libraries
target_link_libraries(openmp-test-sa
                        neptune
#                        ${LIBSPICE}
                        ${LIBSLAM})

# Dynamic linking against our libraries
target_link_libraries(neptune-sa
                        neptune
#                        ${LIBSPICE}
                        ${LIBSLAM})

# Executables will be built with compiler specific libraries necessary to run them statically linked
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
	set_target_properties(neptune-sa openmp-test-sa PROPERTIES CMAKE_Fortran_FLAGS "-static-libgfortran")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    set_target_properties(neptune-sa openmp-test-sa PROPERTIES CMAKE_Fortran_FLAGS "-static-intel")
endif()

# https://conda.io/docs/user-guide/tasks/build-packages/use-shared-libraries.html
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Make sure we can find the library after installation
    SET_TARGET_PROPERTIES (neptune
      PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                 INSTALL_NAME_DIR "@loader_path/../lib"
    )

else()
    # Make sure we can find the library after installation
    SET_TARGET_PROPERTIES (neptune
      PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                 INSTALL_NAME_DIR "\$ORIGIN/../lib"
    )

    # Make the binaries look for the libraries in the right place
    SET_TARGET_PROPERTIES(neptune-sa openmp-test-sa
      PROPERTIES INSTALL_RPATH "\$ORIGIN/../lib")
endif()

### Install Fortran binaries and library
# Install the libraries into the lib directory
 install(TARGETS neptune-sa openmp-test-sa neptune
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

if(ENABLE_OPI_SUPPORT)
  MESSAGE ("Creating OPI plugin")
  # Build the neptune library based on the given source files
  add_library(neptune-opi SHARED propagatorschnittstelle.f90)

  find_library(LIBOPI OPI ../lib)
  # Dynamic linking against our libraries
  target_link_libraries(neptune-opi
                          neptune
                        ${LIBOPI}
                        ${LIBOPIFORTRAN}
#                        ${LIBSPICE}
                        ${LIBSLAM})

  if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # Make sure we can find the library after installation
    SET_TARGET_PROPERTIES (neptune-opi
      PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                 INSTALL_NAME_DIR "@loader_path/../lib"
    )
    #SET_TARGET_PROPERTIES(neptune-opi
    #  PROPERTIES INSTALL_RPATH "@loader_path/../lib")
  else()
    # Make sure we can find the library after installation
    SET_TARGET_PROPERTIES (neptune-opi
      PROPERTIES BUILD_WITH_INSTALL_RPATH 1
                 INSTALL_NAME_DIR "\$ORIGIN/../lib"
    )
  endif()

  ### Install Fortran binaries and library
  # Install the libraries into the lib directory
  install(TARGETS neptune-opi
        RUNTIME DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")

endif(ENABLE_OPI_SUPPORT)


### Install Fortran Mod files
# Install the Fortran module files into the include directory
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include/NEPTUNE)

# Note on OSX using OpenMP together with the intel compiler might require some manual action:
# > install_name_tool -change libiomp5.dylib /opt/intel/lib/libiomp5.dylib neptune-sa
# > install_name_tool -change libgomp.dylib /usr/local/lib/libgomp.dylib neptune-sa
