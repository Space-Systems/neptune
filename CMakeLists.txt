cmake_minimum_required (VERSION 3.12)

# Name the project NEPTUNE
project (NEPTUNE Fortran)

OPTION(BUILD_DOCUMENTATION "Build the doxygen documentation (default=OFF)" OFF)
OPTION(ENABLE_OpenMP_SUPPORT "Enable OpenMP support (default=OFF)" OFF)
OPTION(ENABLE_OPI_SUPPORT "Enable OPI support (default=OFF)" OFF)
OPTION(ENABLE_STATIC_SUPPORT "Enable static library support (default=OFF)" OFF)

SET(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: None, Debug, Release.")

MESSAGE ("Build Type: " ${CMAKE_BUILD_TYPE})
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DDEBUG_BUILD)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")

endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
	set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -O0 -cpp -DDBLE -g -fbounds-check -fbacktrace -ffree-line-length-none -finit-local-zero  -Wall")
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -DDBLE -ffree-line-length-none -cpp")
elseif("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    # Possible solution when segfaults arise with OpenMP is using the following flag: -heap-arrays 32768
    # Reason: some dynamically allocated arrays are replicated when OpenMP is used with the neptuneClass.
    # The (small) stack is used for the replication of the arrays. The flag moves arrays from the stack to the heap.
	#set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -heap-arrays 32768")
	set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -fpp -DDBLE -g -fno-omit-frame-pointer -traceback -check bounds -check output-conversion -check format -check pointers -check uninit -ftrapuv -assume realloc_lhs -fstack-protector -assume protect_parens")
	set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -O3 -fpp -DDBLE")
elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mpreprocess")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -g -O -DDBLE")
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -O3 -DDBLE")
endif()

if (ENABLE_STATIC_SUPPORT)
	if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
	  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -static-libgfortran")
	  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -static-libgfortran")
	elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
	  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -static-intel")
	  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -static-intel")
	elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Flang")
	  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -static-flang-libs")
	  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -static-flang-libs")
	endif()
endif()

if (ENABLE_OpenMP_SUPPORT)
	FIND_PACKAGE(OpenMP)
	MESSAGE ("Using the following OpenMP flags: " ${OpenMP_Fortran_FLAGS})
	if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${OpenMP_Fortran_FLAGS}")
	  MESSAGE ("Fortran Debug flags: " ${CMAKE_Fortran_FLAGS_DEBUG})
	elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
	  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${OpenMP_Fortran_FLAGS}")
	  MESSAGE ("Fortran Release flags: " ${CMAKE_Fortran_FLAGS_RELEASE})
	else()
	  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
	  MESSAGE ("Fortran flags: " ${CMAKE_Fortran_FLAGS})
	endif()
endif()

# Set installation directories
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

# Set path for FORTRAN mod files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod_files)

# SLAM module files are located in this directory:
if(INC_DIR_SLAM)
    include_directories(${INC_DIR_SLAM})
else()
    include_directories(include/SLAM)
endif()
if(INC_DIR_OPI)
    include_directories(${INC_DIR_OPI})
else()
    include_directories(include/OPI)
endif()
if(INC_DIR_OPI_FORTRAN)
    include_directories(${INC_DIR_OPI_FORTRAN})
else()
    include_directories(include/OPI/fortran)
endif()

# Find the slam and spice
find_library(LIBSLAM slam-Fortran lib)
#find_library(LIBSPICE spice lib)
find_library(LIBOPIFORTRAN OPI-Fortran lib)

# Build NEPTUNE
add_subdirectory(src)

# Build validation executables
add_subdirectory(validation)

# Create test scenarios
find_package(PFUNIT REQUIRED)
enable_testing()
#add_test(NAME execute_neptune COMMAND ../bin/neptune-sa WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/work/)
if (PFUNIT_FOUND)
	add_subdirectory(tests)
endif()


if (BUILD_DOCUMENTATION)
	find_package(Doxygen REQUIRED dot)

	if (DOXYGEN_FOUND)
		# note the option ALL which allows to build the docs together with the application
		add_custom_target( doc_doxygen ALL
			COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/documentation/02-Doxygen/nepdoc.cfg
			WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/documentation/02-Doxygen/
			COMMENT "Generating API documentation with Doxygen"
			VERBATIM )
	else (DOXYGEN_FOUND)
	message("Doxygen needs to be installed to generate the source code documentation")
	endif (DOXYGEN_FOUND)
endif (BUILD_DOCUMENTATION)