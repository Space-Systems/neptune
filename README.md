## NEPTUNE
NPI Ephemeris Propagation Tool with Uncertainty Extrapolation is a state of the art numerical orbit propagator. It allows the extrapolation of a state vector and the associated uncertainty forward and backward in time. It regards perturbative forces using
  * gravitational  potential models: EIGEN-GL04C, EGM2008 or EGM96,
  * atmospheric drag based on NRLMSISE-00 model, JB2008 model, or simple power law,
  * 3rd body models (Sun & Moon),
  * solar radiation model,
  * Earth albedo model,
  * solid Earth tides model.

Additionally spacecraft manoeuvres are regarded. It can run as standalone application or interface for integration into your solution via 
  * the [Orbital Propagation Interface (OPI-2015)](https://github.com/ILR/OPI) or 
  * directly using the libraries (*.so/*.dylib).

The propagator is written in Fortran-2008 and HPC ready by enabling OpenMP.

### Building the software
The software can be compiled on 
* Linux, 
* macOS and 
* Windows (using the Windows Subsystem for Linux on Windows 10, cygwin or MinGW)

To build NEPTUNE you will need:
* a fortran compiler (e.g. gfortran >= 8.0 or ifort > 19.0)
* cmake >= 3.12
* the [Orbital Propagation Interface (OPI-2015)](https://github.com/ILR/OPI) (retrieved automatically)
* the support library [libslam](https://github.com/IRAS/libslam) (retrieved automatically)
* doxygen (only for the source code documentation)
* graphviz (only for the source code documentation)
* latex (only for editing the technical documentation)
  * On Ubuntu 18.04 LTS this works: `apt install texlive texlive-lang-german texlive-latex-extra texlive-tubs texlive-science`

To get started execute the build script:
```
$ bash build.sh
```
or do it manually:
* create a build folder within the retrieved directory
* execute cmake to create the make files
* compile
```
$ mkdir build; \
  cd build; \
  cmake -DCMAKE_Fortran_COMPILER=gfortran -DCMAKE_BUILD_TYPE=Debug -DENABLE_OpenMP_SUPPORT=OFF -DENABLE_OPI_SUPPORT=OFF ../; \
  make install
```
This will retrieve the projects [Orbital Propagation Interface (OPI-2015)](https://github.com/ILR/OPI) and [libslam](https://github.com/IRAS/libslam), build and install them and do the same with NEPTUNE itself.
In the bin directory you will find neptune-sa, which is the stand-alone executable.
In addition there are:
* neptune-valsent1b
* openmp-test-sa.

You can run the tests by calling:
```
$ cd build
$ LD_LIBRARY_PATH=../../lib ctest --verbose
```

While the first is a validation executables enabling the comparison of NEPTUNE results against measurements,
the last executable is a test to show that parallel execution using OpenMP works.

#### CMake Options
Enabling OpenMP for parallel propagation:
set  `ENABLE_OpenMP_SUPPORT=ON` 

Enabling the Orbit Propagation Interface:
set  `ENABLE_OPI_SUPPORT=ON` 

Enabling static compilation of dependencies:
set  `ENABLE_STATIC_SUPPORT=ON` 

### Run the standalone version
Execute the stand-alone neptune change into the *work* directory and run
```
$ ../bin/neptune-sa
```
NEPTUNE will execute a default project. The results can be found in the *output* directory of *work*.
The *scripts* folder offers simple python examples for plotting the results. 
The project can be adepted in *input/neptune.inp*. The *data* directory contains necessary data that can be retrieved from external sources.

#### Retrieving up to date EOP and solar and geomagnetic activity data

In the *data* directory you can find the files 
* `eop19620101.txt`,
* `fap_day.dat`,
* `fap_mon.dat` and
* `sw19571001.txt`

The `eop19620101.txt` and `sw19571001.txt` files contain the Earth orientation parameters or solar and geomagnetic activity data, respectively. They can be updated by downloading via [celestrak](http://celestrak.com/SpaceData/).
The files `fap_day.dat` and `fap_mon.dat` are also solar and geomagnetic activity data, but supplied by [ESA](https://static.sdo.esoc.esa.int/SOLMAG/).

### Using the interfaces
When building NEPTUNE a libneptune.so (or libneptune.dylib on macOS) will be created in the *lib* directory, while module files will be copied to the *include* directory. These can be used to configure and run neptune externally.

If you have enabled the OPI, there is also the libneptune-opi.so (or libneptune-opi.dylib on macOS) plugin, which can be called from any OPI host. Example implementations for Java and python hosts can be found in the *work-opi-java* and *work-opi-python* directories, respectively. 

The API of NEPTUNE provides five different functions which allow to set or get
input variables, to perform an initialization, finally run NEPTUNE for a
pre-defined scenario and then, optionally, to retrieve ephemerides for a
pre-defined tabular interval. In order to use the NEPTUNE library, it
has to be compiled and statically linked to the calling process. For
Fortran programs, the NEPTUNE modules *neptuneInput*, *neptuneGet* and
*neptune* have to be made accessible via an *USE* statement. The
individual functions, including the modules which contain them, shall be
described in the following.

#### Function init_neptune

This function performs the initialization of the NEPTUNE library, where
the main tasks are to read all required input and data files. It can be
called after a prior *USE* of the *libneptune* module. An initialization
shall be performed at least once before starting a propagation run.
Another initialization could be required, depending on how input
parameters change after calling the *setNeptuneVar* function. For
example, the geopotential model has been changed with a *setNeptuneVar*
call. This will require to read the appropriate data file, so that an
initialization is obligatory in that case. However, NEPTUNE has been
designed to check internally, whether a re-initialization is required
after some parameter has changed. If this is the case, calling the
*neptune* function will immediately result in an abort and the
appropriate initialization error code will be returned. In an
example is shown for the initialization.

``` {.fortran bgcolor="mintedbg"}
use neptuneClass,           only: Neptune_class
use libneptune,             only: propagate, init_neptune
use orbit_types,            only: state_t
use time,                   only: time_t

integer       :: ierr           ! error flag
type(state_t) :: initialState   ! initial state vector
type(time_t)  :: initialEpoch   ! initial epoch

!** create a neptune instance
neptune = Neptune_class()

!** the initial epoch is only required in AC mode

initialState%mjd = 0.d0

!** the initial state is only required in AC mode,
!   however, NEPTUNE will check the radius vector 
!   for validity (NOT below Earth's surface)

initialState%r = 0.d0
initialState%r = 7.d3

initialState%v = 0.d0

!** now do the initialization..

ierr = init_neptune(neptune, initialState, initialEpoch)

!** alternatively, call without parameters is also possible

ierr = init_neptune(neptune)

! ...
```

The function is called with two arguments, the first being the initial
state vector and the second the initial epoch. Both parameters are
obligatory when called in AC mode, but are optional otherwise.
However, as NEPTUNE always checks the initial state vector for validity
as soon as it is available, the initial radius vector will be checked
for not being below Earth's surface. It is thus recommended to pass the
true initial state vector if it is already available at the
initialization time.

The return value of the function is an error flag, where a value of zero
means that no error occurred, while any other integer value tells about
what kind of error occurred.

#### Function setNeptuneVar

This function allows to set individual input variables in NEPTUNE. It is
a generic function and thus provides different definitions of argument
lists. It can be accessed after an *USE* statement for the
*neptuneInput* module. An example is shown in .

``` {.fortran bgcolor="mintedbg"}
integer :: ierr     ! error flag

ierr = neptune%setNeptuneVar(ARGUMENT-LIST)
```

The return value of the function is an error flag, where a value of zero
means that no error occurred, while any other integer value tells about
what kind of error occurred. A more detailed description of the
argument list and allowed values is given in

#### Function getNeptuneVar

This function allows to get individual variables in NEPTUNE. It is a
generic function and thus provides different definitions of argument
lists. It can be accessed after an *USE* statement for the
*neptuneInput* module. An example is shown in .

``` {.fortran bgcolor="mintedbg"}
<TYPE-SPEC> :: val      ! being assigned the value of the 
                        ! get function

val = neptune%getNeptuneVar(ARGUMENT-LIST)
```

The return value of the function is the value of the requested variable,
which means that it must have the same type and shape as the variable
within NEPTUNE. A more detailed description of the argument list and
returned values is given in

#### Calling propagate

This is the main routine for the propagation of a state vector. 
It needs to have three type definitions available, which are *state\_t*,
*covariance\_t* and *time\_t*. In an example is shown for
performing the propagation, where the type definitions are provided via
three library modules that are also used within NEPTUNE.

``` {.fortran bgcolor="mintedbg"}
use orbit_types, only: state_t, covariance_t
use time,        only: time_t

type(state_t)              :: state_in, state_out
type(covariance_t)         :: covar_in, covar_out
type(time_t), dimension(2) :: epoch
logical                    :: reset

call neptune%propagate(             &    
                        state_in,   &
                        covar_in,   &
                        epoch,      &
                        state_out,  &
                        covar_out,  &
                        reset       &
                    )
```

The argument list contains six variables with the shown types. The
*reset* flag can be used to reset the numerical integration, which is
required, for example, as soon as the satellite configuration changes
(e.g. after a maneuver, where mass has changed). The *reset* flag then
needs to be set to *.TRUE.* and the integrator will initialize for the
new configuration. Thus, a reset will always be necessary, as soon as
the force model experiences a change in its parameters.

#### Function getNeptuneData

This function allows to retrieve ephemerides for intermediate steps
after a call of the propagation routine of NEPTUNE. This functionality
is realized via the *neptuneGet* module, which has to be made accessible
after an *USE* statement for that module. It is then possible, by
providing an array index, to address individual state vectors for
intermediate steps within the propagation span. The index $1$ is for the
initial epoch, the index $2$ would return the state vector for the first
intermediate step, the tabular interval of which can be defined via the
option *OPT\_STORE\_DATA*. If it is set to $300$ seconds, for
example, the index $2$ would provide the state vector after $5$ minutes.
An example is shown in .

``` {.fortran bgcolor="mintedbg"}
type(state_t) :: state 

! perform propagation first
! ...
! ...

! then retrieve data
state = neptune%getNeptuneData(30) ! returns the 30th entry in the table
```

The return value of the function is of derived type *state\_t* and thus
contains the radius and the velocity vector, as well as the epoch for a
given index in the ephemerides array. The reason, why this function is
not part of the module *neptuneGet* but comes along within its own
module, is its complexity. While the *getNeptuneVar* function is quite
simple and returns single values of distinct parameters and variables,
the *getNeptuneData* function provides access to an array being
integrated into the core propagation routine of NEPTUNE.

In this section, the interface shall be described, which allows to set
variables in NEPTUNE. For that purpose, the NEPTUNE API
provides a routine, which is called *setNeptuneVar*. In its argument
list it always receives an identifier *key* first and, depending on the
identifier context, one or more additional variables which allow to set
a specific variable in the NEPTUNE library. Thus, the *setNeptuneVar*
function is generic and NEPTUNE handles which specific function to call
according to the passed arguments of *setNeptuneVar*.

#### Type 1: Set input with key (string) and value (string)

The first possibility is to call *setNeptuneVar* with a *key*/*value*
pair, which are both strings, as shown in . The first line shows the
general definition, while the second line shows exemplarily, how to
switch the atmosphere perturbations on.

``` {.fortran bgcolor="mintedbg"}
integer :: ierr

! ierr = setNeptuneVar(KEY, VAL)
ierr = neptune%setNeptuneVar('ATMOSPHERE', 'ON')
```

The key identifiers, which are allowed, are shown in . The version
column specifies, which NEPTUNE version introduced the key, while in the
*Values* column, the default value for that variable is always followed
by a *(d)*. This means that it is not required to set all variables
explicitly for the propagation, because NEPTUNE will always have a
default configuration for some variables as shown in the table.

 | **Key**                           | **Description**                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | **Values**                            | **Vers.**  |
| :-------------------------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :------------------------------------ | :--------: |
| ALBEDO                            | Switch for the perturbations due to Earth radiation pressure.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | ON (d), OFF                           | 1\.0.0     |
| ATMOSPHERE                        | Switch for the perturbations due to atmospheric drag. If set to OFF, then horizontal wind will not be considered either.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | ON (d), OFF                           | 1\.0.0     |
|   CORRELATION\_ <br> MATRIX       | Switch for considering the auto-correlation matrix within the propagation of the covariance matrix. Only if the propagation of the latter is switched ON, this setting may have an effect.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | ON, OFF (d)                           | 1\.0.0     |
|   COVARIANCE\_ <br> DRAG          | Switch for considering partial derivatives due to the drag perturbations in the covariance matrix propagation.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | ON, OFF (d)                           | 1\.0.0     |
|   COVARIANCE\_ <br> GEOPOTENTIAL  | The degree of the geopotential for the variational equations in the covariance matrix propagation. If set to '0' or '1', the two-body equations will be used. For all integers (*NOTE:* the integers are passed as strings!) higher than '1', the degree of the geopotential for the covariance matrix propagation will be set to that value. E.g. a value '30' would result in a 30\\times30 geopotential. The maximum allowed value is '85' and is further limited by the value which is set by the key GEOPOTENTIAL. The geopotential for the covariance matrix propagation can not have a degree that is higher than the potential used for the state vector propagation.  | 0 (d), ..., 30, ... 85 (max.)  | 1\.0.0     |
|   COVARIANCE\_ <br> MOON          | Switch for considering partial derivatives due to the Moon perturbations in the covariance matrix propagation.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | ON, OFF (d)                           | 1\.0.0     |
|  COVARIANCE\_ <br> PROPAGATION    | Switch for the covariance matrix propagation. If set to OFF, then the covariance matrix will not be propagated.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | ON, OFF (d)                           | 1\.0.0     |
|   COVARIANCE\_ <br> SUN           | Switch for considering partial derivatives due to the Sun perturbations in the covariance matrix propagation.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | ON, OFF (d)                           | 1\.0.0     |
|   COVARIANCE\_ <br> SRP           | Switch for considering partial derivatives due to the solar radiation pressure perturbations in the covariance matrix propagation.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | ON, OFF (d)                           | 1\.0.0     |
|   EPOCH\_END\_ <br> GD            | Propagation end epoch provided in format according to iso 8601 'YYYY-MM-DDThh:mm:ss' with *Y* being the year, *M* being the month, *D* being the day of the month, 'T' being a separator (not substituted!), *h* being the hours, *m* the minutes and *s* the seconds. The time is provided in UTC.                                                                                                                                                                                                                                                                                                                                                                  | 2013-05-01T12:00:00                   | 1\.0.0     |
|   EPOCH\_START\_ <br> GD          | Propagation start epoch provided in format according to iso 8601 'YYYY-MM-DDThh:mm:ss' with *Y* being the year, *M* being the month, *D* being the day of the month, 'T' being a separator (not substituted!), *h* being the hours, *m* the minutes and *s* the seconds. The time is provided in UTC.                                                                                                                                                                                                                                                                                                                                                                | 2013-05-01T12:00:00                   | 1\.0.0     |
| FILE\_ACTGEN                      | File name of the auto-configuration generation configuration file. This file will be searched in the directory as specified via the PATH\_INPUT key. The file will be referenced as soon as a new auto-configuration table shall be generated by NEPTUNE. It contains the orbit region and bin definitions for the discretisation of the regions.                                                                                                                                                                                                                                                                                                                              | dwm07b\_ 104i.dat (d), dwm.dat        | 1\.0.0     |
| FILE\_DWIND                       | File name of the disturbance wind data file (required by hwm07). This file will be searched in the directory as specified via the PATH\_DATA key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | dwm07b\_ 104i.dat (d), dwm.dat        | 1\.0.0     |
| FILE\_EOP                         | File name of the eop data file. This file will be searched in the directory as specified via the PATH\_DATA key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | eop19620101\.txt (d), eop.dat         | 1\.0.0     |
| FILE\_HWIND                       | File name of the horizontal wind data file. This file will be searched in the directory as specified via the PATH\_DATA key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | hwm071308e. dat (d), hwm.dat          | 1\.0.0     |
|   FILE\_ <br> MANEUVERS           | File name of the maneuver specification input file. This file will be searched in the directory as specified via the PATH\_INPUT key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | neptune.mnv (d), manv.inp             | 1\.0.0     |
| FILE\_QDGRID                      | File name of the pre-computed interpolation grid data file for quasi-dipole coordinates (required by hwm07). This file will be searched in the directory as specified via the PATH\_DATA key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | apexgrid.dat (d), qdgrid.dat          | 1\.0.0     |
| FILE\_SOLMAG                      | File name of the daily solar and geomagnetic activity data file. This file will be searched in the directory as specified via the PATH\_DATA key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | fap\_day.dat (d), daily.dat           | 1\.0.0     |
|   FILE\_SOLMAG\_ <br> MONTHLY     | File name of the monthly solar and geomagnetic activity data file. This file will be searched in the directory as specified via the PATH\_DATA key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | fap\_mon.dat (d), monthly.dat         | 1\.0.0     |
|   FILE\_ <br> SURFACES            | File name of the surfaces definition input file. This file will be searched in the directory as specified via the PATH\_INPUT key.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | neptune.srf (d), surf.inp             | 1\.0.0     |
| GEOPOTENTIAL                      | The degree of the geopotential. If set to '0' or '1', a two-body propagation will be performed. For all integers (*NOTE:* the integers are passed as strings!) higher than '1', the degree of the geopotential will be set to that value. E.g. a value '30' would result in a 30\\times30 geopotential. The maximum allowed value is '85'.                                                                                                                                                                                                                                                                                                                                     | 0,..., 30, ... 85 (max.)      | 1\.0.0     |
|   HORIZONTAL\_ <br> WIND          | Switch for the perturbations due to horizontal wind.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | ON (d), OFF                           | 1\.0.0     |
| MANEUVERS                         | Switch for the consideration of orbit maneuvers. This allows also to introduce additional accelerations at any point in the propagation time frame                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | ON, OFF (d)                           | 1\.0.0     |
| MOON                              | Switch for the perturbations due to Moon's gravity.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | ON (d), OFF                           | 1\.0.0     |
| OCEAN\_TIDES                      | Switch for the perturbations due to ocean tides.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | ON (d), OFF                           | 1\.0.0     |
|   OPT\_AP\_ <br> FORECAST         | Constant (integer) value which shall be used for the geomagnetic planetary amplitude A\_p for long-term forecasts.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | 10 (d), 35                            | 1\.0.0     |
| OPT\_EOP                          | Switch for considering eop. If switched OFF, the transformation between inertial and Earth-fixed frame is only based on the GMST.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          | ON (d), OFF                           | 1\.0.0     |
|   OPT\_PN\_ <br> LOOKUP           | Switch for using lookup tables for the precession-nutation theory. It allows to significantly reduce computation time while reducingthe accuracy to the sub-meter regime.   If switched OFF, the relevant quantities X and Y as well as CIO locator s will be determined according to the proper theory.                                                                                                                                                                                                                                                                                                                                                                       | ON, OFF (d)                           | 1\.0.0     |
|   OPT\_GEO\_ <br> MODEL           | Mode switch to select, which geopotential model shall be used. Here, a '1' would result in the egm-96, a '2' in the egm-2008 and a '3' in the eigen-GL04C model.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | 1, 2, 3 (d)                           | 1\.0.0     |
|   OPT\_ <br> HARMONICS            | Switch to enable the analysis of distinct spherical harmonics, the latter being defined individually via the HARMONICS option (seetab:annex:neptune:set:03). If no harmonics are defined via the HARMONICS key, this option will be ignored.                                                                                                                                                                                                                                                                                                                                                                                                                                   | ON, OFF (d)                           | 1\.0.0     |
|   OPT\_INT\_ <br> LOGFILE         | Switch for the generation of the numerical integration logfile. For each call of the integrator a line is dropped into that file containing the current stepsize, the error for that step, etc.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | ON, OFF (d)                           | 1\.0.0     |
|   OPT\_SAT\_ <br> PROPERTIES      | Mode switch for the definition of the satellite's shape. If set to '1', a sphere will be assumed and the cross-section related to drag, the cross-section related to srp as well as the drag and the srp coefficient will have to be provided through subsequent calls of the *setNeptuneVar* routine with the appropriate keys and desired values. If another number is set, the shape definition will be read from an input file (sec:annex-neptune-interfaces-input-srf).                                                                                                                                                                                           | 1 (d), 2,...                      | 1\.0.0     |
|   OPT\_SOL\_ <br> FORECAST        | Constant value which shall be used for the solar activity (10\.7 cm) for long-term forecasts (i.e. if no observed or predicted data is available). Provided in sfu (10^{-22}\\; W/m^2/Hz).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | 80\.0 (d), 123.45                     | 1\.0.0     |
| OPT\_SHADOW                       | Select which Earth shadow model shall be used. In NEPTUNE 1\.0.0 only the conical model is available. Optionally, the shadow can also be switched off.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | CONICAL (d), NONE                     | 1\.0.0     |
|   OPT\_SRP\_ <br> CORRECT         | Correction algorithm based on lundberg1991 for handling srp discontinuities at shadow boundary crossings.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | OFF (d), ON                           | 1\.0.0     |
|   OPT\_STORE\_ <br> DATA          | Allows for the storage of ephemerides in an internal array. The passed integer value gives the step size in the ephemerides table in seconds. If the value is set to zero (default), no ephemerides will be stored.                                                                                                                                                                                                                                                                                                                                                                                                                                                            | 0 (d), 60, 300                        | 1\.0.0     |
| OUTPUT\_FILES                     | General switch for the generation of output files. If switched OFF, no output files will be written, except for the NEPTUNE log and the numerical integration log file the numerical integration logfile.                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACA                       | Switch for the generation of the output file for accelerations due to Earth's radiation pressure.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACC                       | Switch for the generation of the total acceleration output file.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACG                       | Switch for the generation of the output file for accelerations due to Earth's gravity field.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACD                       | Switch for the generation of the output file for accelerations due to atmospheric drag.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACM                       | Switch for the generation of the output file for accelerations due to gravitational perturbations of the Moon.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACS                       | Switch for the generation of the output file for accelerations due to gravitational perturbations of the Sun.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACR                       | Switch for the generation of the output file for accelerations due to solar radiation pressure.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACT                       | Switch for the generation of the output file for accelerations due to solid Earth tides.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ACO                       | Switch for the generation of the output file for accelerations due to ocean tides.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_AMN                       | Switch for the generation of the output file for accelerations due to orbit maneuvers.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_ATM                       | Switch for the generation of the output file for the total atmospheric density. Additionally, also the cross-section (in ) and the   relative velocity (in  for all three components) are provided.                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | ON, OFF (d)                           | 1\.0.0     |
|   OUTPUT\_COV\_ <br> ECI          | Switch for the generation of the output file for the covariances of the state vector (radius and velocity) in the eci frame (gcrf). This file contains the 15 off-diagonal elements making use of the symmetry of the variance/covariance matrix.                                                                                                                                                                                                                                                                                                                                                                                                                      | ON, OFF (d)                           | 1\.0.0     |
|   OUTPUT\_COV\_ <br> UVW          | Switch for the generation of the output file for the covariances of the state vector (radius and velocity) in the UVW frame. This file contains the 15 off-diagonal elements making use of the symmetry of the variance/covariance matrix.                                                                                                                                                                                                                                                                                                                                                                                                                                     | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_CSV                       | Switch for the generation of the output file for the cartesian state vector (radius and velocity).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_GLL                       | Switch for the generation of the output file for geodetic latitude, longitude and altitude.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_OSC                       | Switch for the generation of the output file for osculating kepler elements.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | ON, OFF (d)                           | 1\.0.0     |
| OUTPUT\_STEP                      | Step size for all output files, excluding log files, which has to be provided in seconds.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | 120, 300 (d)                          | 1\.0.0     |
|   OUTPUT\_VAR\_ <br> ECI          | Switch for the generation of the output file for the variances of the state vector (radius and velocity) in the eci frame (gcrf).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      | ON, OFF (d)                           | 1\.0.0     |
|   OUTPUT\_VAR\_ <br> UVW          | Switch for the generation of the output file for the variances of the state vector (radius and velocity) in the UVW frame.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | ON, OFF (d)                           | 1\.0.0     |
| PAR\_CDRAG                        | Set the drag coefficient of the satellite.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     | 2\.2, 2.0                             | 1\.0.0     |
| PAR\_CREFL                        | Set the srp coefficient of the satellite.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  | 1\.3, 1.5                             | 1\.0.0     |
|   PAR\_CROSS\_ <br> SECTION       | Set the cross-section of the satellite in m^2, if a spherical shape is assumed (cannon-ball model).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            | 10\.0, 1.23                           | 1\.0.0     |
|   PAR\_EARTH\_ <br> RADIUS        | Set the Earth's radius (in km) independent of the radius used for the geopotential. While the latter one is selected automatically based on the geopotential model used, it is possible to use this option to set a differing value for the Earth's geocentric radius, which is then used e.g. in the determination of the altitude for the applied atmosphere model. If this parameter is not specified explicitly, the value set by the geopotential will be used throughout the program.                                                                                                                                                                                    | 6378\.0, 6378.135                     | 1\.0.0     |
|   PAR\_INT\_ <br> ABSEPS          | Set the absolute tolerance for the numerical integration. It should be about an order of magnitude lower than the relative tolerance (PAR\_INT\_RELEPS)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | 1E-10, 1E-13 (d)                      | 1\.0.0     |
|   PAR\_INT\_ <br> COV\_STEP       | Set the step size for the numerical integration (rk4) of the covariance matrix in seconds.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 | 60\.0, 120.0                          | 1\.0.0     |
|   PAR\_INT\_ <br> RELEPS          | Set the relative tolerance for the numerical integration. It should be about an order of magnitude higher than the absolute tolerance (PAR\_INT\_ABSEPS)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       | 1E-9, 1E-14 (d)                       | 1\.0.0     |
| PAR\_MASS                         | Set the mass of the satellite in kg.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | 1000\.0, 123.45                       | 1\.0.0     |
| PAR\_REENTRY                      | Set the minimum altitude for decaying satellites. Propagation will be stopped as soon as the altitude is below that value. Has to be provided in km.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           | 50\.0, 10.0 (d)                       | 1\.0.0     |
| PATH\_DATA                        | Path where all data files (the files which are not edited by the user) can be found, e.g. the solar and geomagnetic activity or eop data files. The final delimiter shall be omitted. Maximum length of the path is 512 characters.                                                                                                                                                                                                                                                                                                                                                                                                                                        | data (d), ./test/data                 | 1\.0.0     |
| PATH\_INPUT                       | Path where all input files (the files which are edited by the user) can be found, e.g. the NEPTUNE main input file. The final delimiter shall be omitted. Maximum length of the path is 512 characters.                                                                                                                                                                                                                                                                                                                                                                                                                                                                        | input (d), ./test/input               | 1\.0.0     |
| PATH\_OUTPUT                      | Path where all output files (except for the NEPTUNE and the numerical integration log files, which are written in the same directory where the program is executed) are written to, The final delimiter shall be omitted. Maximum length of the path is 512 characters.                                                                                                                                                                                                                                                                                                                                                                                                        | output (d), ./test/output             | 1\.0.0     |
| RUN\_ID                           | Set the run ID for NEPTUNE. This will result in all output files (except for the numerical integration log file) having the run ID as their file prefix, e.g <RUN-ID>.osc would be the name for the osculating kepler elements output file. The maximum length is 50 characters.                                                                                                                                                                                                                                                                                                                                                                                               | neptune (d), test\_run                | 1\.0.0     |
| SOLID\_TIDES                      | Switch for the perturbations due to solid Earth tides.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         | ON, OFF (d)                           | 1\.0.0     |
| SRP                               | Switch for the perturbations due to solar radiation pressure                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   | ON, OFF (d)                           | 1\.0.0     |
| SUN                               | Switch for the perturbations due to Sun's gravity.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             | ON, OFF (d)                           | 1\.0.0     |

#### Type 2: Set input with key (string) and value (type(kepler\_t))

Another possibility is to call *setNeptuneVar* with a *key* (being a
string, as for Type 1) and a *value*, which is of derived type. In this
section, this type will be *kepler\_t*, which is described in more
detail in . An example is shown in .

``` {.fortran bgcolor="mintedbg"}
use orbit_types, only: kepler_t

type(kepler_t) :: kepel

kepel%sma  = 7000.d0
kepel%ecc  =    0.1d0
kepel%inc  =   98.4d0
kepel%raan =   12.3d0
kepel%aop  =   45.6d0
kepel%tran =   78.9d0

kepel%angles_unit = 'DEG'

ierr = neptune%setNeptuneVar('INITIAL_STATE', kepel)
```

In the initial version (1.0.0), NEPTUNE provides only one key word:
INITIAL\_STATE. Thus, it is possible to set the initial state by
providing osculating Kepler elements. In the example in , six Kepler
elements are defined and then set in NEPTUNE. The angles may be provided
in degrees or radians depending on the *angles\_unit* variable, which
can have the values *DEG* (for degrees) or *RAD* for radians. If that
variable is not explicitly set, NEPTUNE will assume that the angles are
provided in radians.

#### Type 3: Set input with key (string) and value (type(state\_t))

The Type 3 call of *setNeptuneVar* has a value of derived type
*state\_t*. This type is described in more detail in . An example is
shown in .

``` {.fortran bgcolor="mintedbg"}
use orbit_types, only: state_t

type(state_t) :: state_ini

state_ini%r(1) = 1234.56d0
state_ini%r(2) =   78.90d0
state_ini%r(3) = 6789.01d0
state_ini%v(1) =   -6.78d0
state_ini%v(2) =   -2.34d0
state_ini%v(3) =    3.45d0

ierr = neptune%setNeptuneVar('INITIAL_STATE', state_ini)
```

In the initial version (1.0.0), NEPTUNE provides only one key word:
INITIAL\_STATE. Thus, it is possible to set the initial state by
providing a radius and a velocity vector via the derived type. While the
initial state vector is also provided when calling the main routine
*neptune* for the first time with the actual initial state vector,
setting it explicitly is mainly required for the initial state vector to
appear in the headers of the output files, as NEPTUNE does not store the
information, which is passed via the arguments list when called.

The dimension of the radius is in km, while the velocity is in km/s. The
routine will then perform a check, whether the provided state vector
results in the altitude being above the Earth's surface before it sets
the initial state and returns.

#### Type 4: Set input with key (string) and value (type(covariance\_t))

The Type 4 call of *setNeptuneVar* has a value of derived type
*covariance\_t*. This type is described in more detail in . An example
is shown in .

``` {.fortran bgcolor="mintedbg"}
use orbit_types, only: covariance_t

type(covariance_t) :: covariance_ini

covariance_ini%elem = 0.d0          ! initialize all elements 
                                    ! to zero

covariance_ini%elem(1,1) = 1.d-6
covariance_ini%elem(2,2) = 2.d-6
covariance_ini%elem(3,3) = 3.d-6
covariance_ini%elem(4,4) = 4.d-10
covariance_ini%elem(5,5) = 5.d-10
covariance_ini%elem(6,6) = 6.d-10

ierr = neptune%setNeptuneVar('INITIAL_STATE', covariance_ini)
```

In the initial version (1.0.0), NEPTUNE provides only one key word:
INITIAL\_COVARIANCE. Thus, it is possible to set the initial covariance
matrix by providing a $6\times6$ matrix of the derived type. While the
initial covariance matrix is also provided when calling the main routine
*neptune* for the first time with the actual initial covariance matrix,
setting it explicitly is mainly required for the initial covariance
matrix to appear in the headers of the output files, as NEPTUNE does not
store the information, which is passed via the arguments list when
called.

The dimension of the quantities containing radius components is in km,
while the velocity components are in km/s.

#### Type 5: Set input for the geopotential

The Type 5 call of *setNeptuneVar* is a special function, which has been
implemented to provide a means to adapt the geopotential, which, for
example, is useful for the analysis of the properties of the
auto-correlation function of the geopotential (see or ). The argument
list thus has four different values, the first one, again, being the key
(string). An example is shown in .

``` {.fortran bgcolor="mintedbg"}
! Set the J_2 parameter
ierr = neptune%setNeptuneVar('HARMONIC_C', 2, 0, -1.08263d-3)

! Set S_(41,1)
ierr = neptune%setNeptuneVar('HARMONIC_S', 41, 1, -4.13425d-9) 
```

The key identifiers, which are allowed, are shown in , which are the
spherical harmonic coefficients C and S as well as their standard deviations. 
Each coefficient is addressed by two integer indices, the first one being the degree
n and the second one the order 0 <= m <= n of the coefficient.

| **Key**         | **Description**                                                                                                                                                                                                                                                                                                                      | **Value**    | **Vers.**  |
| :-------------- | :----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :----------- | :--------: |
| HARMONIC\_C     | The spherical harmonics C (linked to cosine terms in the geopotential series representation). The first index n may have values from 2 to the set degree of the geopotential, which is checked by NEPTUNE. The second index m shall fulfill the condition 0 <= m <= n.          | -1\.08263d-3 | 1\.0.0     |
| HARMONIC\_S     | The spherical harmonics S (linked to sine terms in the geopotential series representation). The first index n may have values from 2 to the set degree of the geopotential, which is checked by NEPTUNE. The second index m shall fulfill the condition 0 <= m <= n.   | -4\.13425d-9 | 1\.0.0     |
| HARMONIC\_SD\_C | The standard deviations of the spherical harmonics C (linked to cosine terms in the geopotential series representation).                                                                                                                                                                                                    | -5\.4321d-11 | 1\.0.0     |
| HARMONIC\_SD\_S | The standard deviations of the spherical harmonics S (linked to sine terms in the geopotential series representation).                                                                                                                                                                                                      | -9\.8765d-13 | 1\.0.0     |

#### Type 6: Set input for integer arrays

The Type 6 call of *setNeptuneVar* allows to pass integer arrays. The
argument list has two values, the first one, again, being the key
(string), while the second is the integer array. An example is shown in
.

``` {.fortran bgcolor="mintedbg"}
! define array
integer, dimension(10) :: arr

arr(:) = 0
arr(1) = 30
arr(2) = 31

! set the J30 and J31 harmonics only, while all other will be zero...  
ierr = neptune%setNeptuneVar('HARMONICS', arr)
```
| **Key**   | **Description**                                                                                                                                                                                                                                                                                                                                                                                                  | **Value**     | **Vers.**  |
| :-------: | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :-----------: | :--------: |
| HARMONICS | The array contains a list of spherical harmonics to be used exclusively in the geopotential force model. The syntax is as follows: Each entry contains atwo-digit integer, where the first position contains the degree and the second the order of the harmonic under consideration. For instance, an entry '20' wouldcorrespond to using C(2,0) and S(2.0), respectively. | (/20,30,0,0/) | 1\.0.0     |

### Get input parameters {#annex:neptune:get}

In this section, the interface shall be described, which allows to get
variables from . For that purpose, the [acr:api]{acronym-label="acr:api"
acronym-form="singular+short"} provides a routine, which is called
*getNeptuneVar*. In its argument list it always receives an identifier
*key* first and, depending on the identifier context, one or more
additional variables which allow to set a specific variable in the
library. In version 1.0, defines only one function (based on the generic
*getNeptuneVar*, which returns the spherical harmonics of the
geopotential.

This function is similar to Type 5 of the routine *setNeptuneVar* for
setting the spherical harmonics, while now it will return the current
value of each spherical harmonic coefficient. In an example is shown.

``` {.fortran bgcolor="mintedbg"}
use neptuneInput

real*8 :: J_2, S_41_1

! Get the J_2 parameter
J2 = getNeptuneVar('HARMONIC_C', 2, 0)

! Get S_(41,1)
S_41_1 = getNeptuneVar('HARMONIC_S', 41, 1) 
```

The function returns a double precision float value containing the
requested coefficient, as specified via an argument list, as described
in .

### Derived type definitions for NEPTUNE API {#annex:neptune:types}

In order to have the complete functionality of the
[acr:api]{acronym-label="acr:api" acronym-form="singular+short"}
available, one has to provide several type definitions, which are
extensively used within and therefore also in the interface. All derived
types, which were designed for being used in combination with follow a
common naming scheme, where the name of each type is always followed by
an underscore and the letter *t* ($<$NAME$>$\_t). The different types
may require some additional information, which will be explained in
detail in the following paragraphs.

#### Time type

The first derived type, which means that it does not depend on any other
derived types, is a type which provides a means to store a calendar
date. In the definition, as shown in , a gregorian date and the julian
day are combined.

``` {.fortran bgcolor="mintedbg"}
  type time_t        

    integer  :: year
    integer  :: month
    integer  :: day
    integer  :: hour
    integer  :: minute

    real(dp) :: second
    real(dp) :: mjd        ! MJD associated with GD
    real(dp) :: jd         !  JD associated with GD

  end type time_t
```

It is thus possible to provide a gregorian date together with the
associated julian day. For computations within , typically the
[acr:mjd]{acronym-label="acr:mjd" acronym-form="singular+short"} is
used, which was therefore also included in the type definition. The kind
specification for the floats, which is given with *(dp)*, meaning
*double precision*, has to be provided via the *types* module.

#### State vector type

The state vector type provides the radius and velocity vector which, in
combination with the epoch, form the satellite's state vector. The
definition is shown in .

``` {.fortran bgcolor="mintedbg"}
  type state_t

    type(time_t)           :: epoch          ! epoch of state

    integer                :: radius_unit    ! unit of radius
    integer                :: velocity_unit  ! unit of velocity
    real(dp), dimension(3) :: r              ! radius vector
    real(dp), dimension(3) :: v              ! velocity vector

  end type state_t
```

For the epoch, a *time\_t* derived type is required (see previous
paragraph). The units of the radius and velocity vectors are stored as
integer ID's, which can be referenced using the *units* module.

#### Covariance matrix type

The covariance matrix type provides the complete information of the
$6\times6$ covariance matrix for the radius and the velocity vector, as
shown in .

``` {.fortran bgcolor="mintedbg"}
 type covariance_t

    type(time_t) :: epoch

    integer                  :: ref_frame  ! reference frame
    integer,  dimension(6,6) :: unit       ! covar. matrix units
    real(dp), dimension(6,6) :: elem       ! covar. matrix elements

  end type covariance_t
```

As the covariance matrix is symmetric, it would not be necessary to
store 36 elements, however, the implementation with 15 redundant
elements simplifies the computations and source code readability. For
the epoch, a *time\_t* derived type is required (see corresponding
paragraph in this section). The unit of each element is stored as an
integer ID, which can be referenced using the *units* module. An
additional ID is used for the reference frame, e.g. an inertial frame,
or the UVW frame.

#### Kepler elements type

The classical Kepler elements can be provided by using the Kepler
elements type. As shown in , a total of 14 variables is defined within
this type.

``` {.fortran bgcolor="mintedbg"}
  type kepler_t

    type(time_t) :: epoch  ! epoch of kepler elements 

    real(dp) :: sma        ! semi-major axis
    real(dp) :: ecc        ! eccentricity
    real(dp) :: inc        ! inclination
    real(dp) :: raan       ! right ascension of ascending node
    real(dp) :: aop        ! argument of pericenter
    real(dp) :: man        ! mean anomaly
    real(dp) :: ecan       ! eccentric anomaly
    real(dp) :: tran       ! true anomaly

    real(dp) :: arglat     ! argument of true latitude
    real(dp) :: lonper     ! longitude of perigee
    real(dp) :: truelon    ! true longitude

    integer :: sma_unit    ! dimension of semi-major axis  
    integer :: angles_unit ! dimension of all angles

  end type kepler_t
```

The epoch, as for the state and covariance type, is of derived type
*time\_t*. The six classical Kepler elements are augmented by the
eccentric and the mean anomaly. For circular orbits, the argument of
true latitude (positive in flight direction from the ascending node) is
available, as the argument of pericenter is not defined. For circular
equatorial orbits, the true longitude (positive in flight direction from
the Greenwich meridian) is available, due to the line of nodes and the
argument of pericenter not being defined, and, finally, for eccentric
equatorial orbits, the longitude of the pericenter (positive eastward
from the vernal equinox) is defined, because there is no line of nodes.
The unit of the semi-major axis, as well as the units of all angles are
stored in form of two integer ID's, which can be referenced using the
*units* module.
