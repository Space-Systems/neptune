!>-----------------------------------------------------------------------------------------------
!> @brief NEPTUNE specific parameters
!!
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li>VB:  09.10.2013 (initial design) </li>
!!                <li>VB:  12.02.2014 (added parameters and options for auto configuration mode) </li>
!!                <li>VB:  01.02.2015 (added atmospheric density output) </li>
!!                <li>VB:  03.02.2015 (added new option for considering only distinct spherical harmonics)</li>
!!                <li>VB:  08.06.2015 (added new option for an additional integrator)</li>
!!                <li>VB:  06.10.2015 (added other solar system bodies)</li>
!!                <li>VB:  15.05.2016 (updated with implicit none statement)</li>
!!                <li>VB:  14.06.2016 (removed autoconfiguration parameters)</li>
!!              </ul>
!!
!> @details     This module contains all parameters (constants) as required by NEPTUNE, e.g. for
!!              input and output handling.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      neptuneParameters
!!
!!------------------------------------------------------------------------------------------------

module neptuneParameters

  implicit none

  integer, parameter :: n_options       = 12   ! number of options
  integer, parameter :: n_parameters    = 12   ! number of propagation parameters

  integer, parameter :: PAR_MASS          = 1  ! object mass
  integer, parameter :: PAR_CROSS_SECTION = 2  ! cross-section
  integer, parameter :: PAR_CDRAG         = 3  ! drag coefficient
  integer, parameter :: PAR_CREFL         = 4  ! reflectivity coefficient
  integer, parameter :: PAR_REENTRY       = 5  ! re-entry altitude stop condition
  integer, parameter :: PAR_INT_RELEPS    = 6  ! num. integration relative tolerance
  integer, parameter :: PAR_INT_METHOD    = 12 ! num. integration method
  integer, parameter :: PAR_INT_ABSEPS    = 7  ! num. integration absolute tolerance
  integer, parameter :: PAR_INT_COV_STEP  = 8  ! covariance matrix integration step
  integer, parameter :: PAR_MAX_ERROR     = 9  ! max. error in auto configuration mode
  integer, parameter :: PAR_MIN_ERROR     = 10 ! min. error in auto configuration mode
  integer, parameter :: PAR_EARTH_RADIUS  = 11 ! Earth's radius / km

  integer, parameter    :: PROP_CONSTANT = 1   ! constant satellite properties
  integer, parameter    :: PROP_FILE     = 2   ! satellite surface definition read from file

  character(len=*), parameter :: C_PATH_DATA      = "PATH_DATA"
  character(len=*), parameter :: C_PATH_INPUT     = "PATH_INPUT"
  character(len=*), parameter :: C_PATH_OUTPUT    = "PATH_OUTPUT"
  character(len=*), parameter :: C_FILE_EOP       = "FILE_EOP"
  character(len=*), parameter :: C_FILE_MANEUVERS = "FILE_MANEUVERS"
  character(len=*), parameter :: C_FILE_SURFACES  = "FILE_SURFACES"
  character(len=*), parameter :: C_FILE_SOLMAG    = "FILE_SOLMAG"
  character(len=*), parameter :: C_FILE_SOLMAG_MONTHLY  = "FILE_SOLMAG_MONTHLY"
  character(len=*), parameter :: C_FILE_DWIND     = "FILE_DWIND"
  character(len=*), parameter :: C_FILE_HWIND     = "FILE_HWIND"
  character(len=*), parameter :: C_FILE_QDGRID    = "FILE_QDGRID"
  character(len=*), parameter :: C_FILE_GEOP      = "FILE_GEOP"
  character(len=*), parameter :: C_FILE_DE_EPHEM  = "FILE_DE_EPHEM"
  character(len=*), parameter :: C_FILE_LEAP_SPICE = "FILE_SPICE_LEAP_SECOND"
  character(len=*), parameter :: C_FILE_TXYS       = "FILE_TXYS"
  character(len=*), parameter :: C_FILE_INPUT_DUMP = "FILE_INPUT_DUMP"
  character(len=*), parameter :: C_FILE_PROGRESS   = "FILE_PROGRESS"

  character(len=*), parameter :: C_RUN_ID       = "RUN_ID"
  character(len=*), parameter :: C_GEOPOTENTIAL = "GEOPOTENTIAL"
  character(len=*), parameter :: C_ATMOSPHERE   = "ATMOSPHERE"
  character(len=*), parameter :: C_WIND         = "HORIZONTAL_WIND"
  character(len=*), parameter :: C_SUN          = "SUN"
  character(len=*), parameter :: C_MOON         = "MOON"
  character(len=*), parameter :: C_MERCURY      = "MERCURY"
  character(len=*), parameter :: C_VENUS        = "VENUS"
  character(len=*), parameter :: C_MARS         = "MARS"
  character(len=*), parameter :: C_JUPITER      = "JUPITER"
  character(len=*), parameter :: C_SATURN       = "SATURN"
  character(len=*), parameter :: C_URANUS       = "URANUS"
  character(len=*), parameter :: C_NEPTUNE      = "NEPTUNE"
  character(len=*), parameter :: C_SRP          = "SRP"
  character(len=*), parameter :: C_ALBEDO       = "ALBEDO"
  character(len=*), parameter :: C_SOLID_TIDES  = "SOLID_TIDES"
  character(len=*), parameter :: C_OCEAN_TIDES  = "OCEAN_TIDES"
  character(len=*), parameter :: C_MANEUVERS    = "MANEUVERS"
  character(len=*), parameter :: C_COV_PROP     = "COVARIANCE_PROPAGATION"
  character(len=*), parameter :: C_COV_GEOPOTENTIAL = "COVARIANCE_GEOPOTENTIAL"
  character(len=*), parameter :: C_COV_DRAG     = "COVARIANCE_DRAG"
  character(len=*), parameter :: C_COV_MOON     = "COVARIANCE_MOON"
  character(len=*), parameter :: C_COV_SUN      = "COVARIANCE_SUN"
  character(len=*), parameter :: C_COV_SRP      = "COVARIANCE_SRP"
  character(len=*), parameter :: C_CORRELATION  = "CORRELATION_MATRIX"

  character(len=*), parameter :: C_HARMONIC_C         = "HARMONIC_C"
  character(len=*), parameter :: C_HARMONIC_SD_C      = "HARMONIC_SD_C"
  character(len=*), parameter :: C_HARMONIC_SD_S      = "HARMONIC_SD_S"
  character(len=*), parameter :: C_HARMONIC_S         = "HARMONIC_S"
  character(len=*), parameter :: C_HARMONICS          = "HARMONICS"
  character(len=*), parameter :: C_INITIAL_COVARIANCE = "INITIAL_COVARIANCE"
  character(len=*), parameter :: C_INITIAL_STATE      = "INITIAL_STATE"
  character(len=*), parameter :: C_OPT_CANONICAL      = "OPT_CANONICAL"
  character(len=*), parameter :: C_OPT_AP_FORECAST    = "OPT_AP_FORECAST"
  character(len=*), parameter :: C_OPT_SOL_FORECAST   = "OPT_SOL_FORECAST"
  character(len=*), parameter :: C_OPT_GEO_MODEL      = "OPT_GEO_MODEL"
  character(len=*), parameter :: C_OPT_ATMOSPHERE_MODEL = "OPT_ATMOSPHERE_MODEL"
  character(len=*), parameter :: C_OPT_SAT_PROPERTIES = "OPT_SAT_PROPERTIES"
  character(len=*), parameter :: C_OPT_EOP            = "OPT_EOP"
  character(len=*), parameter :: C_OPT_HARMONICS      = "OPT_HARMONICS"
  character(len=*), parameter :: C_OPT_PROGRESS       = "OPT_PROGRESS"
  character(len=*), parameter :: C_OPT_PN_LOOKUP      = "OPT_PN_LOOKUP"
  character(len=*), parameter :: C_OPT_INT_LOG        = "OPT_INT_LOGFILE"
  character(len=*), parameter :: C_OPT_SHADOW         = "OPT_SHADOW"
  character(len=*), parameter :: C_OPT_STORE_DATA     = "OPT_STORE_DATA"
  character(len=*), parameter :: C_OPT_SRP_CORRECT    = "OPT_SRP_CORRECT"
  character(len=*), parameter :: C_PAR_MASS           = "PAR_MASS"
  character(len=*), parameter :: C_PAR_CROSS_SECTION  = "PAR_CROSS_SECTION"
  character(len=*), parameter :: C_PAR_CDRAG          = "PAR_CDRAG"
  character(len=*), parameter :: C_PAR_CREFL          = "PAR_CREFL"
  character(len=*), parameter :: C_PAR_EARTH_RADIUS   = "PAR_EARTH_RADIUS"
  character(len=*), parameter :: C_PAR_REENTRY        = "PAR_REENTRY"
  character(len=*), parameter :: C_PAR_INT_RELEPS     = "PAR_INT_RELEPS"
  character(len=*), parameter :: C_PAR_INT_METHOD     = "PAR_INT_METHOD"
  character(len=*), parameter :: C_PAR_INT_COV_METHOD = "PAR_INT_COV_METHOD"
  character(len=*), parameter :: C_PAR_INT_COV_STEP   = "PAR_INT_COV_STEP"
  character(len=*), parameter :: C_PAR_INT_ABSEPS     = "PAR_INT_ABSEPS"
  character(len=*), parameter :: C_EPOCH_START_GD     = "EPOCH_START_GD"
  character(len=*), parameter :: C_EPOCH_END_GD       = "EPOCH_END_GD"
  character(len=*), parameter :: C_INPUT_CARTESIAN    = "Cartesian (r,v)"
  character(len=*), parameter :: C_INPUT_OSCULATING   = "Osculating Kepler Elements"
  character(len=*), parameter :: C_OUTPUT_FILES       = "OUTPUT_FILES"
  character(len=*), parameter :: C_OUTPUT_ACC         = "OUTPUT_ACC"
  character(len=*), parameter :: C_OUTPUT_ACG         = "OUTPUT_ACG"
  character(len=*), parameter :: C_OUTPUT_ACD         = "OUTPUT_ACD"
  character(len=*), parameter :: C_OUTPUT_ACM         = "OUTPUT_ACM"
  character(len=*), parameter :: C_OUTPUT_ACS         = "OUTPUT_ACS"
  character(len=*), parameter :: C_OUTPUT_AME         = "OUTPUT_AME"
  character(len=*), parameter :: C_OUTPUT_ACV         = "OUTPUT_ACV"
  character(len=*), parameter :: C_OUTPUT_AMA         = "OUTPUT_AMA"
  character(len=*), parameter :: C_OUTPUT_ACJ         = "OUTPUT_ACJ"
  character(len=*), parameter :: C_OUTPUT_ASA         = "OUTPUT_ASA"
  character(len=*), parameter :: C_OUTPUT_ACU         = "OUTPUT_ACU"
  character(len=*), parameter :: C_OUTPUT_ACN         = "OUTPUT_ACN"
  character(len=*), parameter :: C_OUTPUT_ACR         = "OUTPUT_ACR"
  character(len=*), parameter :: C_OUTPUT_ACA         = "OUTPUT_ACA"
  character(len=*), parameter :: C_OUTPUT_ACT         = "OUTPUT_ACT"
  character(len=*), parameter :: C_OUTPUT_ACO         = "OUTPUT_ACO"
  character(len=*), parameter :: C_OUTPUT_AMN         = "OUTPUT_AMN"
  character(len=*), parameter :: C_OUTPUT_ATM         = "OUTPUT_ATM"
  character(len=*), parameter :: C_OUTPUT_OSC         = "OUTPUT_OSC"
  character(len=*), parameter :: C_OUTPUT_CSV         = "OUTPUT_CSV"
  character(len=*), parameter :: C_OUTPUT_GLL         = "OUTPUT_GLL"
  character(len=*), parameter :: C_OUTPUT_VAR_ECI     = "OUTPUT_VAR_ECI"
  character(len=*), parameter :: C_OUTPUT_VAR_UVW     = "OUTPUT_VAR_UVW"
  character(len=*), parameter :: C_OUTPUT_COV_ECI     = "OUTPUT_COV_ECI"
  character(len=*), parameter :: C_OUTPUT_COV_UVW     = "OUTPUT_COV_UVW"
  character(len=*), parameter :: C_OUTPUT_STEP        = "OUTPUT_STEP"

  !** switches/options
  !----------------------------------------------------------------


  integer, parameter    :: OPTION_PROPS       = 1    ! Satellite properties: read from input (=1), read from tabular data file (=2)
  integer, parameter    :: OPTION_EOP         = 2    ! Switch index in options array for EOP
  integer, parameter    :: OPTION_GEO_MODEL   = 3    ! Geopotential model to be used
  integer, parameter    :: OPTION_INT_LOGFILE = 4    ! Generate numerical integration logfile
  integer, parameter    :: OPTION_OUTPUT      = 5    ! Write output to files
  integer, parameter    :: OPTION_COV_METHOD  = 6    ! covariance matrix integration method
  integer, parameter    :: OPTION_CORRELATION = 7    ! correlation matrix processing on/off
  integer, parameter    :: OPTION_STORE_DATA  = 8    ! storing ephemeris data in an internal array, which can be accessed via "getNeptuneData" function
  integer, parameter    :: OPTION_SRP_CORRECT = 9    ! Shadow boundary correction for SRP
  integer, parameter    :: OPTION_SHADOW      = 10   ! Shadow model to be used
  integer, parameter    :: OPTION_PN_LOOKUP   = 11   ! Lookup tables are used for precession-nutation theory
  integer, parameter    :: OPTION_HARMONICS   = 12   ! Switch to use only distinct spherical harmonics

  !----------------------------------------------------------------

  !** input parameters
  !----------------------------------------------------------------
  integer, parameter    :: INPUT_UNDEFINED  = -1  ! index for undefined input
  integer, parameter    :: INPUT_CARTESIAN  =  1  ! index for cartesian state vector input
  integer, parameter    :: INPUT_OSCULATING =  2  ! index for osculating state vector input
  integer, parameter    :: INPUT_TEME       =  3  ! index for TEME frame state vector input
  integer, parameter    :: INPUT_COV_GCRF   =  1  ! index for GCRF covariance matrix input
  integer, parameter    :: INPUT_COV_UVW    =  2  ! index for UVW covariance matrix input

contains

end module neptuneParameters


