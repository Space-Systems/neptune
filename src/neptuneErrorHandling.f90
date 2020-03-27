!>----------------------------------------------------------------------------------------
!!
!> @brief       NEPTUNE-specific error handling subroutines and definition of error codes
!!
!> @anchor      neptune_error_handling
!!
!> @author      Christopher Kebschull
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 13.11.2015 (adopted error codes to be complementary to slam's error handling)</li>
!!                <li> 14.05.2016 (added NEPTUNE-specific error messages - which were missing after libslam implementation)</li>
!!                <li> 14.06.2016 (removed autoconfiguration error messages)</li>
!!              </ul>
!!
!> @details      This module contains specific error handling codes for neptune.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @ancor       neptune_error_handling
!!
!> @todo        Add 'setNeptuneVar("ERROR_HANDLING"...' initialization
!!
!--------------------------------------------------------------------------------------------
module neptune_error_handling

  use slam_error_handling

  implicit none

  !------------------------------------------------------------------
  !
  ! Error codes
  !
  !------------------------------------------------------------------

  !** XML parsing errors (codes 100)
  integer, parameter, public :: E_BODY        = 1101      !< <body> tag missing in XML data
  integer, parameter, public :: E_CODING      = 1102      !< decoding error
  integer, parameter, public :: E_COVARIANCE  = 1103      !< covariance data processing not possible as data is missing or corrupt
  integer, parameter, public :: E_DATA        = 1104      !< <data> tag missing in XML data
  integer, parameter, public :: E_STATE       = 1105      !< state processing not possible as data is missing
  integer, parameter, public :: E_SAT_META    = 1106      !< Satellite meta data element missing in XML file
  integer, parameter, public :: E_SEGMENT         = 1108  !< Data missing within segment struct
  integer, parameter, public :: E_OEM_UNSUPPORTED = 1109  !< Unsupported OEM file
  integer, parameter, public :: E_OMM_UNSUPPORTED = 1110  !< Unsupported OMM file
  integer, parameter, public :: E_OEM_READ    = 1111      !< OEM read error, no tree available
  integer, parameter, public :: E_OEM_TYPE    = 1112      !< OEM type unknown
  integer, parameter, public :: E_OEM_VERSION = 1113      !< Unsupported OEM version
  integer, parameter, public :: E_OMM_READ    = 1114      !< OMM read error, no tree available
  integer, parameter, public :: E_OMM_TYPE    = 1115      !< OMM type unknown
  integer, parameter, public :: E_OMM_VERSION = 1116      !< Unsupported OMM version
  integer, parameter, public :: E_TIME_STAMP  = 1117      !< Time stamp (creation date) missing
  integer, parameter, public :: E_ORIGINATOR  = 1118      !< Originator missing
  integer, parameter, public :: E_METHOD      = 1119      !< Propagation method missing
  integer, parameter, public :: E_OBJECT_NAME = 1120      !< Object name missing
  integer, parameter, public :: E_INVALID_ID  = 1121      !< Invalid object ID
  integer, parameter, public :: E_MISSING_ID  = 1122      !< Object ID missing
  integer, parameter, public :: E_INVALID_CAT_NUM       = 1123  !< Invalid catalog number
  integer, parameter, public :: E_MISSING_CAT_NUM       = 1124  !< Object catalog number missing
  integer, parameter, public :: E_FRAME_CENTER_MISSING  = 1125  !< Reference frame center missing
  integer, parameter, public :: E_REF_FRAME_MISSING     = 1126  !< Reference frame missing
  integer, parameter, public :: E_REF_FRAME_EPOCH_FORMAT = 1127 !< Reference frame epoch format error
  integer, parameter, public :: E_TIME_TYPE_MISSING     = 1128  !< Time type missing (e.g. 'JD' or 'DATETIME')
  integer, parameter, public :: E_MANEUVER_EPOCH        = 1129  !< maneuver epoch missing
  integer, parameter, public :: E_MANEUVER_VECTOR       = 1130  !< maneuver vector component missing
  integer, parameter, public :: E_EPHEM_ENVELOPE        = 1131  !< ephemerides envelope flag missing
  integer, parameter, public :: E_EPHEM_EPOCH           = 1132  !< ephemerides epoch missing
  integer, parameter, public :: E_EPHEM_MISSING         = 1133  !< missing ephemerides data
  integer, parameter, public :: E_COVARIANCE_MISSING    = 1134  !< missing covariance data
  integer, parameter, public :: E_STATE_MISSING         = 1135  !< state vector missing
  integer, parameter, public :: E_SMA_MISSING           = 1136  !< semi-major axis missing
  integer, parameter, public :: E_ECC_MISSING           = 1137  !< eccentricity missing
  integer, parameter, public :: E_INC_MISSING           = 1138  !< inclination missing
  integer, parameter, public :: E_RAAN_MISSING          = 1139  !< RAAN missing
  integer, parameter, public :: E_AOP_MISSING           = 1140  !< AoP missing
  integer, parameter, public :: E_MEAN_ANOMALY_MISSING  = 1141  !< Mean anomaly missing
  integer, parameter, public :: E_OMM_MEAN_ELS_MISSING  = 1142  !< Mean elements missing
  integer, parameter, public :: E_XML_BUFFER            = 1152  !< XML buffer error
  integer, parameter, public :: E_XML_PARSING           = 1153  !< XML parsing error
  integer, parameter, public :: E_OEM_ENCODING_VERSION  = 1154  !< OEM encoding not possible due to version conflict
  integer, parameter, public :: E_OEM_DECODING_VERSION  = 1155  !< OEM decoding not possible due to version conflict
  integer, parameter, public :: E_OEM_DECODING_ID       = 1156  !< OEM decoding not possible due to wrong file format
  integer, parameter, public :: E_OMM_ENCODING_VERSION  = 1157  !< OMM encoding not possible due to version conflict
  integer, parameter, public :: E_OMM_DECODING_VERSION  = 1158  !< OMM decoding not possible due to version conflict
  integer, parameter, public :: E_OMM_DECODING_ID       = 1159  !< OMM decoding not possible due to wrong file format
  integer, parameter, public :: E_UNSUPPORTED_FRAME     = 1160  !< Unsupported reference frame

  !** perturbation methods errors (codes 200)
  integer, parameter, public :: E_GEO_INIT             = 1201  !< geopotential coefficients not initialized
  integer, parameter, public :: E_GEO_MAX_DEGREE       = 1202  !< requested degree of geopotential too high
  integer, parameter, public :: E_GRAVITY_COEFFICIENTS = 1203  !< gravity coefficients are not available

  integer, parameter, public :: E_ATMOSPHERE_INIT      = 1204  !< atmospheric drag routine not initialized
  integer, parameter, public :: E_RADIATION_INIT       = 1205  !< radiation routine not initialized
  integer, parameter, public :: E_MIN_ALTITUDE         = 1206  !< altitude less than lower boundary
  integer, parameter, public :: E_THIRD_BODY_SERIES    = 1207  !< selected planetary ephemerides series not supported
  integer, parameter, public :: E_SOLAR_SYSTEM_INIT    = 1208  !< solar system bodies not initialized
  integer, parameter, public :: E_WRONG_THIRD_BODY     = 1209  !< selected planet/moon is not available
  integer, parameter, public :: E_MANEUVER_INIT        = 1210  !< maneuver handling module not initialized
  integer, parameter, public :: E_SRP_CORRECTION       = 1211  !< something went wrong in the penumbra crossing computation while doing SRP correction
  integer, parameter, public :: E_CHEB_COEFF_MISSING   = 1212  !< Chebyshev polynomial coefficients are missing

  !** data file errors
  integer, parameter, public :: E_SGA_INPUT               = 1401  !< SGA initialization error due to wrong input
  integer, parameter, public :: E_SGA_NO_VERSION          = 1402  !< SGA data file version not supported
  integer, parameter, public :: E_SGA_MISSING             = 1403  !< SGA data missing for given epoch
  integer, parameter, public :: E_SGA_UNSUPPORTED_TYPE    = 1404  !< SGA data file type not supported
  integer, parameter, public :: E_SGA_UNSUPPORTED_VERSION = 1405  !< SGA data file version not supported
  integer, parameter, public :: E_SGA_TIME_TAG            = 1406  !< SGA time tag not found

  !** numerical errors
  integer, parameter, public :: E_REL_TOLERANCE        = 1601  !< relative tolerance less than machine epsilon
  integer, parameter, public :: E_ABS_TOLERANCE        = 1602  !< absolute tolerance less than machine epsilon
  integer, parameter, public :: E_INTEGRATION_STARTUP  = 1603  !< initialisation of integrator failed
  integer, parameter, public :: E_INTEGRATION_METHOD   = 1604  !< integration method undefined
  integer, parameter, public :: E_INTEGRATION_ABORT    = 1605  !< integrator not able to take next step after a number of tries
  integer, parameter, public :: E_INTEGRATION_RESET    = 1606  !< reset of integrator has to be performed first

  !** NEPTUNE specific
  integer, parameter, public :: E_NEPTUNE_INIT         = 1700  !< NEPTUNE has not been initialised prior to call
  integer, parameter, public :: E_MEAN_ELS_MISSING     = 1701  !< Mean elements have not been computed yet
  integer, parameter, public :: E_CORR_INIT            = 1702  !< Correlation matrix update not performed yet
  integer, parameter, public :: E_CORR_INTERPOLATION   = 1703  !< Correlation matrix interpolation error
  integer, parameter, public :: E_CORR_SINGULAR        = 1704  !< Correlation matrix computation not possible for equatorial orbits
  integer, parameter, public :: E_MAX_MANEUVERS        = 1705  !< Maximum number of possible maneuvers exceeded in definition file
  integer, parameter, public :: E_MANEUVER_DEFINITION  = 1706  !< Error in maneuver definition: Starting time of a new maneuver < end time of previous maneuver
  integer, parameter, public :: E_MAX_SURFACES         = 1707  !< Maximum number of possible satellite surfaces exceeded
  integer, parameter, public :: E_NO_MANEUVERS         = 1708  !< No maneuver data given in data file
  integer, parameter, public :: E_NEGATIVE_MASS        = 1709  !< Mass of object becomes less than zero kg
  integer, parameter, public :: E_PERTURBATION_SWITCH  = 1710  !< Perturbation switch access out of array bounds
  integer, parameter, public :: E_INVALID_SATELLITE    = 1711  !< Initialization of satellite geometry failed (e.g. due to surfaces array not allocated, or negative mass)
  integer, parameter, public :: E_INVALID_STEP         = 1712  !< Invalid step size for propagation output
  integer, parameter, public :: E_UNDEFINED_SATELLITE  = 1713  !< Satellite configuration not defined yet.

  !** observations
  integer, parameter, public :: E_UNKNOWN_OBS_TYPE     = 1901  !< unknown observation type

contains

!==============================================================================
!
!>  @brief      Wrapping the setError to pass NEPTUNE-specific messages
!!
!!  @author     Vitali Braun
!!  @date       <ul>
!!                <li>VB: 06.07.2018 (initial implementation)</li>
!!              </ul>
!!
!!  @param[in]  code      Error code
!!  @param[in]  err_type  Error type or level, e.g. FATAL
!!  @param[in]  par       Optional array containing additional message content
!!
!!  @details    This method is required in order to pass NEPTUNE-specific
!!              messages to the setError method of the SLAM library.
!!              If setError is called directly from within NEPTUNE,
!!              it will not be able to determine error messages that are
!!              not in SLAM's namescope.
!!
!--------------------------------------------------------------------------
subroutine setNeptuneError(code, err_type, par)

    implicit none
    integer, intent(in)                                  :: code
    integer, intent(in)                                  :: err_type
    character(len=*), optional, dimension(:), intent(in) :: par

    character(len=512) :: errorMessage

    if(present(par)) then
        call getNeptuneErrorMessage(code, errorMessage, par)
    else
        call getNeptuneErrorMessage(code, errorMessage)
    end if
    call setError(E_SPECIAL, err_type, par=par, errorMessage=errorMessage)
    return

end subroutine setNeptuneError

!==============================================================================
!
!>  @brief      Return error message for a given error code
!!
!!  @author     Vitali Braun
!!  @date       <ul>
!!                <li>VB: 14.05.2016 (initial implementation after separation from libslam)</li>
!!              </ul>
!!
!!  @param[in]   code      Error code
!!  @param[in]   par       Additional parameter array for error description
!!  @param[out]  message   Returned message text
!!
!--------------------------------------------------------------------------
subroutine getNeptuneErrorMessage(code, message, par)

  implicit none
  integer,          intent(in)  :: code
  character(len=*), optional, dimension(:), intent(in) :: par
  character(len=*), intent(out) :: message

  character(len=400)            :: ctemp

  integer :: errorLanguage

  errorLanguage  = getErrorLanguage()

  select case(code)
    case(E_SGA_INPUT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Begin date greater than end date."
      end select

    case(E_SGA_NO_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Version of solar and geomagnetic activity data file could not be determined. The currently "// &
                                               "supported version is '"//trim(par(1))//"'."
      end select

    case(E_SGA_UNSUPPORTED_VERSION)
      select case(errorLanguage)
        case default
          if(index(par(1), "N/A") /= 0) then
            ctemp = ""
          else
            ctemp = "Please provide data file according to version '"//trim(par(1))//"'."
          end if

          write(message(1:len(message)),'(a)') "Unsupported version of solar and geomagnetic data file."//trim(ctemp)

      end select

    case(E_SGA_UNSUPPORTED_TYPE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unsupported type of solar and geomagnetic data file. Please provide supported data file, e.g. "// &
                                               "'"//trim(par(2))//"'."
      end select

    case(E_SGA_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "No solar and geomagnetic data available for given start epoch."
      end select

    case(E_SGA_TIME_TAG)
      select case(errorLanguage)
        case default
          if(index(par(1), "N/A") /= 0) then
            ctemp = ""
          else
            ctemp = " for file '"//trim(par(1))//"'"
          end if

         write(message(1:len(message)),'(a)') "Solar and geomagnetic data file time tag could not be determined"//trim(ctemp)// &
                                               ". Please check for available updates."
      end select

    !-----------------------------
    !
    ! XML parsing errors
    !
    !-----------------------------
    case(E_BODY)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Body tag missing in file '"//trim(par(1))//"'."
      end select

    case(E_CODING)
    case(E_COVARIANCE)
    case(E_DATA)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Data tag missing in file '"//trim(par(1))//"'."
      end select

    case(E_UNSUPPORTED_FRAME)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Unsupported reference frame: '"//trim(par(1))//"'."
      end select

    case(E_STATE)
    case(E_SAT_META)
    case(E_SEGMENT)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "No segment tag found in file '"//trim(par(1))//"'."
      end select


    case(E_OEM_READ)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "XML tree was not initialised... invalid OEM file '"//trim(par(1))//"'."
      end select

    case(E_OMM_READ)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "XML tree was not initialised... invalid OMM file '"//trim(par(1))//"'."
      end select

    case(E_OEM_UNSUPPORTED)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Non-NPI OEM in file '"//trim(par(1))//"' not supported."// &
                                               " Skipping file..."
      end select

   case(E_OMM_UNSUPPORTED)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Non-NPI OMM in file '"//trim(par(1))//"' not supported."// &
                                               " Skipping file..."
      end select

    case(E_OEM_TYPE)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Provided file '"//trim(par(1))//"' is not an OEM"// &
                                               " data file. The tag 'oem' could not be found. File parsing cancelled."
      end select

    case(E_OEM_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unsupported version of OEM: '"//trim(par(1))//"'."

      end select

    case(E_OMM_TYPE)
      select case(errorLanguage)
        case default  ! english as default
          write(message(1:len(message)),'(a)') "Provided file '"//trim(par(1))//"' is not an OMM"// &
                                               " data file. The tag 'omm' could not be found. File parsing cancelled."
      end select

   case(E_OMM_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unsupported version of OMM: '"//trim(par(1))//"'."
      end select

    case(E_TIME_STAMP)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Time stamp ("//trim(par(1))//") missing in file '"//trim(par(2))//"'."
      end select

    case(E_ORIGINATOR)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Originator missing in file '"//trim(par(1))//"'."
      end select

    case(E_METHOD)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Propagation method missing in file '"//trim(par(1))//"'."
      end select

    case(E_OBJECT_NAME)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Object name could not be found in segment no. '"//trim(par(1))// &
                                               "' in file '"//trim(par(2))//"'. Set to '"//trim(par(3))//"'."
      end select

    case(E_INVALID_ID)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Invalid international designator in file '"//trim(par(1))// &
                                               "'. ID set to '"//trim(par(2))//"'."
      end select

    case(E_MISSING_ID)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Missing international designator in segment no. '"//trim(par(1))// &
                                               "' in file '"//trim(par(2))//"'. Set to '"//trim(par(3))//"'."

      end select

    case(E_INVALID_CAT_NUM)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Invalid catalog number in file '"//trim(par(1))// &
                                               "'. Set to '"//trim(par(2))//"'."
      end select

    case(E_MISSING_CAT_NUM)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Catalog number could not be found in segment no. '"//trim(par(1))// &
                                               "' in file '"//trim(par(2))//"'. Set to '"//trim(par(3))//"'."
      end select

    case(E_REF_FRAME_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Reference frame could not be found for segment no. "//trim(par(1))// &
                                               " in file '"//trim(par(2))//"'. Setting to default: '"//trim(par(3))//"'."
      end select

    case(E_REF_FRAME_EPOCH_FORMAT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unknown reference frame epoch format for segment no. "//trim(par(1))// &
                                               " in file '"//trim(par(2))//"'."
      end select

    case(E_TIME_TYPE_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Missing time type in segment no. '"//trim(par(1))//"' in file '"//trim(par(2))//"'."
      end select

    case(E_MANEUVER_EPOCH)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Epoch for maneuver no. "//trim(par(1))//" in segment no. "//trim(par(2))// &
                                               " could not be found in file '"//trim(par(3))//"'."
      end select

    case(E_MANEUVER_VECTOR)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Delta-V vector "//trim(par(1))//" for maneuver no. "//trim(par(2))// &
                                               " in segment no. "//trim(par(3))//" could not be found in file '"//trim(par(4))//"'."
      end select

    case(E_EPHEM_ENVELOPE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Envelope flag in segment no. "//trim(par(1))//" in data set no. "// &
                                               trim(par(2))//" could not be found in attribute list in file '"//trim(par(3))//"'."
      end select

    case(E_EPHEM_EPOCH)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Epoch for data set no. "//trim(par(1))//" in segment no. "//trim(par(2))// &
                                               " could not be found in file '"//trim(par(3))//"'."
      end select

    case(E_STATE_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "State vector for data set no. "//trim(par(1))//" in segment no. "//trim(par(2))// &
                                               " not found in file '"//trim(par(3))//"'."
      end select

    case(E_SMA_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Semi-major axis not found in file '"//trim(par(1))//"'."
      end select

    case(E_ECC_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Eccentricity not found in file '"//trim(par(1))//"'."
      end select

    case(E_INC_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Inclination not found in file '"//trim(par(1))//"'."
      end select

    case(E_RAAN_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "RAAN not found in file '"//trim(par(1))//"'."
      end select

    case(E_AOP_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Argument of perigee not found in file '"//trim(par(1))//"'."
      end select

    case(E_MEAN_ANOMALY_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Mean anomaly not found in file '"//trim(par(1))//"'."
      end select

    case(E_OMM_MEAN_ELS_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Tag 'meanElements' not found in file '"//trim(par(1))// &
                                               "'. Skipping file."
      end select

    case(E_COVARIANCE_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)')  "Covariance matrix element '"//trim(par(1))//"' in data set no. "//trim(par(2))// &
                                                " in segment no. "//trim(par(3))//" could not be found in file '"//trim(par(4))//"'."
      end select

    case(E_EPHEM_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Ephemeris data in segment no. "//trim(par(1))// &
                                               " not found in file '"//trim(par(2))//"'."
      end select

    case(E_XML_BUFFER)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Processing of XML buffer failed for file '"//trim(par(1))//"'."

      end select

    case(E_XML_PARSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Error while parsing the XML file '"//trim(par(1))//"'."

      end select

    case(E_OEM_ENCODING_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "OEM could not be encoded due to unsupported version: '"//trim(par(1))//"' in file '"// &
                                               trim(par(2))//"'."

      end select

    case(E_OMM_ENCODING_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "OMM could not be encoded due to unsupported version: '"//trim(par(1))//"' in file '"// &
                                               trim(par(2))//"'."

      end select

    case(E_OEM_DECODING_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "OEM could not be decoded due to unsupported version: '"//trim(par(1))//"' in file '"// &
                                               trim(par(2))//"'."

      end select

    case(E_OMM_DECODING_VERSION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "OMM could not be decoded due to unsupported version: '"//trim(par(1))//"' in file '"// &
                                               trim(par(2))//"'."

      end select

    case(E_OEM_DECODING_ID)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unknown file type for passed file '"//trim(par(1))//"'. OEM decoding not possible."
      end select

   case(E_OMM_DECODING_ID)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unknown file type for passed file '"//trim(par(1))//"'. OMM decoding not possible."
      end select

    !-----------------------------
    !
    ! Perturbation methods errors
    !
    !-----------------------------
    case(E_GEO_MAX_DEGREE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Maximum degree for gravity potential is "//trim(par(1))// &
                                               ". Degree set to maximum degree."
      end select

    case(E_GEO_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Geopotential coefficients are not available. Initialize first."
      end select

    case(E_GRAVITY_COEFFICIENTS)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') trim(par(1))//" Coefficients of gravity potential could not be read."
      end select

    case(E_ATMOSPHERE_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Atmospheric perturbations not initialized yet. Initialize first."
      end select

    case(E_MANEUVER_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Maneuver handling module not initialized yet. Initialize first."
      end select

    case(E_RADIATION_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Radiation pressure perturbations not initialized yet. Initialize first."
      end select

    case(E_MIN_ALTITUDE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Object re-entered Earth's atmosphere."
      end select

    case(E_THIRD_BODY_SERIES)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Third body ephemerides series '"//trim(par(1))//"' not supported."
      end select

    case(E_SOLAR_SYSTEM_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Solar system bodies not initialized yet. Initialize first."
      end select

    case(E_WRONG_THIRD_BODY)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Selected third body '"//trim(par(1))//"' not supported."
      end select

   case(E_SRP_CORRECTION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Something went wrong in the penumbra crossing computation while doing SRP correction."
      end select

    case(E_CHEB_COEFF_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Chebyshev polynomials coefficients for solarsystem bodies not available."
      end select

    !-----------------------------
    !
    ! Numerical errors
    !
    !-----------------------------
    case(E_ABS_TOLERANCE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a,e10.4e2)') "Absolute tolerance less than machine epsilon. Reduced to: ", 1.d1*epsilon(1.d0)
      end select

    case(E_REL_TOLERANCE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a,e10.4e2)') "Relative tolerance less than machine epsilon. Reduced to: ", 1.d1*epsilon(1.d0)
      end select

    case(E_INTEGRATION_STARTUP)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unable to take first step after 10 tries - giving up."
      end select

    case(E_INTEGRATION_METHOD)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Unknown integration method. Please set one of the following: '"// &
                                               trim(par(1))//" ("//trim(par(2))//")' or '"// &
                                               trim(par(3))//" ("//trim(par(4))//")'."
      end select

    case(E_INTEGRATION_ABORT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Not able to take next step in integration after "//trim(par(1))//" resets. "// &
                                               "Try to reduce relative/absolute tolerance in numerical integration."
      end select

    case(E_INTEGRATION_RESET)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Integrator reset is performed for first call, although reset flag was not set explicitly."
      end select

    case(E_FRAME_CENTER_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Missing reference frame center for segment no. '"//trim(par(1))// &
                                               " in file '"//trim(par(2))//"'. Set to default: '"//trim(par(3))//"'."
      end select

    !-----------------------------
    !
    ! NEPTUNE errors
    !
    !-----------------------------
    case(E_NEPTUNE_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Please initialize first (call initNeptune) before calling NEPTUNE.'
      end select

    case(E_MEAN_ELS_MISSING)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Mean elements not available. Please compute them first before calling this routine.'
      end select

    case(E_CORR_INIT)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Correlation matrix has not been updated prior to interpolation attempt. Update first.'
      end select

    case(E_CORR_INTERPOLATION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Correlation matrix interpolation attempt failed. Illegal values for array index.'
      end select

    case(E_CORR_SINGULAR)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Correlation matrix can not be done for equatorial orbits. Will be ignored in the following.'
      end select

    case(E_MAX_MANEUVERS)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Maximum number of maneuver phases exceeded: nmax = '//trim(par(1))//'.'
      end select

    case(E_MAX_SURFACES)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Maximum number of surfaces exceeded: nmax = '//trim(par(1))//'.'
      end select

    case(E_INVALID_SATELLITE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Invalid satellite definition, initialization not possible. '//trim(par(1))//'.'
      end select

    case(E_INVALID_STEP)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') 'Invalid step size definition, initialization not possible.'
      end select

    case(E_MANEUVER_DEFINITION)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "At least one specified maneuver starts before the previous maneuver ends."
      end select

    case(E_NO_MANEUVERS)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "No maneuvers specified. Propagation is done without considering maneuvers."
      end select

    case(E_NEGATIVE_MASS)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Mass of object can not be negative."
      end select

    case(E_UNDEFINED_SATELLITE)
      select case(errorLanguage)
        case default
          write(message(1:len(message)),'(a)') "Satellite configuration not defined yet."
      end select

    !----------------------------
    !
    ! Observations
    !
    !----------------------------
    case(E_UNKNOWN_OBS_TYPE)
      write(message(1:len(message)),'(a)') "Unknown observation type '"//trim(par(1))//"'."

    case default
        ! try to obtain error message from libslam
        call getErrorMessage(code, message, par)

  end select
  return

end subroutine getNeptuneErrorMessage


end module neptune_error_handling
