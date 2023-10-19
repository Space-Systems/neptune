!>-----------------------------------------------------------------------------------------------
!> @brief NEPTUNE output handling module
!!
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  15.09.2013 (initial design - ported from neptune module)  </li>
!!                <li>VB:  06.06.2014 (added function to close all open output files)</li>
!!                <li>VB:  01.02.2015 (added total atmospheric densitiy output file, fixed some bugs, added 'version' module) </li>
!!                <li>VB:  06.10.2015 (added other planets of the solar system)</li>
!!                <li>CHK: 16.11.2015 (updated to use libslam)</li>
!!                <li>VB:  20.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>VB:  14.03.2017 (updated to new input structure in derivatives module)</li>
!!                <li>CHK: 03.01.2017 (created Output_class)</li>
!!              </ul>
!!
!> @details     This is the NEPTUNE output handling module. It is used by NEPTUNE to generate
!!              each of the selected output files during a program run. Also file headers
!!              are created by this module.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      neptuneOutput
!!
!!------------------------------------------------------------------------------------------------

module neptuneOutput

  use atmosphere,             only: Atmosphere_class
  use slam_astro_conversions, only: ELLIPTICAL_INCLINED, CIRCULAR_EQUATORIAL, ELLIPTICAL_EQUATORIAL, CIRCULAR_INCLINED, rv2coe, &
                                    getGeodeticLatLon
  use correlation,            only: Correlation_class
  use derivatives,            only: Derivatives_class, PERT_MOON, PERT_OCEAN_TIDES, PERT_SOLID_TIDES, &
                                    PERT_SRP, PERT_SUN, PERT_MANEUVERS, PERT_DRAG, PERT_ALBEDO, PERT_GRAVITY, &
                                    PERT_WIND, PERT_VENUS, PERT_URANUS, PERT_MERCURY, PERT_SATURN, PERT_NEPTUNE, PERT_MARS, &
                                    PERT_JUPITER
  use slam_error_handling,    only: hasFailed, hasToReturn, isControlled, checkIn, checkOut, E_ID_MATCH, FATAL, E_OUT_OF_BOUNDS
  use neptune_error_handling, only: setNeptuneError
  use gravity,                only: Gravity_class
  use slam_io,                only: SWITCHED_ON, SWITCHED_OFF, cdelimit, SEQUENTIAL, OUT_FORMATTED_OVERWRITE, openFile, C_ON, C_OFF, closeFile
  use maneuvers,              only: Manoeuvres_class
  use slam_math,              only: rad2deg, mag
  use neptuneParameters,      only: PROP_CONSTANT, &
                                    INPUT_UNDEFINED, INPUT_CARTESIAN, INPUT_TEME, INPUT_ITRF, INPUT_ITRF_TEME, INPUT_COV_GCRF, INPUT_OSCULATING, &
                                    INPUT_COV_UVW, C_VENUS, C_URANUS, C_SUN, C_SRP, C_SATURN, C_NEPTUNE, C_RUN_ID, &
                                    C_OUTPUT_VAR_ECI, C_OUTPUT_VAR_UVW, C_OUTPUT_COV_ECI, C_OUTPUT_COV_UVW, C_OUTPUT_ACC, &
                                    C_OUTPUT_ACG, C_OUTPUT_ACD, C_OUTPUT_ACM, C_OUTPUT_ACS, C_OUTPUT_ACR, C_OUTPUT_ACA, C_OUTPUT_ACT, &
                                    C_OUTPUT_ACO, C_OUTPUT_AMN, C_OUTPUT_ATM, C_OUTPUT_GLL, C_OUTPUT_AME, C_OUTPUT_ACV, &
                                    C_OUTPUT_AMA, C_OUTPUT_ACJ, C_OUTPUT_ASA, C_OUTPUT_ACU, C_OUTPUT_ACN, C_OUTPUT_OSC, &
                                    C_OUTPUT_CSV, C_OUTPUT_STEP, &
                                    C_PAR_MASS, C_PAR_CREFL, C_PAR_CDRAG, C_PAR_CROSS_SECTION, &
                                    C_OPT_SAT_PROPERTIES, C_EPOCH_START_GD, C_EPOCH_END_GD, &
                                    C_ALBEDO, C_ATMOSPHERE, C_GEOPOTENTIAL, &
                                    C_INPUT_CARTESIAN, C_INPUT_OSCULATING, C_JUPITER, C_MANEUVERS, C_MARS, C_MERCURY, C_MOON, &
                                    C_SOLID_TIDES, C_OCEAN_TIDES
  use numint,                only: Numint_class
  use slam_orbit_types,       only: state_t, covariance_t, kepler_t, C_RZ, C_RY, C_RX, C_RU, C_RV, C_RW, C_VX, C_VY, C_VZ, C_VW, C_VU, C_VV, C_TAN, &
                                    C_RAN_LONG, C_SMA, C_AOP_LONG, C_ECC, C_INC
  use radiation,              only: Radiation_class
  use slam_reduction_class,   only: Reduction_type
  use slam_rframes,           only: REF_FRAME_GCRF, REF_FRAME_UVW, C_REF_FRAME_UVW, C_REF_FRAME_GCRF
  use satellite,              only: Satellite_class
  use solarsystem,            only: Solarsystem_class, ID_SUN, ID_MOON
  use slam_strings,           only: toLowercase
  use thirdbody,              only: ThirdBody_class
  use tides,                  only: Tides_class, SOLID_TIDES, OCEAN_TIDES
  use slam_time,              only: time_t, mjd2gd, dayFraction2hms, date2string
  use slam_types,             only: dp
  use slam_units,             only: getUnitString, UNIT_RAD, C_UNIT_DEG
  use version,                only: Version_class

  implicit none

  private

    ! all neptune output files are stored via the following struct
    ! ------------------------------------------------------------------------
    integer, parameter :: LEN_PAR_NAME  = 1000                                  !< maximum length of a parameter's name
    integer, parameter :: LEN_FILE_NAME = 1000                                  !< maximum length of output file name
    integer, parameter :: NUM_OUTPUT_FILES = 17                                 !< number of perturbations

    type, public :: neptune_out_t
        character(len=LEN_PAR_NAME)     :: par_name                             !< parameter name
        character(len=LEN_FILE_NAME)    :: file_name                            !< file name
        integer                         :: id                                   !< ID of output file
        integer                         :: file_unit                            !< unit the file is connected to (if it was opened)
        logical                         :: switch                               !< if true, output is written to par_name
        logical                         :: file_name_set                        !< if true, there is a file name associated with the output
    end type

    ! this array contains the set of output files of NEPTUNE
    integer, parameter :: OUTPUT_ARR_DEFAULT_SIZE = 30

    ! file IDs
    integer, parameter    :: OUTPUT_ACC     =  1                                ! index for total acceleration output file
    integer, parameter    :: OUTPUT_ACG     =  2                                ! index for geopotential acceleration output file
    integer, parameter    :: OUTPUT_ACD     =  3                                ! index for drag acceleration output file
    integer, parameter    :: OUTPUT_ACM     =  4                                ! index for Moon acceleration output file
    integer, parameter    :: OUTPUT_ACS     =  5                                ! index for Sun acceleration output file
    integer, parameter    :: OUTPUT_AME     = 19                                ! index for Mercury acceleration output file
    integer, parameter    :: OUTPUT_ACV     = 20                                ! index for Venus acceleration output file
    integer, parameter    :: OUTPUT_AMA     = 21                                ! index for Mars acceleration output file
    integer, parameter    :: OUTPUT_ACJ     = 22                                ! index for Jupiter acceleration output file
    integer, parameter    :: OUTPUT_ASA     = 23                                ! index for Saturn acceleration output file
    integer, parameter    :: OUTPUT_ACN     = 24                                ! index for Uranus acceleration output file
    integer, parameter    :: OUTPUT_ACU     = 25                                ! index for Neptune acceleration output file
    integer, parameter    :: OUTPUT_ACR     =  6                                ! index for solar radiation pressure acceleration output file
    integer, parameter    :: OUTPUT_ACT     =  7                                ! index for solid Earth tides acceleration output file
    integer, parameter    :: OUTPUT_ACO     =  8                                ! index for ocean tides acceleration output file
    integer, parameter    :: OUTPUT_OSC     =  9                                ! index for osculating kepler elements output file
    integer, parameter    :: OUTPUT_CSV     = 10                                ! index for cartesian state vector output file
    integer, parameter    :: OUTPUT_GLL     = 11                                ! index for geocentric longitude and latitude output file
    integer, parameter    :: OUTPUT_VAR_ECI = 12                                ! index for variances (diagonal elements of var/covar matrix) output file in ECI frame
    integer, parameter    :: OUTPUT_VAR_UVW = 13                                ! index for variances (diagonal elements of var/covar matrix) output file in UVW frame
    integer, parameter    :: OUTPUT_COV_ECI = 14                                ! index for covariances (non-diagonal elements) output file in ECI frame
    integer, parameter    :: OUTPUT_COV_UVW = 15                                ! index for covariances (non-diagonal elements) output file in UVW frame
    integer, parameter    :: OUTPUT_AMN     = 16                                ! index for accelerations due to orbital maneuvers
    integer, parameter    :: OUTPUT_ACA     = 17                                ! index for accelerations due to Earth radiation pressure (albedo)
    integer, parameter    :: OUTPUT_ATM     = 18                                ! index for total atmospheric density

  ! ------------------------------------------------------------------------

    type, public :: Output_class

        type(neptune_out_t), allocatable, dimension(:) :: output_arr

        character(len=LEN_FILE_NAME) :: output_path                             !< output path where all output files go to
        character(len=LEN_FILE_NAME) :: run_id                                  !< run ID

        logical :: flag_output                                                  !< indicating that output is requested, no output as default
        logical :: isInitialised                                                !< initialisation flag for output module
        logical :: isSetInitialState                                            !< true as soon as initial state is set

        real(dp) :: cdrag                                                       !< drag coefficient
        real(dp) :: crefl                                                       !< reflectivity coefficient
        real(dp) :: mass                                                        !< mass
        real(dp) :: cross_section                                               !< cross-section

        type(state_t)      :: initial_orbit_csv                                 !< initial orbit written to output
        type(kepler_t)     :: initial_orbit_osc                                 !< initial kepler elements written to output
        type(covariance_t) :: initial_covariance                                !< initial covariance matrix written to output

        type(time_t) :: start_epoch
        type(time_t) :: end_epoch

        integer :: input_type     = INPUT_UNDEFINED                             !< input type for state vector
        integer :: input_type_cov = INPUT_UNDEFINED                             !< input type for covariance matrix
        integer :: satellite_properties                                         !< indicates where satellite properties come from
        integer :: output_step                                                  !< output step size

        character(len=100), dimension(0:14) :: chead                            ! file header

    contains

        procedure :: destroy
        procedure :: close_open_files
        procedure :: get_output_switch
        procedure :: get_neptune_output_file
        procedure :: get_neptune_output_file_number
        procedure :: output
        procedure :: prepare_output
        procedure :: set_input_type
        procedure :: set_input_type_cov
        procedure :: set_output_file
        procedure :: switch_output
        procedure :: write_header
        procedure :: write_state_to_output
        procedure :: write_kepler_elements_to_output
        procedure :: write_covariance_to_output
        procedure :: write_integer_to_output
        procedure :: write_float_to_output
        procedure :: write_string_to_output
        procedure :: write_date_to_output
        generic   :: write_to_output => write_state_to_output, &
                         write_kepler_elements_to_output, &
                         write_covariance_to_output, write_integer_to_output, &
                         write_float_to_output, write_string_to_output, &
                         write_date_to_output

        procedure, private :: init_output
        procedure, private :: exists_output_file
        procedure, private :: get_output_file_id
        procedure, private :: get_default_output_file_name
        procedure, private :: writeInput2Header
    end type Output_class

      ! Constructor
    interface Output_class
        module procedure constructor
    end interface Output_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !!  @author     Christopher Kebschull
    !!  @date       <ul>
    !!                  <li>ChK: 03.01.2018 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Output_class) function constructor()
        constructor%output_path = 'output'                                      !< output path where all output files go to
        constructor%run_id      = 'neptune'                                     !< run ID

        constructor%flag_output       = .false.                                 !< indicating that output is requested, no output as default
        constructor%isInitialised     = .false.                                 !< initialisation flag for output module
        constructor%isSetInitialState = .false.

        allocate(constructor%output_arr(OUTPUT_ARR_DEFAULT_SIZE))
        constructor%output_arr%par_name         = 'None'
        constructor%output_arr%file_name        = ''
        constructor%output_arr(:)%switch        = .false.
        constructor%output_arr(:)%file_name_set = .false.
        constructor%output_arr(:)%id            = 0
        constructor%output_arr(:)%file_unit     = 0

    end function constructor

    !===========================================================================
    !!
    !>  @anchor     destroy
    !!
    !!  @brief      Destroys all that needs destruction
    !!  @author     Christopher Kebschull
    !!
    !!
    !!  @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Output_class)    :: this

        if (allocated(this%output_arr)) deallocate(this%output_arr)

    end subroutine destroy

!==============================================================================================
!
!> @anchor      init_output
!!
!> @brief       Initialises the output module
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li>VB: 21.03.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    subroutine init_output(this)
        implicit none

        class(Output_class)         :: this

        character(len=*), parameter :: routine_id = 'init_output'
        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(routine_id)
        end if

        if(allocated(this%output_arr)) deallocate(this%output_arr)
        allocate(this%output_arr(OUTPUT_ARR_DEFAULT_SIZE))
        this%output_arr%par_name         = 'None'
        this%output_arr%file_name        = ''
        this%output_arr(:)%switch        = .false.
        this%output_arr(:)%file_name_set = .false.
        this%output_arr(:)%id            = 0
        this%output_arr(:)%file_unit     = 0

        this%output_arr(OUTPUT_ACC)%par_name = C_OUTPUT_ACC
        this%output_arr(OUTPUT_ACG)%par_name = C_OUTPUT_ACG
        this%output_arr(OUTPUT_ACD)%par_name = C_OUTPUT_ACD
        this%output_arr(OUTPUT_ACM)%par_name = C_OUTPUT_ACM
        this%output_arr(OUTPUT_ACS)%par_name = C_OUTPUT_ACS
        this%output_arr(OUTPUT_AME)%par_name = C_OUTPUT_AME
        this%output_arr(OUTPUT_ACV)%par_name = C_OUTPUT_ACV
        this%output_arr(OUTPUT_AMA)%par_name = C_OUTPUT_AMA
        this%output_arr(OUTPUT_ACJ)%par_name = C_OUTPUT_ACJ
        this%output_arr(OUTPUT_ASA)%par_name = C_OUTPUT_ASA
        this%output_arr(OUTPUT_ACU)%par_name = C_OUTPUT_ACU
        this%output_arr(OUTPUT_ACN)%par_name = C_OUTPUT_ACN
        this%output_arr(OUTPUT_ACR)%par_name = C_OUTPUT_ACR
        this%output_arr(OUTPUT_ACA)%par_name = C_OUTPUT_ACA
        this%output_arr(OUTPUT_ACT)%par_name = C_OUTPUT_ACT
        this%output_arr(OUTPUT_ACO)%par_name = C_OUTPUT_ACO
        this%output_arr(OUTPUT_AMN)%par_name = C_OUTPUT_AMN
        this%output_arr(OUTPUT_ATM)%par_name = C_OUTPUT_ATM
        this%output_arr(OUTPUT_OSC)%par_name = C_OUTPUT_OSC
        this%output_arr(OUTPUT_CSV)%par_name = C_OUTPUT_CSV
        this%output_arr(OUTPUT_GLL)%par_name = C_OUTPUT_GLL
        this%output_arr(OUTPUT_VAR_ECI)%par_name = C_OUTPUT_VAR_ECI
        this%output_arr(OUTPUT_VAR_UVW)%par_name = C_OUTPUT_VAR_UVW
        this%output_arr(OUTPUT_COV_ECI)%par_name = C_OUTPUT_COV_ECI
        this%output_arr(OUTPUT_COV_UVW)%par_name = C_OUTPUT_COV_UVW

        this%isInitialised = .true.
        if(isControlled()) then
            call checkOut(routine_id)
        end if
        return
    end subroutine

  !=====================================================================
  !
  !> @brief     Return the number of output files available to NEPTUNE
  !!
  !> @author    Vitali Braun
  !!
  !> @return    Number of output files
  !!
  !> @date      <ul>
  !!              <li>24.03.2017 (initial design) </li>
  !!            </ul>
  !!
  !> @anchor    get_neptune_output_file_number
  !-----------------------------------------------
    pure integer function get_neptune_output_file_number(this) result(res)

        implicit none

        class(Output_class),intent(in)  :: this

        res = NUM_OUTPUT_FILES

        return

    end function

  !=====================================================================
  !
  !> @brief     Return the output file type for a given ID
  !!
  !> @author    Vitali Braun
  !!
  !> @param[in] id      ID of output file
  !!
  !> @return    Output file as neptune_out_t
  !!
  !> @date      <ul>
  !!              <li>24.03.2017 (initial design) </li>
  !!            </ul>
  !!
  !! @anchor   get_neptune_output_file
  !-----------------------------------------
    type(neptune_out_t) function get_neptune_output_file(this,id) result(res)

        implicit none

        class(Output_class) :: this
        integer, intent(in) :: id

        character(len=*), parameter :: csubid = "get_neptune_output_file"

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        ! check initialisation status
        if(.not. allocated(this%output_arr) .or. .not. this%isInitialised) then
            call this%init_output()
            if(hasFailed()) return
        end if

        if(id > 0 .and. id <= NUM_OUTPUT_FILES) then
            res = this%output_arr(id)
            if(isControlled()) then
                call checkOut(csubid)
            end if
        else
            call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
        end if
        return

        return

    end function

!==============================================================================================
!
!> @anchor      set_output_file
!!
!> @brief       Sets the output file name and toggles the output switch
!> @author      Vitali Braun
!!
!> @param[in]   par_name    ! output file identifier
!> @param[in]   file_name   ! name of output file (optional)
!> @param[in]   switch      ! if .true., output will be written (optional)
!!
!> @date        <ul>
!!                <li>VB: 21.03.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    subroutine set_output_file(this, par_name, file_name, switch)

        implicit none

        class(Output_class)                     :: this
        character(len=*), intent(in)            :: par_name
        character(len=*), intent(in), optional  :: file_name
        logical,          intent(in), optional  :: switch

        character(len=*), parameter                     :: csubid = 'set_output'
        type(neptune_out_t), dimension(:), allocatable  :: temp_output_arr
        integer :: i
        integer :: idx
        integer :: id
        logical :: new_entry
        logical :: idx_found

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        if(.not. allocated(this%output_arr) .or. .not. this%isInitialised) then
            call this%init_output()
            if(hasFailed()) return
        end if

        ! check if output file is registered
        if(.not. this%exists_output_file(par_name)) then
            call setNeptuneError(E_ID_MATCH, FATAL)
            return
        else
            id = this%get_output_file_id(par_name)
        end if

        ! get the ID of the file in the output array
        idx_found = .false.
        new_entry = .false.
        do i = 1, size(this%output_arr)
            if(trim(this%output_arr(i)%par_name) == 'None') then
                ! add new element to this index, but continue to search as we could still find the parameter in the array
                idx       = i
                idx_found = .true.
                new_entry = .true.
            else if(index(par_name, trim(this%output_arr(i)%par_name)) > 0) then
                ! here we have found the parameter
                idx       = i
                idx_found = .true.
                new_entry = .false. ! it could be true from the above if
                exit
            end if
        end do

        if(.not. idx_found) then    ! this can only happen if the array is already filled - re-allocate and add to the end
            allocate(temp_output_arr(size(this%output_arr)))
            do i = 1, size(this%output_arr)
                temp_output_arr(i)%par_name  = this%output_arr(i)%par_name
                temp_output_arr(i)%file_name = this%output_arr(i)%file_name
                temp_output_arr(i)%file_name_set = this%output_arr(i)%file_name_set
                temp_output_arr(i)%switch    = this%output_arr(i)%switch
                temp_output_arr(i)%file_unit = this%output_arr(i)%file_unit
                temp_output_arr(i)%id        = this%output_arr(i)%id
            end do
            deallocate(this%output_arr)
            allocate(this%output_arr(size(temp_output_arr) + OUTPUT_ARR_DEFAULT_SIZE))
            do i = 1, size(temp_output_arr)
                this%output_arr(i)%par_name  = temp_output_arr(i)%par_name
                this%output_arr(i)%file_name = temp_output_arr(i)%file_name
                this%output_arr(i)%file_name_set = temp_output_arr(i)%file_name_set
                this%output_arr(i)%switch    = temp_output_arr(i)%switch
                this%output_arr(i)%file_unit = temp_output_arr(i)%file_unit
                this%output_arr(i)%id        = temp_output_arr(i)%id
            end do
            ! add the entry at the end
            idx = size(temp_output_arr) + 1
            deallocate(temp_output_arr)
        end if

        ! check switch
        if(present(switch)) then
            this%output_arr(idx)%switch = switch
        end if

        ! if file name is not present and file is switched ON - get the default file name
        if(present(file_name)) then
            this%output_arr(idx)%file_name = file_name(:min(len_trim(file_name),LEN_FILE_NAME))
        else if(.not. this%output_arr(idx)%file_name_set .and. this%output_arr(idx)%switch) then
            this%output_arr(idx)%file_name     = this%get_default_output_file_name(this%output_arr(idx)%par_name)
            this%output_arr(idx)%file_name_set = .true.
        end if

        this%output_arr(idx)%id = id

        !** done!
        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

    end subroutine

!==============================================================================================
!
!> @anchor      exists_output_file
!!
!> @brief       Check if an output file is registered in NEPTUNE
!> @author      Vitali Braun
!!
!> @param[in]   file_id - output file identifier
!!
!> @return      .true. when the output file exists
!!
!> @date        <ul>
!!                <li>VB: 21.03.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    pure logical function exists_output_file(this,file_id) result(res)

        implicit none

        class(Output_class),intent(in)  :: this
        character(len=*), intent(in)    :: file_id

        select case(trim(file_id))
            case(C_OUTPUT_ACC, C_OUTPUT_ACG, C_OUTPUT_ACD, C_OUTPUT_ACM, C_OUTPUT_ACS, &
                 C_OUTPUT_AME, C_OUTPUT_ACV, C_OUTPUT_AMA, C_OUTPUT_ACJ, C_OUTPUT_ASA, &
                 C_OUTPUT_ACU, C_OUTPUT_ACN, C_OUTPUT_ACR, C_OUTPUT_ACA, C_OUTPUT_ACT, &
                 C_OUTPUT_ACO, C_OUTPUT_AMN, C_OUTPUT_ATM, C_OUTPUT_OSC, C_OUTPUT_CSV, &
                 C_OUTPUT_GLL, C_OUTPUT_VAR_ECI, C_OUTPUT_VAR_UVW, C_OUTPUT_COV_ECI, C_OUTPUT_COV_UVW)
                res = .true.
            case default
                res = .false.
        end select
        return

    end function

!==============================================================================================
!
!> @anchor      get_output_file_id
!!
!> @brief       Returns the integer ID of a file given the string identifier
!> @author      Vitali Braun
!!
!> @param[in]   file_id   output file identifier string
!!
!> @return      identifier of the output file
!!
!> @date        <ul>
!!                <li>VB: 21.03.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------

    pure integer function get_output_file_id(this,file_id) result(res)

        implicit none

        class(Output_class),intent(in)  :: this
        character(len=*), intent(in)    :: file_id

        select case(trim(file_id))
            case(C_OUTPUT_ACC)
                res = OUTPUT_ACC
            case(C_OUTPUT_ACG)
                res = OUTPUT_ACG
            case(C_OUTPUT_ACD)
                res = OUTPUT_ACD
            case(C_OUTPUT_ACM)
                res = OUTPUT_ACM
            case(C_OUTPUT_ACS)
                res = OUTPUT_ACS
            case(C_OUTPUT_AME)
                res = OUTPUT_AME
            case(C_OUTPUT_ACV)
                res = OUTPUT_ACV
            case(C_OUTPUT_AMA)
                res = OUTPUT_AMA
            case(C_OUTPUT_ACJ)
                res = OUTPUT_ACJ
            case(C_OUTPUT_ASA)
                res = OUTPUT_ASA
            case(C_OUTPUT_ACU)
                res = OUTPUT_ACU
            case(C_OUTPUT_ACN)
                res = OUTPUT_ACN
            case(C_OUTPUT_ACR)
                res = OUTPUT_ACR
            case(C_OUTPUT_ACA)
                res = OUTPUT_ACA
            case(C_OUTPUT_ACT)
                res = OUTPUT_ACT
            case(C_OUTPUT_ACO)
                res = OUTPUT_ACO
            case(C_OUTPUT_AMN)
                res = OUTPUT_AMN
            case(C_OUTPUT_ATM)
                res = OUTPUT_ATM
            case(C_OUTPUT_OSC)
                res = OUTPUT_OSC
            case(C_OUTPUT_CSV)
                res = OUTPUT_CSV
            case(C_OUTPUT_GLL)
                res = OUTPUT_GLL
            case(C_OUTPUT_VAR_ECI)
                res = OUTPUT_VAR_ECI
            case(C_OUTPUT_VAR_UVW)
                res = OUTPUT_VAR_UVW
            case(C_OUTPUT_COV_ECI)
                res = OUTPUT_COV_ECI
            case(C_OUTPUT_COV_UVW)
                res = OUTPUT_COV_UVW
            case default
                res = 0
        end select
        return

    end function

!==============================================================================================
!
!> @anchor      get_default_output_file_name
!!
!> @brief       Returns the default name of a file - the run ID and the file extension
!> @author      Vitali Braun
!!
!> @param[in]   file_id - output file identifier
!!
!> @return      suffix as string
!!
!> @date        <ul>
!!                <li>VB: 21.03.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    pure character(len=LEN_FILE_NAME) function get_default_output_file_name(this,file_id) result(res)

        implicit none

        class(Output_class),intent(in)  :: this
        character(len=*), intent(in)    :: file_id

        character(len=4) :: suffix  ! file suffix

        select case(trim(file_id))
            case(C_OUTPUT_ACC)
                suffix = '.acc'
            case(C_OUTPUT_ACG)
                suffix = '.acg'
            case(C_OUTPUT_ACD)
                suffix = '.acd'
            case(C_OUTPUT_ACM)
                suffix = '.acm'
            case(C_OUTPUT_ACS)
                suffix = '.acs'
            case(C_OUTPUT_AME)
                suffix = '.ame'
            case(C_OUTPUT_ACV)
                suffix = '.acv'
            case(C_OUTPUT_AMA)
                suffix = '.ama'
            case(C_OUTPUT_ACJ)
                suffix = '.acj'
            case(C_OUTPUT_ASA)
                suffix = '.asa'
            case(C_OUTPUT_ACU)
                suffix = '.acu'
            case(C_OUTPUT_ACN)
                suffix = '.acn'
            case(C_OUTPUT_ACR)
                suffix = '.acr'
            case(C_OUTPUT_ACA)
                suffix = '.aca'
            case(C_OUTPUT_ACT)
                suffix = '.act'
            case(C_OUTPUT_ACO)
                suffix = '.aco'
            case(C_OUTPUT_AMN)
                suffix = '.amn'
            case(C_OUTPUT_ATM)
                suffix = '.atm'
            case(C_OUTPUT_OSC)
                suffix = '.osc'
            case(C_OUTPUT_CSV)
                suffix = '.csv'
            case(C_OUTPUT_GLL)
                suffix = '.gll'
            case(C_OUTPUT_VAR_ECI)
                suffix = '.vre'
            case(C_OUTPUT_VAR_UVW)
                suffix = '.vru'
            case(C_OUTPUT_COV_ECI)
                suffix = '.cve'
            case(C_OUTPUT_COV_UVW)
                suffix = '.cvu'
            case default
                suffix = ''
        end select

        res = trim(this%output_path)//'/'//trim(this%run_id)//trim(suffix)
        return

    end function
  !------------------------------------------------------------------------------------
  !> @brief       Close all NEPTUNE output files
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 06.06.2014 (first implementation)</li>
  !!                <li> 21.03.2017 (adapted to new structure)</li>
  !!              </ul>
  !!
  !> @details     This routine closes all NEPTUNE output files to free the channels.
  !!
  !> @anchor      closeOutputFiles
  !!
  !------------------------------------------------------------------------------------
  subroutine close_open_files(this, numint)

    class(Output_class)                 :: this
    type(Numint_class),intent(inout)    :: numint
    integer :: i    ! loop counter

    do i = 1, size(this%output_arr)
        if(this%output_arr(i)%file_unit > 0) then
            this%output_arr(i)%file_unit = closeFile(this%output_arr(i)%file_unit)
        end if
    end do

    ! also close integration logfile
    call numint%setIntegrationLogfile('OFF')
    return

  end subroutine close_open_files

  !------------------------------------------------------------------------------------
  !> @brief       Main switch to turn ON/OFF output to files
  !!
  !> @author      Vitali Braun
  !!
  !> @param[in]   switch        Toggle ON/OFF
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      switch_output
  !!
  !------------------------------------------------------------------------------------
  subroutine switch_output(this,switch)
    class(Output_class) :: this
    logical, intent(in) :: switch

    this%flag_output = switch
    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Set input type for state vector
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !------------------------------------------------------------------------------------
  subroutine set_input_type(this,i)
    class(Output_class) :: this
    integer, intent(in) :: i
    this%input_type = i
    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Set input type for covariance matrix
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !------------------------------------------------------------------------------------
  subroutine set_input_type_cov(this,i)
    class(Output_class) :: this
    integer, intent(in) :: i
    this%input_type_cov = i
    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Get main output switch
  !!
  !> @author      Vitali Braun
  !!
  !> @return      .true. when output is requested
  !!
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      get_output_switch
  !!
  !------------------------------------------------------------------------------------

  logical function get_output_switch(this)
    class(Output_class) :: this
    get_output_switch = this%flag_output
    return
  end function

  !------------------------------------------------------------------------------------
  !> @brief       Writing an integer value to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @param[in]     id      variable ID
  !> @param[in]     val     value
  !!
  !> @date        <ul>
  !!                <li> 23.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_integer_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_integer_to_output(this, id, val)

    class(Output_class)             :: this
    character(len=*), intent(in)    :: id
    integer, intent(in)             :: val

    select case(id)
        case(C_OPT_SAT_PROPERTIES)
            this%satellite_properties = val
        case(C_OUTPUT_STEP)
            this%output_step = val
    end select
    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Writing a date value to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @param[in]     id      variable ID
  !> @param[in]     val     value
  !!
  !> @date        <ul>
  !!                <li> 23.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_date_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_date_to_output(this, id, val)

    class(Output_class)             :: this
    character(len=*), intent(in)    :: id
    type(time_t),     intent(in)    :: val

    select case(id)
        case(C_EPOCH_START_GD)
            this%start_epoch = val
        case(C_EPOCH_END_GD)
            this%end_epoch = val
    end select
    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Writing a string value to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @param[in]     id      variable ID
  !> @param[in]     val     value
  !!
  !> @date        <ul>
  !!                <li> 23.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_string_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_string_to_output(this, id, val)

    class(Output_class)             :: this
    character(len=*), intent(in)    :: id
    character(len=*), intent(in)    :: val

    select case(id)
        case(C_RUN_ID)
            this%run_id = val
    end select
    return
  end subroutine


  !------------------------------------------------------------------------------------
  !> @brief       Writing a double precision value to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @param[in]     id      variable ID
  !> @param[in]     val     value
  !!
  !> @date        <ul>
  !!                <li> 23.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_float_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_float_to_output(this,id, val)

    class(Output_class)             :: this
    character(len=*), intent(in)    :: id
    real(dp), intent(in)            :: val

    select case(id)
        case(C_PAR_MASS)
            this%mass = val
        case(C_PAR_CROSS_SECTION)
            this%cross_section = val
        case(C_PAR_CDRAG)
            this%cdrag = val
        case(C_PAR_CREFL)
            this%crefl = val
    end select
    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Writing the initial covariance matrix to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_covariance_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_covariance_to_output(this,cov)
    class(Output_class)             :: this
    type(covariance_t), intent(in)  :: cov
    this%initial_covariance = cov
    return
  end subroutine


  !------------------------------------------------------------------------------------
  !> @brief       Writing the initial Kepler elements to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_kepler_elements_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_kepler_elements_to_output(this,state)

    class(Output_class)             :: this
    type(kepler_t), intent(in)      :: state

    this%initial_orbit_osc = state
    this%isSetInitialState = .true.
    return
  end subroutine


  !------------------------------------------------------------------------------------
  !> @brief       Writing the initial state vector to the output
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 22.03.2017 (first implementation)</li>
  !!              </ul>
  !!
  !> @anchor      write_state_to_output
  !!
  !------------------------------------------------------------------------------------
  subroutine write_state_to_output(this,state)

    class(Output_class)             :: this
    type(state_t), intent(in)       :: state

    this%initial_orbit_csv = state
    this%isSetInitialState = .true.

    return
  end subroutine

  !------------------------------------------------------------------------------------
  !> @brief       Preparation of NEPTUNE output files
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 09.10.2013 (first implementation)</li>
  !!                <li> 21.03.2017 (adapted to new structure)</li>
  !!              </ul>
  !!
  !> @details     This routine initializes the NEPTUNE output files. For example,
  !!              it triggers the writing of the file headers.
  !!
  !> @anchor      prepare_output
  !!
  !------------------------------------------------------------------------------------
  subroutine prepare_output(this,               &
                            gravity_model,      &
                            atmosphere_model,   &
                            manoeuvres_model,   &
                            radiation_model,    &
                            satellite_model,    &
                            solarsystem_model,  &
                            thirdbody_model,    &
                            numint,             &
                            derivatives_model,  &
                            version_model,      &
                            reduction,          &
                            correlation_model   &
                            )

    class(Output_class)                     :: this
    type(Gravity_class),intent(inout)       :: gravity_model                    ! Gravity model
    type(Atmosphere_class),intent(inout)    :: atmosphere_model                 ! Atmosphere model
    type(Manoeuvres_class),intent(inout)    :: manoeuvres_model                 ! Manoeuvres model
    type(Radiation_class),intent(inout)     :: radiation_model                  ! Radiation model
    type(Satellite_class),intent(inout)     :: satellite_model                  ! Satellite model
    type(Solarsystem_class),intent(inout)   :: solarsystem_model                ! Satellite model
    type(ThirdBody_class),intent(inout)     :: thirdbody_model                  ! Third body model
    type(Numint_class),intent(inout)        :: numint                           ! Numerical integrator
    type(Derivatives_class),intent(inout)   :: derivatives_model                ! Derivatives model
    type(Version_class),intent(inout)       :: version_model                    ! Version model
    type(Reduction_type),intent(inout)      :: reduction
    type(Correlation_class),intent(inout)   :: correlation_model                ! Correlation model

    integer :: i    ! loop counter

    this%flag_output = .false.

    do i = 1, size(this%output_arr)
        if(this%output_arr(i)%switch) then
            this%flag_output = .true.
            this%output_arr(i)%file_unit = openFile(this%output_arr(i)%file_name, SEQUENTIAL, OUT_FORMATTED_OVERWRITE)
            !** write header
            call this%write_header(gravity_model,                   &
                                    atmosphere_model,               &
                                    manoeuvres_model,               &
                                    radiation_model,                &
                                    satellite_model,                &
                                    thirdbody_model,                &
                                    numint,                         &
                                    derivatives_model,              &
                                    version_model,                  &
                                    reduction,                      &
                                    correlation_model,              &
                                    this%output_arr(i)%file_unit,   &
                                    this%output_arr(i)%file_name)
        end if
    end do
    return

  end subroutine prepare_output

  !------------------------------------------------------------------------------------
  !> @brief       NEPTUNE file output routine
  !!
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 01.02.2015 (added Doxygen header and atmospheric density output file handler)</li>
  !!              </ul>
  !!
  !> @details     This routine writes NEPTUNE output to dedicated output files. F
  !!
  !> @param[in] gravity_model       the gravity model
  !> @param[in] atmosphere_model    the atmosphere model
  !> @param[in] manoeuvres_model    the manoeuvres model
  !> @param[in] radiation_model     the gravity model
  !> @param[in] satellite_model     the satellite model
  !> @param[in] solarsystem_model   the solar system model
  !> @param[in] thirdbody_model     the third body model
  !> @param[in] tides_model         the tides model
  !> @param[in] numint              the numerical integrator
  !> @param[in] derivatives_model   the derivatives model
  !> @param[in] reduction           the reduction model
  !> @param[in] time                MJD of output epoch
  !> @param[in] state_out           Output state vector
  !> @param[in] covar_out           Output covariance matrix
  !!
  !> @anchor      output
  !!
  !------------------------------------------------------------------------------------
  subroutine output( this,              &
                     gravity_model,     &  ! <-> TYPE gravity model
                     atmosphere_model,  &  ! <-> TYPE atmospheric model
                     manoeuvres_model,  &  ! <-> TYPE manoeuvres model
                     radiation_model,   &  ! <-> TYPE radiation model
                     satellite_model,   &  ! <-> TYPE satellite model
                     solarsystem_model, &  ! <-> TYPE solarsystem model
                     thirdbody_model,   &  ! <-> TYPE thirdbody model
                     tides_model,       &  ! <-> TYPE tides model
                     numint,            &  ! <-> TYPE numerical intergator
                     derivatives_model, &  ! <-> TYPE derivatives model
                     reduction,         &  ! <-> TYPE reduction
                     time_offset,       &  ! <-- DBL time offset wrt. start epoch / s
                     state_out,         &  ! <-- TYP output state vector
                     covar_out          &  ! <-- TYP covariance matrix output
                   )

    implicit none

    !** interface
    !-------------------------------------------
    class(Output_class)                             :: this
    type(Gravity_class),intent(inout)               :: gravity_model
    type(Atmosphere_class),intent(inout)            :: atmosphere_model
    type(Manoeuvres_class),intent(inout)            :: manoeuvres_model
    type(Radiation_class),intent(inout)             :: radiation_model
    type(Satellite_class),intent(inout)             :: satellite_model
    type(Solarsystem_class),intent(inout)           :: solarsystem_model
    type(ThirdBody_class),intent(inout)             :: thirdbody_model
    type(Tides_class),intent(inout)                 :: tides_model
    type(Numint_class),intent(inout)                :: numint
    type(Derivatives_class),intent(inout)           :: derivatives_model
    type(Reduction_type),intent(inout)              :: reduction
    real(dp),           intent(in)                  :: time_offset
    type(state_t),      intent(in)                  :: state_out
    type(covariance_t), intent(in)                  :: covar_out
    !-------------------------------------------

    character(len=*), parameter :: csubid = "output"
    integer                 :: i,j,k                                            ! loop counter
    integer                 :: orbit_type                                       ! orbit type

    logical                 :: flag_cov_uvw                                     ! indicating that covariance matrix already transformed to UVW frame
    logical                 :: flag_deriv                                       ! flag indicating that force model has been already called
    logical                 :: flag_itrf                                        ! GCRF->ITRF transformation already performed

    real(dp), dimension(3)   :: acc_out                                         ! acceleration output vector in ECI system (km/s**2)
    real(dp), dimension(3)   :: acc_out_uvw                                     ! acceleration output vector in UVW system (km/s**2)
    real(dp)                 :: alt                                             ! current altitude (km)
    real(dp)                 :: crossSection                                    ! cross-section / m**2
    real(dp)                 :: fracday                                         ! fraction of day (0<=fracday<1)
    real(dp)                 :: dsec                                            ! second 
    real(dp), dimension(6,6) :: jacobi                                          ! jacobian matrix to convert covariance matrix from ECI to UVW system
    real(dp)                 :: lat                                             ! current latitude (degrees)
    real(dp)                 :: lon                                             ! current longitude (degrees)
    real(dp), dimension(3)   :: r_itrf                                          ! radius vector in ITRF
    real(dp)                 :: rho                                             ! total atmospheric density / g/cm**3
    real(dp), dimension(3)   :: v_itrf                                          ! velocity vector in ITRF
    real(dp), dimension(3)   :: v_rel                                           ! relative velocity wrt. drag / km/s

    type(covariance_t)   :: covar_out_uvw                                       ! covariance matrix in UVW frame
    type(kepler_t)       :: kepler                                              ! kepler elements
    type(time_t)         :: epoch                                               ! output epoch

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    flag_cov_uvw = .false.
    flag_deriv   = .false.
    flag_itrf    = .false.

    epoch%mjd = this%start_epoch%mjd + time_offset/86400.d0

    call mjd2gd(epoch%mjd,epoch%year,epoch%month,epoch%day,fracday)
    call dayFraction2HMS(fracday,epoch%hour,epoch%minute,dsec)

    do i = 1, size(this%output_arr)

      if(this%output_arr(i)%file_unit > 0) then

        select case(this%output_arr(i)%id)

          case(OUTPUT_OSC) !** osculating kepler elements vs. time

            !** convert state vector to kepler elements
            call rv2coe      (                &
                               state_out,     &
                               kepler,        &
                               orbit_type     &
                             )

            select case(orbit_type)

              case(ELLIPTICAL_INCLINED)

                write(this%output_arr(i)%file_unit,110) epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                    kepler%sma, kepler%ecc, kepler%inc*rad2deg, kepler%raan*rad2deg,               &
                                    kepler%aop*rad2deg, kepler%tran*rad2deg, kepler%man*rad2deg

              case(CIRCULAR_INCLINED)

                write(this%output_arr(i)%file_unit,110) epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                    kepler%sma, kepler%ecc, kepler%inc*rad2deg, kepler%raan*rad2deg,               &
                                    0.d0, kepler%arglat*rad2deg, kepler%man*rad2deg

              case(ELLIPTICAL_EQUATORIAL)

                write(this%output_arr(i)%file_unit,110) epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                    kepler%sma, kepler%ecc, kepler%inc*rad2deg, 0.d0, kepler%lonper*rad2deg,       &
                                    kepler%tran*rad2deg, kepler%man*rad2deg

              case(CIRCULAR_EQUATORIAL)

                write(this%output_arr(i)%file_unit,110) epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                    kepler%sma, kepler%ecc, kepler%inc*rad2deg, 0.d0, 0.d0, kepler%truelon*rad2deg, &
                                    kepler%man*rad2deg

            end select

          case(OUTPUT_CSV)  !** cartesian state vs. time

            write(this%output_arr(i)%file_unit,120)   epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd,   &
                                  state_out%r, state_out%v

          case(OUTPUT_VAR_ECI)  !** variances (diagonal elements of covariance matrix) in ECI frame vs. time

            if(numint%getCovariancePropagationFlag()) then

              write(this%output_arr(i)%file_unit,150)   epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd,   &
                                    (covar_out%elem(j,j), j=1,6)

            end if

          case(OUTPUT_VAR_UVW)  !** variances (diagonal elements of covariance matrix) in UVW frame vs. time

            if(numint%getCovariancePropagationFlag()) then

              if(.not. flag_cov_uvw) then

                !** convert covariance matrix from ECI to UVW using current state vector
                call reduction%getJacobianEci2uvw(state_out%r, state_out%v, jacobi)

                covar_out_uvw%elem = matmul(matmul(jacobi,covar_out%elem), transpose(jacobi))
                flag_cov_uvw       = .true.

              end if

              write(this%output_arr(i)%file_unit,150)   epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd,   &
                                    (covar_out_uvw%elem(j,j), j=1,6)
                                   ! (sqrt(max(0.d0,covar_out_uvw%elem(j,j)))*1.d3, j=1,6)

            end if

          case(OUTPUT_COV_ECI)  !** covariances (non-diagonal) in ECI frame vs. time

            if(numint%getCovariancePropagationFlag()) then

              write(this%output_arr(i)%file_unit,160)   epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd,   &
                                    ((covar_out%elem(j,k), k=1,j-1), j=2,6)

            end if

          case(OUTPUT_COV_UVW)  !** covariances (non-diagonal) in UVW frame vs. time

            if(numint%getCovariancePropagationFlag()) then

              if(.not. flag_cov_uvw) then

                !** convert covariance matrix from ECI to UVW using current state vector (if not done before..)
                call reduction%getJacobianEci2uvw(state_out%r, state_out%v, jacobi)

                covar_out_uvw%elem = matmul(matmul(jacobi,covar_out%elem), transpose(jacobi))
                flag_cov_uvw       = .true.

              end if

              write(this%output_arr(i)%file_unit,160)   epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd,   &
                                    ((covar_out_uvw%elem(j,k), k=1,j-1), j=2,6)

            end if

          case(OUTPUT_ACC)  !** accelerations vs. time

            if(.not. flag_deriv) then

              !** get accelerations for current state vector
              call derivatives_model%deriv(                     &
                                            gravity_model,      &
                                            atmosphere_model,   &
                                            manoeuvres_model,   &
                                            radiation_model,    &
                                            satellite_model,    &
                                            solarsystem_model,  &
                                            thirdbody_model,    &
                                            tides_model,        &
                                            reduction,          &
                                            state_out%r,        &
                                            state_out%v,        &
                                            time_offset,        &
                                            acc_out)
              if(hasFailed()) return

              flag_deriv = .true.

            end if

            !** convert acceleration to UVW coordinates
            call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

            write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                 (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

          case(OUTPUT_ACG) !** acceleration due to geopotential vs. time

            !** convert state to ITRF
            if(.not. flag_itrf) then

              call reduction%inertial2earthFixed(state_out%r, state_out%v, epoch%mjd, r_itrf, v_itrf)
              if(hasFailed()) return
              flag_itrf = .true.

            end if

            !** get acceleration due to geopotential
            !------------------------------------------------------
            call gravity_model%getGravityAcceleration(&
                                         reduction,   &
                                         r_itrf,      &   ! <-- DBL(3) radius vector in ITRF / km
                                         v_itrf,      &   ! <-- DBL(3) velocity vector in ITRF / km/s
                                         epoch%mjd,   &   ! <-- DBL    epoch in MJD
                                         acc_out      &   ! --> DBL(3) acceleration vector in GCRF / km/s**2
                                       )
            if(hasFailed()) return

            !write(*,*) "------------------------------------------"
            !write(*,*) "<OUTPUT>"
            !write(*,*)
            !write(*,*) "r = ", state_out%r
            !write(*,*) "mjd = ", epoch%mjd
            !write(*,*)
            !write(*,*) "------------------------------------------"
            !** convert acceleration to UVW coordinates
            call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

            write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                 (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

          case(OUTPUT_ACD)

            if(derivatives_model%getPertSwitch(PERT_DRAG)) then !** acceleration due to drag vs. time

              !** convert state to ITRF
              if(.not. flag_itrf) then

                call reduction%inertial2EarthFixed(state_out%r, state_out%v, epoch%mjd, r_itrf, v_itrf)
                if(hasFailed()) return
                flag_itrf = .true.

              end if

              !** get acceleration due to atmosphere
              !------------------------------------------------------
              call atmosphere_model%getAtmosphereAcceleration(  &
                                              gravity_model,    & ! <-> TYPE   gravity model
                                              satellite_model,  & ! <-> TYPE   satellite model
                                              solarsystem_model,& ! <-> TYPE   solarsystem model
                                              reduction,        & ! <-> TYPE   reduction
                                              state_out%r,      & ! <-- DBL(3) radius vector in GCRF
                                              state_out%v,      & ! <-- DBL(3) velocity vector in GCRF
                                              r_itrf,           & ! <-- DBL(3) radius vector in ITRF
                                              v_itrf,           & ! <-- DBL(3) velocity vector in ITRF
                                              epoch%mjd,        & ! <-- DBL    time MJD
                                              acc_out           & ! --> DBL(3) acceleration vector in inertial frame
                                            )
              if(hasFailed()) return

              !** convert acceleration to UVW coordinates
              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ACM)

            if(derivatives_model%getPertSwitch(PERT_MOON)) then !** accelerations due to moon perturbations

              call thirdbody_model%getThirdBodyAcceleration(    &
                                              solarsystem_model,&   ! <-> TYPE   solarsystem model
                                              state_out%r,      &   ! <-- DBL(3) radius vector in GCRF
                                              epoch%mjd,        &   ! <-- DBL    current time MJD
                                              ID_MOON,          &   ! <-- CHR()  third body identifier (e.g. 'SUN' or 'MOON')
                                              acc_out           &   ! --> DBL(3) acceleration vector in inertial frame
                                           )
              if(hasFailed()) return

              !** convert acceleration to UVW coordinates
              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ACS)

            if(derivatives_model%getPertSwitch(PERT_SUN)) then !** accelerations due to sun perturbations

              call thirdbody_model%getThirdBodyAcceleration(    &
                                            solarsystem_model,  & ! <-> TYPE   solarsystem model
                                            state_out%r,        &  ! <-- DBL(3) radius vector in GCRF
                                            epoch%mjd,          &  ! <-- DBL    current time MJD
                                            ID_SUN,             &  ! <-- CHR()  third body identifier (e.g. 'SUN' or 'MOON')
                                            acc_out             &  ! --> DBL(3) acceleration vector in inertial frame
                                           )
              if(hasFailed()) return

              !** convert acceleration to UVW coordinates
              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ACR)

            if(derivatives_model%getPertSwitch(PERT_SRP)) then !** accelerations due to solar radiation pressure

              call radiation_model%getSrpAcceleration(      &
                                        satellite_model,    &  ! <-> TYPE   satellite model
                                        solarsystem_model,  &  ! <-> TYPE   solarsystem model
                                        reduction,          &  ! <-> TYPE   reduction
                                        state_out%r,        &  ! <-- DBL(3) GCRF radius of spacecraft
                                        state_out%v,        &  ! <-- DBL(3) GCRF velocity of spacecraft
                                        epoch%mjd,          &  ! <-- DBL    current MJD
                                        acc_out             &  ! --> DBL(3) acceleration in inertial frame
                                     )
              if(hasFailed()) return

              !** convert acceleration to UVW coordinates
              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ACA)

            if(derivatives_model%getPertSwitch(PERT_ALBEDO)) then !** accelerations due to solar radiation pressure

              call radiation_model%getAlbedoAcceleration(       &
                                           satellite_model,     &  ! <-> TYPE   satellite model
                                           solarsystem_model,   &  ! <-> TYPE   solarsystem model
                                           reduction,           &  ! <-> TYPE   reduction
                                           state_out%r,         &  ! <-- DBL(3) GCRF radius of spacecraft
                                           state_out%v,         &  ! <-- DBL(3) GCRF velocity of spacecraft
                                           epoch%mjd,           &  ! <-- DBL    current MJD
                                           acc_out              &  ! --> DBL(3) acceleration in inertial frame
                                        )
              if(hasFailed()) return

              !** convert acceleration to UVW coordinates
              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ACT)

            if(derivatives_model%getPertSwitch(PERT_SOLID_TIDES)) then !** accelerations due to solid Earth tides

              !** convert state to ITRF
              if(.not. flag_itrf) then

                call reduction%inertial2earthFixed(state_out%r, state_out%v, epoch%mjd, r_itrf, v_itrf)
                if(hasFailed()) return
                flag_itrf = .true.

              end if

              call tides_model%getTidesAcceleration(            &
                                           solarsystem_model,   &   ! <-> TYPE   solarsystem model
                                           reduction,           &   ! <-> TYPE   reduction
                                           state_out%r,         &   ! <-- DBL(3) radius vector in GCRF
                                           r_itrf,              &   ! <-- DBL(3) radius vector in ITRF
                                           v_itrf,              &   ! <-- DBL(3) velocity vector in ITRF
                                           epoch%mjd,           &   ! <-- DBL    current epoch
                                           SOLID_TIDES,         &   ! <-- INT    solid earth tides
                                           acc_out              &   ! --> DBL(3) acceleration vector in inertial frame
                                       )
              if(hasFailed()) return

              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ACO)

            if(derivatives_model%getPertSwitch(PERT_OCEAN_TIDES)) then !** accelerations due to solid Earth tides

              !** convert state to ITRF
              if(.not. flag_itrf) then

                call reduction%inertial2earthFixed(state_out%r, state_out%v, epoch%mjd, r_itrf, v_itrf)
                if(hasFailed()) return
                flag_itrf = .true.

              end if

              call tides_model%getTidesAcceleration(            &
                                           solarsystem_model,   &   ! <-> TYPE   solarsystem model
                                           reduction,           &   ! <-> TYPE   reduction
                                           state_out%r,         &   ! <-- DBL(3) radius vector in GCRF
                                           r_itrf,              &   ! <-- DBL(3) radius vector in ITRF
                                           v_itrf,              &   ! <-- DBL(3) velocity vector in ITRF
                                           epoch%mjd,           &   ! <-- DBL    current epoch
                                           OCEAN_TIDES,         &   ! <-- INT    ocean tides
                                           acc_out              &   ! --> DBL(3) acceleration vector in inertial frame
                                       )
              if(hasFailed()) return

              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_AMN)    ! orbital maneuvers accelerations

            if(derivatives_model%getPertSwitch(PERT_MANEUVERS)) then !** accelerations due to orbital maneuvers

              call manoeuvres_model%get_maneuver_acceleration(    &
                                             reduction,         &   ! <-> TYPE reduction
                                             state_out%r,       &   ! <-- DBL(3) radius vector in GCRF
                                             state_out%v,       &   ! <-- DBL(3) velocity vector in GCRF
                                             epoch%mjd,         &   ! <-- DBL    current epoch
                                             acc_out            &   ! --> DBL(3) acceleration vector in inertial frame
                                          !   .true.       &   ! --> LOG    output flag set
                                          )
              if(hasFailed()) return

              call reduction%eci2uvw(state_out%r, state_out%v, acc_out, acc_out_uvw)

              write(this%output_arr(i)%file_unit,130)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   (acc_out(j), j=1,3), (acc_out_uvw(j), j=1,3), mag(acc_out)

            end if

          case(OUTPUT_ATM)

            if(derivatives_model%getPertSwitch(PERT_DRAG)) then !** atmospheric density only available if drag is propagated...

              !** convert state to ITRF
              if(.not. flag_itrf) then

                call reduction%inertial2EarthFixed(state_out%r, state_out%v, epoch%mjd, r_itrf, v_itrf)
                if(hasFailed()) return
                flag_itrf = .true.

              end if

              !** get atmospheric density
              rho = atmosphere_model%getAtmosphericDensity(gravity_model, state_out%r, state_out%v, r_itrf, v_itrf, epoch%mjd)*1.d-12 ! in g/cm**3

              !** get relative velocity (considering wind IF selected)
              v_rel = atmosphere_model%getRelativeVelocity(reduction, state_out%r, state_out%v, epoch%mjd)
              if(hasFailed()) return

              !** get cross-section
              crossSection = satellite_model%getObjectCrossSection(solarsystem_model, reduction, epoch%mjd, state_out%r, v_rel)
              if(hasFailed()) return

              write(this%output_arr(i)%file_unit,170)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, &
                                   rho, crossSection*1.d6, (v_rel(j), j=1,3)    ! '1.d6' for cross-section to convert to m**2

            end if

          case(OUTPUT_GLL) !** geodetic longitude, latitude and altitude vs. time

            if(.not. flag_deriv) then

              call derivatives_model%deriv(                     &
                                            gravity_model,      &
                                            atmosphere_model,   &
                                            manoeuvres_model,   &
                                            radiation_model,    &
                                            satellite_model,    &
                                            solarsystem_model,  &
                                            thirdbody_model,    &
                                            tides_model,        &
                                            reduction,          &
                                            state_out%r,        &
                                            state_out%v,        &
                                            time_offset,        &
                                            acc_out)
              if(hasFailed()) return
              flag_deriv = .true.

            end if

            !** convert state to ITRF
            if(.not. flag_itrf) then

              call reduction%inertial2earthFixed(state_out%r, state_out%v, epoch%mjd, r_itrf, v_itrf)
              if(hasFailed()) return
              flag_itrf = .true.

            end if

            call getGeodeticLatLon(r_itrf, alt, lat, lon)

            write(this%output_arr(i)%file_unit,140)  epoch%year, epoch%month, epoch%day, epoch%hour, epoch%minute, dsec, epoch%mjd, alt, lat*rad2deg, lon*rad2deg

        end select

      end if

    end do

    !call flush()

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

110   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,2x,f12.5,2x,f10.8,5(1x,f12.6))
120   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,3(1x,f14.6),3(1x,f12.8))
130   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,7(1x,e16.9e2))
140   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,1x,f9.2,2(1x,f8.3))
150   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,6(1x,e15.8e2))
160   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,15(1x,e15.8e2))
170   format(4x,i4.4,2('-',i2.2),1x,2(i2.2,':'),f9.6,5x,f16.9,1x,e15.8e2,x,f10.4,3(1x,f12.8))

  end subroutine output

  !==================================================================================================
  !
  !> @anchor      writeInput2Header
  !!
  !> @brief       Writes the settings from the NEPTUNE input file into the header of an output file
  !> @author      Vitali Braun
  !!
  !> @param[in] gravity_model       the gravity model
  !> @param[in] atmosphere_model    the atmosphere model
  !> @param[in] manoeuvres_model    the manoeuvres model
  !> @param[in] radiation_model     the gravity model
  !> @param[in] satellite_model     the satellite model
  !> @param[in] solarsystem_model   the solar system model
  !> @param[in] thirdbody_model     the third body model
  !> @param[in] numint              the numerical integrator
  !> @param[in] derivatives_model   the derivatives model
  !> @param[in] reduction           the reduction model
  !> @param[in] correlation_model   the correlation model
  !!
  !> @date        <ul>
  !!                <li> 09.10.2013 (initial design) </li>
  !!              </ul>
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine writeInput2Header( this,               &
                                gravity_model,      &
                                atmosphere_model,   &
                                manoeuvres_model,   &
                                radiation_model,    &
                                satellite_model,    &
                                thirdbody_model,    &
                                numint,             &
                                derivatives_model,  &
                                reduction,          &
                                correlation_model,  &
                                ich)

    implicit none

    class(Output_class)                             :: this
    type(Gravity_class),intent(inout)               :: gravity_model            ! Gravity model
    type(Atmosphere_class),intent(inout)            :: atmosphere_model         ! Atmosphere model
    type(Manoeuvres_class),intent(inout)            :: manoeuvres_model         ! Manoeuvres model
    type(Radiation_class),intent(inout)             :: radiation_model          ! Radiation model
    type(Satellite_class),intent(inout)             :: satellite_model          ! Satellite model
    type(ThirdBody_class),intent(inout)             :: thirdbody_model          ! Third body model
    type(Numint_class),intent(inout)                :: numint                   ! Numerical integrator
    type(Derivatives_class),intent(inout)           :: derivatives_model        ! Derivatives model
    type(Reduction_type),intent(inout)              :: reduction                ! Reduction
    type(Correlation_class),intent(inout)           :: correlation_model        ! Correlation model
    integer, intent(in)                             :: ich                      ! output channel

    character(len=30)   :: ctemp    ! temporary storage
    character(len=90)   :: ctemp2   ! temporary storage
    integer  :: i,j                 ! loop counter
    real(dp) :: conv                ! temporary used for conversion purposes

    write(ich,'(A)')         '#  Input Data (Run ID: '//trim(this%run_id)//'): '
    write(ich,'(A)')         '#'

    if(this%input_type == INPUT_CARTESIAN     &
      .or. this%input_type == INPUT_TEME      &
      .or. this%input_type == INPUT_ITRF      &
      .or. this%input_type == INPUT_ITRF_TEME &
      .or. this%input_type == INPUT_UNDEFINED) then   ! assuming that an undefined input only obtains a cartesian state vector through the call of neptune
      ctemp = C_INPUT_CARTESIAN
    else if(this%input_type == INPUT_OSCULATING) then
      ctemp = C_INPUT_OSCULATING
    else
      ctemp = '*ERROR*'
    end if

    write(ich,'(A)')         '#  Input type:               '//trim(ctemp)
    write(ich,'(A)')         '#  Reference Frame:          '//C_REF_FRAME_GCRF
    write(ich,'(A)')         '#'
    write(ich,'(a16,a20)')   '#  Begin epoch: ', trim(adjustl(date2string(this%start_epoch)))
    write(ich,'(a16,a20)')   '#  End epoch:   ', trim(adjustl(date2string(this%end_epoch)))
    write(ich,'(A)')         '#'

    if(this%isSetInitialState) then

      write(ich,'(A)')         '#  Initial state vector:'
      if(this%input_type == INPUT_CARTESIAN     &
        .or. this%input_type == INPUT_TEME      &
        .or. this%input_type == INPUT_ITRF      &
        .or. this%input_type == INPUT_ITRF_TEME &
        .or. this%input_type == INPUT_UNDEFINED) then

          ctemp = getUnitString(this%initial_orbit_csv%radius_unit)

          write(ich,'(A,f16.8)')         '#    '//trim(C_RX)//' / '//trim(ctemp)//' :   ',   this%initial_orbit_csv%r(1)
          write(ich,'(A,f16.8)')         '#    '//trim(C_RY)//' / '//trim(ctemp)//' :   ',   this%initial_orbit_csv%r(2)
          write(ich,'(A,f16.8)')         '#    '//trim(C_RZ)//' / '//trim(ctemp)//' :   ',   this%initial_orbit_csv%r(3)

          ctemp = getUnitString(this%initial_orbit_csv%velocity_unit)

          write(ich,'(A,f16.8)')         '#    '//trim(C_VX)//' / '//trim(ctemp)//' : ',   this%initial_orbit_csv%v(1)
          write(ich,'(A,f16.8)')         '#    '//trim(C_VY)//' / '//trim(ctemp)//' : ',   this%initial_orbit_csv%v(2)
          write(ich,'(A,f16.8)')         '#    '//trim(C_VZ)//' / '//trim(ctemp)//' : ',   this%initial_orbit_csv%v(3)

      else if(this%input_type == INPUT_OSCULATING) then

        if(this%initial_orbit_osc%angles_unit == UNIT_RAD) then
          conv = rad2deg
        else
          conv = 1.d0
        end if

        ctemp = getUnitString(this%initial_orbit_osc%sma_unit)

        write(ich,'(A,f16.8)')         '#    '//trim(C_SMA)//' / '//trim(ctemp)//' :                    ', this%initial_orbit_osc%sma
        write(ich,'(A,f16.8)')         '#    '//trim(C_ECC)//' / - :                        ', this%initial_orbit_osc%ecc
        write(ich,'(A,f16.8)')         '#    '//trim(C_INC)//' / '//C_UNIT_DEG//     ' :                       ', this%initial_orbit_osc%inc*conv
        write(ich,'(A,f16.8)')         '#    '//trim(C_RAN_LONG)//' / '//C_UNIT_DEG//' : ', this%initial_orbit_osc%raan*conv
        write(ich,'(A,f16.8)')         '#    '//trim(C_AOP_LONG)//' / '//C_UNIT_DEG//' :               ', this%initial_orbit_osc%aop*conv
        write(ich,'(A,f16.8)')         '#    '//trim(C_TAN)//' / '//C_UNIT_DEG//     ' :                  <    ', this%initial_orbit_osc%tran*conv

      end if

    else

      write(ich,'(A)')         '#  <Initial state not set in initialization>'

    end if

    write(ich,'(A)') '#'

    if(numint%getCovariancePropagationFlag()) then
      if(this%input_type_cov == INPUT_COV_UVW) then
        ctemp  = C_REF_FRAME_UVW
        ctemp2 = '#'//repeat(' ', 10)//C_RU//repeat(' ',11)//C_RV//repeat(' ',11)//C_RW//repeat(' ',11)//C_VU//repeat(' ',11)//C_VV//repeat(' ',11)//C_VW
      else if(this%input_type_cov == INPUT_COV_GCRF) then
        ctemp = C_REF_FRAME_GCRF
        ctemp2 = '#'//repeat(' ',10)//C_RX//repeat(' ',11)//C_RY//repeat(' ',11)//C_RZ//repeat(' ',11)//C_VX//repeat(' ',11)//C_VY//repeat(' ',11)//C_VZ
      else ! undefined(?) - obtain from covariance reference frame
        if(this%initial_covariance%frame == REF_FRAME_GCRF) then
          ctemp2 = '#'//repeat(' ',10)//C_RX//repeat(' ',11)//C_RY//repeat(' ',11)//C_RZ//repeat(' ',11)//C_VX//repeat(' ',11)//C_VY//repeat(' ',11)//C_VZ
        else if(this%initial_covariance%frame == REF_FRAME_UVW) then
          ctemp2 = '#'//repeat(' ', 10)//C_RU//repeat(' ',11)//C_RV//repeat(' ',11)//C_RW//repeat(' ',11)//C_VU//repeat(' ',11)//C_VV//repeat(' ',11)//C_VW
        else
          ctemp2 = '# (NOTE: unknown reference frame for covariance matrix. not specified via input.)'
        end if
      end if

      write(ich,'(A)') '#  Initial covariance matrix ('//trim(adjustl(ctemp))//' / km / km/s):'
      write(ich,'(A)') adjustl(ctemp2)
      do i=1,6
        write(ich,'(a5,6(e13.6e2,x))') '#    ', (this%initial_covariance%elem(i,j), j=1,6)
      end do
      write(ich,'(A)') '#'

    end if

    !** satellite characteristics
    !----------------------------------------------------------
    write(ich,'(A)') '#  Satellite properties:'
    write(ich,'(A,f10.3)') '#    - Mass / kg :                              ', this%mass

    if(this%satellite_properties == PROP_CONSTANT) then
      write(ctemp,'(f10.4)') this%cross_section
    else
      ctemp = '<from file: '//trim(satellite_model%getSurfaceDefinitionFileName())//'>'
    end if

    write(ich,'(A,f10.4)') '#    - Cross-section drag / m**2 :              '//trim(ctemp)

    if(this%satellite_properties == PROP_CONSTANT) then
      write(ctemp,'(f10.4)') this%cdrag
    else
      ctemp = '<from file: '//trim(satellite_model%getSurfaceDefinitionFileName())//'>'
    end if

    write(ich,'(A,f10.4)') '#    - Drag coefficient / - :                   '//trim(ctemp)

    if(this%satellite_properties == PROP_CONSTANT) then
      write(ctemp,'(f10.4)') this%crefl
    else
      ctemp = '<from file: '//trim(satellite_model%getSurfaceDefinitionFileName())//'>'
    end if

    write(ich,'(A,f10.3)') '#    - Reflectivity coefficient / - :           '//trim(ctemp)

    write(ich,'(A)') '#'

    !** perturbations
    !----------------------------------------------------------
    write(ich,'(A)') '#  Considered perturbations:'

    if(gravity_model%getGeopotentialDegree() == 0) then

      write(ich,'(A)') '#    - '//toLowercase(C_GEOPOTENTIAL,1)//' (spherical Earth)'

    else if(gravity_model%usedDistinctHarmonics()) then

      write(ich,'(A)') '#    - '//toLowercase(C_GEOPOTENTIAL,1)//' (Only distinct spherical harmonics were used, see below)'
      write(ich,'(A)') '#         Harmonics used (degree,order): '//trim(adjustl(gravity_model%getDistinctHarmonicsString()))

    else if(gravity_model%getGeopotentialDegree() > 0) then

      write(ctemp,'(i2)') gravity_model%getGeopotentialDegree()
      write(ich,'(A)') '#    - '//toLowercase(C_GEOPOTENTIAL,1)//' ('//trim(adjustl(gravity_model%getGeopotentialModelName()))//', '//trim(adjustl(ctemp))//'x'//trim(adjustl(ctemp))//')'

    end if

    if(derivatives_model%getPertSwitch(PERT_DRAG)) then

      if(derivatives_model%getPertSwitch(PERT_WIND)) then
        ctemp = '   ('//trim(adjustl(atmosphere_model%getAtmosphereModelName()))//' + '//trim(adjustl(atmosphere_model%getWindModelName()))//')'
      else
        ctemp = '   ('//trim(adjustl(atmosphere_model%getAtmosphereModelName()))//')'
      end if

      write(ich,'(A)') '#    - '//toLowercase(C_ATMOSPHERE,1)//trim(ctemp)

    end if

    if(derivatives_model%getPertSwitch(PERT_SUN)) then
      write(ich,'(A)') '#    - '//toLowercase(C_SUN,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_MOON)) then
      write(ich,'(A)') '#    - '//toLowercase(C_MOON,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_MERCURY)) then
      write(ich,'(A)') '#    - '//toLowercase(C_MERCURY,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_VENUS)) then
      write(ich,'(A)') '#    - '//toLowercase(C_VENUS,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_MARS)) then
      write(ich,'(A)') '#    - '//toLowercase(C_MARS,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_JUPITER)) then
      write(ich,'(A)') '#    - '//toLowercase(C_JUPITER,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_SATURN)) then
      write(ich,'(A)') '#    - '//toLowercase(C_SATURN,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_URANUS)) then
      write(ich,'(A)') '#    - '//toLowercase(C_URANUS,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_NEPTUNE)) then
      write(ich,'(A)') '#    - '//toLowercase(C_NEPTUNE,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_SRP)) then
      write(ich,'(A)') '#    - '//C_SRP
    end if

    if(derivatives_model%getPertSwitch(PERT_ALBEDO)) then
      write(ich,'(A)') '#    - '//toLowercase(C_ALBEDO,1)
    end if

    if(derivatives_model%getPertSwitch(PERT_SOLID_TIDES)) then
      write(ich,'(A)') '#    - '//toLowercase(C_SOLID_TIDES(1:5),1)//' '//toLowercase(C_SOLID_TIDES(7:11),2)
    end if

    if(derivatives_model%getPertSwitch(PERT_OCEAN_TIDES)) then
      write(ich,'(A)') '#    - '//toLowercase(C_OCEAN_TIDES(1:5),1)//' '//toLowercase(C_OCEAN_TIDES(7:11),2)
    end if

    if(derivatives_model%getPertSwitch(PERT_MANEUVERS)) then
      write(ich,'(A)') '#    - '//toLowercase(C_MANEUVERS,1)//'    (from file: '//trim(adjustl(manoeuvres_model%get_man_file_name()))//')'
    end if

    if(reduction%getEopFlag()) then
      ctemp = toLowercase(C_ON,1)//" ("//reduction%getEopMethod()//")"
    else
      ctemp = toLowercase(C_OFF,1)
    end if

    write(ich,'(A)')         '#'
    write(ich,'(A)')         '#  Earth orientation parameters  : '//trim(ctemp)

    if(numint%getCovariancePropagationFlag()) then
      ctemp = toLowercase(C_ON,1)
    else
      ctemp = toLowercase(C_OFF,1)
    end if

    write(ich,'(A)')         '#  Covariance  matrix propagation: '//trim(ctemp)

    if(gravity_model%getGeoCovDegree() >= 2) then
      write(ctemp,'(i2)') gravity_model%getGeoCovDegree()
      write(ich, '(A)'),     '#    - Geopotential degree: '//trim(ctemp)
    end if

    if(atmosphere_model%getDragCovFlag()) then
      write(ich,'(A)')       '#    - Drag variational equations'
    end if

    if(thirdbody_model%getThirdBodyCovFlag(ID_SUN)) then
      write(ich,'(A)')       '#    - Sun variational equations'
    end if

    if(thirdbody_model%getThirdBodyCovFlag(ID_MOON)) then
      write(ich,'(A)')       '#    - Moon variational equations'
    end if

    if(radiation_model%getSrpCovFlag()) then
      write(ich,'(A)')       '#    - SRP variational equations'
    end if

    if(correlation_model%getNoisePropagationFlag()) then
      ctemp = toLowercase(C_ON,1)
    else
      ctemp = toLowercase(C_OFF,1)
    end if

    write(ich,'(A)')         '#  Correlation matrix propagation: '//trim(ctemp)

    write(ich,'(A)')         '#'
    write(ich,'(A,f6.1)')    '#  Re-entry stop condition / km :                ', atmosphere_model%getMinAltitude()
    write(ich,'(A,i6)')      '#  Const. geomagnetic activity Ap (in forecast): ', atmosphere_model%getApForecast()
    write(ich,'(A,f6.1)')    '#  Const. solar flux F10.7 (in forecast) / sfu:  ', atmosphere_model%getSolarFluxForecast()
    write(ich,'(A)')         '#'

    write(ich,'(A,i6)')      '#  Output time step / s :                        ', this%output_step
    if(numint%getCovariancePropagationFlag()) then
      write(ich,'(A,i6)')    '#  Covariance matrix integration step / s :      ', int(numint%getCovarianceIntegrationStep())
    end if

    write(ich,'(A)')         '#'
    write(ich,'(A,e13.6e2)') '#  Relative error tolerance: ', numint%getRelativeTolerance()
    write(ich,'(A,e13.6e2)') '#  Absolute error tolerance: ', numint%getAbsoluteTolerance()
    write(ich,'(A)')         '#'


  end subroutine writeInput2Header

  !==================================================================================================
  !
  !> @anchor      write_header
  !!
  !> @brief       Writes the header of an output file
  !> @author      Vitali Braun
  !!
  !> @param[in] gravity_model       the gravity model
  !> @param[in] atmosphere_model    the atmosphere model
  !> @param[in] manoeuvres_model    the manoeuvres model
  !> @param[in] radiation_model     the gravity model
  !> @param[in] satellite_model     the satellite model
  !> @param[in] solarsystem_model   the solar system model
  !> @param[in] thirdbody_model     the third body model
  !> @param[in] numint              the numerical integrator
  !> @param[in] derivatives_model   the derivatives model
  !> @param[in] version_model       the version model
  !> @param[in] reduction           the reduction model
  !> @param[in] correlation_model   the correlation model
  !!
  !> @date        <ul>
  !!                <li> 09.10.2013 (initial design) </li>
  !!              </ul>
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine write_header(                  &
                         this,              &
                         gravity_model,     & ! <-> Gravity model
                         atmosphere_model,  & ! <-> Atmosphere model
                         manoeuvres_model,  & ! <-> Manoeuvres Model
                         radiation_model,   & ! <-> Radiation Model
                         satellite_model,   & ! <-> Satellite Model
                         thirdbody_model,   & ! <-> Third body Model
                         numint,            & ! <-> Numerical integrator
                         derivatives_model, & ! <-> Derivatives model
                         version_model,     & ! <-> Version model
                         reduction,         & ! <-> Reduction
                         correlation_model, & ! <-> Correlation model
                         ich,               & ! <-- INT   channel number
                         cfile              & ! <-- CHR() file path and name
                        )

    implicit none

    !** interface
    !-----------------------------------------
    class(Output_class)                     :: this
    type(Gravity_class),intent(inout)       :: gravity_model                    ! Gravity model
    type(Atmosphere_class),intent(inout)    :: atmosphere_model                 ! Atmosphere model
    type(Manoeuvres_class),intent(inout)    :: manoeuvres_model                 ! Manoeuvres model
    type(Radiation_class),intent(inout)     :: radiation_model                  ! Radiation model
    type(Satellite_class),intent(inout)     :: satellite_model                  ! Satellite model
    type(ThirdBody_class),intent(inout)     :: thirdbody_model                  ! Third body model
    type(Numint_class),intent(inout)        :: numint                           ! Numerical integrator
    type(Derivatives_class),intent(inout)   :: derivatives_model                ! Derivatives model
    type(Version_class),intent(inout)       :: version_model                    ! Version model
    type(Reduction_type),intent(inout)      :: reduction                        ! Reduction
    type(Correlation_class),intent(inout)   :: correlation_model                ! Correlation model
    integer,          intent(in)            :: ich
    character(len=*), intent(in)            :: cfile
    !-----------------------------------------

    character(len=*), parameter :: csubid = "write_header"
    !** locals
    integer i1
    integer i2

    character(len=8)   :: cdate           ! system time
    character(len=10)  :: ctime           ! system time
    character(len=80)  :: csep            ! separator
    character(len=3)   :: suffix          ! file suffix, e.g. ".log"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    csep = '#'//repeat('-',79)

    if(ich == 0) then

      !** Get current date and time
      call date_and_time(cdate,ctime)

      write(this%chead(1),'(A)')  "#"
      write(this%chead(2),'(A)')  "#     _/     _/  _/_/_/_/  _/_/_/   _/_/_/_/_/  _/   _/  _/      _/  _/_/_/_/"
      write(this%chead(3),'(A)')  "#    _/_/   _/  _/        _/    _/     _/      _/   _/  _/_/    _/  _/"
      write(this%chead(4),'(A)')  "#   _/ _/  _/  _/_/_/    _/_/_/       _/      _/   _/  _/  _/  _/  _/_/_/"
      write(this%chead(5),'(A)')  "#  _/   _/_/  _/        _/           _/      _/   _/  _/    _/_/  _/"
      write(this%chead(6),'(A)')  "# _/     _/  _/_/_/_/  _/           _/        _/_/   _/      _/  _/_/_/_/"
      write(this%chead(7),'(A)')  "#"
      write(this%chead(8),'(A)')  "#                Networking/Partnering Initiative ESA/TUBS"
      write(this%chead(9),'(A)')  "#       "//version_model%getNeptuneAbbreviation()
      write(this%chead(10),'(A)') '#                   (NEPTUNE '//trim(version_model%to_string(version_model%get_neptune_version()))//' :: '//trim(version_model%getVersionDate(version_model%to_string(version_model%get_neptune_version())))//')'
      write(this%chead(11),'(A)') trim(csep)
      write(this%chead(12),'(A)') "# This file was generated by the program NEPTUNE on"
      write(this%chead(13),'(A)') "# "//cdate(1:4)//"-"//cdate(5:6)//"-"//cdate(7:8)//"T"//ctime(1:2)//":"//ctime(3:4)//":"//ctime(5:6)
      write(this%chead(14),'(A)') trim(csep)

    else

      !---------------------------------------------------------------------
      !** Combine FIRST HEADER LINE contains file name
      this%chead(0) = ''
      this%chead(0) = '#--<NPI/NEPTUNE/'//trim(cfile)//'>'
      this%chead(0)(54:80) = '-<ESA-TUBS/VB>------------'

      do i1=1,80
        if (this%chead(0)(i1:i1) == ' ') this%chead(0)(i1:i1) = '-'
      end do
      !---------------------------------------------------------------------

      !---------------------------------------------------------------------
      !** DUMP general head lines to output file
      do i1=0,14
        write(ich,'(A)') trim(this%chead(i1))
      end do
      !---------------------------------------------------------------------

      !** find type of file by looking at suffix
      !---------------------------------------------------------------------
      do i1=len_trim(cfile),1,-1
        if(cfile(i1:i1) == ".") then
          suffix = cfile(i1+1:i1+3)
        end if
      end do

      if(index(cfile,numint%getIntegrationLogfileName()).ne.0) then    ! integration logfile

        write(ich,'(A)')'#'

        !** write title
        write(ich,'(A)') '#  Title  : Numerical Integration Logfile'
        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'
        write(ich,'(A)') "#  Explanation of 'Flag': "
        write(ich,'(A)') '#   INI = The integrator has been initialized in that call'
        write(ich,'(A)') '#   INT = An interpolation was performed (no integration)'
        write(ich,'(A)') '#   EXT = An integration was performed (AND possible interpolation, if required)'
        write(ich,'(A)') '#   NST = Unlikely case that nothing was done, as requested time is exactly the integration time (No step in this case!)'
        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'

        !** write input parameters
        call this%writeInput2Header(                    &
                                    gravity_model,      &
                                    atmosphere_model,   &
                                    manoeuvres_model,   &
                                    radiation_model,    &
                                    satellite_model,    &
                                    thirdbody_model,    &
                                    numint,             &
                                    derivatives_model,  &
                                    reduction,          &
                                    correlation_model,  &
                                    ich)

        !** write long separation line
        write(ich,'(A)') "#"//repeat('-', 182)

        write(ich,'(A,$)') "#          Date           "
        write(ich,'(A,$)') "    Modified Julian Day   "
        write(ich,'(A,$)') "         Stepsize         "
        write(ich,'(A,$)') "     Integration Order    "
        write(ich,'(A,$)') " Error Double Integration "
        write(ich,'(A,$)') " Error Single Integration "
        write(ich,'(A)')   "     Call No. [Flag]      "

        write(ich,'(A,$)') "#_ [YYYY-MM-DD hh:mm:ss]__"
        write(ich,'(A,$)') "_____[JD - 2400000.5]_____"
        write(ich,'(A,$)') "____________[s]___________"
        do i2 = 1,4

          write(ich,'(A,$)') "____________[-]___________"

        end do

        write(ich,'(A)') "_"
        write(ich,'(A)') "#"

      else if(suffix == 'csv') then   ! cartesian state vector output file

        write(ich,'(A)')'#'

        !** write title
        write(ich,'(A)') '#  Title  : Cartesian State Vector Output File'
        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'

        !** write input parameters
        call this%writeInput2Header(gravity_model,      &
                                    atmosphere_model,   &
                                    manoeuvres_model,   &
                                    radiation_model,    &
                                    satellite_model,    &
                                    thirdbody_model,    &
                                    numint,             &
                                    derivatives_model,  &
                                    reduction,          &
                                    correlation_model,  &
                                    ich)

        !** write long separation line
        write(ich,'(A)') "#"//repeat("-",131)

        write(ich,'(A,$)') "#          Date         "
        write(ich,'(A,$)') "  Modified Julian Day "
        write(ich,'(A,$)') "  Position-X   "
        write(ich,'(A,$)') "  Position-Y   "
        write(ich,'(A,$)') "  Position-Z   "
        write(ich,'(A,$)') " Velocity-X  "
        write(ich,'(A,$)') " Velocity-Y  "
        write(ich,'(A)')   " Velocity-Z "

        write(ich,'(A,$)') "#_ [YYYY-MM-DD hh:mm:ss]"
        write(ich,'(A,$)') "___[JD - 2400000.5]___"

        do i2 = 1,3

          write(ich,'(A,$)') "_____[km]______"

        end do

        do i2 = 1,3

          write(ich,'(A,$)') "____[km/s]___"

        end do

        write(ich,'(A)') "_"
        write(ich,'(A)') "#"

      else if(suffix == 'osc') then   ! osculating kepler elements output file

        write(ich,'(A)')'#'

        !** write title
        write(ich,'(A)') '#  Title  : Osculating Kepler Elements Output File'
        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'

        !** write input parameters
        call this%writeInput2Header(                    &
                                    gravity_model,      &
                                    atmosphere_model,   &
                                    manoeuvres_model,   &
                                    radiation_model,    &
                                    satellite_model,    &
                                    thirdbody_model,    &
                                    numint,             &
                                    derivatives_model,  &
                                    reduction,          &
                                    correlation_model,  &
                                    ich)

        !** write long separation line
        write(ich,'(A)') "#"//repeat('-',137)

        write(ich,'(A,$)') "#          Date         "
        write(ich,'(A,$)') "  Modified Julian Day "
        write(ich,'(A,$)') "     SMA     "
        write(ich,'(A,$)') "     ECC     "
        write(ich,'(A,$)') "     INC     "
        write(ich,'(A,$)') "    RAAN     "
        write(ich,'(A,$)') "     AOP     "
        write(ich,'(A,$)') "    TRAN     "
        write(ich,'(A)')   "    MEAN     "

        write(ich,'(A,$)') "#_ [YYYY-MM-DD hh:mm:ss]"
        write(ich,'(A,$)') "___[JD - 2400000.5]___"

        write(ich,'(A,$)') "_____[km]____"
        write(ich,'(A,$)') "_____[-]_____"

        do i2 = 1,5

          write(ich,'(A,$)') "____[deg]____"

        end do

        write(ich,'(A)') "_"
        write(ich,'(A)') "#"

      else if(suffix == 'acc' .or. suffix == 'acg' .or. suffix == 'acd' .or. &
              suffix == 'acm' .or. suffix == 'acs' .or. suffix == 'acr' .or. &
              suffix == 'ame' .or. suffix == 'acv' .or. suffix == 'ama' .or. &
              suffix == 'acj' .or. suffix == 'asa' .or. suffix == 'acn' .or. &
              suffix == 'acu' .or. suffix == 'aca' .or. suffix == 'amn') then   ! accelerations output file

        write(ich,'(A)')'#'

        !** write title
        select case(suffix)
          case('acc')
            write(ich,'(A)') '#  Title  : Accelerations output file'
          case('acg')
            write(ich,'(A)') '#  Title  : Geopotential accelerations output file'
          case('acd')
            write(ich,'(A)') '#  Title  : Drag accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_DRAG)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Drag perturbations not considered (no output)!'
            end if
          case('acm')
            write(ich,'(A)') '#  Title  : Moon accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_MOON)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Moon perturbations not considered (no output)!'
            end if
          case('acs')
            write(ich,'(A)') '#  Title  : Sun accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_SUN)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Sun perturbations not considered (no output)!'
            end if
          case('ame')
            write(ich,'(A)') '#  Title  : Mercury accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_MERCURY)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Mercury perturbations not considered (no output)!'
            end if
          case('acv')
            write(ich,'(A)') '#  Title  : Venus accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_VENUS)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Venus perturbations not considered (no output)!'
            end if
          case('ama')
            write(ich,'(A)') '#  Title  : Mars accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_MARS)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Mars perturbations not considered (no output)!'
            end if
          case('acj')
            write(ich,'(A)') '#  Title  : Jupiter accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_JUPITER)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Jupiter perturbations not considered (no output)!'
            end if
          case('asa')
            write(ich,'(A)') '#  Title  : Saturn accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_SATURN)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Saturn perturbations not considered (no output)!'
            end if
          case('acn')
            write(ich,'(A)') '#  Title  : Neptune accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_NEPTUNE)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Neptune perturbations not considered (no output)!'
            end if
          case('acu')
            write(ich,'(A)') '#  Title  : Uranus accelerations output file'
            if(.not. derivatives_model%getPertSwitch(PERT_URANUS)) then
              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Uranus perturbations not considered (no output)!'
            end if
          case('acr')
            write(ich,'(A)') '#  Title  : Solar radiation pressure accelerations output file'

            if(.not. derivatives_model%getPertSwitch(PERT_SRP)) then

              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: SRP perturbations not considered (no output)!'

            end if

          case('aca')
            write(ich,'(A)') '#  Title  : Earth radiation pressure accelerations output file'

            if(.not. derivatives_model%getPertSwitch(PERT_ALBEDO)) then

              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Earth radiation perturbations not considered (no output)!'

            end if

          case('act')
            write(ich,'(A)') '#  Title  : Solid tides accelerations output file'

            if(.not. derivatives_model%getPertSwitch(PERT_SOLID_TIDES)) then

              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Solid tides perturbations not considered (no output)!'

            end if

          case('aco')
            write(ich,'(A)') '#  Title  : Ocean tides accelerations output file'

            if(.not. derivatives_model%getPertSwitch(PERT_OCEAN_TIDES)) then

              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Ocean tides perturbations not considered (no output)!'

            end if

          case('amn')
            write(ich,'(A)') '#  Title  : Maneuver accelerations output file'

            if(.not. derivatives_model%getPertSwitch(PERT_MANEUVERS)) then

              write(ich,'(A)') '#'
              write(ich,'(A)') '#  NOTE: Maneuvers not considered (no output)!'

            end if

        end select

        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'

        !** write input parameters
        call this%writeInput2Header(                    &
                                    gravity_model,      &
                                    atmosphere_model,   &
                                    manoeuvres_model,   &
                                    radiation_model,    &
                                    satellite_model,    &
                                    thirdbody_model,    &
                                    numint,             &
                                    derivatives_model,  &
                                    reduction,          &
                                    correlation_model,  &
                                    ich)

        !** write long separation line
        write(ich,'(A)') "#"//repeat('-', 165)

        write(ich,'(A,$)') "#          Date         "
        write(ich,'(A,$)') "  Modified Julian Day "
        write(ich,'(A,$)') "     Acc.-X      "
        write(ich,'(A,$)') "     Acc.-Y      "
        write(ich,'(A,$)') "     Acc.-Z      "
        write(ich,'(A,$)') "     Acc.-R      "
        write(ich,'(A,$)') "     Acc.-S      "
        write(ich,'(A,$)') "     Acc.-W      "
        write(ich,'(A)')   "   Acc. Total    "

        write(ich,'(A,$)') "#_ [YYYY-MM-DD hh:mm:ss]"
        write(ich,'(A,$)') "___[JD - 2400000.5]___"

        do i2 = 1,7

          write(ich,'(A,$)') "____[km/s**2]____"

        end do

        write(ich,'(A)') "_"
        write(ich,'(A)') "#"

      else if(suffix == 'gll') then   ! geodetic latitude, longitude and altitude output file

        write(ich,'(A)')'#'

        !** write title
        write(ich,'(A)') '#  Title  : Altitude, Geodetic Latitude and Longitude Output File'
        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'

        !** write input parameters
        call this%writeInput2Header(                    &
                                    gravity_model,      &
                                    atmosphere_model,   &
                                    manoeuvres_model,   &
                                    radiation_model,    &
                                    satellite_model,    &
                                    thirdbody_model,    &
                                    numint,             &
                                    derivatives_model,  &
                                    reduction,          &
                                    correlation_model,  &
                                    ich)

        !** write long separation line
        write(ich,'(A)') "#"//repeat('-',137)

        write(ich,'(A,$)') "#          Date         "
        write(ich,'(A,$)') "    Mod. Julian Day  "
        write(ich,'(A,$)') " Altitude "
        write(ich,'(A,$)') "    Lat. "
        write(ich,'(A)')   "    Lon. "

        write(ich,'(A,$)') "#_ [YYYY-MM-DD hh:mm:ss]"
        write(ich,'(A,$)') "___[JD - 2400000.5]____"
        write(ich,'(A,$)') "_[km]___"
        write(ich,'(A,$)') "___[deg]_"
        write(ich,'(A,$)') "___[deg]__"

        write(ich,'(A)') "_"
        write(ich,'(A)') "#"

      else if(suffix == 'atm') then   ! total atmospheric density

        write(ich,'(A)')'#'

        !** write title
        write(ich,'(A)') '#  Title  : Total Atmospheric Density Output File'
        write(ich,'(A)') '#'
        write(ich,'(A)') csep
        write(ich,'(A)') '#'

        !** write input parameters
        call this%writeInput2Header(                    &
                                    gravity_model,      &
                                    atmosphere_model,   &
                                    manoeuvres_model,   &
                                    radiation_model,    &
                                    satellite_model,    &
                                    thirdbody_model,    &
                                    numint,             &
                                    derivatives_model,  &
                                    reduction,          &
                                    correlation_model,  &
                                    ich)

        !** write long separation line
        write(ich,'(A)') "#"//repeat('-',137)

        write(ich,'(A,$)') "#          Date         "
        write(ich,'(A,$)') "    Mod. Julian Day  "
        write(ich,'(A,$)') " Tot. Density "
        write(ich,'(A,$)') " Cross-section "
        write(ich,'(A)')   "          Rel. velocity"

        write(ich,'(A,$)') "#_ [YYYY-MM-DD hh:mm:ss]"
        write(ich,'(A,$)') "___[JD - 2400000.5]____"
        write(ich,'(A,$)') "__[g/cm**3]___"
        write(ich,'(A,$)') "___[m**2]___"
        write(ich,'(A,$)') "__[x / km/s]__[y / km/s]__[z / km/s]__"

        write(ich,'(A)') "_"
        write(ich,'(A)') "#"

      end if

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine write_header


end module neptuneOutput
