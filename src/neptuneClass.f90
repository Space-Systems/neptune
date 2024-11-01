!>-----------------------------------------------------------------------------------------------
!> @brief       NEPTUNE - NPI Ephemeris Propagation Tool with Uncertainty Extrapolation
!!
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>CHK:  26.12.2017 (Adding atmospheric model class)</li>
!!                <li>CHK:  27.12.2017 (Incorporating neptuneInput module)</li>
!!              </ul>
!!
!> @details     This is the NEPTUNE class. It holds all initialized models.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      neptuneClass
!!
!!------------------------------------------------------------------------------------------------
module neptuneClass
    use slam_astro,             only: getEarthRadius, setEarthRadius
    use atmosphere,             only: Atmosphere_class, MONTHLY, DAILY
    use correlation,            only: Correlation_class
    use derivatives,            only: Derivatives_class, perturbation_t, PERT_SUN, PERT_MOON, PERT_MERCURY, PERT_VENUS,  &
                                    PERT_MARS, PERT_SATURN, PERT_JUPITER, PERT_URANUS, PERT_NEPTUNE, PERT_OCEAN_TIDES, PERT_SOLID_TIDES,  &
                                    PERT_SRP, PERT_ALBEDO, PERT_MANEUVERS, PERT_GRAVITY, PERT_WIND, PERT_DRAG
    use neptune_error_handling, only: E_SEMI_MAJOR_AXIS, E_ECCENTRICITY, E_INCLINATION, E_ALTITUDE, E_INVALID_STEP, setNeptuneError
    use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, getLatestError, &
                                    E_UNKNOWN_PARAMETER, E_SPECIAL, E_OUT_OF_BOUNDS, E_INPUT_PATH, E_OUTPUT_PATH, E_DATA_PATH, E_FILE_LENGTH, E_RUN_ID, &
                                    E_STRING_LENGTH, checkIn, checkOut, E_FILE_READ_ERROR
    use slam_io,                only: C_ON, C_OFF, SWITCHED_ON, SWITCHED_OFF, openFile, closeFile, OUT_FORMATTED_OVERWRITE, SEQUENTIAL, IN_FORMATTED, LOG_AND_STDOUT, &
                                    message
    use gravity,                only: Gravity_class, COEFF_C, COEFF_S, COEFF_DC, COEFF_DS
    use maneuvers,              only: Manoeuvres_class
    use slam_math,              only: mag, infinite, undefined
    use neptuneParameters,      only: n_parameters, n_options, C_HARMONIC_S, C_HARMONIC_C, C_HARMONICS, C_HARMONIC_SD_S, C_HARMONIC_SD_C,   &
                                    C_INITIAL_COVARIANCE, C_INITIAL_STATE, C_OPT_AP_FORECAST, C_OPT_SAT_PROPERTIES, C_OPT_GEO_MODEL,      &
                                    OPTION_COV_METHOD, OPTION_GEO_MODEL, OPTION_PROPS, C_FILE_QDGRID, C_FILE_HWIND, C_FILE_GEOP, &
                                    C_FILE_DWIND, C_FILE_SOLMAG_MONTHLY, C_FILE_SOLMAG, C_FILE_EOP, C_FILE_SURFACES, C_FILE_MANEUVERS,    &
                                    C_FILE_INPUT_DUMP, &
                                    C_PAR_INT_RELEPS, C_PAR_INT_ABSEPS, C_PAR_EARTH_RADIUS, C_PAR_REENTRY, PAR_MAX_ERROR, &
                                    PAR_MIN_ERROR, PAR_CREFL, C_PAR_CREFL, C_SUN, C_MOON, C_MERCURY, C_VENUS, C_MARS,     &
                                    C_JUPITER, C_SATURN, C_NEPTUNE, C_URANUS, C_SRP, PAR_INT_RELEPS, PAR_INT_ABSEPS, PAR_INT_COV_STEP, OPTION_OUTPUT, &
                                    OPTION_PN_LOOKUP, OPTION_SRP_CORRECT, OPTION_INT_LOGFILE, OPTION_HARMONICS, OPTION_EOP, OPTION_CORRELATION, &
                                    C_EPOCH_END_GD, C_EPOCH_START_GD, C_OPT_SOL_FORECAST, C_PAR_INT_COV_STEP, C_PAR_CDRAG, &
                                    PAR_CDRAG, C_PAR_CROSS_SECTION, PAR_CROSS_SECTION, C_PAR_MASS, PAR_MASS, C_COV_GEOPOTENTIAL,          &
                                    C_PAR_INT_COV_METHOD, C_OPT_STORE_DATA, C_OPT_ATMOSPHERE_MODEL, C_PAR_INT_METHOD, C_OUTPUT_STEP,      &
                                    C_GEOPOTENTIAL, C_OUTPUT_COV_UVW, C_OUTPUT_COV_ECI, C_OUTPUT_VAR_ECI, C_OUTPUT_VAR_UVW, &
                                    C_OUTPUT_AMA, C_WIND, C_ATMOSPHERE, &
                                    C_OPT_SHADOW, C_SOLID_TIDES, C_PATH_DATA, C_PATH_INPUT, C_PATH_OUTPUT, C_OUTPUT_CSV, C_OUTPUT_ACO,    &
                                    C_OUTPUT_AMN, C_OUTPUT_OSC, C_OUTPUT_GLL, C_OUTPUT_ATM, &
                                    C_OUTPUT_ACT, C_OUTPUT_ACR, C_OUTPUT_ACA, C_OUTPUT_ACU, &
                                    C_OUTPUT_ASA, C_OUTPUT_ACD, C_OUTPUT_ACG, C_OUTPUT_ACN, &
                                    C_OUTPUT_ACM, C_OUTPUT_ACS, C_OUTPUT_ACJ, C_OUTPUT_ACV, &
                                    C_OUTPUT_AME, C_OUTPUT_ACC, C_OUTPUT_FILES, C_OPT_HARMONICS, C_OPT_SRP_CORRECT, C_OPT_INT_LOG, &
                                    C_OPT_PN_LOOKUP, C_OPT_EOP, C_CORRELATION, C_COV_MOON, C_COV_SUN, C_COV_SRP, C_COV_DRAG, C_COV_PROP, &
                                    C_MANEUVERS, C_OCEAN_TIDES, C_ALBEDO, C_RUN_ID, INPUT_UNDEFINED, &
                                    C_FILE_DE_EPHEM, C_FILE_LEAP_SPICE, C_FILE_TXYS, C_FILE_PROGRESS, C_OPT_PROGRESS, C_BOUNDARY_CHECK
    use numint,                 only: Numint_class
    use neptuneOutput,          only: Output_class, neptune_out_t
    use slam_orbit_types,       only: covariance_t, state_t, kepler_t, convertToRadians, toString, parse_state_from_string, parse_covariance_from_string, assignment(=)
    use radiation,              only: Radiation_class
    use slam_reduction_class,   only: Reduction_type
    use satellite,              only: Satellite_class
    use solarsystem,            only: Solarsystem_class, ID_SUN, ID_MOON
    use slam_strings,           only: toLowercase, toString
    use thirdbody,              only: ThirdBody_class
    use tides,                  only: Tides_class
    use slam_time,              only: time_t, tokenizeDate, mjd2gd, gd2mjd, date2string, getDateTimeNow, assignment(=), jd245
    use slam_types,             only: dp
    use slam_txys,              only: getTxysFileName
    use slam_units,             only: UNIT_DEG, UNIT_KM, UNIT_KMPS
    use version,                only: Version_class

    ! all neptune inputs are stored via the following struct
    integer, parameter :: LEN_PARAMETER_NAME = 1000                             !< maximum length of a parameter's name
    integer, parameter :: LEN_PARAMETER_TYPE = 8                                !< maximum length of a parameter's type description

    type :: neptuneIn
        character(len=LEN_PARAMETER_NAME) :: parName = 'None'
        character(len=LEN_PARAMETER_NAME) :: val = 'None'
        character(len=LEN_PARAMETER_TYPE) :: valType = 'None'
        character(len=LEN_PARAMETER_NAME) :: info = ''
        logical :: parameterSet  = .false.                                                 !< this one allows to identify if the input parameter is defined (but it still could have an unset value (see below))
        logical :: valueSet  = .false.                                                   !< is true if val is set
    end type

    ! this array contains the entire input set of NEPTUNE
    integer, parameter :: INPUT_ARR_DEFAULT_SIZE = 30
    integer, parameter :: isize = 1000000                                        ! max. number of array elements

    type(state_t),      dimension(:),allocatable :: ephem                               !< ephemerides
    type(covariance_t), dimension(:),allocatable :: covMatrix                           !< covariance matrix
    type(covariance_t), dimension(:),allocatable :: setMatrix                           !< state error transition matrix

! Using threadprivate as a memory optimization here. ephem, covMatrix and setMatrix
! should be part of the Neptune_class. But because the stack might be too small when
! using OpenMP, the heap is used instead, thus moving these arrays from the class to the
! module and declaring them threadprivate so that they are not shared among all threads.
!$omp threadprivate(ephem,covMatrix,setMatrix)

    ! Neptune class (type defintion)
    type, public :: Neptune_class

        !=======================================================================
        !
        ! Models
        !
        !-----------------------------------------------------------------------
        type(Gravity_class)     :: gravity_model                                ! Gravity model
        type(Atmosphere_class)  :: atmosphere_model                             ! Atmospheric model
        type(Manoeuvres_class)  :: manoeuvres_model                             ! Manoeuvres model
        type(Radiation_class)   :: radiation_model                              ! Radiation model
        type(Solarsystem_class) :: solarsystem_model                            ! Solarsystem model
        type(ThirdBody_class)   :: thirdBody_model                              ! Third body model
        type(Tides_class)       :: tides_model                                  ! Tides model
        type(Satellite_class)   :: satellite_model                              ! Satellite model
        type(Numint_class)      :: numerical_integrator                         ! Numint class
        type(Derivatives_class) :: derivatives_model                            ! Force model evaluation
        type(Output_class)      :: output                                       ! Output
        type(Version_class)     :: version_model                                ! Version model
        type(Reduction_type)    :: reduction                                    ! Reduction model
        type(Correlation_class) :: correlation_model                            ! Correlation model

        !=======================================================================
        !
        ! Variables
        !
        !-----------------------------------------------------------------------
        integer                         :: output_step                          ! output time step in seconds, 5 mins as default
        integer, public                 :: input_type_cov                       ! input type for covariance matrix

        logical, dimension(2), public   :: flag_init_tolerances                 ! rel. and abs. tolerance initialization for numerical integration
        logical, public                 :: isSetMass                            ! becomes TRUE as soon as mass has been set
        logical, public                 :: isSetArea                            ! becomes TRUE as soon as area has been set
        logical, public                 :: isSetCdrag                           ! becomes TRUE as soon as drag coefficient has been set
        logical, public                 :: isSetCsrp                            ! becomes TRUE as soon as SRP coefficient has been set

        logical                         :: isSetStartEpoch                      ! becomes TRUE as soon as start epoch has been set
        logical                         :: isSetEndEpoch                        ! becomes TRUE as soon as end epoch has been set

        type(time_t)                    :: start_epoch                          ! start epoch
        type(time_t)                    :: end_epoch                            ! end epoch

        type(state_t)                   :: initial_orbit_csv                    !< initial orbit written to output
        type(kepler_t)                  :: initial_orbit_osc                    !< initial kepler elements written to output
        type(covariance_t)              :: initial_covariance                   !< initial covariance matrix written to output

        real(dp)                        :: progress_step                        !< step threshold for progress output to file

        !=======================================================================
        !
        ! Initialization
        !
        !-----------------------------------------------------------------------

        character(len=50)               :: runid                                ! Run ID
        character(len=512)              :: cdatapath                            ! data path
        character(len=512)              :: cinpath                              ! input path
        character(len=512)              :: coutpath                             ! output path

        logical, public                 :: neptuneInitialized                   ! NEPTUNE initialization flag

        character(len=512)              :: dump_file_name                       !< default file the neptune input dump goes to
        character(len=512)              :: progress_file_name                   !< default file the neptune progress goes to
        logical                         :: hasToReadInputFromDump               !< if true, input will be read from a dump file
        logical                         :: has_to_dump_input                    !< if true, input will be dumped to stdout
        logical                         :: write_progress                       !< if true, the progress will be written to a file

        !================================================================
        !
        ! Input options
        !
        !---------------------------------------------
        integer,dimension(n_options)    :: iopt                                 ! further options
        real(dp),dimension(n_parameters):: ipar                                 ! parameters

        integer :: current_index                                                !< current index in ephem array, while data is stored
        integer :: current_index_set                                            !< current index in setMatrix array, while data is stored
        integer :: current_index_cov                                            !< current index in covariance matrix array, while data is stored
        integer :: ephem_step                                                   !< step size (in seconds) between subsequent ephemerides. If set to '0', no data will be available

        logical :: storeEphemData                                               ! default: not stored
        logical :: warned                                                       ! default: not warned (if ephem array is full)

        !================================================
        !
        ! Data array(s)
        !
        !------------------------------------------
        type(neptuneIn),    dimension(:), allocatable :: input_arr

    contains
        ! ** setters
        procedure :: setNeptuneVar_char                                         ! character array passed as value
        procedure :: setNeptuneVar_kepl                                         ! kepler elements type passed as value
        procedure :: setNeptuneVar_stat                                         ! state type passed as value
        procedure :: setNeptuneVar_covr                                         ! covariance type passed as value
        procedure :: setNeptuneVar_int_arr                                      ! integer array passed as value
        procedure :: setNeptuneVar_geop                                         ! harmonic coefficients of the geopotential
        procedure :: setNeptuneVar_geop_toggle
        procedure :: setNeptuneVar_logical
        generic   :: setNeptuneVar =>           &
                     setNeptuneVar_char,        &                                   ! character array passed as value
                     setNeptuneVar_kepl,        &                                   ! kepler elements type passed as value
                     setNeptuneVar_stat,        &                                   ! state type passed as value
                     setNeptuneVar_covr,        &                                   ! covariance type passed as value
                     setNeptuneVar_int_arr,     &                                   ! integer array passed as value
                     setNeptuneVar_geop,        &                                   ! harmonic coefficients of the geopotential
                     setNeptuneVar_geop_toggle, &
                     setNeptuneVar_logical
        procedure :: setEndEpoch
        procedure :: setStartEpoch
        procedure :: set_input
        procedure :: setStep
        procedure :: setDataPath
        procedure :: storeData
        procedure :: storeCovData
        procedure :: storeSetData

        ! ** getters
        procedure :: getRunID
        procedure :: getInputPath
        procedure :: getDataPath
        procedure :: getEndEpoch
        procedure :: getStartEpoch
        procedure :: getOutputPath
        procedure :: get_progress_file_name
        procedure :: get_progress_step
        procedure :: get_output_step
        procedure :: has_to_write_progress
        procedure :: getDataArraySize
        procedure :: getNeptuneCovarianceData
        procedure :: getNeptuneData
        procedure :: getNeptuneSetMatrixData
        procedure :: getNumberOfEphemerides
        procedure :: getStep
        procedure :: getStoreDataFlag

        !** other
        procedure :: resetStoreCounter
        procedure :: read_input_from_dump
        procedure :: dump_input
        procedure :: write_input_to_dump
        procedure :: initialize_input_array
        procedure :: destroy
        procedure :: reallocate
        procedure :: switchStoreDataFlag

    end type Neptune_class

    ! Constructor
    interface Neptune_class
        module procedure constructor
    end interface Neptune_class

contains

    !===========================================================================
    !!
    !>  @anchor     constructor
    !!
    !>  @brief      Creates all models
    !>  @author     Christopher Kebschull
    !!
    !!
    !>  @date       <ul>
    !!                <li>26.12.2017 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    type(Neptune_class) function constructor()

        character(len=255) :: cmess !
        

        ! Instantiate model classes
        constructor%gravity_model = Gravity_class()
        constructor%atmosphere_model = Atmosphere_class()
        constructor%manoeuvres_model = Manoeuvres_class()
        constructor%radiation_model = Radiation_class()
        constructor%solarsystem_model = Solarsystem_class()
        constructor%thirdbody_model = ThirdBody_class()
        constructor%tides_model = Tides_class()
        constructor%satellite_model = Satellite_class()
        constructor%numerical_integrator = Numint_class()
        constructor%derivatives_model = Derivatives_class()
        constructor%output = Output_class()
        constructor%version_model = Version_class()
        constructor%reduction = Reduction_type()
        constructor%correlation_model = Correlation_class()

        constructor%output_step    = 300                                        ! output time step in seconds, 5 mins as default
        constructor%input_type_cov = INPUT_UNDEFINED                            ! input type for covariance matrix

        constructor%flag_init_tolerances = (/.false.,.false./)                  ! rel. and abs. tolerance initialization for numerical integration
        constructor%isSetMass  = .false.                                        ! becomes TRUE as soon as mass has been set
        constructor%isSetArea  = .false.                                        ! becomes TRUE as soon as area has been set
        constructor%isSetCdrag = .false.                                        ! becomes TRUE as soon as drag coefficient has been set
        constructor%isSetCsrp  = .false.                                        ! becomes TRUE as soon as SRP coefficient has been set

        constructor%isSetStartEpoch = .false.                                   ! becomes TRUE as soon as start epoch has been set
        constructor%isSetEndEpoch   = .false.                                   ! becomes TRUE as soon as end epoch has been set

        constructor%progress_step = 1.d-4                                       !< step threshold for progress output to file

        !=================================================================
        !
        ! Initialization
        !
        !--------------------------------------------

        constructor%runid     = 'neptune'                                       ! Run ID
        constructor%cdatapath = 'data'                                          ! data path
        constructor%cinpath   = 'input'                                         ! input path
        constructor%coutpath  = 'output'                                        ! output path

        constructor%neptuneInitialized   = .false.                              ! NEPTUNE initialization flag

        constructor%dump_file_name     = 'neptune.dump'                         !< default file the neptune input dump goes to
        constructor%progress_file_name = 'neptune.tmp'                          !< default file the neptune progress goes to
        constructor%hasToReadInputFromDump = .false.                            !< if true, input will be read from a dump file
        constructor%has_to_dump_input      = .false.                            !< if true, input will be dumped to stdout
        constructor%write_progress         = .false.                            !< if true, the progress will be written to a file

        constructor%current_index     = 0                                       !< current index in ephem array, while data is stored
        constructor%current_index_set = 0                                       !< current index in setMatrix array, while data is stored
        constructor%current_index_cov = 0                                       !< current index in covariance matrix array, while data is stored
        constructor%ephem_step        = 0                                       !< step size (in seconds) between subsequent ephemerides. If set to '0', no data will be available

        constructor%storeEphemData = .false.                                    ! default: not stored
        constructor%warned         = .false.                                    ! default: not warned (if ephem array is full)

        if (.not. allocated(ephem)) allocate(ephem(isize))
        if (.not. allocated(covMatrix)) allocate(covMatrix(isize))
        if (.not. allocated(setMatrix)) allocate(setMatrix(isize))

        ! write(cmess,'(a)') 'Initialised index = '//toString(constructor%current_index)//' of '//toString(isize)
        ! call message(cmess, LOG_AND_STDOUT)

    end function constructor


    !===========================================================================
    !!
    !>  @anchor     reallocate
    !!
    !>  @brief      Re-allocates the ephem, covMatrix and setMatrix arrays.
    !!
    !>  @details    This is only needed when openmp creates a new thread of this class.
    !!
    !>  @author     Christopher Kebschull
    !!
    !!
    !>  @date       <ul>
    !!                <li>04.04.2021 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine reallocate(this)
      class(Neptune_class)    :: this

      if (.not. allocated(ephem)) allocate(ephem(isize))
      if (.not. allocated(covMatrix)) allocate(covMatrix(isize))
      if (.not. allocated(setMatrix)) allocate(setMatrix(isize))

    end subroutine reallocate

    !===========================================================================
    !!
    !>  @anchor     destroy
    !!
    !>  @brief      Destroys all that needs destruction
    !>  @author     Christopher Kebschull
    !!
    !!
    !>  @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Neptune_class)    :: this

        ! Cleanup
        !call this%gravity_model%destroy()
        !call this%atmosphere_model%destroy()
        !call this%manoeuvres_model%destroy()
        !call this%solarsystem_model%destroy()
        !call this%output%destroy()
        !call this%reduction%destroy()

        if(allocated(this%input_arr)) deallocate(this%input_arr)
        if(allocated(ephem)) deallocate(ephem)
        if(allocated(covMatrix)) deallocate(covMatrix)
        if(allocated(setMatrix)) deallocate(setMatrix)

    end subroutine destroy

    !==============================================================================================
    !
    !> @anchor      initialize_input_array
    !!
    !> @brief       Initialize the input array with all input parameters available
    !> @author      Vitali Braun
    !!
    !> @date        <ul>
    !!                <li>VB: 16.02.2017 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    subroutine initialize_input_array(this)
        implicit none

        class(Neptune_class)                :: this

        character(len=*), parameter         :: csubid = 'initialize_input_array'
        type(perturbation_t)                :: pert
        type(neptune_out_t)                 :: outFile

        character(len=100) :: ctemp !< auxiliary
        integer :: i

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        ! FILE/PATH parameters
        call this%set_input(parName=C_PATH_DATA,   valType='dir', initFlag=.true.)
        call this%set_input(parName=C_PATH_INPUT,  valType='dir', initFlag=.true.)
        call this%set_input(parName=C_PATH_OUTPUT, valType='dir', initFlag=.true.)
        call this%set_input(parName=C_FILE_MANEUVERS, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_SURFACES,  valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_EOP, valType='file', initFlag=.true.)
        ctemp = this%reduction%getEopFileVersion(this%getDataPath())
        if(hasFailed()) return
        call this%set_input(parName=C_FILE_EOP, val=this%reduction%getEopFileName(), set=.true., info=ctemp)
        call this%set_input(parName=C_FILE_SOLMAG, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_SOLMAG, val=this%atmosphere_model%getSolmagFileName(DAILY), set=.true., info=this%atmosphere_model%getSolmagFileVersion(DAILY))
        call this%set_input(parName=C_FILE_SOLMAG_MONTHLY, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_SOLMAG_MONTHLY, val=this%atmosphere_model%getSolmagFileName(MONTHLY), set=.true., info=this%atmosphere_model%getSolmagFileVersion(MONTHLY))
        call this%set_input(parName=C_FILE_DWIND,  valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_HWIND,  valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_QDGRID, valType='file', initFlag=.true.)

        ! dump and progress files also have defaults
        call this%set_input(parName=C_FILE_INPUT_DUMP, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_INPUT_DUMP, val = this%dump_file_name, set = .true.)
        call this%set_input(parName=C_FILE_PROGRESS,   valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_PROGRESS,   val = this%progress_file_name, set = .true.)

        ! COVARIANCE
        call this%set_input(parName=C_INITIAL_COVARIANCE, valType='matrix', initFlag=.true.)

        ! STATE
        call this%set_input(parName=C_INITIAL_STATE, valType='state', initFlag=.true.)

        ! DATETIME parameters
        call this%set_input(parName=C_EPOCH_START_GD, valType='datetime', initFlag=.true.)
        call this%set_input(parName=C_EPOCH_END_GD, valType='datetime', initFlag=.true.)

        ! STRING parameters
        call this%set_input(parName=C_RUN_ID, valType='string', initFlag=.true.)
        call this%set_input(parName=C_OPT_SHADOW, valType='string', initFlag=.true.)

        ! INTEGER parameters
        call this%set_input(parName=C_GEOPOTENTIAL, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_SAT_PROPERTIES, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OUTPUT_STEP, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_COV_METHOD, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_GEO_MODEL, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_AP_FORECAST, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_STORE_DATA, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_COV_GEOPOTENTIAL, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_METHOD, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_ATMOSPHERE_MODEL, valType='integer', initFlag=.true.)

        ! REAL parameters
        call this%set_input(parName=C_PAR_MASS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_CROSS_SECTION, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_CDRAG, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_CREFL, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_REENTRY, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_RELEPS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_ABSEPS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_COV_STEP, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_EARTH_RADIUS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_OPT_SOL_FORECAST, valType='double', initFlag=.true.)

        ! ON/OFF parameters (also set default values, if available)
        do i = 1, this%derivatives_model%get_neptune_perturbation_number()
            pert = this%derivatives_model%get_neptune_perturbation(this%gravity_model,i)
            if(pert%id == PERT_GRAVITY) cycle    ! skip geopotential (as we set its degree)
            call this%set_input(parName = pert%name, valType = 'boolean', initFlag = .true.)
            call this%set_input(parName = pert%name, val = toString(pert%switch), set = .true.)
        end do

        do i = 1, this%output%get_neptune_output_file_number()
            outFile = this%output%get_neptune_output_file(i)
            call this%set_input(parName = outFile%par_name, valType = 'boolean', initFlag = .true.)
            call this%set_input(parName = outFile%par_name, val = toString(outFile%switch), set = .true.)
        end do

        call this%set_input(parName=C_CORRELATION, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_EOP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_PN_LOOKUP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_INT_LOG, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_SRP_CORRECT, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OUTPUT_FILES, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_PROP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_SUN, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_MOON, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_DRAG, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_SRP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OUTPUT_VAR_ECI, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_HARMONICS, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_PROGRESS, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_PROGRESS, val = toString(this%write_progress), set = .true.)

        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

    end subroutine

    !==============================================================================================
    !
    !> @anchor      has_to_write_progress
    !!
    !> @brief       Tells if progress is written to a file
    !> @author      Vitali Braun
    !!
    !> @return      .true. when the progress should be written to a file
    !!
    !> @date        <ul>
    !!                <li>VB: 26.03.2017 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    pure logical function has_to_write_progress(this) result(h)
        implicit none
        class(Neptune_class),intent(in)     :: this
        h = this%write_progress
        return
    end function

    !==============================================================================================
    !
    !> @anchor      get_progress_file_name
    !!
    !> @brief       Returns the name of the file the progress is written to
    !< @author      Vitali Braun
    !!
    !> @return      file name as string
    !!
    !> @date        <ul>
    !!                <li>VB: 26.03.2017 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    pure character(len=LEN_PARAMETER_NAME) function get_progress_file_name(this) result(g)
        implicit none
        class(Neptune_class),intent(in)     :: this
        ! at the moment, there is only the neptune.tmp - which is not configurable, but this can be done at some later point
        g = this%progress_file_name
        return
    end function

    !==============================================================================================
    !
    !> @anchor      get_output_step
    !!
    !> @brief       Returns the output step size
    !> @author      Vitali Braun
    !!
    !> @return      output step size
    !!
    !> @date        <ul>
    !!                <li>VB: 02.05.2017 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    pure integer function get_output_step(this) result(s)
        implicit none
        class(Neptune_class),intent(in)     :: this
        s = this%output_step
        return
    end function


    !==============================================================================================
    !
    !> @anchor      get_progress_step
    !!
    !> @brief       Returns the step threshold for the progress writing
    !> @author      Vitali Braun
    !!
    !> @return      progress step
    !!
    !> @date        <ul>
    !!                <li>VB: 26.03.2017 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    pure real(dp) function get_progress_step(this) result(s)
        implicit none
        class(Neptune_class),intent(in)     :: this
        ! at the moment, there is only one fixed value - which is not configurable, but this can be done at some later point
        s = this%progress_step
        return
    end function

    !!==============================================================================================
    !!
    !> @anchor      set_input
    !!
    !> @brief       Sets an entry in the input array
    !> @author      Vitali Braun
    !!
    !> @param[in]   parName     name of input parameter
    !> @param[in]   val         value of parameter (optional)
    !> @param[in]   valType     data type of value of parameter (optional)
    !> @param[in]   set         either true or false to toggle on or off (optional)
    !> @param[in]   info        additional information about the parameter (optional)
    !> @param[in]   initFlag    if true, it allows to define inputs - while it is just checking for the available ones otherwise
    !!
    !> @date        <ul>
    !!                <li>VB: 18.11.2016 (initial implementation)</li>
    !!                <li>VB: 16.02.2017 (changed logic, now everything goes via this routine)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    subroutine set_input(this, parName, val, valType, set, info, initFlag)
        implicit none
        class(Neptune_class)                   :: this
        character(len=*), intent(in)           :: parName
        character(len=*), intent(in), optional :: val
        character(len=*), intent(in), optional :: valType
        character(len=*), intent(in), optional :: info
        logical,          intent(in), optional :: set
        logical,          intent(in), optional :: initFlag

        character(len=*), parameter       :: csubid = 'set_input'
        character(len=LEN_PARAMETER_NAME) :: paramValue, paramName
        type(neptuneIn), dimension(:), allocatable :: tempInputArr
        logical :: found
        logical :: newInput
        logical :: paramSet
        logical :: initMode
        integer :: i

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        if(present(initFlag)) then
            initMode = initFlag
        else
            initMode = .false.
        end if

        ! if val is not present, this essentially means to reserve a spot in the input_arr (= uninitialised input)
        if(present(val)) then
            paramValue = val(:min(len(val),LEN_PARAMETER_NAME))
            newInput   = .false.
        else
            paramValue = 'None'
            newInput   = .true.
        end if

        if(present(set)) then
            paramSet = set
        else
            paramSet = .false.
        end if

        paramName = parName(:min(len(parName),LEN_PARAMETER_NAME))

        if(.not. allocated(this%input_arr)) then
            allocate(this%input_arr(INPUT_ARR_DEFAULT_SIZE))
!            this%input_arr(:)%parName = 'None'
!            this%input_arr(:)%val     = 'None'
!            this%input_arr(:)%valType = 'None'
!            this%input_arr(:)%info    = ''
!            this%input_arr(:)%parameterSet = .false.                         !< this one allows to identify if the input parameter is defined (but it still could have an unset value (see below))
!            this%input_arr(:)%valueSet     = .false.                         !< is true if val is set
        end if

        ! Run through the array and search for the name of the parameter
        ! - if it exists, replace its settings accordingly. If not, create a new entry
        ! ------------------------------------------------------------------------------
        found = .false.
        do i = 1,size(this%input_arr)
            if(this%input_arr(i)%parName == trim(paramName)) then
                if(.not. newInput) then
                    this%input_arr(i)%val      = trim(paramValue)
                    this%input_arr(i)%valueSet = paramSet
                    if(present(valType)) then
                        this%input_arr(i)%valType  = valType(:min(len(valType),LEN_PARAMETER_TYPE))
                    end if

                    if(present(info)) then
                        this%input_arr(i)%info = info(:min(len(info),LEN_PARAMETER_NAME))
                    end if
                else  ! essentially resetting
                    this%input_arr(i)%val      = 'None'
                    this%input_arr(i)%valueSet = .false.
                    this%input_arr(i)%info     = ''
                end if
                found = .true.
                exit
            end if
        end do

        ! already return here, unless we run in init mode - the next part will arrange
        ! for a new input parameter
        if(.not. initMode) then
            if(isControlled()) then
                call checkOut(csubid)
            end if
            return
        end if

        ! If not found, loop array again to find a free spot
        ! ------------------------------------------------------------------------------
        if(.not. found) then
            do i = 1, size(this%input_arr)
                if(.not. this%input_arr(i)%parameterSet) then
                    this%input_arr(i)%parName      = trim(paramName)
                    if(present(valType)) then
                        this%input_arr(i)%valType  = valType(:min(len(valType),LEN_PARAMETER_TYPE))
                    else
                        this%input_arr(i)%valType  = 'None'
                    end if
                    this%input_arr(i)%parameterSet = .true.
                    if(.not. newInput) then
                        this%input_arr(i)%val      = trim(paramValue)
                        this%input_arr(i)%valueSet = paramSet
                        if(present(info)) then
                            this%input_arr(i)%info = info(:min(len(info),LEN_PARAMETER_NAME))
                        else
                            this%input_arr(i)%info = ''
                        end if
                    else
                        this%input_arr(i)%val      = 'None'
                        this%input_arr(i)%valueSet = .false.
                        this%input_arr(i)%info     = ''
                    end if
                    found = .true.
                    exit
                end if
            end do

            ! if there was no new spot available, increase the size of the array
            if(.not. found) then
                allocate(tempInputArr(size(this%input_arr)))
                do i = 1,size(tempInputArr)
                    tempInputArr(i)%parName      = this%input_arr(i)%parName
                    tempInputArr(i)%val          = this%input_arr(i)%val
                    tempInputArr(i)%valType      = this%input_arr(i)%valType
                    tempInputArr(i)%info         = this%input_arr(i)%info
                    tempInputArr(i)%valueSet     = this%input_arr(i)%valueSet
                    tempInputArr(i)%parameterSet = this%input_arr(i)%parameterSet
                end do

                ! allocate new array
                deallocate(this%input_arr)
                allocate(this%input_arr(size(tempInputArr)+INPUT_ARR_DEFAULT_SIZE))
!                this%input_arr(:)%parName = 'None'
!                this%input_arr(:)%val     = 'None'
!                this%input_arr(:)%valType = 'None'
!                this%input_arr(:)%info    = ''
!                this%input_arr(:)%parameterSet = .false.                         !< this one allows to identify if the input parameter is defined (but it still could have an unset value (see below))
!                this%input_arr(:)%valueSet     = .false.                         !< is true if val is set
                do i = 1,size(tempInputArr)
                    this%input_arr(i)%parName      = tempInputArr(i)%parName
                    this%input_arr(i)%val          = tempInputArr(i)%val
                    this%input_arr(i)%valType      = tempInputArr(i)%valType
                    this%input_arr(i)%info         = tempInputArr(i)%info
                    this%input_arr(i)%valueSet     = tempInputArr(i)%valueSet
                    this%input_arr(i)%parameterSet = tempInputArr(i)%parameterSet
                end do

                ! now put new entry
                this%input_arr(size(tempInputArr)+1)%parName      = trim(paramName)
                if(present(valType)) then
                    this%input_arr(size(tempInputArr)+1)%valType      = valType(:min(len(valType),LEN_PARAMETER_TYPE))
                else
                    this%input_arr(size(tempInputArr)+1)%valType      = 'None'
                end if
                this%input_arr(size(tempInputArr)+1)%parameterSet = .true.
                if(.not. newInput) then
                    this%input_arr(size(tempInputArr)+1)%val      = trim(paramValue)
                    this%input_arr(size(tempInputArr)+1)%valueSet = .true.
                    if(present(info)) then
                        this%input_arr(size(tempInputArr)+1)%info = info(:min(len(info),LEN_PARAMETER_NAME))
                    else
                        this%input_arr(size(tempInputArr)+1)%info = ''
                    end if
                else
                    this%input_arr(size(tempInputArr)+1)%val      = 'None'
                    this%input_arr(size(tempInputArr)+1)%valueSet = .false.
                    this%input_arr(size(tempInputArr)+1)%info     = ''
                end if
                deallocate(tempInputArr)
            end if

        end if

        !** done!
        if(isControlled()) then
            call checkOut(csubid)
        end if
        return
    end subroutine

    !-----------------------------------------------------------------------------------------------
    !
    !> @anchor      setStartEpoch
    !!
    !> @brief       Set the start epoch of NEPTUNE
    !> @author      Vitali Braun
    !!
    !> @param[in]   sepoch start epoch as time_t
    !!
    !> @date        <ul>
    !!                <li>VB 23.04.2014 (initial design)</li>
    !!                <li>VB 14.07.2017 (added routines to set start epoch in downstream modules)</li>
    !!              </ul>
    !!
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine setStartEpoch(this,sepoch)

        class(Neptune_class)                    :: this
        type(time_t), intent(in)                :: sepoch

        this%start_epoch     = sepoch

        ! pass the message also to all the other dependent modules
        call this%output%write_to_output(C_EPOCH_START_GD, this%start_epoch)
        call this%derivatives_model%set_start_epoch(this%start_epoch)
        call this%numerical_integrator%set_start_epoch(this%start_epoch)

        this%isSetStartEpoch = .true.
        return

    end subroutine setStartEpoch

    !-----------------------------------------------------------------------------------------------
    !
    !> @anchor      setEndEpoch
    !!
    !> @brief       Set the end epoch of NEPTUNE
    !> @author      Vitali Braun
    !!
    !> @param[in]   eepoch end epoch as time_t
    !!
    !> @date        <ul>
    !!                <li> 23.04.2014: initial design</li>
    !!              </ul>
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine setEndEpoch(this,eepoch)

        class(Neptune_class)                    :: this
        type(time_t), intent(in)                :: eepoch

        this%end_epoch     = eepoch
        call this%output%write_to_output(C_EPOCH_START_GD, this%start_epoch)
        this%isSetEndEpoch = .true.

    end subroutine setEndEpoch


    !-----------------------------------------------------------------------------------------------
    !
    !> @anchor      getStartEpoch
    !!
    !> @brief       Get the start epoch, which is set in the initialization of NEPTUNE
    !> @author      Vitali Braun
    !!
    !> @param[out]   sepoch start epoch as time_t
    !> @param[out]   ierr -1 when no epoch is available
    !!
    !> @date        <ul>
    !!                <li> 23.04.2014: initial design</li>
    !!              </ul>
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine getStartEpoch(this,sepoch,ierr)

        class(Neptune_class)                    :: this
        type(time_t), intent(out)               :: sepoch
        integer,      intent(out)               :: ierr

        if(this%isSetStartEpoch) then
          sepoch = this%start_epoch
          ierr   = 0
        else
          sepoch = 0.d0
          ierr   = -1
        end if

        return

    end subroutine getStartEpoch

    !-----------------------------------------------------------------------------------------------
    !
    !> @anchor      getEndEpoch
    !!
    !> @brief       Get the end epoch, which is set in the initialization of NEPTUNE
    !> @author      Vitali Braun
    !!
    !> @param[out]   eepoch end epoch as time_t
    !> @param[out]   ierr -1 when no epoch is available
    !!
    !> @date        <ul>
    !!                <li> 23.04.2014: initial design</li>
    !!              </ul>
    !!
    !-----------------------------------------------------------------------------------------------
    subroutine getEndEpoch(this,eepoch,ierr)

        class(Neptune_class)                    :: this
        type(time_t), intent(out)               :: eepoch
        integer,      intent(out)               :: ierr

        if(this%isSetEndEpoch) then
          eepoch = this%end_epoch
          ierr   = 0
        else
          eepoch = 0.d0
          ierr   = -1
        end if

        return

    end subroutine getEndEpoch

    !==============================================================================================
    !
    !> @anchor      getRunID
    !!
    !> @brief       Get the NEPTUNE run ID
    !> @author      Vitali Braun
    !!
    !> @return      The run id as string with a length of 50 characters
    !!
    !> @date        <ul>
    !!                <li> 13.02.2014 (initial design)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    character(len=50) pure function getRunID(this)

        class(Neptune_class),intent(in)     :: this
        getRunID = this%runid
        return

    end function getRunID

    !==============================================================================================
    !
    !> @anchor      getDataPath
    !!
    !> @brief       Get the NEPTUNE data path
    !< @author      Vitali Braun
    !!
    !> @return      The data path as string with a length of 512 characters
    !!
    !< @date        <ul>
    !!                <li> 13.02.2014 (initial design)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    character(len=512) pure function getDataPath(this)

        class(Neptune_class),intent(in)     :: this
        getDataPath = this%cdatapath
        return

    end function getDataPath

    !==============================================================================================
    !
    !> @anchor      setDataPath
    !!
    !! @brief       Set the NEPTUNE data path
    !! @author      Christopher Kebschull
    !!
    !! @date        <ul>
    !!                <li> 24.10.2019 (initial design)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    subroutine setDataPath(this,data_path)

        class(Neptune_class),intent(inout)  :: this
        character(len=*),intent(in)         :: data_path

        this%cdatapath = trim(data_path)

    end subroutine setDataPath

    !==============================================================================================
    !
    !> @anchor      getInputPath
    !!
    !> @brief       Get the NEPTUNE input path
    !> @author      Vitali Braun
    !!
    !> @return      The input path as string with a length of 512 characters
    !!
    !> @date        <ul>
    !!                <li> 13.02.2014 (initial design)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    character(len=512) pure function getInputPath(this)

        class(Neptune_class),intent(in)     :: this
        getInputPath = this%cinpath
        return

    end function getInputPath

    !==============================================================================================
    !
    !> @anchor      getOutputPath
    !!
    !> @brief       Get the NEPTUNE output path
    !> @author      Vitali Braun
    !!
    !> @return      The output path as string with a length of 512 characters
    !!
    !> @date        <ul>
    !!                <li> 13.02.2014 (initial design)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    character(len=512) pure function getOutputPath(this)

        class(Neptune_class),intent(in)     :: this
        getOutputPath = this%coutpath
        return

    end function getOutputPath

    !==============================================================================================
    !
    !> @anchor      dump_input
    !!
    !> @brief       If switched via this subroutine, all input will be dumped to specified file
    !> @author      Vitali Braun
    !!
    !> @param[in]   dumpFile
    !!
    !> @date        <ul>
    !!                <li> 16.10.2016 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    subroutine dump_input(this,dumpFile)

        implicit none

        class(Neptune_class)                    :: this
        character(len=*), intent(in)            :: dumpFile

        character(len=*), parameter             :: csubid = 'dump_input'

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        if(len_trim(dumpFile) > len(this%dump_file_name)) then
            call setNeptuneError(E_STRING_LENGTH, FATAL)
            return
        else
            this%dump_file_name = trim(dumpFile)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        ! register in input array
        call this%set_input(parName=C_FILE_INPUT_DUMP, val=this%dump_file_name, set=.true.)

        ! toggle this one ON
        this%has_to_dump_input = .true.

        !** done!
        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

    end subroutine dump_input

    !==============================================================================================
    !
    !> @anchor      write_input_to_dump
    !!
    !> @brief       Write all input parameters to dump file
    !> @author      Vitali Braun
    !!
    !> @date        <ul>
    !!                <li> 18.10.2016 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    subroutine write_input_to_dump(this)

        implicit none

        class(Neptune_class)                    :: this

        !character(len=32)                       :: md5Sum
        character(len=LEN_PARAMETER_NAME)       :: cout
        integer                                 :: iunit                        ! unit the output is written to
        integer                                 :: i

        if(.not. this%has_to_dump_input) return

        iunit = openFile(this%dump_file_name, SEQUENTIAL, OUT_FORMATTED_OVERWRITE)

        write(iunit,'(a)') 'NEPTUNE ('//trim(this%version_model%to_string(this%version_model%get_neptune_version()))//') input dump on '//date2string(getDateTimeNow())
        ! now dump the input array
        do i = 1, size(this%input_arr)
            if(.not. this%input_arr(i)%parameterSet) cycle
            cout = trim(this%input_arr(i)%parName)//' = '//trim(this%input_arr(i)%val)
            if(len_trim(this%input_arr(i)%info) > 0) then
                cout = trim(cout)//' ('//trim(this%input_arr(i)%info)//')'
            end if
            write(iunit,'(a)') trim(cout)
        end do

        iunit = closeFile(iunit)
        return

    end subroutine write_input_to_dump

!==============================================================================================
!
!> @anchor      getHashMD5
!!
!> @brief       Computes the MD5 checksum for a given file
!> @author      Vitali Braun
!!
!> @param[in]   fname   ! file name
!!
!> @date        <ul>
!!                <li>VB: 23.10.2016 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
!    character(len=32) function getHashMD5(this,fname) result(hash)
!
!        implicit none
!
!        class(Neptune_class)                    :: this
!        character(len=*), intent(in)            :: fname
!
!        integer                                 :: iexit
!        integer                                 :: icmd
!        integer                                 :: ihash
!        integer                                 :: ios
!        character                               :: msg
!
!        character(len=len_trim(fname)+4)        :: md5TempFile
!        md5TempFile = trim(fname)//'.md5'
!
!        call execute_command_line('md5 -q '//trim(fname)//' > '//md5TempFile, exitstat=iexit, cmdstat=icmd, cmdmsg=msg)
!        if(iexit /= 0) then
!            hash = 'no hash'
!        else if(icmd /= 0) then
!            hash = 'no hash'
!        else
!            ihash = openFile(md5TempFile, SEQUENTIAL, IN_FORMATTED)
!            read(ihash,'(a)',iostat=ios) hash
!            if(ios /= 0) then
!                hash = 'no hash'
!            end if
!            ihash = closeFile(ihash)
!            call execute_command_line('rm '//md5TempFile, exitstat=iexit, cmdstat=icmd, cmdmsg=msg)
!        end if
!
!        return
!
!    end function

    !==============================================================================================
    !
    !> @anchor      read_input_from_dump
    !!
    !> @brief       Read all input from a dump file
    !> @author      Vitali Braun
    !!
    !> @param[in]   fdump   ! dump file the input is read from
    !!
    !> @date        <ul>
    !!                <li> 16.10.2016 (initial implementation)</li>
    !!              </ul>
    !!
    !!------------------------------------------------------------------------------------------------
    subroutine read_input_from_dump(this,fdump)

        implicit none

        class(Neptune_class)                    :: this
        character(len=*), intent(in)            :: fdump

        character(len=*), parameter             :: csubid = 'read_input_from_dump'
        character(len=3)                        :: c3
        character(len=100)                      :: cmatch
        character(len=100)                      :: key
        character(len=1900)                     :: val
        character(len=2000)                     :: line
        integer                                 :: ios, ierr
        integer                                 :: iunit

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        call message(' - Reading input from dump file: '//trim(fdump), LOG_AND_STDOUT)

        this%hasToReadInputFromDump = .true.
        iunit = openFile(fdump, SEQUENTIAL, IN_FORMATTED)

        do
            read(iunit,'(a)',iostat=ios) line
            if(IS_IOSTAT_END(ios)) then
                exit
            else if(ios /= 0) then
                write(c3,'(i3)') iunit
                call setNeptuneError(E_FILE_READ_ERROR, FATAL, (/c3/))
                return
            end if
            ! read version information
            if(index(line, 'NEPTUNE (') /= 0) then
                cmatch = line(index(line,'NEPTUNE (')+9:index(line,')')-1)
                if(hasFailed()) return
                if(trim(cmatch) /= trim(this%version_model%to_string(this%version_model%get_neptune_version()))) then
                    call message('   * WARNING: different NEPTUNE version in dump file: '//trim(cmatch), LOG_AND_STDOUT)
                end if
            else
                ! all other inputs are key value pairs, separated by '='
                key = line(1:index(line, '=') - 1)
                val = adjustl(line(index(line, '=') + 1:))
                ! now just set the input - but there are some special cases
                select case(trim(key))
                    case(C_FILE_EOP)
                        ! get the file version
                        cmatch = val(index(val, '(')+1:index(val,')')-1)
                        if(trim(cmatch) /= trim(this%reduction%getEopFileVersion(this%getDataPath()))) then
                            call message('   * WARNING: different EOP file version in dump file: '//trim(cmatch), LOG_AND_STDOUT)
                        end if
                    case(C_FILE_SOLMAG)
                        ! get the file version for the daily file
                        cmatch = val(index(val, '(')+1:index(val,')')-1)
                        if(trim(cmatch) /= trim(this%atmosphere_model%getSolmagFileVersion(DAILY))) then
                            call message('   * WARNING: different SOLMAG daily file version in dump file: '//trim(cmatch), LOG_AND_STDOUT)
                        end if
                    case(C_FILE_SOLMAG_MONTHLY)
                        ! get the file version for the monthly file
                        cmatch = val(index(val, '(')+1:index(val,')')-1)
                        if(trim(cmatch) /= trim(this%atmosphere_model%getSolmagFileVersion(MONTHLY))) then
                            call message('   * WARNING: different SOLMAG monthly file version in dump file: '//trim(cmatch), LOG_AND_STDOUT)
                        end if
                    case(C_INITIAL_STATE)
                        ierr = this%setNeptuneVar(trim(key), parse_state_from_string(trim(val)))
                    case(C_INITIAL_COVARIANCE)
                        ierr = this%setNeptuneVar(trim(key), parse_covariance_from_string(trim(val)))
                    case default
                        ierr = this%setNeptuneVar(trim(key), trim(val))
                        if(hasFailed()) return
                end select
            end if

        end do

        iunit = closeFile(iunit)

        !** done!
        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

    end subroutine read_input_from_dump

    !========================================================================
    !
    ! Set NEPTUNE parameters required for initialization
    !
    !-----------------------------------------------------
    !> @brief Set NEPTUNE initialization parameters
    !!
    !> @author      Vitali Braun
    !!
    !> @param[in]   key   key as string
    !> @param[in]   val   value as string
    !!
    !> @date        <ul>
    !!                <li> 20.04.2013 (first implementation)</li>
    !!                <li> 16.10.2013 (providing an interface to call with derived types as value)</li>
    !!              </ul>
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting parameters which are required for the initialization
    !!              procedure. Parameters are passed via two strings, the first
    !!              one being an identifier (key), the second one the value.
    !!
    !!              \par Overview
    !!
    !!              <ol>
    !!                <li> Parse input string </li>
    !!                <li> Set NEPTUNE parameter to given value  </li>
    !!                <li> Generate error message, if unknown parameter </li>
    !!                <li> Finish. </li>
    !!              </ol>
    !!
    !> @anchor   setNeptuneVar_char
    !!
    !> @todo     Using the CORRELATION_MATRIX option via the API does not work (array error?!?)
    !> @todo     Re-setting the start and end epoch does not affect a reset of the EOP and space weather (?) arrays
    !!
    !!------------------------------------------------------------------------------------------------
    integer function setNeptuneVar_char(        &
                                        this,   &
                                        key,    & ! <-- CHR() identifier string
                                        val     & ! <-- CHR() value
                                     )

        class(Neptune_class)                    :: this
        character(len=*), intent(in)            :: key
        character(len=*), intent(in)            :: val

        character(len=*), parameter  :: csubid = "setNeptuneVar"  ! ID
        character(len=255)           :: cmess     ! message string
        character(len=3)             :: ctemp     ! temporary string

        integer :: dy                             ! day
        integer :: hr                             ! hour
        integer :: ios                            ! I/O status
        integer :: itemp                          ! temporary
        integer :: mi                             ! minute
        integer :: mo                             ! month
        integer :: temp_len                       ! auxiliary holding the length of a string
        integer :: yr                             ! year
        logical :: ltemp                          ! temporary
        logical :: flag_init_status               ! saving initialization status right after call of this routine (see note below)
        real(dp)  :: dtemp                          ! temporary
        real(dp)  :: sc                             ! second

        setNeptuneVar_char = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call initializeInputArray(this)
        end if

        !=========================================================================================
        !
        ! Handle initialization flag, as some action requires a re-initialization of NEPTUNE,
        ! which is not done automatically. So, if a parameter can be set without requiring a
        ! re-initialization, the 'neptuneInitialized' flag can be set to its PRIOR status
        ! (it may be 'F' upon call of this routine!!). That is why the status is saved
        ! before the initialization flag is set to .false. (default behaviour shall be a
        ! re-initialization)
        !
        !-----------------------------------------------------------------------------------------
        flag_init_status   = this%neptuneInitialized
        this%neptuneInitialized = .false.

        select case(trim(adjustl(key)))

            !========================================
            !
            ! CASE for STRING parameters
            !
            !-------------------------
            case(C_RUN_ID)

                if(len_trim(val) > len(this%runid)) then ! too long...
                    call setNeptuneError(E_RUN_ID, WARNING)
                    setNeptuneVar_char = E_RUN_ID
                else
                    call this%set_input(parName=C_RUN_ID, val=val, set=.true.)
                    call this%output%write_to_output(C_RUN_ID, val)
                end if

            case(C_OPT_SHADOW)
                call this%set_input(parName=C_OPT_SHADOW, val=val, set=.true.)
                call this%radiation_model%setShadowModel(trim(val))

            !======================================
            !
            !  CASE for ON/OFF parameters
            !
            !---------------------------------
            case(C_ATMOSPHERE,     C_SUN, C_MOON,    C_SRP,            C_SOLID_TIDES,     &
                 C_MERCURY, C_VENUS, C_JUPITER, C_MARS, C_SATURN, C_URANUS, C_NEPTUNE,    &
                 C_OCEAN_TIDES,    C_MANEUVERS,      C_CORRELATION,    C_WIND,            &
                 C_ALBEDO,         C_OPT_EOP,        C_OPT_PN_LOOKUP,  C_OPT_INT_LOG,     &
                 C_OPT_SRP_CORRECT,  &
                 C_OUTPUT_FILES,   C_OUTPUT_ACC,     C_OUTPUT_ACG,     C_OUTPUT_ACD,      &
                 C_OUTPUT_ACM,     C_OUTPUT_ACS,     C_OUTPUT_ACR,     C_OUTPUT_ACA,      &
                 C_OUTPUT_AME,     C_OUTPUT_ACV,     C_OUTPUT_AMA,     C_OUTPUT_ACJ,      &
                 C_OUTPUT_ASA,     C_OUTPUT_ACU,     C_OUTPUT_ACN,                        &
                 C_OUTPUT_ACT,     C_OUTPUT_ATM,                                          &
                 C_OUTPUT_ACO,     C_OUTPUT_OSC,     C_OUTPUT_CSV,     C_OUTPUT_GLL,      &
                 C_COV_PROP,       C_COV_SUN,        C_COV_MOON,       C_COV_DRAG,        &
                 C_COV_SRP,        C_OUTPUT_VAR_ECI, C_OPT_HARMONICS, C_OPT_PROGRESS,     &
                 C_OUTPUT_VAR_UVW, C_OUTPUT_COV_ECI, C_OUTPUT_COV_UVW, C_OUTPUT_AMN)

                if(trim(adjustl(val)) == C_ON .or. trim(adjustl(toLowercase(val))) == 'true') then
                    itemp = SWITCHED_ON
                    ltemp = .true.
                else if(trim(adjustl(val)) == C_OFF .or. trim(adjustl(toLowercase(val))) == 'false') then
                    itemp = SWITCHED_OFF
                    ltemp = .false.
                else
                    cmess = "Input not accepted: '"//val//" for "//(trim(adjustl(key)))//"'."
                    call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
                    setNeptuneVar_char = E_SPECIAL
                    return
                end if

                select case(trim(adjustl(key)))

                    !** Perturbations (Geopotential is handled differently...)
                    case(C_ATMOSPHERE)
                        !** check whether solar & geomagnetic data are already available...if so then NEPTUNE does not
                        !   need to be re-initialized, if switched off or on again.
                        if(this%atmosphere_model%getAtmosphereInitFlag()) then
                            this%neptuneInitialized = flag_init_status
                        end if
                        call this%set_input(parName = C_ATMOSPHERE, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_DRAG, ltemp)

                    case(C_WIND)
                        call this%set_input(parName = C_WIND, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_WIND, ltemp)

                    case(C_SUN)
                        call this%set_input(parName = C_SUN, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_SUN, ltemp)

                    case(C_MOON)
                        call this%set_input(parName = C_MOON, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_MOON, ltemp)

                    case(C_MERCURY)
                        call this%set_input(parName = C_MERCURY, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_MERCURY, ltemp)

                    case(C_VENUS)
                        call this%set_input(parName = C_VENUS, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_VENUS, ltemp)

                    case(C_MARS)
                        call this%set_input(parName = C_MARS, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_MARS, ltemp)

                    case(C_JUPITER)
                        call this%set_input(parName = C_JUPITER, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_JUPITER, ltemp)

                    case(C_SATURN)
                        call this%set_input(parName = C_SATURN, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_SATURN, ltemp)

                    case(C_URANUS)
                        call this%set_input(parName = C_URANUS, val= toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_URANUS, ltemp)

                    case(C_NEPTUNE)
                        call this%set_input(parName = C_NEPTUNE, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_NEPTUNE, ltemp)

                    case(C_SRP)
                        call this%set_input(parName = C_SRP, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_SRP, ltemp)

                    case(C_ALBEDO)
                        call this%set_input(parName = C_ALBEDO, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_ALBEDO, ltemp)

                    case(C_SOLID_TIDES)
                        call this%set_input(parName=C_SOLID_TIDES, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_SOLID_TIDES, ltemp)

                    case(C_OCEAN_TIDES)
                        call this%set_input(parName=C_OCEAN_TIDES, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_OCEAN_TIDES, ltemp)

                    case(C_MANEUVERS)
                        call this%set_input(parName=C_MANEUVERS, val = toString(ltemp), set=.true.)
                        call this%derivatives_model%setPertSwitch(this%gravity_model, PERT_MANEUVERS, ltemp)

                  !** covariance matrix propagation
                  case(C_COV_PROP)
                        call this%set_input(parName=C_COV_PROP, val=val, set=.true.)

                  if(itemp == SWITCHED_ON) then
                    call this%numerical_integrator%setCovariancePropagationFlag(.true.)
                  else
                    call this%numerical_integrator%setCovariancePropagationFlag(.false.)
                  end if

                  !** set the initialization status prior to routine call,
                  !   as switching on/off covariance matrix propagation does
                  !   not require a re-initialization in any case.
                  this%neptuneInitialized = flag_init_status

                case(C_COV_DRAG)
                  call this%set_input(parName=C_COV_DRAG, val=val, set=.true.)

                  if(itemp == SWITCHED_ON) then
                    call this%atmosphere_model%setDragCovFlag(.true.)
                  else
                    call this%atmosphere_model%setDragCovFlag(.false.)
                  end if

                case(C_COV_SRP)
                  call this%set_input(parName=C_COV_SRP, val=val, set=.true.)

                  if(itemp == SWITCHED_ON) then
                    call this%radiation_model%setSrpCovFlag(.true.)
                  else
                    call this%radiation_model%setSrpCovFlag(.false.)
                  end if

                case(C_COV_SUN)
                  call this%set_input(parName=C_COV_SUN, val=val, set=.true.)

                  if(itemp == SWITCHED_ON) then
                    call this%thirdbody_model%setThirdBodyCovFlag(ID_SUN, .true.)
                  else
                    call this%thirdbody_model%setThirdBodyCovFlag(ID_SUN, .false.)
                  end if

                case(C_COV_MOON)
                  call this%set_input(parName=C_COV_MOON, val=val, set=.true.)

                  if(itemp == SWITCHED_ON) then
                    call this%thirdbody_model%setThirdBodyCovFlag(ID_MOON, .true.)
                  else
                    call this%thirdbody_model%setThirdBodyCovFlag(ID_MOON, .false.)
                  end if

                !** correlation matrix
                case(C_CORRELATION)
                  call this%set_input(parName=C_CORRELATION, val=val, set=.true.)

                  this%iopt(OPTION_CORRELATION) = itemp
                  if(itemp == SWITCHED_ON) then
                    call this%correlation_model%setNoisePropagationFlag(.true.)
                  else
                    call this%correlation_model%setNoisePropagationFlag(.false.)
                  end if

                case(C_OPT_EOP)
                  call this%set_input(parName=C_OPT_EOP, val=val, set=.true.)

                  !** check whether EOP are already available...if so then NEPTUNE does not
                  !   need to be re-initialized, if switched off or on again.
                  if(this%reduction%getEopInitFlag()) then
                    this%neptuneInitialized = flag_init_status
                  end if

                  this%iopt(OPTION_EOP) = itemp

                  if(itemp == SWITCHED_ON) then
                    call this%reduction%setEopFlag(.true.)
                  else
                    call this%reduction%setEopFlag(.false.)
                  end if

                case(C_OPT_PN_LOOKUP)
                  call this%set_input(parName=C_OPT_PN_LOOKUP, val=val, set=.true.)

                  this%iopt(OPTION_PN_LOOKUP) = itemp

                  if(itemp == SWITCHED_ON) then
                    call this%reduction%setPNLookup(.true.)
                  else
                    call this%reduction%setPNLookup(.false.)
                  end if

                case(C_OPT_INT_LOG)
                  call this%set_input(parName=C_OPT_INT_LOG, val=val, set=.true.)
                  this%iopt(OPTION_INT_LOGFILE) = itemp

                case(C_OPT_SRP_CORRECT)
                  call this%set_input(parName=C_OPT_SRP_CORRECT, val=val, set=.true.)
                  this%iopt(OPTION_SRP_CORRECT) = itemp

                  if(itemp == SWITCHED_ON) then
                    call this%numerical_integrator%setSrpCorrect(.true.)
                  else
                    call this%numerical_integrator%setSrpCorrect(.false.)
                  end if

                case(C_OPT_HARMONICS)
                  call this%set_input(parName=C_OPT_HARMONICS, val=val, set=.true.)
                  this%iopt(OPTION_HARMONICS) = itemp

                  if(itemp == SWITCHED_ON) then
                    call this%gravity_model%setDistinctHarmonics(.true.)
                  else
                    call this%gravity_model%setDistinctHarmonics(.false.)
                  end if

                case(C_OPT_PROGRESS)
                  call this%set_input(parName=C_OPT_PROGRESS, val=val, set=.true.)
                  this%write_progress = ltemp

                !** Output files
                case(C_OUTPUT_FILES)
                  call this%set_input(parName=C_OUTPUT_FILES, val=toString(ltemp), set=.true.)
                  call this%output%switch_output(ltemp)
                case(C_OUTPUT_ACC)
                  call this%set_input(parName=C_OUTPUT_ACC, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACC, switch = ltemp)
                case(C_OUTPUT_ACG)
                  call this%set_input(parName=C_OUTPUT_ACG, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACG, switch = ltemp)
                case(C_OUTPUT_ACD)
                  call this%set_input(parName=C_OUTPUT_ACD, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACD, switch = ltemp)
                case(C_OUTPUT_ACM)
                  call this%set_input(parName=C_OUTPUT_ACM, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACM, switch = ltemp)
                case(C_OUTPUT_ACS)
                  call this%set_input(parName=C_OUTPUT_ACS, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACS, switch = ltemp)
                case(C_OUTPUT_AME)
                  call this%set_input(parName=C_OUTPUT_AME, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_AME, switch = ltemp)
                case(C_OUTPUT_ACV)
                  call this%set_input(parName=C_OUTPUT_ACV, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACV, switch = ltemp)
                case(C_OUTPUT_AMA)
                  call this%set_input(parName=C_OUTPUT_AMA, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_AMA, switch = ltemp)
                case(C_OUTPUT_ACJ)
                  call this%set_input(parName=C_OUTPUT_ACJ, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACJ, switch = ltemp)
                case(C_OUTPUT_ASA)
                  call this%set_input(parName=C_OUTPUT_ASA, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ASA, switch = ltemp)
                case(C_OUTPUT_ACN)
                  call this%set_input(parName=C_OUTPUT_ACN, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACN, switch = ltemp)
                case(C_OUTPUT_ACU)
                  call this%set_input(parName=C_OUTPUT_ACU, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACU, switch = ltemp)
                case(C_OUTPUT_ACR)
                  call this%set_input(parName=C_OUTPUT_ACR, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACR, switch = ltemp)
                case(C_OUTPUT_ACA)
                  call this%set_input(parName=C_OUTPUT_ACA, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACA, switch = ltemp)
                case(C_OUTPUT_ACT)
                  call this%set_input(parName=C_OUTPUT_ACT, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACT, switch = ltemp)
                case(C_OUTPUT_ACO)
                  call this%set_input(parName=C_OUTPUT_ACO, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ACO, switch = ltemp)
                case(C_OUTPUT_AMN)
                  call this%set_input(parName=C_OUTPUT_AMN, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_AMN, switch = ltemp)
                case(C_OUTPUT_ATM)
                  call this%set_input(parName=C_OUTPUT_ATM, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_ATM, switch = ltemp)
                case(C_OUTPUT_OSC)
                  call this%set_input(parName=C_OUTPUT_OSC, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_OSC, switch = ltemp)
                case(C_OUTPUT_CSV)
                  call this%set_input(parName=C_OUTPUT_CSV, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_CSV, switch = ltemp)
                case(C_OUTPUT_GLL)
                  call this%set_input(parName=C_OUTPUT_GLL, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_GLL, switch = ltemp)
                case(C_OUTPUT_VAR_ECI)
                  call this%set_input(parName=C_OUTPUT_VAR_ECI, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_VAR_ECI, switch = ltemp)
                case(C_OUTPUT_VAR_UVW)
                  call this%set_input(parName=C_OUTPUT_VAR_UVW, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_VAR_UVW, switch = ltemp)
                case(C_OUTPUT_COV_ECI)
                  call this%set_input(parName=C_OUTPUT_COV_ECI, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_COV_ECI, switch = ltemp)
                case(C_OUTPUT_COV_UVW)
                  call this%set_input(parName=C_OUTPUT_COV_UVW, val=toString(ltemp), set=.true.)
                  call this%output%set_output_file(par_name=C_OUTPUT_COV_UVW, switch = ltemp)
              end select

          !========================================
          !
          ! CASE for INTEGER parameters
          !
          !------------------------
          case(C_GEOPOTENTIAL,  C_OPT_SAT_PROPERTIES, C_OUTPUT_STEP,    C_PAR_INT_COV_METHOD, &
               C_OPT_GEO_MODEL, C_OPT_AP_FORECAST,    C_OPT_STORE_DATA, C_COV_GEOPOTENTIAL,   &
               C_PAR_INT_METHOD, C_OPT_ATMOSPHERE_MODEL)

            read(val,*,iostat=ios) itemp

            if(ios /= 0) then

              cmess = "Format error for key '"//key//"'. Integer expected but got "//toString(itemp)
              call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
              setNeptuneVar_char = E_SPECIAL
              return

            else

              select case(trim(adjustl(key)))

                !** Perturbations
                case(C_GEOPOTENTIAL)
                  call this%set_input(parName = C_GEOPOTENTIAL, val=val, set=.true.)
                  call this%gravity_model%setGeoDegree(itemp)
                  call this%derivatives_model%setPertSwitch(this%gravity_model,PERT_GRAVITY, .true.)  ! geopotential is never switched OFF!
                  this%neptuneInitialized = flag_init_status

                !** Options
                case(C_OPT_AP_FORECAST)
                  call this%set_input(parName=C_OPT_AP_FORECAST, val=val, set=.true.)
                  call this%atmosphere_model%setApForecast(itemp)

                case(C_OPT_SAT_PROPERTIES)
                  call this%set_input(parName=C_OPT_SAT_PROPERTIES, val=val, set=.true.)
                  call this%output%write_to_output(C_OPT_SAT_PROPERTIES, itemp)
                  this%iopt(OPTION_PROPS) = itemp

                case(C_OPT_GEO_MODEL)
                  itemp = max(1,min(3,itemp))
                  write(ctemp,'(i1)') itemp
                  call this%set_input(parName=C_OPT_GEO_MODEL, val=trim(ctemp), set=.true.)
                  this%iopt(OPTION_GEO_MODEL) = itemp

                case(C_OPT_ATMOSPHERE_MODEL)
                  call this%set_input(parName=C_OPT_ATMOSPHERE_MODEL, val=val, set=.true.)
                  call this%atmosphere_model%setAtmosphereModel(itemp)

                case(C_OPT_STORE_DATA)
                  call this%set_input(parName=C_OPT_STORE_DATA, val=val, set=.true.)
                  call this%setStep(itemp)

                !** Output
                case(C_OUTPUT_STEP)
                  call this%set_input(parName=C_OUTPUT_STEP, val=val, set=.true.)
                  call this%output%write_to_output(C_OUTPUT_STEP, itemp)
                  this%output_step = itemp

                !** state vector integration method
                case(C_PAR_INT_METHOD)
                  call this%set_input(parName=C_PAR_INT_METHOD, val=val, set=.true.)
                  call this%numerical_integrator%setIntegrationMethod(itemp)

                !** Covariance matrix integration method
                case(C_PAR_INT_COV_METHOD)
                  call this%set_input(parName=C_PAR_INT_COV_METHOD, val=val, set=.true.)
                  this%iopt(OPTION_COV_METHOD) = itemp

                !** Geopotential covariance propagation degree
                case(C_COV_GEOPOTENTIAL)
                  call this%set_input(parName=C_COV_GEOPOTENTIAL, val=val, set=.true.)
                  call this%gravity_model%setGeoCovDegree(itemp)

              end select

            end if

          !========================================
          !
          ! CASE for REAL parameters
          !
          !-------------------------
          case(C_PAR_MASS, C_PAR_CROSS_SECTION, C_PAR_CDRAG,         &
               C_PAR_CREFL, C_PAR_REENTRY, C_PAR_INT_RELEPS, C_PAR_INT_ABSEPS, &
               C_PAR_INT_COV_STEP, C_PAR_EARTH_RADIUS, &
               C_OPT_SOL_FORECAST)

            read(val,*,iostat=ios) dtemp

            if(ios /= 0) then

              cmess = "Format error for key '"//key//"'. Float expected but got "//toString(dtemp)
              call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
              setNeptuneVar_char = E_SPECIAL
              return

            else    !** set parameters and force re-initialization of NEPTUNE

              select case(trim(adjustl(key)))
                case(C_PAR_MASS)
                  call this%set_input(parName=C_PAR_MASS, val=val, set=.true.)
                  call this%output%write_to_output(C_PAR_MASS, dtemp)
                  this%ipar(PAR_MASS) = dtemp
                  this%isSetMass      = .true.

                case(C_PAR_CROSS_SECTION)
                  call this%set_input(parName=C_PAR_CROSS_SECTION, val=val, set=.true.)
                  call this%output%write_to_output(C_PAR_CROSS_SECTION, dtemp)
                  this%ipar(PAR_CROSS_SECTION) = dtemp
                  this%isSetArea = .true.

                case(C_PAR_CDRAG)
                  call this%set_input(parName=C_PAR_CDRAG, val=val, set=.true.)
                  call this%output%write_to_output(C_PAR_CDRAG, dtemp)
                  this%ipar(PAR_CDRAG) = dtemp
                  this%isSetCdrag      = .true.

                case(C_PAR_CREFL)
                  call this%set_input(parName=C_PAR_CREFL, val=val, set=.true.)
                  call this%output%write_to_output(C_PAR_CREFL, dtemp)
                  this%ipar(PAR_CREFL) = dtemp
                  this%isSetCsrp       = .true.

                case(C_PAR_REENTRY)
                  call this%set_input(parName=C_PAR_REENTRY, val=val, set=.true.)
                  call this%atmosphere_model%setMinAltitude(dtemp)
                  this%neptuneInitialized = flag_init_status ! no re-initialization required after changing min. altitude
                  !this%ipar(PAR_REENTRY)    = dtemp

                case(C_PAR_EARTH_RADIUS)
                  call this%set_input(parName=C_PAR_EARTH_RADIUS, val=val, set=.true.)
                  call setEarthRadius(dtemp)

                case(C_PAR_INT_RELEPS)
                  call this%set_input(parName=C_PAR_INT_RELEPS, val=val, set=.true.)
                  this%ipar(PAR_INT_RELEPS)    = dtemp
                  this%flag_init_tolerances(1) = .true.

                case(C_PAR_INT_ABSEPS)
                  call this%set_input(parName=C_PAR_INT_ABSEPS, val=val, set=.true.)
                  this%ipar(PAR_INT_ABSEPS)    = dtemp
                  this%flag_init_tolerances(2) = .true.

                case(C_PAR_INT_COV_STEP)
                  call this%set_input(parName=C_PAR_INT_COV_STEP, val=val, set=.true.)
                  this%ipar(PAR_INT_COV_STEP) = dtemp

                case(C_OPT_SOL_FORECAST)
                  call this%set_input(parName=C_OPT_SOL_FORECAST, val=val, set=.true.)
                  call this%atmosphere_model%setSolarFluxForecast(dtemp)

              end select

            end if

          !========================================
          !
          ! CASE for DATE parameters
          !
          !-------------------------
          case(C_EPOCH_START_GD, C_EPOCH_END_GD)
            call tokenizeDate(val, yr, mo, dy, hr, mi, sc)
            if(hasFailed()) then
              setNeptuneVar_char = getLatestError()
              return
            end if

            select case(trim(adjustl(key)))
              case(C_EPOCH_START_GD)
                this%start_epoch%year   = yr
                this%start_epoch%month  = mo
                this%start_epoch%day    = dy
                this%start_epoch%hour   = hr
                this%start_epoch%minute = mi
                this%start_epoch%second = sc
                call gd2mjd(this%start_epoch)
                call this%set_input(parName=C_EPOCH_START_GD, val=val, set=.true.)
                call this%output%write_to_output(C_EPOCH_START_GD, this%start_epoch)
                this%isSetStartEpoch    = .true.
              case(C_EPOCH_END_GD)
                this%end_epoch%year   = yr
                this%end_epoch%month  = mo
                this%end_epoch%day    = dy
                this%end_epoch%hour   = hr
                this%end_epoch%minute = mi
                this%end_epoch%second = sc
                call gd2mjd(this%end_epoch)
                call this%set_input(parName=C_EPOCH_END_GD, val=val, set=.true.)
                call this%output%write_to_output(C_EPOCH_END_GD, this%end_epoch)
                this%isSetEndEpoch    = .true.
            end select

          !========================================
          !
          ! CASE for data path
          !
          !-------------------------
          case(C_PATH_DATA)

            if(len_trim(val) > len(this%cdatapath)) then ! too long...
              call setNeptuneError(E_DATA_PATH, FATAL)
              setNeptuneVar_char = E_DATA_PATH
              return
            else
              read(val(1:min(len(val),len(this%cdatapath))),*) this%cdatapath
              call this%set_input(parName=C_PATH_DATA, val=val, set=.true.)
            end if

          !========================================
          !
          ! CASE for input path
          !
          !-------------------------
          case(C_PATH_INPUT)

            if(len_trim(val) > len(this%cinpath)) then ! too long...
              call setNeptuneError(E_INPUT_PATH, FATAL)
              setNeptuneVar_char = E_INPUT_PATH
              return
            else
              read(val(1:min(len(val),len(this%cinpath))),*) this%cinpath
              call this%set_input(parName=C_PATH_INPUT, val=val, set=.true.)
            end if

          !========================================
          !
          ! CASE for output path
          !
          !-------------------------
          case(C_PATH_OUTPUT)

            if(len_trim(val) > len(this%coutpath)) then ! too long...
              call setNeptuneError(E_OUTPUT_PATH, FATAL)
              setNeptuneVar_char = E_OUTPUT_PATH
              return
            else
              read(val(1:min(len(val),len(this%coutpath))),*) this%coutpath
              call this%set_input(parName=C_PATH_OUTPUT, val=val, set=.true.)
            end if

          !========================================
          !
          ! CASE for NEPTUNE input dump file
          !
          !-------------------------
          case(C_FILE_INPUT_DUMP)
            ! file is checked in dumpInput
            call this%dump_input(val)

          !========================================
          !
          ! CASE for maneuvers definition file
          !
          !-------------------------
          case(C_FILE_MANEUVERS)

            temp_len = len(this%manoeuvres_model%get_man_file_name())

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%manoeuvres_model%set_man_file_name(val(1:min(temp_len,len_trim(val))))
              call this%set_input(parName=C_FILE_MANEUVERS, val=val, set=.true.)
            end if

          !========================================
          !
          ! CASE for surfaces definition file
          !
          !-------------------------
          case(C_FILE_SURFACES)

            temp_len = len(this%satellite_model%getSurfaceDefinitionFileName())

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%satellite_model%setSurfaceDefinitionFileName(val(1:min(temp_len,len_trim(val))))
              call this%set_input(parName=C_FILE_SURFACES, val=val, set=.true.)
            end if

          !=====================================================
          !
          ! CASE for Earth orientation parameters
          !
          !--------------------------------------------------
          case(C_FILE_EOP)

            temp_len = len(this%reduction%getEopFileName())

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%reduction%setEopFileName(val(1:min(temp_len,len_trim(val))))
              call this%set_input(parName=C_FILE_EOP, val=val, set=.true., info=this%reduction%getEopFileVersion(this%getDataPath()))
            end if

          !=====================================================
          !
          ! CASE for solar and geomagnetic activity data file
          !
          !--------------------------------------------------
          case(C_FILE_SOLMAG)

            temp_len = len(this%atmosphere_model%getSolMagFileName(DAILY))

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%atmosphere_model%setSolMagFileName(val(1:min(temp_len,len_trim(val))),DAILY)
              call this%set_input(parName=C_FILE_SOLMAG, val=val, set=.true., info=this%atmosphere_model%getSolmagFileVersion(DAILY))
            end if

          !=====================================================
          !
          ! CASE for monthly solar and geomagnetic activity data file
          !
          !--------------------------------------------------
          case(C_FILE_SOLMAG_MONTHLY)

            temp_len = len(this%atmosphere_model%getSolMagFileName(MONTHLY))

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%atmosphere_model%setSolMagFileName(val(1:min(temp_len,len_trim(val))),MONTHLY)
              call this%set_input(parName=C_FILE_SOLMAG_MONTHLY, val=val, set=.true., info=this%atmosphere_model%getSolmagFileVersion(MONTHLY))
            end if

          !========================================
          !
          ! CASE for disturbance winds file (HWM07)
          !
          !-------------------------
          case(C_FILE_DWIND)

            temp_len = len(this%atmosphere_model%getDisturbanceWindFileName())

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%atmosphere_model%setDisturbanceWindFileName(val(1:min(temp_len,len_trim(val))))
              call this%set_input(parName=C_FILE_DWIND, val=val, set=.true.)
            end if

          !===================================================
          !
          ! CASE for horizontal wind model data file (HWM07)
          !
          !-------------------------------------------------
          case(C_FILE_HWIND)

            temp_len = len(this%atmosphere_model%getHorizontalWindFileName())

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%atmosphere_model%setHorizontalWindFileName(val(1:min(temp_len,len_trim(val))))
              call this%set_input(parName=C_FILE_HWIND, val=val, set=.true.)
            end if

          !=================================================
          !
          ! CASE for QD interpolation grid file (HWM07)
          !
          !-----------------------------------------------
          case(C_FILE_QDGRID)

            temp_len = len(this%atmosphere_model%getQDGridFileName())

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              call this%atmosphere_model%setQDGridFileName(val(1:min(temp_len,len_trim(val))))
              call this%set_input(parName=C_FILE_QDGRID, val=val, set=.true.)
            end if

          !=================================================
          !
          ! CASE for progress file name
          !
          !-----------------------------------------------
          case(C_FILE_PROGRESS)

            temp_len = len(this%progress_file_name)

            if(len_trim(val) > temp_len) then ! too long...
              write(ctemp,'(i3)') temp_len
              call setNeptuneError(E_FILE_LENGTH, FATAL, (/val, ctemp/))
              setNeptuneVar_char = E_FILE_LENGTH
              return
            else
              this%progress_file_name = val(1:min(temp_len,len_trim(val)))
              call this%set_input(parName=C_FILE_PROGRESS, val=val, set=.true.)
            end if

          !========================================
          !
          ! Default CASE
          !
          !-------------------------

          case default
            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_char = E_UNKNOWN_PARAMETER
            return
        end select

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if
        return

    end function setNeptuneVar_char

    !==============================================================================================
    !
    !> @brief Set NEPTUNE initialization parameters provided as derived type of kepler elements
    !!
    !> @author      Vitali Braun
    !!
    !> @param[in]   key   key as string
    !> @param[in]   val   value as type \ref kepler_t
    !!
    !> @date        <ul>
    !!                <li> 16.10.2013 (first implementation)</li>
    !!              </ul>
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting the initial orbit of an object during the initialization
    !!              procedure. The orbit is passed via a string being the key and a kepler element
    !!              set as the second argument (derived type 'kepler_t').
    !!
    !!              \par Overview
    !!
    !!              <ol>
    !!                <li> Parse input string </li>
    !!                <li> Set NEPTUNE initial orbit to given value  </li>
    !!                <li> Generate error message, if invalid data </li>
    !!                <li> Finish. </li>
    !!              </ol>
    !!
    !> @anchor   setnepvar_kepl
    !!
    !> @todo     check dimensions of sma and angles!
    !!
    !!------------------------------------------------------------------------------------------------
    integer function setNeptuneVar_kepl(        &
                                        this,   &
                                        key,    & ! <-- CHR() identifier string
                                        val     & ! <-- TYP() value
                                     )

        class(Neptune_class)                    :: this
        character(len=*), intent(in)            :: key
        type(kepler_t),   intent(in)            :: val

        character(len=*), parameter  :: csubid = "setNeptuneVar_kepl"  ! ID

        type(kepler_t) :: val_temp                  ! temporary array to hold kepler state (conversions possible, as opposed to val, which has intent(in) attribute.)

        setNeptuneVar_kepl = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        select case(key)

          case(C_INITIAL_STATE)

            !** check eccentricity
            if(val%ecc < 0.d0) then

              call setNeptuneError(E_ECCENTRICITY, FATAL)
              setNeptuneVar_kepl = E_ECCENTRICITY
              return

            end if

            !** check semi-major axis
            if(val%sma < 0.d0) then

              call setNeptuneError(E_SEMI_MAJOR_AXIS, FATAL)
              setNeptuneVar_kepl = E_SEMI_MAJOR_AXIS
              return

            end if

            !** check for altitude being above Earth's surface
            if(val%sma*(1.d0-val%ecc**2.d0)/(1.d0+val%ecc*cos(val%tran)) < getEarthRadius()) then

              call setNeptuneError(E_ALTITUDE, FATAL)
              setNeptuneVar_kepl = E_ALTITUDE
              return

            end if

            !** check inclination
            if(val%inc < 0.d0) then

              call setNeptuneError(E_INCLINATION, FATAL)
              setNeptuneVar_kepl = E_INCLINATION
              return

            end if

            val_temp = val

            !** convert angles to radians if required
            if(val%angles_unit == UNIT_DEG) then

              call convertToRadians(val_temp)

            end if

            this%initial_orbit_osc = val_temp
            call this%output%write_to_output(val_temp)

            !** check dimensions
            !---------------------------------
            !if(sma_unit

          case default

            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_kepl = E_UNKNOWN_PARAMETER
            return

        end select

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if

        return

    end function setNeptuneVar_kepl

    !!==============================================================================================
    !> @brief Set NEPTUNE initialization parameters provided as derived type of state vector
    !!
    !> @author      Vitali Braun
    !!
    !> @param[in]   key   key as string
    !> @param[in]   val   value as type \ref state_t
    !!
    !> @date        <ul>
    !!                <li> 16.10.2013 (first implementation)</li>
    !!              </ul>
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting the initial orbit of an object during the initialization
    !!              procedure. The orbit is passed via a string being the key and a cartesian state vector
    !!              set as the second argument (derived type 'state_t').
    !!
    !!              \par Overview
    !!
    !!              <ol>
    !!                <li> Parse input string </li>
    !!                <li> Set NEPTUNE initial orbit to given value  </li>
    !!                <li> Generate error message, if invalid data </li>
    !!                <li> Finish. </li>
    !!              </ol>
    !!
    !> @anchor   setNeptuneVar_stat
    !!
    !!------------------------------------------------------------------------------------------------
    integer function setNeptuneVar_stat(        &
                                        this,   &
                                        key,    & ! <-- CHR() identifier string
                                        val     & ! <-- TYP() value
                                     )

        class(Neptune_class)                    :: this
        character(len=*), intent(in)            :: key
        type(state_t),   intent(in)             :: val

        character(len=*), parameter  :: csubid = "setNeptuneVar_stat"  ! ID

        setNeptuneVar_stat = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        select case(key)

          case(C_INITIAL_STATE)

            !** check for altitude being above Earth's surface
            if(mag(val%r) < getEarthRadius()) then

              call setNeptuneError(E_ALTITUDE, FATAL)
              setNeptuneVar_stat = E_ALTITUDE
              return

            end if

            this%initial_orbit_csv = val
            call this%output%write_to_output(val)
            call this%set_input(parName=C_INITIAL_STATE, val=toString(val), set=.true.)

          case default

            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_stat = E_UNKNOWN_PARAMETER
            return

        end select

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if

        return

    end function setNeptuneVar_stat

    !!==============================================================================================
    !> @brief Set NEPTUNE initialization parameters provided as derived type of covariance type
    !!
    !> @author      Vitali Braun
    !!
    !> @param[in]   key   key as string
    !> @param[in]   val   value as type \ref covariance_t
    !!
    !> @date        <ul>
    !!                <li> 21.10.2013 (first implementation)</li>
    !!              </ul>
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting the initial covariance matrix of an object during the initialization
    !!              procedure. The orbit is passed via a string being the key and a covariance matrix
    !!              as the second argument (derived type 'covariance_t').
    !!
    !!              \par Overview
    !!
    !!              <ol>
    !!                <li> Parse input string </li>
    !!                <li> Set NEPTUNE initial covariance to given value  </li>
    !!                <li> Generate error message, if invalid data </li>
    !!                <li> Finish. </li>
    !!              </ol>
    !!
    !> @anchor   setnepvar_covr
    !!
    !!------------------------------------------------------------------------------------------------
    integer function setNeptuneVar_covr(        &
                                        this,   &
                                        key,    & ! <-- CHR() identifier string
                                        val     & ! <-- TYP() value
                                     )

        class(Neptune_class)                    :: this
        character(len=*),   intent(in)          :: key
        type(covariance_t), intent(in)          :: val

        character(len=*), parameter  :: csubid = "setNeptuneVar_covr"  ! ID

        setNeptuneVar_covr = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        select case(key)
          case(C_INITIAL_COVARIANCE)
            call this%set_input(parName=C_INITIAL_COVARIANCE, val=toString(val), set=.true.)
            this%initial_covariance = val
            call this%output%write_to_output(val)
          case default
            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_covr = E_UNKNOWN_PARAMETER
            return
        end select

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if
        return

    end function setNeptuneVar_covr

    !==============================================================================================
    !> @brief       Toggle harmonic coefficients ON or OFF
    !!
    !< @author      Vitali Braun
    !!
    !> @date        <ul>
    !!                <li> 17.01.2016 (first implementation)</li>
    !!              </ul>
    !!
    !> @param[in]   key     identifies whether coefficient C or the coefficient S are to be set
    !> @param[in]   n       index n of coefficient (degree)
    !< @param[in]   m       index m of coefficient (order)
    !> @param[in]   switch  value which shall be set
    !!
    !> @returns     error code
    !!
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              toggling the harmonic coefficients of the geopotential.
    !!              As argument it receives the degree (n) and
    !!              order (m) as well as either a '.true.' or '.false.' to set the
    !!              coefficient. If the optional switch is omitted, the currently set value is toggled.
    !!              In principle, this is a wrapper to the toggler in the gravity module.
    !!
    !> @anchor      setnepvar_geop_toggle
    !!
    !!----------------------------------------------------------------------------------
    integer function setNeptuneVar_geop_toggle( &
                                        this,   &
                                        key,    &
                                        n,      &
                                        m,      &
                                        switch)

        class(Neptune_class)                    :: this
        character(len=*),  intent(in)           :: key
        integer,           intent(in)           :: n
        integer,           intent(in)           :: m
        logical, optional, intent(in)           :: switch

        character(len=*), parameter             :: csubid = "setNeptuneVar_geop_toggle"  ! ID

        setNeptuneVar_geop_toggle = 0

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        ! (the input is checked in the toggle routine)
        select case(key)
            case(C_HARMONIC_C)
                if(present(switch)) then
                    call this%gravity_model%toggleHarmonicCoefficient(COEFF_C, n, m, switch)
                else
                    call this%gravity_model%toggleHarmonicCoefficient(COEFF_C, n, m)
                end if
            case(C_HARMONIC_S)
                if(present(switch)) then
                    call this%gravity_model%toggleHarmonicCoefficient(COEFF_S, n, m, switch)
                else
                    call this%gravity_model%toggleHarmonicCoefficient(COEFF_S, n, m)
                end if
            case default
                call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
                setNeptuneVar_geop_toggle = E_UNKNOWN_PARAMETER
                return
        end select

        if(hasFailed()) then
            setNeptuneVar_geop_toggle = getLatestError()
            return
        end if

        !** done!
        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

    end function setNeptuneVar_geop_toggle

    !==============================================================================================
    !> @brief       Set harmonic coefficients or their standard deviations of the geopotential in NEPTUNE
    !!
    !> @author      Vitali Braun
    !!
    !> @date        <ul>
    !!                <li> 10.01.2014 (first implementation)</li>
    !!                <li> 17.01.2016 (added exception handler for set function)</li>
    !!              </ul>
    !!
    !> @param[in]   key     identifies whether coefficient C, its st.dev. dC, the coefficient S, or its st. dev. dS is to be set
    !> @param[in]   n       index n of coefficient (degree)
    !> @param[in]   m       index m of coefficient (order)
    !> @param[in]   val     value which shall be set
    !!
    !> @returns     Error code
    !!
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting the harmonic coefficients or their respective standard deviations
    !!              of the geopotential. As argument it receives the degree (n) and
    !!              order (m) as well as the associated coefficient Cnm, Snm, dCnm or dSnm and sets it accordingly.
    !!
    !> @anchor      setNeptuneVar_geop
    !!
    !!----------------------------------------------------------------------------------
    integer function setNeptuneVar_geop( this,&
                                       key, & ! <-- CHR() identifier string
                                       n,   & ! <-- INT   index n of geopotential series
                                       m,   & ! <-- INT   index m of geopotential series
                                       val  & ! <-- DBL   value to set the variable to
                                     )

        class(Neptune_class)                    :: this
        character(len=*), intent(in)            :: key
        integer,          intent(in)            :: n
        integer,          intent(in)            :: m
        real(dp),         intent(in)            :: val

        character(len=*), parameter             :: csubid = "setNeptuneVar_geop"  ! ID

        setNeptuneVar_geop = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        select case(key)
          case(C_HARMONIC_C)
            call this%gravity_model%setHarmonicCoefficient(COEFF_C, n, m, val)
          case(C_HARMONIC_SD_C)
            call this%gravity_model%setHarmonicCoefficient(COEFF_DC, n, m, val)
          case(C_HARMONIC_S)
            call this%gravity_model%setHarmonicCoefficient(COEFF_S, n, m, val)
          case(C_HARMONIC_SD_S)
            call this%gravity_model%setHarmonicCoefficient(COEFF_DS, n, m, val)
          case default
            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_geop = E_UNKNOWN_PARAMETER
            return
        end select
        if(hasFailed()) then
          setNeptuneVar_geop = getLatestError()
          return
        end if

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if
        return

    end function setNeptuneVar_geop

    !==============================================================================================
    !!
    !> @brief       Set parameters provided by a whole integer array, e.g. solar and geomagnetic activity or distinct harmonics
    !!
    !> @author      Vitali Braun
    !!
    !> @date        <ul>
    !!                <li> 04.02.2015 (first implementation)</li>
    !!              </ul>
    !!
    !> @param[in]   key
    !> @param[in]   val     array containing values which shall be set
    !!
    !> @returns     Error code
    !!
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting input parameters provided in an integer array, for example distinct harmonic coefficients
    !!              to be analysed in the geopotential.
    !!
    !> @anchor      setNeptuneVar_int_arr
    !!
    !!----------------------------------------------------------------------------------
    integer function setNeptuneVar_int_arr(     &
                                        this,   &
                                        key,    & ! <-- CHR() identifier string
                                        val     & ! <-- INT() array of values to set the variables to
                                        )

        class(Neptune_class)                    :: this
        character(len=*),      intent(in)       :: key
        integer, dimension(:), intent(in)       :: val

        character(len=*), parameter  :: csubid = "setNeptuneVar_int_arr"  ! ID

        setNeptuneVar_int_arr = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        select case(key)

          case(C_HARMONICS)
            call this%gravity_model%setDistinctHarmonicsArray(val)

          case default

            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_int_arr = E_UNKNOWN_PARAMETER
            return

        end select

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if

        return

    end function setNeptuneVar_int_arr

    !==============================================================================================
    !!
    !> @brief       Set parameters provided by a logical, right now only the flag whether or not to check for c_d and c_r boundaries
    !!
    !> @author      Daniel Lück
    !!
    !> @date        <ul>
    !!                <li> 27.03.2024 (first implementation)</li>
    !!              </ul>
    !!
    !> @param[in]   key
    !> @param[in]   val     logical containing values which shall be set
    !!
    !> @returns     Error code
    !!
    !> @details     This function is part of the NEPTUNE API and serves for
    !!              setting input parameters provided in a logical, for example check boundaries for cr and cd
    !!
    !> @anchor      setNeptuneVar_logical
    !!
    !!----------------------------------------------------------------------------------
    integer function setNeptuneVar_logical(     &
                                        this,   &
                                        key,    & ! <-- CHR() identifier string
                                        val     & ! <-- INT() array of values to set the variables to
                                        )

        class(Neptune_class)                    :: this
        character(len=*),      intent(in)       :: key
        logical, intent(in)                     :: val

        character(len=*), parameter  :: csubid = "setNeptuneVar_logical"  ! ID

        setNeptuneVar_logical = 0

        if(isControlled()) then
          if(hasToReturn()) return
          call checkIn(csubid)
        end if

        ! check if input array is already set up
        if(.not. allocated(this%input_arr)) then
            call this%initialize_input_array()
        end if

        select case(key)

          case(C_BOUNDARY_CHECK)
            call this%satellite_model%setBoundaryCheck(val)

          case default

            call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/key/))
            setNeptuneVar_logical = E_UNKNOWN_PARAMETER
            return

        end select

        !** done!
        if(isControlled()) then
          call checkOut(csubid)
        end if

        return

    end function setNeptuneVar_logical

!==============================================================================================
!
!> @anchor      initializeInputArray
!!
!> @brief       Initialize the input array with all input parameters available
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li>VB: 16.02.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    subroutine initializeInputArray(this)
        implicit none

        class(Neptune_class)                    :: this

        character(len=*), parameter :: csubid = 'initializeInputArray'
        type(perturbation_t) :: pert
        type(neptune_out_t)  :: outFile

        integer :: i

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        ! FILE/PATH parameters
        call this%set_input(parName=C_PATH_DATA,   valType='dir', initFlag=.true.)
        call this%set_input(parName=C_PATH_INPUT,  valType='dir', initFlag=.true.)
        call this%set_input(parName=C_PATH_OUTPUT, valType='dir', initFlag=.true.)
        call this%set_input(parName=C_FILE_MANEUVERS, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_SURFACES,  valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_EOP,    valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_EOP,    val=this%reduction%getEopFileName(), set=.true., info=this%reduction%getEopFileVersion(this%getDataPath()))
        call this%set_input(parName=C_FILE_SOLMAG, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_SOLMAG, val=this%atmosphere_model%getSolmagFileName(DAILY), set=.true., info=this%atmosphere_model%getSolmagFileVersion(DAILY))
        call this%set_input(parName=C_FILE_SOLMAG_MONTHLY, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_SOLMAG_MONTHLY, val=this%atmosphere_model%getSolmagFileName(MONTHLY), set=.true., info=this%atmosphere_model%getSolmagFileVersion(MONTHLY))
        call this%set_input(parName=C_FILE_DWIND,  valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_HWIND,  valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_QDGRID, valType='file', initFlag=.true.)

        ! dump and progress files also have defaults
        call this%set_input(parName=C_FILE_INPUT_DUMP, valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_INPUT_DUMP, val = this%dump_file_name, set = .true.)
        call this%set_input(parName=C_FILE_PROGRESS,   valType='file', initFlag=.true.)
        call this%set_input(parName=C_FILE_PROGRESS,   val = this%progress_file_name, set = .true.)

        ! COVARIANCE
        call this%set_input(parName=C_INITIAL_COVARIANCE, valType='matrix', initFlag=.true.)

        ! STATE
        call this%set_input(parName=C_INITIAL_STATE, valType='state', initFlag=.true.)

        ! DATETIME parameters
        call this%set_input(parName=C_EPOCH_START_GD, valType='datetime', initFlag=.true.)
        call this%set_input(parName=C_EPOCH_END_GD, valType='datetime', initFlag=.true.)

        ! STRING parameters
        call this%set_input(parName=C_RUN_ID, valType='string', initFlag=.true.)
        call this%set_input(parName=C_OPT_SHADOW, valType='string', initFlag=.true.)

        ! INTEGER parameters
        call this%set_input(parName=C_GEOPOTENTIAL, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_SAT_PROPERTIES, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OUTPUT_STEP, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_COV_METHOD, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_GEO_MODEL, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_AP_FORECAST, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_STORE_DATA, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_COV_GEOPOTENTIAL, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_METHOD, valType='integer', initFlag=.true.)
        call this%set_input(parName=C_OPT_ATMOSPHERE_MODEL, valType='integer', initFlag=.true.)

        ! REAL parameters
        call this%set_input(parName=C_PAR_MASS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_CROSS_SECTION, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_CDRAG, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_CREFL, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_REENTRY, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_RELEPS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_ABSEPS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_INT_COV_STEP, valType='double', initFlag=.true.)
        call this%set_input(parName=C_PAR_EARTH_RADIUS, valType='double', initFlag=.true.)
        call this%set_input(parName=C_OPT_SOL_FORECAST, valType='double', initFlag=.true.)

        ! ON/OFF parameters (also set default values, if available)
        do i = 1, this%derivatives_model%get_neptune_perturbation_number()
            pert = this%derivatives_model%get_neptune_perturbation(this%gravity_model,i)
            if(pert%id == PERT_GRAVITY) cycle    ! skip geopotential (as we set its degree)
            call this%set_input(parName = pert%name, valType = 'boolean', initFlag = .true.)
            call this%set_input(parName = pert%name, val = toString(pert%switch), set = .true.)
        end do

        do i = 1, this%output%get_neptune_output_file_number()
            outFile = this%output%get_neptune_output_file(i)
            call this%set_input(parName = outFile%par_name, valType = 'boolean', initFlag = .true.)
            call this%set_input(parName = outFile%par_name, val = toString(outFile%switch), set = .true.)
        end do

        call this%set_input(parName=C_CORRELATION, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_EOP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_PN_LOOKUP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_INT_LOG, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_SRP_CORRECT, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OUTPUT_FILES, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_PROP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_SUN, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_MOON, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_DRAG, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_COV_SRP, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OUTPUT_VAR_ECI, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_HARMONICS, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_PROGRESS, valType='boolean', initFlag=.true.)
        call this%set_input(parName=C_OPT_PROGRESS, val = toString(this%write_progress), set = .true.)

        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

    end subroutine

!==================================================================
!
!> @brief     Resets the current_index variable
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  26.04.2014 (initial design) </li>
!!              <li>VB:  20.04.2015 (added cov counter) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @anchor   resetStoreCounter
!!
!------------------------------------
  subroutine resetStoreCounter(this)

    class(Neptune_class)    :: this

    this%current_index     = 0
    this%current_index_set = 0
    this%current_index_cov = 0

  end subroutine resetStoreCounter

!==================================================================
!
!> @brief     Returns one data set from the state error transition matrix array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  04.12.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @param[in] idx   Requested index in setMatrix array
!!
!> @details   This routine returns the stored components of the state error
!!            transition matrix from a NEPTUNE propagation run.
!!            If the index is greater than the
!!            max. index for available data, an empty set will be
!!            returned and a warning will be generated.
!!
!> @anchor    getNeptuneSetMatrixData
!!
!------------------------------------
  type(covariance_t) function getNeptuneSetMatrixData(this,idx) result(set)

    class(Neptune_class)    :: this
    integer, intent(in)     :: idx

    character(len=*), parameter :: csubid = 'getNeptuneSetMatrixData'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(idx < 1 .or. idx > this%current_index_set) then
      call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
      set%elem      = 0.d0
      set%epoch%mjd = 0.d0
      return
    else
      set = setMatrix(idx)
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getNeptuneSetMatrixData

!==================================================================
!
!> @brief     Returns one data set from the covariance matrix array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  20.04.2015 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!< @param[in] idx   Requested index in covMatrix array
!!
!< @details   This routine returns stored covariance matrix
!!            information from a NEPTUNE propagation run.
!!            If the index is greater than the
!!            max. index for available data, an empty set will be
!!            returned and a warning will be generated.
!!
!> @anchor    getNeptuneCovarianceData
!!
!------------------------------------
  type(covariance_t) function getNeptuneCovarianceData(this,idx) result(covOut)

    class(Neptune_class)    :: this
    integer, intent(in)     :: idx

    character(len=*), parameter :: csubid = 'getNeptuneCovarianceData'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(idx < 1 .or. idx > this%current_index_cov) then
      call setNeptuneError(E_OUT_OF_BOUNDS, WARNING)
      covOut%elem      = 0.d0
      covOut%epoch%mjd = 0.d0
      return
    else
      covOut = covMatrix(idx)
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end function getNeptuneCovarianceData

!==================================================================
!
!> @brief     Returns one data set from the ephemerides array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  06.02.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @param[in] idx   Requested index in ephem array
!!
!> @details   This routine returns stored ephemerides from a NEPTUNE
!!            propagation run. If the index is greater than the
!!            max. index for available data, an empty set will be
!!            returned and a warning will be generated.
!!
!> @anchor    getNeptuneData
!!
!------------------------------------
  type(state_t) function getNeptuneData(this,idx) result(state)

    class(Neptune_class)    :: this
    integer, intent(in)     :: idx

    character(len=*), parameter :: csubid = 'getNeptuneData'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(idx < 1 .or. idx > this%current_index) then

      call setNeptuneError(E_OUT_OF_BOUNDS, WARNING)
      state = 0.d0
      return

    else

      state = ephem(idx)

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getNeptuneData

!==================================================================
!
!> @brief     Returns the storeEphemData flag
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  06.02.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @anchor    getStoreDataFlag
!!
!------------------------------------
  logical function getStoreDataFlag(this)

    class(Neptune_class)    :: this

    getStoreDataFlag = this%storeEphemData
    return

  end function getStoreDataFlag

!==================================================================
!
!> @brief     Switches the storeEphemData flag
!!
!> @author    Daniel Lubián Arenillas
!!
!> @date      <ul>
!!              <li>DLA:  21.10.2022 (initial design) </li>
!!            </ul>
!!
!> @anchor    switchStoreDataFlag
!!
!------------------------------------
  subroutine switchStoreDataFlag(this)

    class(Neptune_class)    :: this

    this%storeEphemData = .not. this%storeEphemData
    return

  end subroutine switchStoreDataFlag

!==================================================================
!
!> @brief     Returns the number of the ephemerides in the data array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  07.05.2014 (initial design) </li>
!!              <li>VB:  20.04.2015 (added optional flag discriminating between state, covariance and set matrix index output)</li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @param[in] typ   optional string telling which index to return ('COV' -> covariance matrix, 'SET' -> state error transition matrix, <NONE> -> state vector)
!> @anchor    getNumberOfEphemerides
!!
!------------------------------------
  integer function getNumberOfEphemerides(this,typ)

    class(Neptune_class)                    :: this
    character(len=*), intent(in), optional  :: typ

    if(present(typ)) then
      if(typ == 'COV') then
        getNumberOfEphemerides = this%current_index_cov
      else if(typ == 'SET') then
        getNumberOfEphemerides = this%current_index_set
      else
        getNumberOfEphemerides = this%current_index
      end if
    else
      getNumberOfEphemerides = this%current_index
    end if

    return

  end function getNumberOfEphemerides

!==================================================================
!
!> @brief     Returns the size of the ephemerides data array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  06.02.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @anchor    getDataArraySize
!!
!------------------------------------
  integer function getDataArraySize(this)

    class(Neptune_class)    :: this

    getDataArraySize = isize
    return

  end function getDataArraySize

!==================================================================
!
!> @brief     Gets the time step in seconds between subsequent ephemerides
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  06.02.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @anchor    getStep
!!
!------------------------------------
  integer function getStep(this)

    class(Neptune_class)    :: this

    getStep = this%ephem_step
    return

  end function getStep
!==================================================================
!
!> @brief     Sets the time step in seconds between subsequent ephemerides
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  06.02.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @param[in] time_step   Step size in seconds used for ephemerides data storage
!!
!> @anchor    setStep
!!
!------------------------------------
  subroutine setStep(this,time_step)

    class(Neptune_class)    :: this
    integer, intent(in)     :: time_step

    character(len=*), parameter :: csubid = 'setStep'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(time_step < 0) then
      call setNeptuneError(E_INVALID_STEP, FATAL)
      return
    else if(time_step == 0) then
      this%storeEphemData  = .false.
    else
      this%ephem_step      = time_step
      this%storeEphemData  = .true.
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine setStep

!==================================================================
!
!> @brief     Stores data in the covariance matrix array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  20.04.2015 (initial design) </li>
!!              <li>VB:  15.05.2016 (added hasFailed check for date2string conversion) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!              <li>CHK: 04.01.2018 (current_index_set -> current_index_cov) </li>
!!            </ul>
!!
!> @param[in] cov   Covariance matrix
!> @param[in] mjd   MJD of epoch
!!
!> @anchor    storeCovData
!!
!------------------------------------
  subroutine storeCovData(this, cov, mjd)

    class(Neptune_class)                 :: this
    real(dp), dimension(6,6), intent(in) :: cov
    real(dp),                 intent(in) :: mjd

    character(len=*), parameter :: csubid = 'storeCovData'
    character(len=20)           :: cepoch
    character(len=255)          :: cmess

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if
!
!    if (.not. allocated(covMatrix)) then
!        allocate(covMatrix(isize))
!    end if

    this%current_index_cov = this%current_index_cov + 1

    if(this%current_index_cov > isize) then         ! array is full
      if(.not. this%warned) then                    ! warn that array is full
        cepoch = date2string(ephem(isize)%epoch)
        if(hasFailed()) return
        write(cmess,'(a)') 'No more data can be added to cov-array. Last element at epoch '//cepoch
        call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
        this%warned = .true.
      end if
    else
        ! Note there was a possible bug here:
        !  current_index_cset has been replaced by current_index_cov
      covMatrix(this%current_index_cov)%elem      = cov
      covMatrix(this%current_index_cov)%epoch%mjd = mjd
      call mjd2gd(covMatrix(this%current_index_cov)%epoch)
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine storeCovData

!==================================================================
!
!> @brief     Stores data in the setMatrix array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  04.12.2014 (initial design) </li>
!!              <li>VB:  20.04.2015 (changed routine name, to harmonise with new covariance setter)</li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @param[in] set   State error transition matrix
!> @param[in] mjd   MJD of epoch
!!
!> @anchor    storeSetData
!!
!------------------------------------
  subroutine storeSetData(this, set, mjd)

    class(Neptune_class)                 :: this
    real(dp), dimension(6,6), intent(in) :: set
    real(dp),                 intent(in) :: mjd

    character(len=*), parameter :: csubid = 'storeSetData'
    character(len=20)  :: cepoch
    character(len=255) :: cmess

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if
!
!    if (.not. allocated(this%setMatrix)) then
!        allocate(this%setMatrix(isize))
!    end if

    this%current_index_set = this%current_index_set + 1
    
    if(this%current_index_set > isize) then  ! array is full

      if(.not. this%warned) then ! warn that array is full

        cepoch = date2string(ephem(isize)%epoch)
        write(cmess,'(a)') 'No more data can be added to set-array. Last element at epoch '//cepoch

        call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))

        this%warned = .true.

      end if

    else

      setMatrix(this%current_index_set)%elem      = set
      setMatrix(this%current_index_set)%epoch%mjd = mjd

      call mjd2gd(setMatrix(this%current_index_set)%epoch)

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine storeSetData

!==================================================================
!
!> @brief     Stores data in the ephemerides array
!!
!> @author    Vitali Braun
!> @author    Christopher Kebschull
!!
!> @date      <ul>
!!              <li>VB:  06.02.2014 (initial design) </li>
!!              <li>CHK: 04.01.2018 (Moved to neptuneClass) </li>
!!            </ul>
!!
!> @param[in] r     radius vector (km) in GCRF
!> @param[in] v     velocity vector (km/s) in GCRF
!> @param[in] mjd   MJD of epoch
!!
!> @anchor    storeData
!!
!------------------------------------
  subroutine storeData(this, r, v, mjd)

    use slam_strings,       only: toString
        
    class(Neptune_class)               :: this
    real(dp), dimension(3), intent(in) :: r
    real(dp), dimension(3), intent(in) :: v
    real(dp), intent(in)               :: mjd

    character(len=*), parameter        :: csubid = 'storeData'
    character(len=20)                  :: cepoch
    character(len=255)                 :: cmess

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

!    if (.not. allocated(this%ephem)) then
!        allocate(this%ephem(isize))
!    end if

    this%current_index = this%current_index + 1

    if(this%current_index > isize) then  ! array is full

      if(.not. this%warned) then ! warn that array is full

        cepoch = date2string(ephem(isize)%epoch)
        write(cmess,'(a)') 'No more data can be added to eph-array ('//toString(this%current_index)//'). Last element at epoch: '//cepoch

        call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))

        this%warned = .true.

      end if

    else

      ephem(this%current_index)%r = r
      ephem(this%current_index)%v = v

      ephem(this%current_index)%radius_unit   = UNIT_KM
      ephem(this%current_index)%velocity_unit = UNIT_KMPS

      ephem(this%current_index)%epoch%mjd = mjd
      ephem(this%current_index)%epoch%jd  = mjd + jd245

      call mjd2gd(ephem(this%current_index)%epoch)
      
      ! cepoch = date2string(ephem(this%current_index)%epoch)
      ! write(cmess,'(a)') 'Added date '//cepoch//' at index '//toString(this%current_index)
      ! call message(cmess, LOG_AND_STDOUT)

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine storeData


end module NeptuneClass
