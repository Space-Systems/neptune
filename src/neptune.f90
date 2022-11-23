!>-----------------------------------------------------------------------------------------------
!> @brief       NEPTUNE - NPI Ephemeris Propagation Tool with Uncertainty Extrapolation
!!
!> @author      Vitali Braun (VB)
!< @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  18.12.2012 (Geopotential + EOP)</li>
!!                <li>VB:  24.01.2013 (Atmosphere) </li>
!!                <li>VB:  28.02.2013 (Third bodies) </li>
!!                <li>VB:  14.03.2013 (Solar radiation pressure) </li>
!!                <li>VB:  20.04.2013 (Interface/API redesign) </li>
!!                <li>VB:  21.04.2013 (Logfile initialization) </li>
!!                <li>VB:  09.10.2013 (Interface/API redesign) </li>
!!                <li>CHK: 13.11.2013 (updated to use libslam) </li>
!!                <li>VB:  24.01.2016 (Introduced reference frame check in NEPTUNE propagator)</li>
!!                <li>VB:  20.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>VB:  14.06.2016 (removed autoconfiguration, as it is not used)</li>
!!                <li>ChK: 24.12.2017 (using a neptune class now)</li>
!!              </ul>
!!
!> @details     This is the NEPTUNE (<b>N</b>PI <b>E</b>phemeris
!!              <b>P</b>ropagation <b>T</b>ool with <b>Un</b>certainty <b>E</b>xtrapolation) module.
!!              It is initialized via the routine 'setNeptuneVar' and a subsequent call
!!              to 'init_neptune'. Propagation is done via the routine 'neptune'.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      libneptune
!!
!> @todo        create output file for correlation matrix vs. time
!!
!!------------------------------------------------------------------------------------------------
module libneptune

    use slam_astro,             only: initAstroConstants, isAstroInitialized, getEarthRadius, getEarthGeopotentialRadius
    use averaging,              only: setMeanElements
    use slam_astro_conversions, only: rv2coe
    use derivatives,            only: PERT_ALBEDO, PERT_DRAG, PERT_JUPITER, PERT_MANEUVERS, PERT_MARS, &
                                    PERT_MERCURY, PERT_MOON, PERT_NEPTUNE, PERT_OCEAN_TIDES, PERT_SATURN, PERT_SOLID_TIDES, PERT_SRP, PERT_SUN, &
                                    PERT_URANUS, PERT_VENUS, PERT_WIND
    use neptune_error_handling, only: E_INTEGRATION_ABORT, E_MIN_ALTITUDE, E_NEPTUNE_INIT, E_UNSUPPORTED_FRAME, &
                                    E_INVALID_SATELLITE, setNeptuneError
    use slam_error_handling,    only: hasToReturn, isControlled, hasFailed, isSetErrorHandling, E_FILE_NOT_FOUND, WARNING, FATAL, &
                                    E_MISSING_PARAMETER, ERRORS, initErrorHandler, ALL_MSG, getLogfileChannel, &
                                    setLogFileName, getLogFileName, setLogFileChannel, setLogVerbosity, setCliVerbosity, checkIn, checkOut
    use slam_io,                only: SWITCHED_ON, SWITCHED_OFF, SEQUENTIAL, FT_LOG, openFile, cdelimit, LOGFILE, LOG_AND_STDOUT, message, write_progress
    use slam_math,              only: eps3, eps9, eps15, identity_matrix, mag, undefined, rad2deg, cross
    use neptuneClass,           only: Neptune_class
    use neptuneClock,           only: Clock_class
    use neptuneParameters,      only: OPTION_CORRELATION, OPTION_COV_METHOD, OPTION_GEO_MODEL, OPTION_INT_LOGFILE, &
                                    OPTION_OUTPUT, OPTION_PROPS, PAR_CDRAG, PAR_CREFL, PAR_CROSS_SECTION, PAR_INT_ABSEPS, PAR_INT_RELEPS, PAR_INT_COV_STEP, &
                                    PAR_MASS, PROP_CONSTANT, PROP_FILE, C_INITIAL_STATE, C_PAR_EARTH_RADIUS
    use numint,                 only: MAX_RESETS
    use satellite,              only: ORIENT_SPHERE
    use maneuvers,              only: neptune_maneuver_t => maneuver_t
    use slam_strings,           only: toString
    use slam_time,              only: time_t, gd2mjd, date2longstring, operator(>)
    use slam_types,             only: dp
    use slam_rframes,           only: FRAME_NAME_LENGTH => MAX_ID_LENGTH, getFrameName, REF_FRAME_GCRF
    use slam_orbit_types,       only: state_t, kepler_t, covariance_t

    implicit none

    private

    !========================================================================
    !
    ! PUBLIC elements
    !
    !-------------------------------
    public :: init_neptune                                                      ! initialization
    public :: propagate                                                         ! propagation of state vector and covariance
    public :: propagate_set                                                     ! propagation of state vector, covariance and state error transition matrix
    public :: neptune_maneuver_t                                                ! Maneuver type

contains


  !========================================================================
  !!
  !>  @anchor     init_neptune
  !!
  !>  @brief      Performs the NEPTUNE initialization
  !>  @author     Vitali Braun
  !!
  !>  @param[in]  neptune         The NEPTUNE class instance
  !>  @param[in]  initial_state   (optional) Initial state vector (state_t)
  !!
  !>  @date       <ul>
  !!                <li>13.02.2014 (added doxygen header)</li>
  !!                <li>18.02.2014 (added epoch as input parameter)</li>
  !!                <li>23.04.2014 (changed 'initial_state' and 'epoch' to optional parameters)</li>
  !!                <li>14.06.2016 (removed autoconfiguration)</li>
  !!                <li>14.03.2017 (removed epoch input parameter, as no longer used)</li>
  !!              </ul>
  !!
  !----------------------------------------------------
  integer function init_neptune(neptune,initial_state_in)

    type(Neptune_class)                 :: neptune
    type(state_t),optional,intent(in)   :: initial_state_in

    type(state_t)                       :: initial_state
    character(len=255)                  :: cmess                                ! message string
    character(len=*), parameter :: csubid = "init_neptune"                      ! subroutine id

    integer  :: ierr                                                            ! error flag

    type(time_t)            :: startEpoch, endEpoch, tempStartEpoch, tempEndEpoch

    init_neptune = 0

    !** check whether error handling has been initialized before - if not, NEPTUNE will apply its own scheme...
    if(.not. isSetErrorHandling()) then
      call initErrorHandler(control = 'YES', errAction = 'RETURN', traceback = 'YES')
    end if

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !-------------------------------------------------------
    !
    !   Initialization procedure
    !
    !--------------------------------
    !
    !   1) Initialize logfile
    !   2) Initialize astrodynamics constants
    !   3) Satellite configuration
    !   4) Initialize perturbation methods and integration parameters
    !   5) Initialize correlation routines
    !   6) Prepare output files
    !
    !--------------------------------------------------------


    !==============================================================
    !
    ! 1) Initialize logfile
    !
    !--------------------------------------------------------
    if(getLogfileChannel() == 0) then !** open logfile
      ierr = setLogfileName("neptune.log")
      ierr = setLogfileChannel(openFile(getLogfileName(), SEQUENTIAL, FT_LOG))
      ierr = setLogVerbosity(ALL_MSG)
      ierr = setCliVerbosity(ERRORS)
    end if

    !==============================================================
    !
    ! 2) Initialize astrodynamical constants
    !
    !----------------------------------------------
    if(.not. isAstroInitialized()) then
      call initAstroConstants(neptune%iopt(OPTION_GEO_MODEL))  ! using selected model
    end if
    !----------------------------------------------

    !==============================================================
    !
    ! 3) Satellite configuration
    !
    !----------------------------------------------
    if(neptune%iopt(OPTION_PROPS) == PROP_CONSTANT) then                        ! cannon ball model

      ! Reseting the satellite at this point as the init_neptune procedure could be called multiple times.
      ! Making sure the new satellite has an initialized surfaces array.
      call neptune%satellite_model%resetObject()

      !** Mass
      if(neptune%isSetMass) then
        neptune%satellite_model%mass  = neptune%ipar(PAR_MASS)
      else
        cmess = "Satellite mass not defined."
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return
      end if

      !** Drag coefficient
      if(neptune%isSetCdrag) then
        neptune%satellite_model%cdrag = neptune%ipar(PAR_CDRAG)
      else
        cmess = "Drag coefficient not defined."
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return
      end if

      !** SRP coefficient
      if(neptune%isSetCsrp) then

        neptune%satellite_model%surface(1)%reflSpec = neptune%ipar(PAR_CREFL) - 1.d0             ! the specular reflectivity coefficient is saved,
                                                                                ! while the SRP coefficient is used in a cannon-ball model
        neptune%satellite_model%surface(1)%reflDiff = 0.d0                                       ! no diffuse reflectivity in this case

      else
        cmess = "SRP coefficient not defined."
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return
      end if

      !** Area
      if(neptune%isSetArea) then
        neptune%satellite_model%surface(1)%area = neptune%ipar(PAR_CROSS_SECTION)                ! in km**2
      else
        cmess = "Cross-section area not defined."
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return
      end if

      !** normal angle (not required for sphere, but has to be initialised...)
      neptune%satellite_model%surface(1)%normal_angle(:) = 0.d0

      neptune%satellite_model%surface(1)%orientation = ORIENT_SPHERE
      neptune%satellite_model%surface(1)%id = 1
      neptune%satellite_model%nsurfaces = 1

      call neptune%satellite_model%initObject()

      if(hasFailed()) return

    else if(neptune%iopt(OPTION_PROPS) == PROP_FILE) then  ! read surfaces definition file

      call neptune%satellite_model%initObject()
      call neptune%satellite_model%readSurfaceDefinitionFile(neptune%getInputPath())

    end if

    !** initial orbit
    if(present(initial_state_in)) then
      initial_state = initial_state_in
      ierr = neptune%setNeptuneVar("INITIAL_STATE", initial_state)
    end if

    !** check for start and end epoch being available
    call neptune%getStartEpoch(tempStartEpoch, ierr)
    if(ierr /= 0) then  ! epoch missing
      call setNeptuneError(E_MISSING_PARAMETER, FATAL, (/'startEpoch (provide via setNeptuneVar)'/))
      return
    end if

    call neptune%getEndEpoch(tempEndEpoch, ierr)
    if(ierr /= 0) then  ! epoch missing
      call setNeptuneError(E_MISSING_PARAMETER, FATAL, (/'endEpoch - (provide via setNeptuneVar)'/))
      return
    end if

    !** For the following routines, it is required to have startEpoch < endEpoch.
    !   Therefore, check if there is backward propagation...
    if(tempEndEpoch > tempStartEpoch) then
        endEpoch   = tempEndEpoch
        startEpoch = tempStartEpoch
    else
        endEpoch   = tempStartEpoch
        startEpoch = tempEndEpoch
    end if

    !============================================================
    !
    ! 4) Initialize perturbation methods
    !----------------------------------------------

    ! 4-1) Geopotential

    call neptune%gravity_model%initGravityPotential(                &
                               neptune%getDataPath(),               & ! <-- CHR path of coefficients data file
                               neptune%iopt(OPTION_GEO_MODEL)       & ! <-- INT model to be used
                             )

    ! 4-2) Earth orientation parameters
    if(neptune%reduction%getEopFlag()) then
      call neptune%reduction%initEOP(       &
                    neptune%getDataPath())
    end if

    ! 4-3) Earth's atmosphere

    if(neptune%derivatives_model%getPertSwitch(PERT_DRAG)) then
      if(neptune%derivatives_model%getPertSwitch(PERT_WIND)) then  ! wind only switched on, if drag perturbations are considered
        call neptune%atmosphere_model%setHorizontalWind(.true.)
      else
        call neptune%atmosphere_model%setHorizontalWind(.false.)
      end if

      call neptune%atmosphere_model%initAtmosphere( &
                          neptune%getDataPath())
      if(hasFailed()) return

    end if

    ! 4-4) Perturbations requiring Sun's or Moon's position

    if(neptune%derivatives_model%getPertSwitch(PERT_MOON)        .or. &
       neptune%derivatives_model%getPertSwitch(PERT_SUN)         .or. &
       neptune%derivatives_model%getPertSwitch(PERT_MERCURY)     .or. &
       neptune%derivatives_model%getPertSwitch(PERT_VENUS)       .or. &
       neptune%derivatives_model%getPertSwitch(PERT_JUPITER)     .or. &
       neptune%derivatives_model%getPertSwitch(PERT_SATURN)      .or. &
       neptune%derivatives_model%getPertSwitch(PERT_URANUS)      .or. &
       neptune%derivatives_model%getPertSwitch(PERT_NEPTUNE)     .or. &
       neptune%derivatives_model%getPertSwitch(PERT_MARS)        .or. &
       neptune%derivatives_model%getPertSwitch(PERT_ALBEDO)      .or. &
       neptune%derivatives_model%getPertSwitch(PERT_SRP)         .or. &
       neptune%derivatives_model%getPertSwitch(PERT_SOLID_TIDES) .or. &
       neptune%derivatives_model%getPertSwitch(PERT_OCEAN_TIDES) .or. &
       neptune%iopt(OPTION_PROPS) == PROP_FILE ) then     ! also if surface definition is used, as Sun-oriented surfaces may be used

      call neptune%solarsystem_model%initSolarSystem(   &
                                neptune%getDataPath(),  &
                                'DE-421',               &
                                endEpoch                &
                                )
      if(hasFailed()) return
    end if

    ! 4-5) Albedo
    if(neptune%derivatives_model%getPertSwitch(PERT_ALBEDO)) then
      call neptune%radiation_model%initAlbedo()
      if(hasFailed()) return
    end if


    ! 4-6) Tides
    if(neptune%derivatives_model%getPertSwitch(PERT_SOLID_TIDES) .or. neptune%derivatives_model%getPertSwitch(PERT_OCEAN_TIDES)) then
      call neptune%tides_model%initTides()
      if(hasFailed()) return
    end if

    ! 4-7) Maneuvers
    ! This might be a bit un-nice: IF maneuvers have been initialized using the API, this call is performed but actually "ignored"
    if(neptune%derivatives_model%getPertSwitch(PERT_MANEUVERS)) then
      call neptune%manoeuvres_model%init_maneuvers(neptune%getInputPath(), startEpoch, endEpoch)
      if(hasFailed()) return
    end if

    ! 4-8) numerical integration parameters
    if(neptune%flag_init_tolerances(1)) then    !** relative tolerance
      call neptune%numerical_integrator%setIntegrationTolerance( trel = neptune%ipar(PAR_INT_RELEPS))
      if(hasFailed()) return
    end if

    if(neptune%flag_init_tolerances(2)) then !** absolute tolerance
      call neptune%numerical_integrator%setIntegrationTolerance( tabs = neptune%ipar(PAR_INT_ABSEPS))
      if(hasFailed()) return
    end if

    !===================================================================
    !
    ! 5) Covariance/correlation matrix initialization
    !
    !---------------------------------------------------------------
    if(neptune%iopt(OPTION_CORRELATION) == SWITCHED_ON) then
      call neptune%correlation_model%initCorrelation(neptune%gravity_model)
      if(hasFailed()) return
    end if

    ! 5.1) set covariance matrix integration step size
    call neptune%numerical_integrator%setCovarianceIntegrationStep(neptune%ipar(PAR_INT_COV_STEP))
    if(hasFailed()) return

    ! 5.2) set covariance matrix integration method
    call neptune%numerical_integrator%setCovarianceIntegrationMethod(neptune%iopt(OPTION_COV_METHOD))
    if(hasFailed()) return

    ! 5.3) write general header (ALSO FOR OTHER OUTPUT FILES)
    call neptune%output%write_header(neptune%gravity_model,             &
                                        neptune%atmosphere_model,       &
                                        neptune%manoeuvres_model,       &
                                        neptune%radiation_model,        &
                                        neptune%satellite_model,        &
                                        neptune%thirdbody_model,        &
                                        neptune%numerical_integrator,   &
                                        neptune%derivatives_model,      &
                                        neptune%version_model,          &
                                        neptune%reduction,              &
                                        neptune%correlation_model,      &
                                        0,                              &
                                        "")
    if(hasFailed()) return

    ! 5.4) Open integration logfile (if required)
    if(neptune%iopt(OPTION_INT_LOGFILE) == SWITCHED_ON) then

      call neptune%numerical_integrator%setIntegrationLogfile('ON')
      if(hasFailed()) return

      !** write specific information
      call neptune%output%write_header(neptune%gravity_model,                           &
                        neptune%atmosphere_model,                                       &
                        neptune%manoeuvres_model,                                       &
                        neptune%radiation_model,                                        &
                        neptune%satellite_model,                                        &
                        neptune%thirdbody_model,                                        &
                        neptune%numerical_integrator,                                   &
                        neptune%derivatives_model,                                      &
                        neptune%version_model,                                          &
                        neptune%reduction,                                              &
                        neptune%correlation_model,                                      &
                        neptune%numerical_integrator%getIntegrationLogfileChannel(),    &
                        neptune%numerical_integrator%getIntegrationLogfileName())
      if(hasFailed()) return

    end if

    neptune%neptuneInitialized = .true.

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

  end function init_neptune

  !========================================================================
  !
  !>  @anchor     propagate
  !!
  !>  @brief      NEPTUNE propagation routine
  !<  @author     Vitali Braun
  !!
  !>  @date       <ul>
  !!                <li> 04.07.2013 (documentation started)</li>
  !!                <li> 30.07.2013 (added correlation matrix handling)</li>
  !!                <li> 06.02.2014 (added getData functionality)</li>
  !!                <li> 25.11.2014 (changed 'reset' handling, so that NEPTUNE really stops, if a step can not be passed)</li>
  !!                <li> 24.01.2016 (added reference frame check for input)</li>
  !!                <li> 14.06.2016 (backward propagation now working - but testing still needed)</li>
  !!              </ul>
  !!
  !! @details     This routine is called for the propagation of the state
  !!              vector and the covariance matrix, after NEPTUNE has been
  !!              initialized by calling init_neptune.
  !!
  !> @param[in]   neptune     the neptune class
  !> @param[in]   state_in    input state vector (GCRF expected, but checked)
  !> @param[in]   covar_in    input covariance matrix (GCRF expected, but checked)
  !> @param[in]   epoch       propagation start(1) and end(2) epoch
  !> @param[out]  state_out   output state vector (GCRF)
  !> @param[out]  covar_out   covariance output matrix (GCRF)
  !> @param[in]   flag_reset  resetting the numerical integrator manually
  !!
  !> @todo        Get averaging right
  !!
  !!----------------------------------------------------
  subroutine propagate(          &
                      neptune,   &
                      state_in,  &   ! <-- TYP    input state vector (GCRF)
                      covar_in,  &   ! <-- TYP    input covariance matrix (GCRF)
                      epoch,     &   ! <-- TYP    propagation start and end epoch [epoch(1) : start epoch, epoch(2) : end epoch]
                      state_out, &   ! --> TYP    output state vector
                      covar_out, &   ! --> TYP    output covariance matrix
                      flag_reset &   ! <-- LOG    reset flag
                    )
      implicit none

      !** interface
      !------------------------------------------------------
      type(Neptune_class)        ,intent(inout)  :: neptune
      type(state_t),              intent(in)      :: state_in
      type(covariance_t),         intent(in)      :: covar_in
      type(time_t), dimension(2), intent(in)      :: epoch
      logical,                    intent(in)      :: flag_reset

      type(state_t),              intent(out)     :: state_out
      type(covariance_t),         intent(out)     :: covar_out
      !------------------------------------------------------
      type(covariance_t)                          :: set_in
      type(covariance_t)                          :: set_out

      call identity_matrix(set_in%elem)

      call propagate_set(neptune,  &
                        state_in,  &   ! <-- TYP    input state vector (GCRF)
                        covar_in,  &   ! <-- TYP    input covariance matrix (GCRF)
                        set_in,    &   ! <-- TYP    input set matrix
                        epoch,     &   ! <-- TYP    propagation start and end epoch [epoch(1) : start epoch, epoch(2) : end epoch]
                        state_out, &   ! --> TYP    output state vector
                        covar_out, &   ! --> TYP    output covariance matrix
                        set_out,   &   ! --> TYP    output set matrix
                        flag_reset &   ! <-- LOG    reset flag)
                      )

  end subroutine propagate

  !========================================================================
  !
  !>  @anchor     propagate
  !!
  !>  @brief      NEPTUNE propagation routine
  !<  @author     Vitali Braun
  !!
  !>  @date       <ul>
  !!                <li> 04.07.2013 (documentation started)</li>
  !!                <li> 30.07.2013 (added correlation matrix handling)</li>
  !!                <li> 06.02.2014 (added getData functionality)</li>
  !!                <li> 25.11.2014 (changed 'reset' handling, so that NEPTUNE really stops, if a step can not be passed)</li>
  !!                <li> 24.01.2016 (added reference frame check for input)</li>
  !!                <li> 14.06.2016 (backward propagation now working - but testing still needed)</li>
  !!              </ul>
  !!
  !! @details     This routine is called for the propagation of the state
  !!              vector and the covariance matrix, after NEPTUNE has been
  !!              initialized by calling init_neptune.
  !!
  !> @param[in]   neptune     the neptune class
  !> @param[in]   state_in    input state vector (GCRF expected, but checked)
  !> @param[in]   covar_in    input covariance matrix (GCRF expected, but checked)
  !> @param[in]   epochs      propagation start(1) and end(last) epoch
  !> @param[out]  state_out   output state vector (GCRF)
  !> @param[out]  covar_out   covariance output matrix (GCRF)
  !> @param[in]   flag_reset  resetting the numerical integrator manually
  !!
  !> @todo        Get averaging right
  !!
  !!----------------------------------------------------
  subroutine propagate_set(      &
                      neptune,   &
                      state_in,  &   ! <-- TYP    input state vector (GCRF)
                      covar_in,  &   ! <-- TYP    input covariance matrix (GCRF)
                      set_in,    &   ! <-- TYP    input set matrix
                      epochs,    &   ! <-- TYP    propagation start and end epoch [epoch(1) : start epoch, epoch(last) : end epoch]
                      state_out, &   ! --> TYP    output state vector
                      covar_out, &   ! --> TYP    output covariance matrix
                      set_out,   &   ! --> TYO    output set matrix
                      flag_reset &   ! <-- LOG    reset flag
                    )

    implicit none

    !** interface
    !------------------------------------------------------
    type(Neptune_class),        intent(inout)   :: neptune
    type(state_t),              intent(in)      :: state_in
    type(covariance_t),         intent(in)      :: covar_in
    type(covariance_t),         intent(in)      :: set_in
    type(time_t), dimension(:), intent(in)      :: epochs
    logical,                    intent(in)      :: flag_reset

    type(state_t),              intent(out)     :: state_out
    type(covariance_t),         intent(out)     :: covar_out
    type(covariance_t),         intent(out)     :: set_out
    !------------------------------------------------------

    character(len=*), parameter         :: csubid = "neptune"
    character(len=2)                    :: cresets                              ! number of resets in integration routine
    character(len=FRAME_NAME_LENGTH)    :: frameName                            ! reference frame name
    integer                 :: ierr                                             ! error flag
    integer                 :: nresets                                          ! number of subsequent resets for one integration step
    integer                 :: reset                                            ! initialisation flag
    integer                 :: otype                                            ! characterising orbit type

    logical                 :: flag_backward                                    ! backward propagation flag
    logical                 :: flag_exit                                        ! exit flag

    real(dp),dimension(6,6) :: corrMat                                          ! correlation matrix
    real(dp),dimension(6,6) :: cumSet                                           ! cumulated state transition matrix
    real(dp)                :: dtmp                                             ! auxiliary
    real(dp)                :: delta_time                                       ! step size actually used for each propagation step
    real(dp)                :: end_epoch_sec                                    ! end epoch in seconds (MJD)
    real(dp),dimension(:),allocatable :: step_epochs_sec                        ! intermediate epochs in seconds (MJD)
    real(dp)                :: lastPropCounter                                  ! saving the last prop_counter for which output has taken place
    real(dp)                :: lastPropCounterSuccess                           ! saving the last prop_counter for which a step was successful
    real(dp)                :: prop_counter                                     ! propagation time counter in seconds
    real(dp)                :: request_time                                     ! requested time in numerical integration loop
    real(dp)                :: propCounterAtReset                               ! propCounter at last reset - required to prevent infinite loops
    real(dp),dimension(6,6) :: set                                              ! state error transition matrix
    real(dp)                :: start_epoch_sec                                  ! start epoch in seconds (MJD)
    type(kepler_t)          :: kep                                              ! mean kepler elements for correlation matrix computation
    type(state_t)           :: last_state_out                                   ! saving the last state vector which has been written to output
    type(Clock_class)       :: nep_clock                                        ! scheduling/clock unit to control the stepsize-related actions
    integer                 :: step_counter
    integer                 :: i_epoch
    logical                 :: suppressed_output                                ! Indicates that output and storage should not be triggered due to an intermediate call to the numerical integrator
    logical                 :: intermediate_integrator_call                     ! Indicates that an intermediate call to the integrator is performed
    real(dp)                :: manoeuvre_change_counter                           ! Is the start or end epoch of a manoeuvre
    real(dp)                :: restored_request_time                            ! Requested time before the intermediate step got necessary
    logical                 :: force_no_interpolation
    real(dp)                :: upcoming_maneuver_epoch_mjd
    character(len=255)      :: cmess
    real(dp)                :: diff
    integer                 :: i_next_calls         ! counting runs through getNext while-loop 
    integer                 :: max_next_calls = 100

    restored_request_time = 0.d0
    prop_counter = 0.0d0
    propCounterAtReset = 0.0d0
    lastPropCounterSuccess = 0.0d0
    diff = 0.d0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check if initialization has been performed
    if(.not. neptune%neptuneInitialized) then
      call setNeptuneError(E_NEPTUNE_INIT, FATAL)
      return
    end if

    !============================================================
    !
    !   Input check
    !
    !--------------------------------------------------

    ! check state vector reference frame (should be GCRF, if not: add conversions here...)
    if(state_in%frame /= REF_FRAME_GCRF) then
      ! if necessary, introduce conversions here for other frames - at the moment it only works for GCRF, therefore an exception is thrown
      ! try first, if frame is a known one:
      frameName = getFrameName(state_in%frame)
      if(hasFailed()) return ! unknown frame
        call setNeptuneError(E_UNSUPPORTED_FRAME, FATAL, (/getFrameName(state_in%frame)/))
        return
      end if

    ! check covariance matrix reference frame (should be GCRF, if not: add conversions here...)
    if(covar_in%frame /= REF_FRAME_GCRF) then
      ! if necessary, introduce conversions here for other frames - at the moment it only works for GCRF, therefore an exception is thrown
      ! try first, if frame is a known one
      frameName = getFrameName(covar_in%frame)
      if(hasFailed()) return ! unknown frame
        call setNeptuneError(E_UNSUPPORTED_FRAME, FATAL, (/getFrameName(covar_in%frame)/))
        return
      end if

    !** check for state vector being an orbit actually
    if(mag(state_in%r) < getEarthRadius()) then
      call setNeptuneError(E_MIN_ALTITUDE, FATAL)
      return
    end if

    ! now set eventually missing input parameters
    ierr = neptune%setNeptuneVar(C_INITIAL_STATE, state_in)
    ierr = neptune%setNeptuneVar(C_PAR_EARTH_RADIUS, toString(getEarthGeopotentialRadius()))

    !** dump all inputs
    !call neptune%dump_input("neptune_dump_file.out")

    ! this only happens, if called with dump flag (checked in the routine)
    !call neptune%write_input_to_dump()

    ! prepare output files if requested (checked in the routine)
    call neptune%output%prepare_output(neptune%gravity_model,   &
                            neptune%atmosphere_model,           &
                            neptune%manoeuvres_model,           &
                            neptune%radiation_model,            &
                            neptune%satellite_model,            &
                            neptune%solarsystem_model,          &
                            neptune%thirdbody_model,            &
                            neptune%numerical_integrator,       &
                            neptune%derivatives_model,          &
                            neptune%version_model,              &
                            neptune%reduction,                  &
                            neptune%correlation_model)

    ! init progress writing (if requested)
    if(neptune%has_to_write_progress()) call write_progress(neptune%get_progress_file_name(), -1.d0, 0.d0)  ! step parameter doesn't matter here

    nresets = 0

    ! set start and end epoch - we work with an offset here, where
    ! the start epoch is at 0 seconds and the end epoch reflects the propagation span
    call neptune%setStartEpoch(epochs(1))
    start_epoch_sec = 0.d0 ! epoch(1)%mjd*86400.d0 ! in seconds

    call neptune%setEndEpoch(epochs(size(epochs)))
    end_epoch_sec   = (epochs(size(epochs))%mjd-epochs(1)%mjd)*86400.d0 ! in seconds

    !============================================================
    !
    !   Initialise propagation time counter
    !
    !------------------------------------------------------------
    
    if (size(epochs) > 2) then
      if (allocated(step_epochs_sec)) deallocate(step_epochs_sec)
      allocate(step_epochs_sec(size(epochs)-1))
      do i_epoch = 2, size(epochs)

        ! Handle intermediate steps
        step_epochs_sec(i_epoch-1) = (epochs(i_epoch)%mjd-epochs(i_epoch-1)%mjd)*86400.d0 ! in seconds
      
        ! write(cmess, '(a)') 'diff = '//date2longstring(epochs(i_epoch))//', step_epochs_sec = '//toString(step_epochs_sec(i_epoch-1))
        ! call message(cmess, LOG_AND_STDOUT)
        ! call message('Initial request: '//toString(step_epochs_sec(i_epoch-1)), LOG_AND_STDOUT)
      end do
      !write (*,*) step_epochs_sec
      nep_clock = Clock_class(start_epoch_sec, end_epoch_sec, step_epochs_sec)
    else
      nep_clock = Clock_class(start_epoch_sec, end_epoch_sec)
    end if
    ! Initialize the counters
    call nep_clock%init_counter(neptune)

    prop_counter = start_epoch_sec

    !=====================================================
    !
    ! Consider backward propagation
    !
    !---------------------------------------------
    if(epochs(2)%mjd < epochs(1)%mjd) then
      flag_backward = .true.
      call message(" - Set to backward propagation", LOG_AND_STDOUT)
    else
      flag_backward = .false.
      call message(" - Set to forward propagation", LOG_AND_STDOUT)
    end if

    !============================================================
    !
    !   Initialise covariance matrix propagation
    !
    !--------------------------------------------------
    if(neptune%numerical_integrator%getCovariancePropagationFlag()) then

      cumSet = set_in%elem
      !call identity_matrix(cumSet)  ! initial state error transition matrix is the unity matrix
      call neptune%numerical_integrator%resetCountSetMatrix()    ! the counter for the number of calls to the getStateTransitionMatrix routine is being reset

      !** correlation matrix
      if(neptune%correlation_model%getNoisePropagationFlag()) then
        !** provide mean elements for correlation computations (HAS TO BE ADAPTED LATER ON!!)
        !--------------------------------------------------------------------------------------
        call rv2coe(state_in, kep, otype) ! <-- ADAPT here
        call setMeanElements(kep)
        !-------------------------------------------------------
        call message(' - Computing correlation matrix...', LOG_AND_STDOUT)
        call neptune%correlation_model%updateCorrelationMatrix(     &
                                        neptune%gravity_model,      &
                                        state_in%r,                 &
                                        abs(end_epoch_sec-start_epoch_sec))
      end if

    end if

    !------------------------------------------------------------

    !** set integrator initialisation flag
    if(flag_reset) then
      reset = 1

      !** in case of reset - also reset stored array index in neptuneGet module!
      if(neptune%getStoreDataFlag()) then
        call neptune%resetStoreCounter()
      end if

      ! reset also the count of subroutine calls to the integrator
      call neptune%numerical_integrator%resetCountIntegrator()

    else
      reset = 0
    end if

    !** initialise output state
    state_out = state_in
    covar_out = covar_in
    set_out = set_in

    if(neptune%has_to_write_progress()) call write_progress(neptune%get_progress_file_name(), 0.d0, neptune%get_progress_step())
    step_counter = 0
    suppressed_output = .false.
    intermediate_integrator_call = .false.
    restored_request_time = request_time

    do
      step_counter = step_counter + 1
      ! Check exit clause
      flag_exit = nep_clock%has_finished(prop_counter)

      !===============================================
      !
      ! Output if requested
      !
      !---------------------------------------------
      ! The output is suppressed when an intermediate step is needed
      ! (e.g. for start or end of a manoeuvre, which needs an exact start and end epoch and a reset of the integrator)
      if((nep_clock%has_to_write_output() .and. .not. suppressed_output) .or. flag_exit) then
        ! make sure, that the final result is also written to output
        if(neptune%output%get_output_switch()) then
          call neptune%output%output(neptune%gravity_model, &
                            neptune%atmosphere_model,       &
                            neptune%manoeuvres_model,       &
                            neptune%radiation_model,        &
                            neptune%satellite_model,        &
                            neptune%solarsystem_model,      &
                            neptune%thirdbody_model,        &
                            neptune%tides_model,            &
                            neptune%numerical_integrator,   &
                            neptune%derivatives_model,      &
                            neptune%reduction,              &
                            prop_counter,                   &
                            state_out,                      &
                            covar_out)
          if(hasFailed()) then    ! error occurred during output - abort propagation
            call neptune%output%close_open_files(neptune%numerical_integrator)
            return
          end if
          last_state_out  = state_out
          lastPropCounter = prop_counter
          call nep_clock%toggle_output_flag()
        end if

        if(neptune%getStoreDataFlag()) then
          dtmp = epochs(1)%mjd + prop_counter/86400.d0
          ! call message("output epoch: "//toString(prop_counter), LOG_AND_STDOUT)
          call neptune%storeData(state_out%r, state_out%v, dtmp)
          ! store covariance matrix data if requested
          if(neptune%numerical_integrator%getCovariancePropagationFlag()) then
            call neptune%storeSetData(cumSet, dtmp)
            call neptune%storeCovData(covar_out%elem, dtmp)
          end if
          last_state_out   = state_out
          lastPropCounter  = prop_counter
        end if
        call nep_clock%toggle_output_flag()
      end if
      !---------------------------------------------

      if(flag_exit) exit

      ! Check if the previous call was an intermediate integrator call (due to a manoeuvre start or end)
      if (intermediate_integrator_call) then
        ! ... if so we will restore the original request time
        request_time = restored_request_time
      else
        ! No previous intermediate call - we can just proceed as planned
        request_time = nep_clock%get_next_step(neptune,prop_counter)
        
        ! Save the current requested time in case we determine that an intermediate call must be done due to an immenant manoeuvre start or end
        restored_request_time = request_time
        ! Create a point to restore to when the propagation fails
        last_state_out   = state_out
        lastPropCounter  = prop_counter
      end if

      ! Check whether there is a manoeuvre change (start or end)
      force_no_interpolation = .false.
      neptune%manoeuvres_model%ignore_maneuver_start = .false.
      suppressed_output = .false.
      intermediate_integrator_call = .false.
      if (neptune%derivatives_model%getPertSwitch(PERT_MANEUVERS)) then

        !call message(' Maneuver change at '//toString(epochs(1)%mjd + manoeuvre_change_counter/86400.d0)//' ('//toString(manoeuvre_change_counter)//')', LOG_AND_STDOUT)

        ! Reset when we are starting into the new maneuver or no-maneuver interval
        if (abs(prop_counter - manoeuvre_change_counter) < epsilon(1.d0)) then
          reset = 1
          !force_no_interpolation = .true.
          ! reset also the count of subroutine calls to the integrator
          call neptune%numerical_integrator%resetCountIntegrator()
          call message(' - Performing reset of integrator now '//toString(epochs(1)%mjd + prop_counter/86400.d0)//'('//toString(prop_counter)//') for upcoming maneuver '//toString(epochs(1)%mjd + manoeuvre_change_counter/86400.d0)//' ('//toString(manoeuvre_change_counter)//')', LOGFILE)
        end if

        ! Get the next manoeuvre change epoch
        upcoming_maneuver_epoch_mjd = neptune%manoeuvres_model%get_upcoming_manoeuvre_change_epoch(epochs(1)%mjd + prop_counter/86400.d0, using_backwards_propagation=flag_backward)
        if (upcoming_maneuver_epoch_mjd > 0.d0) then
          manoeuvre_change_counter = (upcoming_maneuver_epoch_mjd - epochs(1)%mjd) * 86400.d0
          !call message(' New maneuver change at: '//toString(epochs(1)%mjd + manoeuvre_change_counter/86400.d0)//' ('//toString(manoeuvre_change_counter)//')', LOG_AND_STDOUT)
          ! Check whether the manoeuvre change is more immenent than the step proposed
          if (.not. flag_backward .and. (manoeuvre_change_counter < request_time .and. manoeuvre_change_counter > prop_counter) &
            .or. flag_backward .and. (manoeuvre_change_counter > request_time .and. manoeuvre_change_counter < prop_counter)) then
            if (abs(manoeuvre_change_counter - request_time) > epsilon(1.d0)) then
              ! The manoeuvre_change_counter does not coincide with a storage and output epoch - surpress the output!
              suppressed_output = .true.
              call message(' - Suppressing intermediate integrator output ', LOGFILE)
            end if
            ! Override the request time with the up coming manoeuvre start or end epoch
            request_time = manoeuvre_change_counter
            intermediate_integrator_call = .true.
            force_no_interpolation = .true.
            neptune%manoeuvres_model%ignore_maneuver_start = .true.
            call message(' - Adding intermediate integrator call now '//toString(epochs(1)%mjd + prop_counter/86400.d0)//'('//toString(prop_counter)//') for planned maneuver at '//toString(epochs(1)%mjd + request_time/86400.d0)//' ('//toString(request_time)//')', LOGFILE)
          end if
        end if
      end if

      !call message(' Going into integrator loop at '//toString(epochs(1)%mjd + prop_counter/86400.d0)//'  ('//toString(prop_counter)//')', LOG_AND_STDOUT)
      !call message(' prop counter: '//toString(prop_counter)//' request_time: '//toString(request_time)//' difference: '//toString(prop_counter - request_time), LOG_AND_STDOUT)


      !====================================================================================
      !
      ! Integration of position, velocity and state error transition matrix
      !
      !--------------------------------------------------------------------------------
      do

        if (neptune%derivatives_model%getPertSwitch(PERT_MANEUVERS)) then
          if (upcoming_maneuver_epoch_mjd > 0.d0) then
            ! Reset when the timestep of the integrator is greater than
            if (.not. flag_backward .and. (manoeuvre_change_counter < (prop_counter + neptune%numerical_integrator%get_current_step_size()) .and. manoeuvre_change_counter > prop_counter) &
            .or. flag_backward .and. (manoeuvre_change_counter > (prop_counter - neptune%numerical_integrator%get_current_step_size()) .and. manoeuvre_change_counter < prop_counter)) then
              !reset = 1
              force_no_interpolation = .true.
              ! reset also the count of subroutine calls to the integrator
              !call neptune%numerical_integrator%resetCountIntegrator()
              call message(' - Forcing no interpolation now '//toString(epochs(1)%mjd + prop_counter/86400.d0)//'('//toString(prop_counter)//') for upcoming maneuver due to time step size '//toString(epochs(1)%mjd + manoeuvre_change_counter/86400.d0)//' ('//toString(manoeuvre_change_counter)//') int step size: '//toString(neptune%numerical_integrator%get_current_step_size()), LOGFILE)
            end if
          end if
        end if

        ! write(cmess,'(a)') ' pre: request_time '//toString(request_time)//', prop_counter '//toString(prop_counter)//', delta_time '//toString(delta_time)//', reset: '//toString(reset)
        ! call message(cmess, LOG_AND_STDOUT)

        if (abs(request_time-prop_counter) > epsilon(1.0d0)) then
          call neptune%numerical_integrator%integrateStep(         &
                              neptune%gravity_model,               &              ! <->  TYPE  Gravity model
                              neptune%atmosphere_model,            &              ! <->  TYPE  Atmosphere model
                              neptune%manoeuvres_model,            &              ! <->  TYPE  Manoeuvres model
                              neptune%radiation_model,             &              ! <->  TYPE  Radiation model
                              neptune%satellite_model,             &              ! <->  TYPE  Satellite  model
                              neptune%solarsystem_model,           &              ! <->  TYPE  Solarsysten model
                              neptune%thirdbody_model,             &              ! <->  TYPE  Third body model
                              neptune%tides_model,                 &              ! <->  TYPE  Tides model
                              neptune%derivatives_model,           &              ! <->  TYPE  Derivatives model
                              neptune%reduction,                   &              ! <->  TYPE  Reduction
                              request_time,                        &              ! <--  DBL   requested time (s)
                              force_no_interpolation,              &              ! <--  LOGI  indicates no interpolation => perfect last step for big changes in acceleration (i.e. maneuvers)
                              prop_counter,                        &              ! <--> DBL   current time (s)
                              state_out%r,                         &              ! <--> DBL() radius vector (km)
                              state_out%v,                         &              ! <--> DBL() velocity vector (km/s)
                              reset,                               &              ! <--> INT   reset flag
                              delta_time                           &              ! -->  DBL   propagated time (s)
                              )
        end if
        
        ! write(cmess,'(a)') 'post: request_time '//toString(request_time)//', prop_counter '//toString(prop_counter)//', delta_time '//toString(delta_time)//', reset: '//toString(reset)
        ! call message(cmess, LOG_AND_STDOUT)
                  
        ! progress output immediately after integration step
        if(neptune%has_to_write_progress()) then
            dtmp = abs(prop_counter - start_epoch_sec)/abs(end_epoch_sec - start_epoch_sec)
            call write_progress(neptune%get_progress_file_name(), dtmp, neptune%get_progress_step())
            ! write(*,*) dtmp, prop_counter
        end if
        
        ! dtmp = abs(prop_counter - start_epoch_sec)/abs(end_epoch_sec - start_epoch_sec)
        ! write(*,"(A,F7.2,A,F15.4,A,F15.6)") "progress =",100*dtmp, " // prop_counter =", prop_counter," // mjd =",epochs(1)%mjd + prop_counter/86400.d0

        if(hasFailed()) then
          call neptune%output%close_open_files(neptune%numerical_integrator)
          return
        end if

        ! =================

        if(reset == 1) then     ! in case of a reset the last output state has to be re-established,
                                ! as the original state_out is overwritten after each call
          nresets            = nresets + 1
          propCounterAtReset = prop_counter

          if(nresets == MAX_RESETS) then    ! cancel execution if 'caught' within one step
            write(cresets,'(i2)') MAX_RESETS
            call setNeptuneError(E_INTEGRATION_ABORT, FATAL, (/cresets/))
            return
          end if

          if(neptune%output%get_output_switch() .or. neptune%getStoreDataFlag()) then    ! this has only to be done in a mode, where output is requested...
            if((.not. flag_backward .and. (lastPropCounterSuccess < lastPropCounter)) .or. &
               (      flag_backward .and. (lastPropCounterSuccess > lastPropCounter))) then
              ! This is the usual case when the intergator cannot go on and needs a valid state to restart from.
              state_out    = last_state_out
              prop_counter = lastPropCounter
            else if ((.not. flag_backward .and. (prop_counter > request_time)) .or. &
                    (       flag_backward .and. (prop_counter < request_time))) then
              ! This case may happen when the integrator oversteps the requested time
              !  and usually would start interpolation. Though, due to integration issues
              !  he never gets to the point of interpolation but sets the prop_counter anyway.
              !  This leads to an invalid state, which is also written to output at a time,
              !  which was not even requested.
              state_out    = last_state_out
              prop_counter = lastPropCounter
            end if
          end if
        else if((.not. flag_backward .and. (prop_counter > propCounterAtReset)) .or. &
                (      flag_backward .and. (prop_counter < propCounterAtReset))) then  ! reset 'nresets' flag
          nresets = 0
          lastPropCounterSuccess = prop_counter
        end if
      
        
        ! diff = prop_counter - request_time
        ! write(cmess, '(a, D15.6, a, a)') 'diff = ', diff, " // intermediate_integrator_call=", toString(intermediate_integrator_call)
        ! call message(cmess, LOG_AND_STDOUT)

        ! Exit the integration loop when the desired time is reached
        if (intermediate_integrator_call) then
          if (abs(prop_counter - request_time) < epsilon(1.0d0)) exit
        else
          if(nep_clock%has_finished_step(prop_counter)) exit
        end if
        ! =================

      end do
      !---------------------------------------------------------------------------------------------------------

      nresets = 0   ! reset the counter as step was successful

      !===========================================================================
      !
      ! Update of state error transition matrix and computation of covariance
      ! matrix for output
      !
      !---------------------------------------------------------------------------
      if(neptune%numerical_integrator%getCovariancePropagationFlag() &
        .and. nep_clock%has_to_update_covariance()) then
          call neptune%numerical_integrator%getStateTransitionMatrix(           &
                                       neptune%gravity_model,                   & ! <--> TYPE Gravity model
                                       neptune%atmosphere_model,                & ! <--> TYPE Atmosphere model
                                       neptune%manoeuvres_model,                & ! <--> TYPE Manoeuvres model
                                       neptune%radiation_model,                 & ! <--> TYPE Radiation model
                                       neptune%satellite_model,                 & ! <--> TYPE Satellite model
                                       neptune%solarsystem_model,               & ! <--> TYPE Solarsystem model
                                       neptune%thirdbody_model,                 & ! <--> TYPE Third body model
                                       neptune%tides_model,                     & ! <--> TYPE Tides model
                                       neptune%derivatives_model,               & ! <->  TYPE Atmosphere model
                                       neptune%reduction,                       & ! <->  TYPE Reduction
                                       state_out%r,                             & ! <--  DBL() radius vector (km)
                                       state_out%v,                             & ! <--  DBL() velocity vector (km/s)
                                       request_time,                            & ! <--  DBL   requested time
                                       set                                      & ! <--> DBL() state error transition matrix
                                     )
          if(hasFailed()) return

          !** cumulate state transition matrix
          cumSet = matmul(set,cumSet)
          set_out%elem = cumSet

          !** compute new covariance matrix for given time
          covar_out%elem = matmul(matmul(cumSet,covar_in%elem),transpose(cumSet))
          if(neptune%correlation_model%getNoisePropagationFlag()) then
              corrMat                 = neptune%correlation_model%getCorrelationMatrix(request_time)
              covar_out%elem(1:6,1:6) = covar_out%elem(1:6,1:6) + corrMat
          end if
          call nep_clock%toggle_cov_update_flag()

      else if(neptune%numerical_integrator%getCovariancePropagationFlag() &
        .and. nep_clock%has_to_save_covariance()) then
          call neptune%numerical_integrator%save_covariance_state(state_out, request_time)
      end if
    end do

    !** close open files
    call neptune%output%close_open_files(neptune%numerical_integrator)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine propagate_set

end module libneptune
