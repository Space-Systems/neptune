!>-----------------------------------------------------------------------------------------------
!> @anchor      numint
!!
!> @brief       Numerical integration based on Cowell's method
!> @author      Vitali Braun (VB)
!> @brief       Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  01.2013 (initial design)</li>
!!                <li>VB:  04.06.2013 (code optimization) </li>
!!                <li>VB:  04.07.2013 (covariance propagation) </li>
!!                <li>VB:  08.06.2015 (added infrastructure for additional integrators) </li>
!!                <li>CHK: 16.11.2015 (updated to use with libslam) </li>
!!                <li>VB:  24.01.2016 (changed the way counters are used: now controlled via module variables, which can be reset by NEPTUNE)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>VB:  14.07.2017 (added method to return covariance integration method)</li>
!!                <li>CHK: 02.01.2018 (updated to use NEPTUNE class) </li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for the
!!              numerical integration of the equations of motion.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      numint
!!
!!------------------------------------------------------------------------------------------------
module numint

  use slam_types,             only: dp
  use slam_io,                only: openFile, SEQUENTIAL, OUT_FORMATTED_OVERWRITE, closeFile
  use slam_math,              only: identity_matrix
  use neptune_error_handling, only: E_INTEGRATION_METHOD, E_INTEGRATION_STARTUP, E_ABS_TOLERANCE, E_REL_TOLERANCE, E_INTEGRATION_RESET, &
                                    E_SRP_CORRECTION, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, &
                                    E_UNKNOWN_PARAMETER, E_SPECIAL, checkIn, checkOut
  use slam_orbit_types,       only: state_t
  use slam_reduction_class,   only: Reduction_type
  use gravity,                only: Gravity_class
  use atmosphere,             only: Atmosphere_class
  use maneuvers,              only: Manoeuvres_class
  use derivatives,            only: Derivatives_class, PERT_SRP
  use radiation,              only: Radiation_class
  use satellite,              only: Satellite_class
  use solarsystem,            only: Solarsystem_class
  use thirdbody,              only: Thirdbody_class
  use tides,                  only: Tides_class
  use slam_time,              only: time_t, sec_per_day, mjd2gd, dayfraction2hms
  use ieee_arithmetic,        only: ieee_is_nan

  implicit none

  private

  character(len=*), parameter   :: IntLogfileName = "numint.log"                ! logfile name
  character(len=*), parameter   :: C_INT_METHOD_VSSC = 'VarStep Stoermer-Cowell'
  character(len=*), parameter   :: C_INT_METHOD_FGJ  = 'FixStep Gauss-Jackson'

  integer, parameter, public    :: TAYLOR = 1
  integer, parameter, public    :: RK4    = 2

  integer, parameter            :: INT_METHOD_VSSC = 1                          !< Variable step Störmer-Cowell according to M. Berry (2004)
  integer, parameter            :: INT_METHOD_FGJ  = 2                          !< Fixed step Gauß-Jackson
  integer, parameter, public    :: MAX_RESETS = 3                               !< maximum number of integrator resets for one single step, after which integration is cancelled.
  integer, parameter            :: setdim = 6                                   ! dimension of state error transition matrix
  integer, parameter            :: kmax   = 9                                   ! maximum value of k

  ! Numerical intergrator class (type defintion)
    type, public :: Numint_class

      integer :: cov_int_method                                                 !< covariance matrix integration method
      integer :: intlog                                                         !< output channel of integration logfile
      integer :: intMethod = INT_METHOD_VSSC                                    !< state vector integration method, default is variable step Störmer-Cowell
      integer :: countCallsVarstormcow                                          !< counts the number of calls to the routine varstormcow - needs a call to resetCountIntegrator(), to be again at zero, which is done by neptune.
      integer :: countCallsGaussJackson                                         !< counts the number of calls to the routine fixedGaussJackson - needs a call to resetCountIntegrator(), to be again at zero, which is done by neptune.
      integer :: countCallsSetMatrix                                            !< counts the number of calls to the routine getStateTransitionMatrix - needs a call to resetCountSetMatrix(), to be again at zero, which is done by neptune.

      logical :: correctSRP                                                     !< flag for SRP Lundberg correction algorithm
      logical :: flog                                                           !< flag for logfile output
      logical :: propagateCovariance                                            !< covariance matrix propagation on/off
      logical :: setReset                                                       !< reset state transition matrix integration
      logical :: is_start_epoch_set                                             !< to check if start epoch has been set
      logical :: is_start_epoch_reset                                           !< to check if start epoch has changed

      real(dp)  :: eps_abs                                                      !< absolute tolerance
      real(dp)  :: eps_rel                                                      !< relative tolerance

      real(dp)  :: covIntegrationStep                                           !< covariance matrix integration step size
      real(dp)  :: saved_time_offset                                            !< saved time offset, e.g. at half-step

      type(state_t) :: initialState                                             !< initial state vector
      type(state_t) :: saved_state                                              !< saved state, e.g. at half-step
      type(time_t)  :: start_epoch                                              !< initial epoch corresponding to initial state

      type(state_t) :: lastState                                                !< stored state of last step

      ! variables with formerly used save statement
      real(dp), dimension(3,kmax+1), private :: diff                            ! modified divided differences
      real(dp), private                      :: r                               ! amount to change the step by on next step
      real(dp), dimension(kmax), private     :: steps                           ! step size history
      real(dp), dimension(kmax), private     :: psi                             ! sums of steps
      real(dp), dimension(kmax-1), private   :: psin                            ! sums of steps,
                                                                                ! starting 1 step back
      real(dp), dimension(3), private        :: lastpos                         ! position from last step
      real(dp), dimension(3), private        :: lastvel                         ! velocity from last step
      real(dp), private                      :: errs                            ! error estimate on velocity
      real(dp), private                      :: errd                            ! error estimate on position
      real(dp), private                      :: eps                             ! tolerance
      real(dp), private                      :: releps                          ! relative error / eps
      real(dp), private                      :: abseps                          ! absolute error / eps
      real(dp), dimension(3), private        :: intpos                          ! last position integrated to
      real(dp), dimension(3), private        :: intvel                          ! last velocity integrated to
      real(dp), private                      :: inttime                         ! last time integrated to
                                                                                ! inttime is the same as currtime
                                                                                ! unless the last step had an
                                                                                ! interpolation
      real(dp), private                      :: stepsize                        ! step size
      integer, private                       :: intdir                          ! integration direction
                                                                                !  1 = forward
                                                                                ! -1 = backward
      integer, private                       :: k                               ! number of backpoints
      integer, private                       :: lk                              ! k in last integration step for interpolator
      integer, private                       :: fail                            ! number of consecutive failures
      logical, private                       :: flag_srp                        ! becomes true as soon as SRP perturbation has been activated and varstormcow is initialized



    contains

        !** get
        procedure,public :: getIntegrationLogfileChannel
        procedure,public :: getIntegrationLogfileName
        procedure,public :: getRelativeTolerance
        procedure,public :: getAbsoluteTolerance
        procedure,public :: getCovarianceIntegrationStep
        procedure,public :: getCovariancePropagationFlag
        procedure,public :: get_covariance_integration_method
        procedure,public :: getStateTransitionMatrix
        procedure,public :: get_current_step_size

        !** set
        procedure,public :: resetCountIntegrator
        procedure,public :: resetCountSetMatrix
        procedure,public :: save_covariance_state
        procedure,public :: setIntegrationTolerance
        procedure,public :: setIntegrationLogfile
        procedure,public :: setIntegrationMethod
        procedure,public :: setCovarianceIntegrationMethod
        procedure,public :: setCovarianceIntegrationStep
        procedure,public :: setCovariancePropagationFlag
        procedure,public :: setSrpCorrect
        procedure,public :: set_start_epoch
        !** others
        procedure,public  :: integrateStep
        procedure,private :: varstormcow
        procedure,private :: srpCorrection

    end type Numint_class

        ! Constructor
    interface Numint_class
        module procedure constructor
    end interface Numint_class

contains

    !===========================================================================
    !!
    !>  @anchor     constructor
    !!
    !!  @brief      Initializes all variables with default values
    !!  @author     Christopher Kebschull
    !!
    !!
    !!  @date       <ul>
    !!                <li>02.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    type(Numint_class) function constructor()
        constructor%countCallsVarstormcow  = 0                                  !< counts the number of calls to the routine varstormcow - needs a call to resetCountIntegrator(), to be again at zero, which is done by neptune.
        constructor%countCallsGaussJackson = 0                                  !< counts the number of calls to the routine fixedGaussJackson - needs a call to resetCountIntegrator(), to be again at zero, which is done by neptune.
        constructor%countCallsSetMatrix    = 0                                  !< counts the number of calls to the routine getStateTransitionMatrix - needs a call to resetCountSetMatrix(), to be again at zero, which is done by neptune.

        constructor%correctSRP          = .false.                               !< flag for SRP Lundberg correction algorithm
        constructor%flog                = .false.                               !< flag for logfile output
        constructor%propagateCovariance = .false.                               !< covariance matrix propagation on/off
        constructor%setReset            = .true.                                !< reset state transition matrix integration
        constructor%is_start_epoch_set  = .false.                               !< to check if start epoch has been set
        constructor%is_start_epoch_reset= .false.                               !< to check if start epoch has changed

        constructor%eps_abs = 1.d-14                                            !< absolute tolerance
        constructor%eps_rel = 1.d-13                                            !< relative tolerance

        constructor%covIntegrationStep = 30.d0                                  !< covariance matrix integration step size
        constructor%flag_srp = .false.
        constructor%cov_int_method = RK4
    end function constructor

    !=====================================================================
    !
    !> @brief     Saves a state vector to be used during the integration, e.g. during half-step
    !!
    !> @brief    Vitali Braun
    !!
    !> @param[in]  s state vector of type state_t
    !> @param[in]  t_offset time offset
    !!
    !> @date      <ul>
    !!              <li>VB: 14.07.2017 (initial design) </li>
    !!            </ul>
    !!
    !! @anchor    save_covariance_state
    !---------------------------------------------------------------------
    subroutine save_covariance_state(this, s, t_offset)
        implicit none
        class(Numint_class),intent(inout)           :: this
        type(state_t), intent(in)                   :: s
        real(dp), intent(in)                        :: t_offset
        this%saved_state = s
        this%saved_time_offset = t_offset
        return
    end subroutine

    !=====================================================================
    !
    !> @brief     Retrieves the covariance integration method
    !!
    !> @brief    Vitali Braun
    !!
    !! @return    Covariance integration method
    !!
    !> @date      <ul>
    !!              <li>2017-07-14 (initial design) </li>
    !!            </ul>
    !!
    !! @anchor    get_covariance_integration_method
    !---------------------------------------------------------------------
    pure function get_covariance_integration_method(this) result(i)
        implicit none
        class(Numint_class),intent(in)              :: this
        integer                                     :: i
        i = this%cov_int_method
        return
    end function

  !=====================================================================
  !
  !> @brief     Set the start epoch which is added to the offset during propagation
  !!
  !> @brief    Vitali Braun
  !!
  !> @param[in] epoch   Start epoch
  !!
  !> @date      <ul>
  !!              <li>2017-07-14 (initial design) </li>
  !!            </ul>
  !!
  !! @anchor    set_start_epoch
  !-----------------------------------------------------------------
    subroutine set_start_epoch(this,epoch)
        implicit none
        class(Numint_class),intent(inout)           :: this
        type(time_t), intent(in)                    :: epoch
        this%start_epoch = epoch
        this%is_start_epoch_reset = .true.
        this%is_start_epoch_set = .true.
        return
    end subroutine

!========================================================================
!
!> @anchor      setSrpCorrect
!!
!> @brief       Set the flag for the SRP shadow boundary correction
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 14.03.2014 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setSrpCorrect(this,val)
    class(Numint_class),intent(inout)           :: this
    logical, intent(in)                         :: val
    this%correctSrp = val
    return
  end subroutine

!========================================================================
!
!> @anchor      getIntegrationLogfileChannel
!!
!> @brief       Returns the device number of the integration logfile
!> @author      Vitali Braun
!!
!> @return      device number
!!
!> @date        <ul>
!!                <li> 01.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!------------------------------------------------------------------------
  integer function getIntegrationLogfileChannel(this)
    class(Numint_class),intent(inout)       :: this
    getIntegrationLogfileChannel = this%intlog
    return
  end function

!========================================================================
!
!> @anchor      getIntegrationLogfileName
!!
!> @brief       Returns the name of the integration logfile
!> @author      Vitali Braun
!!
!> @return      name of the integration logfile as string
!!
!> @date        <ul>
!!                <li> 01.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!------------------------------------------------------------------------
  character(len=len(IntLogFileName)) function getIntegrationLogfileName(this)
    class(Numint_class),intent(inout)       :: this
    getIntegrationLogfileName = IntLogfileName
    return
  end function
!========================================================================
!
!> @anchor      getAbsoluteTolerance
!!
!> @brief       Returns the absolute integration tolerance
!> @author      Vitali Braun
!!
!> @return      absolute integration tolerance
!!
!> @date        <ul>
!!                <li> 01.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!------------------------------------------------------------------------
  real(dp) function getAbsoluteTolerance(this)
    class(Numint_class),intent(inout)       :: this
    getAbsoluteTolerance = this%eps_abs
    return
  end function

!========================================================================
!
!> @anchor      getRelativeTolerance
!!
!> @brief       Returns the relative integration tolerance
!> @author      Vitali Braun
!!
!> @return      relative integration tolerance
!!
!> @date        <ul>
!!                <li> 01.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!------------------------------------------------------------------------
  real(dp) function getRelativeTolerance(this)
    class(Numint_class),intent(inout)       :: this
    getRelativeTolerance = this%eps_rel
    return
  end function

!========================================================================
!
!> @anchor      setIntegrationMethod
!!
!> @brief       Set the integration method
!> @author      Vitali Braun
!!
!> @param[in]   method name as string
!!
!> @date        <ul>
!!                <li> 08.06.2015 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setIntegrationMethod(this,method)

    class(Numint_class),intent(inout)           :: this
    integer, intent(in)                         :: method
    character(len=*), parameter                 :: csubid = "setIntegrationMethod"
    character(len=len(C_INT_METHOD_VSSC))       :: ctemp1, ctemp2

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(method == INT_METHOD_VSSC .or. method == INT_METHOD_FGJ) then
      this%intMethod = method
    else
      write(ctemp1,'(i1)') INT_METHOD_VSSC
      write(ctemp2,'(i1)') INT_METHOD_FGJ
      call setNeptuneError(E_INTEGRATION_METHOD, FATAL, (/C_INT_METHOD_VSSC,ctemp1,C_INT_METHOD_FGJ,ctemp2/))
      return
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine

!========================================================================
!
!> @anchor      integrateStep
!!
!> @brief       Wrapper to perform one integration step with the selected integration method
!> @author      Vitali Braun
!!
!> @param[inout]  atmosphere_model   Atmosphere model instance
!> @param[in]     rqtime    Requested time step
!> @param[inout]  currtime  Current time
!> @param[inout]  pos       Position vector
!> @param[inout]  vel       Velocity vector
!> @param[inout]  reset     Reset flag: 1 if restarting
!> @param[out]    delt      Change in currtime
!!
!> @date        <ul>
!!                <li> 08.06.2015 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine integrateStep(this,                &
                            gravity_model,      &
                            atmosphere_model,   &
                            manoeuvres_model,   &
                            radiation_model,    &
                            satellite_model,    &
                            solarsystem_model,  &
                            thirdbody_model,    &
                            tides_model,        &
                            derivatives_model,  &
                            reduction,          &
                            rqtime,             &
                            force_no_interpolation, &
                            currtime,           &
                            pos,                &
                            vel,                &
                            reset,              &
                            delt)

    class(Numint_class),intent(inout)               :: this
    type(Gravity_class),intent(inout)               :: gravity_model            ! <-> Gravity model
    type(Atmosphere_class),intent(inout)            :: atmosphere_model         ! <-> Atmosphere model
    type(Manoeuvres_class),intent(inout)            :: manoeuvres_model         ! <-> Manoeuvres model
    type(Radiation_class),intent(inout)             :: radiation_model          ! <-> Radiation model
    type(Satellite_class),intent(inout)             :: satellite_model          ! <-> Satellite model
    type(Solarsystem_class),intent(inout)           :: solarsystem_model        ! <-> Solarsystem model
    type(Thirdbody_class),intent(inout)             :: thirdbody_model          ! <-> Third body model
    type(Tides_class),intent(inout)                 :: tides_model              ! <-> Tides model
    type(Derivatives_class),intent(inout)           :: derivatives_model        ! <-> Derivatives model
    type(Reduction_type),intent(inout)              :: reduction                ! <-> Reduction
    real(dp),                 intent(in)            :: rqtime                   ! <-- requested time / s
    logical,                  intent(in)            :: force_no_interpolation   ! <- indicates no interpolation => perfect last step for big changes in acceleration (i.e. maneuvers)
    real(dp),                 intent(inout)         :: currtime                 ! <-> current time
    real(dp), dimension(3),   intent(inout)         :: pos                      ! <-> position vector
    real(dp), dimension(3),   intent(inout)         :: vel                      ! <-> velocity vector
    integer,                  intent(inout)         :: reset                    ! <-> 1 if restarting
    real(dp),                 intent(out)           :: delt                     ! --> change in currtime

    character(len=*), parameter :: csubid = "integrateStep"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%intMethod == INT_METHOD_VSSC) then

      call this%varstormcow(                &
                          gravity_model,    &  ! <->  TYPE   Gravity model
                          atmosphere_model, &  ! <->  TYPE   Atmosphere model
                          manoeuvres_model, &  ! <->  TYPE   Manoeuvres model
                          radiation_model,  &  ! <->  TYPE   Radiation model
                          satellite_model,  &  ! <->  TYPE   Satelllite model
                          solarsystem_model,&  ! <->  TYPE   Solarsystem model
                          thirdbody_model,  &  ! <->  TYPE   Third body model
                          tides_model,      &  ! <->  TYPE   Tides model
                          derivatives_model,&  ! <->  TYPE   Derivatives model
                          reduction,        &  ! <->  TYPE   Reduction
                          rqtime,           &  ! <--  DBL    requested time
                          force_no_interpolation, &
                          currtime,         &  ! <--> DBL    current time
                          pos,              &  ! <--> DBL()  position vector
                          vel,              &  ! <--> DBL()  velocity vector
                          reset,            &  ! <--> INT    1 = restarting
                                       !             0 = normal return
                          delt              &  ! -->  DBL    change in current time
                      )
      if(hasFailed()) return

    else if(this%intMethod == INT_METHOD_FGJ) then

      call fixedgaussjackson(       &
                          rqtime,   &  ! <--  DBL    requested time
                          currtime, &  ! <--> DBL    current time
                          pos,      &  ! <--> DBL()  position vector
                          vel,      &  ! <--> DBL()  velocity vector
                          reset,    &  ! <--> INT    1 = restarting
                                       !             0 = normal return
                          delt      &  ! -->  DBL    change in current time
                      )
      write(*,*) "Gauss-Jackson integrator Not implemented yet. AnH"
      write(*,*) "17.06.2015"
      stop
      if(hasFailed()) return
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine

!========================================================================
!
!> @anchor      setIntegrationLogfile
!!
!> @brief       Opens the integration logfile
!> @author      Vitali Braun
!!
!> @param[in]   cswitch a string that should be either 'ON' or 'OFF'
!!
!> @date        <ul>
!!                <li> 01.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!                <li> 11.01.2015 (added cswitch input to also close the file again) </li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setIntegrationLogfile(this,cswitch)

    class(Numint_class),intent(inout)               :: this
    character(len=*), intent(in)                    :: cswitch

    character(len=*), parameter  :: csubid = "setIntegrationLogfile"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(cswitch == 'ON') then
      this%intlog = openFile(IntLogfileName,SEQUENTIAL,OUT_FORMATTED_OVERWRITE)
      this%flog   = .true.
    else if(cswitch == 'OFF' .and. this%flog) then
      this%intlog = closeFile(this%intlog)
      this%flog   = .false.
    else if(.not. cswitch == 'OFF') then
      call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/cswitch/))
      return
    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine

!========================================================================
!
!> @anchor      setIntegrationTolerance
!!
!> @brief       Set the relative and absolute integration tolerances
!> @author      Vitali Braun
!!
!> @param[in]   tabs - absolute tolerance
!> @param[in]   trel - relative tolerance
!!
!> @date        <ul>
!!                <li> 01.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setIntegrationTolerance( this, &
                                      tabs, &  ! <-- DBL absolute tolerance
                                      trel  &  ! <-- DBL relative tolerance
                                    )

    !** interface
    !------------------------------
    class(Numint_class),intent(inout)           :: this
    real(dp), optional, intent(in)              :: tabs
    real(dp), optional, intent(in)              :: trel
    !------------------------------

    character(len=255) :: cmess   ! message string
    character(len=*), parameter :: csubid = "setIntegrationTolerance"! subroutine id

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(present(tabs)) then
      this%eps_abs = abs(tabs)
    end if

    if(present(trel)) then
      this%eps_rel = abs(trel)
    end if

    if(this%eps_rel < epsilon(1.d0)) then   !** tolerance smaller than machine epsilon
      call setNeptuneError(E_REL_TOLERANCE, WARNING)
      this%eps_rel = 1.d1*epsilon(1.d0)
    end if


    if(this%eps_abs < epsilon(1.d0)) then   !** tolerance smaller than machine epsilon
      call setNeptuneError(E_ABS_TOLERANCE, WARNING)
      this%eps_abs = 1.d1*epsilon(1.d0)
    end if

    if(this%eps_abs/this%eps_rel >= 1.d0) then  !** relative tolerance should not
                                                !   be smaller than absolute tolerance. if so, rel tolerance is
                                                !   increased to be 10 times the abs tolerance

      this%eps_rel = this%eps_abs*10.d0

      write(cmess,'(a,e11.4e2)') "Relative tolerance should not be smaller than absolute tolerance."// &
                                 " Increasing relative tolerance to: ", this%eps_rel

      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine

!========================================================================
!
!> @anchor      getCovariancePropagationFlag
!!
!> @brief       Get the covariance matrix integration flag
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 09.07.2013 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  logical function getCovariancePropagationFlag(this)
    class(Numint_class),intent(inout)       :: this
    getCovariancePropagationFlag = this%propagateCovariance
    return
  end function getCovariancePropagationFlag

!========================================================================
!
!> @anchor      setCovarianceIntegrationMethod
!!
!> @brief       Set the covariance matrix integration method
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 16.07.2013 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setCovarianceIntegrationMethod(this,method)

    class(Numint_class),intent(inout)               :: this
    integer, intent(in)                             :: method

    integer :: last_method

    last_method = this%cov_int_method

    select case(method)
      case(TAYLOR)
        this%cov_int_method = TAYLOR
      case(RK4)
        this%cov_int_method = RK4
      case default
        this%cov_int_method = RK4
    end select

    !** Reset SET matrix integration
    if(this%cov_int_method /= last_method) this%setReset = .true.
    return

  end subroutine setCovarianceIntegrationMethod


!========================================================================
!
!> @anchor      setCovarianceIntegrationFlag
!!
!> @brief       Set the covariance matrix integration flag
!> @author      Vitali Braun
!!
!> @param[in]   flag - should be .true. to enable
!!
!> @date        <ul>
!!                <li> 04.07.2013 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setCovariancePropagationFlag(this,flag)
    class(Numint_class),intent(inout)               :: this
    logical, intent(in)                             :: flag
    this%propagateCovariance = flag
    return
  end subroutine setCovariancePropagationFlag

!========================================================================
!
!> @anchor      setCovarianceIntegrationStep
!!
!> @brief       Set the covariance matrix integration step size
!> @author      Vitali Braun
!!
!> @param[in]   step - integration step size
!!
!> @date        <ul>
!!                <li> 11.07.2013 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine setCovarianceIntegrationStep(this,step)
    class(Numint_class),intent(inout)       :: this
    real(dp), intent(in) :: step
    this%covIntegrationStep = step
    return
  end subroutine setCovarianceIntegrationStep

!========================================================================
!
!> @anchor      getCovarianceIntegrationStep
!!
!> @brief       Get the covariance matrix integration step size
!> @author      Vitali Braun
!!
!> @return      integration step size
!!
!> @date        <ul>
!!                <li> 11.07.2013 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
  real(dp) function getCovarianceIntegrationStep(this)
    class(Numint_class),intent(inout)       :: this
    getCovarianceIntegrationStep = this%covIntegrationStep
    return
  end function getCovarianceIntegrationStep

!========================================================================
!
!> @anchor      resetCountSetMatrix
!!
!> @brief       Resets the counter for the number of calls to the state error transition matrix getter
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 24.01.2016 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
subroutine resetCountSetMatrix(this)
  class(Numint_class),intent(inout)       :: this
  this%countCallsSetMatrix = 0
  return
end subroutine

!========================================================================
!
!> @anchor      getStateTransitionMatrix
!!
!> @brief       Compute the state error transition matrix from (t) to (reqt-t)
!> @author      Vitali Braun
!!
!!  @param[in]  r     Radius vector in GCRF / km
!!  @param[in]  v     Velocity vector in GCRF / km/s
!!  @param[in]  reqt  Time offset wrt. the epoch of the initial state / s
!!  @param[out] set   State error transition matrix
!!
!!  @date       <ul>
!!                <li>VB 10.07.2013 (initial design)</li>
!!                <li>VB 07.08.2013 (adapted to 'getStateTransitionMatrix')</li>
!!                <li>VB 24.01.2016 (fixed issue with call counter: now via modul variable and reset function, which is called by NEPTUNE)<li>
!!                <li>VB 14.07.2017 (changed to time offset in parameter list and adding to start epoch)<li>
!!              </ul>
!!
!------------------------------------------------------------------------
  subroutine getStateTransitionMatrix(  this,               &
                                        gravity_model,      &
                                        atmosphere_model,   &
                                        manoeuvres_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        tides_model,        &
                                        derivatives_model,  &
                                        reduction,          &
                                        r,                  &
                                        v,                  &
                                        reqt,               &
                                        set                 &
                                     )

    class(Numint_class),intent(inout)               :: this
    type(Gravity_class),intent(inout)               :: gravity_model
    type(Atmosphere_class),intent(inout)            :: atmosphere_model
    type(Manoeuvres_class),intent(inout)            :: manoeuvres_model
    type(Radiation_class),intent(inout)             :: radiation_model
    type(Satellite_class),intent(inout)             :: satellite_model
    type(Solarsystem_class),intent(inout)           :: solarsystem_model
    type(Thirdbody_class),intent(inout)             :: thirdbody_model
    type(Tides_class),intent(inout)                 :: tides_model
    type(Derivatives_class),intent(inout)           :: derivatives_model
    type(Reduction_type),intent(inout)             :: reduction
    real(dp), dimension(3), intent(in)              :: r
    real(dp), dimension(3), intent(in)              :: v
    real(dp),               intent(in)              :: reqt
    real(dp), dimension(setdim,setdim), intent(out) :: set

    real(dp), dimension(setdim,setdim) :: k1,k2,k3,k4
    real(dp), dimension(setdim,setdim) :: pdm

    integer :: reset

    real(dp) :: dt
    real(dp) :: reqt_loc
    real(dp) :: time_mjd        ! MJD of current point in time
    real(dp) :: t1,t2,t3        ! auxiliaries
    real(dp) :: t1d, t2d, t3d   ! auxiliaries
    real(dp) :: eps8 = 1.d-8    !< auxiliary


    type(state_t), dimension(3) :: state    ! different states for RK4 method (3 points in time)

    ! NOTE: the time_mjd is computed from the passed offset with the start epoch - the latter is set by NEPTUNE's input routine
    !       Be careful when calling this from any other module.
    reqt_loc = reqt
    time_mjd = this%start_epoch%mjd + reqt_loc/sec_per_day
    reset    = 0

    this%countCallsSetMatrix = this%countCallsSetMatrix + 1

    if(this%cov_int_method == TAYLOR) then

      !** reset state error transition matrix
      call identity_matrix(set)

      !** get partial derivatives matrix
      call derivatives_model%deriv_cov(                     &
                                        gravity_model,      &
                                        atmosphere_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        reduction,          &
                                        r,                  &
                                        v,                  &
                                        set,                &
                                        time_mjd,           &
                                        pdm                 &
                                        )
      if(hasFailed()) return

      set = set + this%covIntegrationStep*pdm &
                            + 0.5d0*this%covIntegrationStep**2.d0*matmul(pdm,pdm) &
                            + 6.d0*this%covIntegrationStep**3.d0*matmul(matmul(pdm,pdm),pdm)

      !** update covariance matrix
      !cov = matmul(matmul(set,cov),transpose(set))

      !** save state, if method is changed to RK4
      if(this%countCallsSetMatrix == 1) then
        state(1) = this%initialState
      else
        state(1)%r = r
        state(1)%v = v
      end if

    else if(this%cov_int_method == RK4) then

      !** reset state error transition matrix
      call identity_matrix(set)

      !** if called for the first time, use initial state vector
      if(this%countCallsSetMatrix == 1) then
        state(1) = this%initialState
      else
        state(1) = this%lastState
      end if

      state(3)%r = r
      state(3)%v = v

      t1 = reqt - this%covIntegrationStep
      t2 = reqt - 0.5d0*this%covIntegrationStep
      t3 = reqt

      t1d = this%start_epoch%mjd + t1/sec_per_day
      t2d = this%start_epoch%mjd + t2/sec_per_day
      t3d = this%start_epoch%mjd + t3/sec_per_day

      if(abs(this%saved_time_offset - t2) < eps8) then
          state(2) = this%saved_state
      else
          !read(*,*)
          !** find states 2 and 3 through interpolation
          call this%varstormcow(                             &
                        gravity_model,                       &
                        atmosphere_model,                    &
                        manoeuvres_model,                    &
                        radiation_model,                     &
                        satellite_model,                     &
                        solarsystem_model,                   &
                        thirdbody_model,                     &
                        tides_model,                         &
                        derivatives_model,                   &
                        reduction,                           &
                        t2,                                  &      ! <--  DBL   requested time (s)
                        .false.,                             &
                        reqt_loc,                            &      ! <--  DBL   current time (s)
                        state(2)%r,                          &      ! <--> DBL() radius vector (km)
                        state(2)%v,                          &      ! <--> DBL() velocity vector (km/s)
                        reset,                               &      ! <--  INT   reset flag
                        dt                                   &      ! -->  DBL   propagated time (s)
                      )
          if(hasFailed()) return
      end if

      !** evaluate derivative of SET matrix for RK states
      call derivatives_model%deriv_cov(                     &
                                        gravity_model,      &
                                        atmosphere_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        reduction,          &
                                        state(1)%r,         &
                                        state(1)%v,         &
                                        set,                &
                                        t1d,                &
                                        k1)
      if(hasFailed()) return
      call derivatives_model%deriv_cov(gravity_model,       &
                                        atmosphere_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        reduction,          &
                                        state(2)%r,         &
                                        state(2)%v,         &
                                        set + 0.5d0*this%covIntegrationStep*k1, &
                                        t2d,                &
                                        k2)
      if(hasFailed()) return
      call derivatives_model%deriv_cov(                     &
                                        gravity_model,      &
                                        atmosphere_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        reduction,          &
                                        state(2)%r,         &
                                        state(2)%v,         &
                                        set + 0.5d0*this%covIntegrationStep*k2,&
                                        t2d,                &
                                        k3)
      if(hasFailed()) return
      call derivatives_model%deriv_cov(gravity_model,       &
                                        atmosphere_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        reduction,          &
                                        state(3)%r,         &
                                        state(3)%v,         &
                                        set + this%covIntegrationStep*k3,&
                                        t3d,                &
                                        k4)
      if(hasFailed()) return

      set = set + this%covIntegrationStep/6.d0*(k1 + 2.0*k2 + 2.d0*k3 + k4)
      !write(52,'(36(e14.7e2,x))') set

    !write(*,*)
    !write(*,*) "set with drag"
    !do i=1,6
    !  write(*,'(6(e22.15e2))') (set(i,j), j=1,6)
    !end do
    !read(*,*)

      !** save last state
      this%lastState = state(3)

    end if

    return
! write(*,*) "-----------------------------------"
! write(*,*) "PDM:"
! do i=1,6
!   write(*,'(6(e15.6e2,x))') (pdm(i,j), j=1,6)
! end do
! write(*,*) "----------------------------------"
!  read(*,*)
! write(*,*) "-----------------------------------"
! write(*,*) "set:"
! do i=1,6
!   write(*,'(6(e15.6e2,x))') (set(i,j), j=1,6)
! end do
! write(*,*) "----------------------------------"
!  read(*,*)
!  write(*,*) "step = ", covIntegrationStep

! write(*,*) "-----------------------------------"
! write(*,*) "cov:"
! do i=1,6
!   write(*,'(6(e15.6e2,x))') (cov(i,j), j=1,6)
! end do
! write(*,*) "----------------------------------"
!  read(*,*)

!   if(called == 1) then

!     open(unit=111, file="set.out")
!     open(unit=112, file="rv.out")

!   end if

!   write(111,'(36(1x,e10.4e2))') ((set(i,j), i=1,6), j=1,6)
!   write(112,'(3(1x,f9.3),3(1x,f9.5))') r, v


  end subroutine getStateTransitionMatrix

!========================================================================
!
!> @anchor      resetCountIntegrator
!!
!> @brief       Resets the counter for the number of calls to the integrator(s)
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 24.01.2016 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
subroutine resetCountIntegrator(this)
  class(Numint_class),intent(inout)       :: this
  this%countCallsVarstormcow  = 0
  this%countCallsGaussJackson = 0
  return
end subroutine

!========================================================================
!
!> @anchor      get_current_step_size
!!
!> @brief       Returns the step size of the integrator
!> @author      Christopher Kebschull
!!
!> @date        <ul>
!!                <li> 29.10.2020 (initial design)</li>
!!              </ul>
!!
!------------------------------------------------------------------------
real(dp) function get_current_step_size(this)
  class(Numint_class),intent(inout)       :: this

  get_current_step_size = this%stepsize

  return
end function


!========================================================================
!
!> @anchor      varstormcow
!!
!> @brief       Störmer-Cowell integrator based on Berry (2004)
!> @author      Vitali Braun
!!
!> @param[in] gravity_model       the gravity model
!> @param[in] atmosphere_model    the atmosphere model
!> @param[in] manoeuvres_model    the manoeuvres model
!> @param[in] radiation_model     the gravity model
!> @param[in] satellite_model     the satellite model
!> @param[in] solarsystem_model   the solar system model
!> @param[in] thirdbody_model     the third body model
!> @param[in] tides_model         the tides model
!> @param[in] reduction           the reduction model
!> @param[in] rqtime              requested time
!> @param[in] currtime            current time
!> @param[in] pos                 position vector
!> @param[in] vel                 velocity vector
!> @param[in] reset               1 = restarting
!> @param[in] delt                change in current time
!!
!> @date        <ul>
!!                <li> 09.2012 (initial design)             </li>
!!                <li> 04.06.2013 (code optimization)       </li>
!!                <li> 04.07.2013 (added covariance matrix -- but did not work, so removed immediately) </li>
!!                <li> 10.04.2015 (fixed a bug in the stepsize control: 'lastpos/lastvel --> pos/vel' in wts/wtd computation </li>
!!                <li> 24.01.2016 (fixed: for rqtime == inttime, the integrated values are returned - was indefinite before! In addition, the number of calls is now
!!                                 counted via a dedicated module counter, which can be reset.</li>
!!                <li>CHK: 02.01.2017 (updated to use NEPTUNE class) </li>
!!              </ul>
!!
!> @detail      The subroutine varstormcow performs one step of the variable-step
!!              Stormer-Cowell algorithm. If the integration step takes the
!!              current time beyond the request time, the routine will
!!              interpolate to the request time. The routine changes the value
!!              of current time, position, and velocity during the call. The
!!              routine sets the value of delt, which is the amount the
!!              current time changed.
!!              Call with reset = 1 to start the integrator. The routine will
!!              normally return with reset = 0. If the routine returns
!!              reset = 1, the method must restart because a step has failed
!!              after 3 tries.
!!              The algorithms used in this code are explained in
!!              "A Variable-Step Double-Integration Multi-Step Integrator" by
!!              Matthew M. Berry. PhD Dissertation, Virginia Polytechnic
!!              Institute and State University, Blacksburg, VA. April 2004.
!!
!------------------------------------------------------------------------
  subroutine varstormcow( this,             &
                          gravity_model,    &                                   ! <->  TYPE   Gravity model
                          atmosphere_model, &                                   ! <->  TYPE   Atmosphere model
                          manoeuvres_model, &                                   ! <->  TYPE   Manoeuvres model
                          radiation_model,  &                                   ! <->  TYPE   Radiation model
                          satellite_model,  &                                   ! <->  TYPE   Satellite model
                          solarsystem_model,&                                   ! <->  TYPE   Solarsystem model
                          thirdbody_model,  &                                   ! <->  TYPE   Third body model
                          tides_model,      &                                   ! <->  TYPE   Tides model
                          derivatives_model,&                                   ! <->  TYPE   Derivatives model
                          reduction,        &                                   ! <->  TYPE   Reduction
                          rqtime,           &                                   ! <--  DBL    requested time
                          force_no_interpolation, &                             ! <--  LOGI   indicates no interpolation => perfect last step for big changes in acceleration (i.e. maneuvers)
                          currtime,         &                                   ! <--> DBL    current time
                          pos,              &                                   ! <--> DBL()  position vector
                          vel,              &                                   ! <--> DBL()  velocity vector
                          reset,            &                                   ! <--> INT    1 = restarting
                                                                                !             0 = normal return
                          delt              &                                   ! -->  DBL    change in current time
                        )

    !** interface
    !------------------------------------------------------------------------
    class(Numint_class),intent(inout)               :: this
    type(Gravity_class),intent(inout)               :: gravity_model            ! Gravity model
    type(Atmosphere_class),intent(inout)            :: atmosphere_model         ! Atmosphere model
    type(Manoeuvres_class),intent(inout)            :: manoeuvres_model         ! Manoeuvres model
    type(Radiation_class),intent(inout)             :: radiation_model          ! Radiation model
    type(Satellite_class),intent(inout)             :: satellite_model          ! Satellite model
    type(Solarsystem_class),intent(inout)           :: solarsystem_model        ! Solarsystem model
    type(Thirdbody_class),intent(inout)             :: thirdbody_model          ! Third body model
    type(Tides_class),intent(inout)                 :: tides_model              ! Tides model
    type(Derivatives_class),intent(inout)           :: derivatives_model        ! Derivatives model
    type(Reduction_type),intent(inout)              :: reduction                ! Reduction
    real(dp),                 intent(in)            :: rqtime                   ! requested time
    logical,                intent(in)              :: force_no_interpolation   !

    real(dp),                 intent(inout)         :: currtime                 ! current time
    real(dp), dimension(3),   intent(inout)         :: pos                      ! position vector
    real(dp), dimension(3),   intent(inout)         :: vel                      ! velocity vector
    integer,                  intent(inout)         :: reset                    ! 1 if restarting

    real(dp),                 intent(out)           :: delt                     ! change in currtime
    !------------------------------------------------------------------------

    ! during call
    ! note: the units of rqtime, currtime, pos, vel, and accel must
    ! all be compatible

    !------------------------------------------------------------------------
    ! parameters
    character(len=*), parameter        :: csubid = "varstormcow"
    !------------------------------------------------------------------------
    ! local variables
    character(len=3)                   :: logInfo
    real(dp), dimension(3)             :: a_srp         ! current SRP acceleration
    real(dp), dimension(3)             :: a_srp_wo      ! SRP acceleration without shadow! (f=1)
    real(dp)                           :: stepsize2     ! step size squared
    real(dp), dimension(3)             :: accel         ! acceleration vector
    real(dp), dimension(3)             :: alpha_e       ! angle satellite -> earth
    real(dp), dimension(3)             :: alpha_s       ! angle satellite -> sun
    real(dp), dimension(3,kmax)        :: diffsave      ! diff from last step
    real(dp), dimension(kmax)          :: alpha         ! ratio of current step to psi
    real(dp), dimension(0:kmax-2)      :: psinm1        ! sums of steps,
                                                        ! starting 2 steps back
    real(dp), dimension(kmax+2,kmax+1) :: g             ! integration coefficients
    real(dp), dimension(kmax+2,kmax+1) :: gp            ! integration coefficients
    real(dp), dimension(kmax)          :: beta          ! ratio of phi’s
    real(dp), dimension(3)             :: pos2back      ! position from 2 steps ago
    real(dp), dimension(3,kmax+1)      :: newdiff       ! differences after evaluation
    real(dp)                           :: savecurrtime  ! last value of current time
    real(dp)                           :: erkd          ! position error estimate
    real(dp)                           :: erks          ! velocity error estimate
    real(dp)                           :: terrd         ! temp value for errd
    real(dp)                           :: terrs         ! temp value for errs
    real(dp)                           :: rs            ! r value for single integration
    real(dp)                           :: rd            ! r value for double integration
    real(dp)                           :: sigma         ! parameter used in error estimate
    real(dp)                           :: savestep      ! saved step size of
                                                  ! furthest backpoint
    real(dp)                           :: ratio         ! ratio of latest step sizes
    real(dp)                           :: invrat        ! inverse of ratio
    real(dp)                           :: sum1          ! integration sum for velocity
    real(dp)                           :: sum2          ! integration sum for position
    real(dp), dimension(kmax)          :: gammastars = (/-0.5d0,                 -0.08333333333333333d0,      &
                                                         -0.04166666666666667d0, -0.02638888888888891d0,      &
                                                         -0.01874999999999999d0, -0.01426917989417986d0,      &
                                                         -0.01136739417989419d0, -0.009356536596119958d0,     &
                                                         -0.007892554012345676d0/)
                                                        ! difference of single fixed-step coefficients
    real(dp), dimension(kmax)          :: gammastard = (/-1.d0, 0.08333333333333333d0, 0.0d0, &
                                                         -0.004166666666666666d0, -0.004166666666666666d0,    &
                                                         -0.003654100529100521d0, -0.003141534391534404d0,    &
                                                         -0.002708608906525564d0, -0.002355324074074072d0/)
                                                        ! difference of double fixed-step coefficients
    real(dp)                           :: macheps       ! machine epsilon
    real(dp), dimension(3)             :: wts           ! error weighting for velocity
    real(dp), dimension(3)             :: wtd           ! error weighting for position
    real(dp)                           :: h1s           ! initial step for single integration
    real(dp)                           :: h1d           ! initial step for double integration
    real(dp)                           :: fracday
    real(dp), dimension(3)             :: fshadow       ! shadow factors [0:1] for three subsequent points
    real(dp), dimension(kmax+2,kmax+1) :: gint          ! interp. coefficients
    real(dp), dimension(kmax+2,kmax+1) :: gpint         ! interp. coefficients
    real(dp), dimension(kmax)          :: gamma1        ! gamma associated with gint
    real(dp), dimension(kmax)          :: gammap        ! gamma associated with gpint
    real(dp)                           :: hI            ! interpolation step
    real(dp), dimension(3,kmax+1)      :: srpdiff       ! SRP acceleration difference
    real(dp), dimension(3,kmax+1)      :: newsrpdiff    ! SRP acceleration difference new
    real(dp), dimension(3,3)           :: posl2         ! matrix containing three position vectors: pos, lastpos, pos2back
    real(dp), dimension(3,kmax+1)      :: srpdiffsave   ! SRP acceleration saved
    real(dp), dimension(2,3)           :: T             ! T-function for SRP correction
    real(dp), dimension(3)             :: theta         ! angle theta

    integer                            :: betaSRP       ! telling whether object moves from shadow to sun (=1) or from sun to shadow (=-1)
    integer                            :: border        ! telling how many boundaries where crossed since last step (0=one boundary, 1=two boundaries, -1=no crossings)
    integer                            :: q             ! coefficient index
    integer                            :: i,m           ! loop control
    integer                            :: num
    integer                            :: year, month, day, hour, minute, second

    logical                            :: again         ! true if repeating initial step
    logical                            :: isBoundary    ! true if shadow boundary crossing happened

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    logInfo = 'UNK'  ! unknown log info for the moment, gets updated later. Essentially required to know, whether prior to extrapolation an integration took place

    !** count number of subroutine calls
    this%countCallsVarstormcow = this%countCallsVarstormcow + 1

    !** correct for reset flag if called for first time
    if(this%countCallsVarstormcow == 1 .and. reset == 0) then
      call setNeptuneError(E_INTEGRATION_RESET, WARNING)
      reset = 1
      this%is_start_epoch_reset = .false.
    else if (reset == 1) then
        ! Nothing to do when we reset anyway
        this%is_start_epoch_reset = .false.
    else if (this%countCallsVarstormcow > 1 .and. this%is_start_epoch_reset) then
        ! Override integration time
        this%inttime = this%inttime - rqtime
        this%is_start_epoch_reset = .false.
    end if


    ! save the current time
    savecurrtime = currtime

    !========================================================================
    !
    !   1) Initialization
    !
    !------------------------------------------------------------------------
    if (reset == 1) then
      logInfo = 'INI'

      ! 1.1) set EPS, this%releps, this%abseps, initial state error transition matrix
      !
      this%eps    = max(this%eps_rel,this%eps_abs)
      this%releps = this%eps_rel / this%eps
      this%abseps = this%eps_abs / this%eps

      ! get machine epsilon
      macheps = epsilon(1.d0)

      this%k     = 1
      reset = 0

      !
      ! 1.2) get initial value of acceleration and save initial state vector in module
      !
      if(this%countCallsVarstormcow == 1) then
        this%initialState%r = pos
        this%initialState%v = vel
      end if

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
                                    pos,                &
                                    vel,                &
                                    currtime,           &
                                    accel)
      if(hasFailed()) return

      terrs = 0.0d0
      terrd = 0.0d0

      this%lastpos = pos
      this%lastvel = vel

      this%diff(1:3,1) = accel     ! initial global acceleration

      !** check whether solar radiation pressure is activated or not - in order to perform Lundberg correction if requested...
      if(derivatives_model%getPertSwitch(PERT_SRP)) this%flag_srp = .true.
      if(this%correctSRP .and. this%flag_srp) then ! get initial SRP acceleration if correction will be required
        call radiation_model%getSrpAcceleration(    satellite_model,    &
                                                    solarsystem_model,  &
                                                    reduction,          &
                                                    pos,                &
                                                    vel,                &
                                                    this%start_epoch%mjd + currtime/sec_per_day, &
                                                    srpdiff(:,1))
      end if

      wts         = abs(vel)*this%releps + this%abseps
      wtd         = abs(pos)*this%releps + this%abseps

      do i=1,3

        terrs       = terrs + (accel(i)/wts(i))**2
        terrd       = terrd + (vel(i)/wtd(i))**2

      end do

      !
      ! 1.3) get initial step estimate for single and double integration
      !
      h1s = 0.25d0 * sqrt(this%eps / sqrt(terrs))
      h1d = 0.25d0 * sqrt(this%eps / sqrt(terrd))

      !
      ! 1.4) get initial step size
      !
      this%stepsize = min(h1s,h1d,abs(rqtime-currtime))
      this%stepsize = max(this%stepsize,4*macheps * currtime)

      ! set step size in the right direction
      if ( (rqtime - currtime) > epsilon(1.0d0) ) then
        this%intdir = 1
      else
        this%intdir = -1
      endif

      this%stepsize = sign(this%stepsize,dble(this%intdir))
      this%fail     = 0
      again    = .true.
      num      = 0

      do while (again .and. this%fail < 10)

        again     = .false.
        num       = num + 1
        stepsize2 = this%stepsize * this%stepsize
        this%steps(1)  = this%stepsize

        ! do first step with general formulation
        ! predict

        vel(1:3) = this%lastvel(1:3) + this%stepsize*this%diff(1:3,1)
        pos(1:3) = this%lastpos(1:3) + this%stepsize*this%lastvel(1:3) + stepsize2*this%diff(1:3,1)*0.5d0

        currtime = currtime + this%stepsize

        ! evaluate
        call derivatives_model%deriv(                       &
                                        gravity_model,      &
                                        atmosphere_model,   &
                                        manoeuvres_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        tides_model,        &
                                        reduction,          &
                                        pos,                &
                                        vel,                &
                                        currtime,           &
                                        accel)
        if(hasFailed()) return

        terrs = 0.0d0
        terrd = 0.0d0

        newdiff(:,1) = accel(:)
        newdiff(:,2) = newdiff(:,1) - this%diff(:,1)

        if(this%correctSrp .and. this%flag_srp) then
          call radiation_model%getSrpAcceleration(  satellite_model,    &
                                                    solarsystem_model,  &
                                                    reduction,          &
                                                    pos,                &
                                                    vel,                &
                                                    this%start_epoch%mjd + currtime/sec_per_day, &
                                                    newsrpdiff(:,1))
        end if

        ! correct
        vel = vel + this%stepsize  * newdiff(:,2)*0.5d0
        pos = pos + stepsize2 * newdiff(:,2)/6.0d0

        ! get error estimate
        wts = abs(this%lastvel) * this%releps + this%abseps
        wtd = abs(this%lastpos) * this%releps + this%abseps

        do i = 1, 3

          terrs  = terrs + (newdiff(i,this%k+1) / wts(i))**2
          terrd  = terrd + (newdiff(i,this%k+1) / wtd(i))**2

        end do

        terrs = sqrt(terrs)
        terrd = sqrt(terrd)
        this%errd  = abs( stepsize2 * ( 1.d0/3.d0 ) ) * terrd
        this%errs  = abs( this%stepsize * ( 1.d0/2.d0 ) ) * terrs

        if (this%errd > this%eps .or. this%errs > this%eps) then
          ! set everything back and try again with half step
          again = .true.
          this%fail = this%fail + 1
          currtime = currtime - this%stepsize
          this%stepsize = 0.5d0 * this%stepsize

        elseif(this%fail == 0) then

          ! step succeeded, try again with a larger step,
          ! to find maximum initial step
          again = .true.
          currtime = currtime - this%stepsize
          this%stepsize = 2.0d0 * this%stepsize

        endif

      enddo

      ! end of initial step loop
      if (again) then
        call setNeptuneError(E_INTEGRATION_STARTUP, FATAL)
        return
      endif

      ! reset fail count
      this%fail = 0

      ! re-evaluate
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
                                    pos,                &
                                    vel,                &
                                    currtime,           &
                                    accel)
      if(hasFailed()) return

      newdiff(:,1) = accel

      if(this%correctSrp .and. this%flag_srp) then
        call radiation_model%getSrpAcceleration(    satellite_model,    &
                                                    solarsystem_model,  &
                                                    reduction,          &
                                                    pos,                &
                                                    vel,                &
                                                    this%start_epoch%mjd + currtime/sec_per_day, &
                                                    newsrpdiff(:,1))

        newsrpdiff(:,2) = newsrpdiff(:,1) - srpdiff(:,1)
        srpdiff(:,1:2)  = newsrpdiff(:,1:2)
      end if

      newdiff(:,2) = newdiff(:,1) - this%diff(:,1)
      this%diff(:,1:2)  = newdiff(:,1:2)

      ! set integration time and state to current values
      this%inttime = currtime
      this%intpos  = pos
      this%intvel  = vel

      ! increment order, double step size on next step
      this%lk = 1
      this%k  = 2
      this%r  = 2.0d0

      !** log this%stepsize, number of backpoints, local error
      if(this%flog) then
        call mjd2gd(currtime/86400.d0,year,month,day,fracday)
        call dayFraction2HMS(fracday,hour,minute,second)
        write(this%intlog,100) year, month, day, hour, minute,second, currtime/86400.d0, &
                          this%stepsize, this%k, this%errd, this%errs, this%countCallsVarstormcow, logInfo
      end if

      if(isControlled()) then
        call checkOut(csubid)
      end if
      return

    endif ! end of initialization

    !** return integrated state for the very unlikely (but still possible!!) case, that this%inttime == rqtime!
    if(this%inttime == rqtime) then
      if(logInfo == 'UNK') logInfo = 'NST' ! no step
      pos      = this%intpos
      vel      = this%intvel
      currtime = rqtime
    end if

    !========================================================================
    !
    ! 2)  Numerical Integration
    !
    !------------------------------------------------------------------------
    !
    ! 2.1) take a step if requested time is past integration time
    !
    if ( (rqtime > this%inttime .and. this%intdir == 1) .or. (rqtime < this%inttime .and. this%intdir == -1) ) then

      logInfo  = 'EXT'

      pos2back = this%lastpos
      this%lastpos  = this%intpos
      this%lastvel  = this%intvel

      ! change the step size
      this%stepsize  = this%stepsize * this%r

      ! Force a step size to meet e.g. changes in acceleration => maneuvers
      if ((rqtime < (this%inttime + this%stepsize)) .and. force_no_interpolation) then
        this%stepsize = rqtime - this%inttime
      end if

      stepsize2 = this%stepsize * this%stepsize

      if (this%lk == kmax) then

        ! If already using maximum value of k, cycle this value into
        ! steps array. Otherwise, it gets appended to the end.
        savestep =  this%steps(1)

        do i=1,this%k-1
           this%steps(i) =  this%steps(i+1)
        enddo

      endif

      this%steps(this%k) = this%stepsize

      ! compute ratio of current step size to last
      ratio  = this%stepsize /  this%steps(this%k-1)

      invrat = 1.d0/ratio

      ! calculate alpha and this%psi(n+1)
      do i = 1,this%k
        if (i > 1) then
          this%psi(i) = this%psi(i-1) +  this%steps(this%k+1-i)
        else
          this%psi(1) = this%stepsize
        endif
        alpha(i) = this%stepsize / this%psi(i)
      enddo

      ! get this%psi(n-1) (need 0 through k-2)
      psinm1(0) = 0.d0

      do i=1,this%k-2
        psinm1(i) = psinm1(i-1) +  this%steps(this%k-1-i)
      enddo

      ! get this%psi(n)
      this%psin(1) =  this%steps(this%k-1)

      do i=2,this%k-1
        this%psin(i) = this%psin(i-1) +  this%steps(this%k-i)
      enddo

      ! calculate coefficients
      do i = 1,this%k+1
        do q = 1,this%k+3-i
          if (i==1) then
            g(q,1) = 1.d0 / dble(q)
            gp(q,1) = g(q,1) * (-invrat)**q
          elseif (i==2) then
            g(q,2) = 1.d0 / dble(q) / ( dble(q) + 1.d0 )
            gp(q,2) = g(q,2) * (-invrat)**(q+1)
          else
            g(q,i) = g(q,i-1) - alpha(i-1) * g(q+1,i-1)
            gp(q,i) = psinm1(i-3) / this%psi(i-1) * gp(q,i-1) - alpha(i-1) * gp(q+1,i-1)
          endif
        enddo
      enddo

      ! calculate beta
      beta(1) = 1.d0

      do m = 2,this%k
        beta(m) = beta(m-1) * this%psi(m-1) / this%psin(m-1)
      enddo

      ! now integrate
      do i=1,3

        sum2 = 0.d0
        sum1 = 0.d0

        do m = 1,this%k

          ! save differences before multiplying by beta
          diffsave(i,m) = this%diff(i,m)
          if(this%correctSrp .and. this%flag_srp) then
            srpdiffsave(i,m) = srpdiff(i,m)
          end if
          ! calculate phi* = phi*beta
          this%diff(i,m) = beta(m) * this%diff(i,m)
          ! calculate sums for integration
          sum2 = sum2 + ( g(2,m) + ratio * gp(2,m) ) * this%diff(i,m)
          sum1 = sum1 + g(1,m) * this%diff(i,m)

        enddo

        ! predict velocity and position
        pos(i) = (1.d0 + ratio) * this%lastpos(i) - ratio*pos2back(i) + stepsize2 * sum2
        vel(i) = this%lastvel(i) + this%stepsize * sum1

      enddo

      ! update the time
      currtime = this%inttime + this%stepsize

      ! evaluate force model
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
                                    pos,                &
                                    vel,                &
                                    currtime,           &
                                    accel)
      if(hasFailed()) return

      if(this%correctSrp .and. this%flag_srp) then
        call radiation_model%getSrpAcceleration(    satellite_model,    &
                                                    solarsystem_model,  &
                                                    reduction,          &
                                                    pos,                &
                                                    vel,                &
                                                    this%start_epoch%mjd + currtime/sec_per_day, &
                                                    newsrpdiff(:,1))
      end if

      terrs = 0.0d0
      terrd = 0.0d0

      do i=1,3
        ! get new differences
        newdiff(i,1) = accel(i)
        do m=2,this%k+1
          newdiff(i,m) = newdiff(i,m-1) - this%diff(i,m-1)
          if(this%correctSrp .and. this%flag_srp) then
            newsrpdiff(i,m) = newdiff(i,m-1) - srpdiff(i,m-1)
            if(ieee_is_nan(newsrpdiff(i,m))) then
              newsrpdiff(i,m) = 0.d0
            end if
          end if
        enddo

        ! corrector
        pos(i) = pos(i) + stepsize2 * ( g(2,this%k+1) + ratio*gp(2,this%k+1) ) * newdiff(i,this%k+1)
        vel(i) = vel(i) + this%stepsize * g(1,this%k+1) * newdiff(i,this%k+1)

        ! get error estimate
        wts(i) = abs(vel(i)) * this%releps + this%abseps
        wtd(i) = abs(pos(i)) * this%releps + this%abseps
        terrs  = terrs + (newdiff(i,this%k+1) / wts(i))**2
        terrd  = terrd + (newdiff(i,this%k+1) / wtd(i))**2

      enddo

      terrs = sqrt(terrs)
      terrd = sqrt(terrd)
      this%errd  = abs( stepsize2 * ( g(2,this%k+1) - g(2,this%k) + ratio * ( gp(2,this%k+1) - gp(2,this%k) ) ) ) * terrd
      this%errs  = abs( this%stepsize * ( g(1,this%k+1) - g(1,this%k) ) ) * terrs

      if (this%errd > this%eps .or. this%errs > this%eps) then

        ! step failed,
        ! reset everything and use half this step next time
        this%fail = this%fail + 1
        this%r    = 0.5d0

        if (this%lk == kmax) then
          do m = this%k,2,-1
             this%steps(m) =  this%steps(m-1)
          enddo
           this%steps(1) = savestep
        endif

        !pos     = this%lastpos
        !vel     = this%lastvel
        !this%lastpos = pos2back

        !diff(1:3,1:k) = diffsave(1:3,1:k)

        !currtime = savecurrtime
        !delt     = 0.d0

        ! if 3rd failure, start over OR do SRP correction in case of shadow boundary crossing...
        if (this%fail >= 3) then

          if(this%correctSrp .and. this%flag_srp) then

            ! check if there was a shadow boundary crossing...
            !-----------------------------------------------------------

            !** prepare arrays
            posl2(:,1) = pos(:)
            posl2(:,2) = this%lastpos(:)
            posl2(:,3) = pos2back(:)

            !** call checking routine in radiation module
            call radiation_model%checkBoundaryCrossing(       &
                                        solarsystem_model,    &
                                        posl2,                & ! <-- DBL() current + 2 prior position vectors
                                        this%start_epoch%mjd + currtime/sec_per_day, & ! <-- DBL   current time / MJD
                                        this%steps,           & ! <-- DBL() step history
                                        isBoundary,           & ! --> LOG   flag indicating that boundary crossing occurred
                                        T,                    & ! --> DBL() values of T-function
                                        fshadow,              & ! --> DBL() shadow function
                                        theta,                & ! --> DBL() Angle Theta
                                        alpha_e,              & ! --> DBL() Angle Alpha earth
                                        alpha_s               & ! --> DBL() Angle Alpha
                                      )

            !-----------------------------------------------------------

            if(isBoundary) then

              !** perform update of position and velocity vector
              call this%srpCorrection(                                          &
                                  radiation_model,                              &
                                  satellite_model,                              &
                                  solarsystem_model,                            &
                                  reduction,                                    &
                                  posl2,                                        & ! <-- DBL() current + 2 prior position vectors
                                  vel,                                          & ! <-- DBL() current velocity vector
                                  this%start_epoch%mjd*sec_per_day+ currtime,   & ! <-- DBL   current time / sec
                                  this%stepsize,                                & ! <-- DBL   current stepsize / sec
                                  T,                                            & ! <-- DBL() values of T-function
                                  theta,                                        & ! <-- DBL() angles theta
                                  alpha_e,                                      & ! <-- DBL() angles alpha earth
                                  alpha_s,                                      & ! <-- DBL() angles alpha sun
                                  a_srp_wo,                                     & ! --> DBL() acceleration as if satellite was in full sunlight
                                  betaSRP,                                      & ! --> INT   decides which crossing there is: 1 = shadow->sun, -1 = sun->shadow
                                  border                                        & ! --> INT   tells where the satellite is: 0 = only one boundary crossed,
                                                                                  !                                         1 = two boundaries crossed,
                                                                                  !                                        -1 = no crossings
                                )

              if (border == 0) then ! Only one boundary crossed (node in penumbra region), no correction applied yet

                ! The integration is carried out through the penumbra as if the shadow factor has not changed.
                ! This means: the SRP acceleration has 'virtually' not changed within this step.

                ! calculating pos and vel as if shadow function has not changed:

                ! 1) going back to last valid pos and vel ("delete" pos and vel in penumbra)

                pos     = this%lastpos                        ! resetting position to last node
                vel     = this%lastvel                        ! resetting velocity to last node
                this%lastpos = pos2back                       ! resetting past position to pos2back node

                ! predict pos and vel for n+1
                forall(i = 1:3)
                  pos(i) = (1.d0 + ratio) * this%lastpos(i) - ratio*pos2back(i) + stepsize2 * sum2
                  vel(i) = this%lastvel(i) + this%stepsize * sum1
                end forall

                currtime = this%inttime + this%stepsize

                !getting overall acceleration and SRP acceleration
                call derivatives_model%deriv(                       &
                                                gravity_model,      &
                                                atmosphere_model,   &
                                                manoeuvres_model,   &
                                                radiation_model,    &
                                                satellite_model,    &
                                                solarsystem_model,  &
                                                thirdbody_model,    &
                                                tides_model,        &
                                                reduction,          &
                                                pos,                &
                                                vel,                &
                                                currtime,           &
                                                accel)                       ! overall acceleration vector
                if(hasFailed()) return
                call radiation_model%getSrpAcceleration(    satellite_model,    &
                                                            solarsystem_model,  &
                                                            reduction,          &
                                                            pos,                &
                                                            vel,                &
                                                            this%start_epoch%mjd + currtime/sec_per_day, &
                                                            a_srp)  ! current srp-acceleration vector
                if(hasFailed()) return

                ! 2) modify the accelerations and the acceleration differences
                if (betaSRP == 1) then   ! shadow -> penumbra (shadow factor still 0)

                  ! modifying overall acceleration vector [accel] with srp-acceleration vector [a_srp]
                  ! "virtual" acceleration for umbra (f = 0)
                  accel = accel - a_srp

                  ! get new differences
                  do i=1,3

                    ! first column is the modified acceleration itself
                    newdiff(i,1) = accel(i)

                    do m=2,this%k+1
                      ! the following columns are the back differences of the overall accelerations.
                      ! They do not need to be modified, because their real position is still in the shadow +
                      ! the virtual acceleration at the leading node is matched to no SRP (as the states before)
                      newdiff(i,m) = newdiff(i,m-1) - this%diff(i,m-1)

                    end do!m=2,k+1

                  end do!i=1,3

                else ! sun -> penumbra (shadow factor still 1)

                  ! modifying overall acceleration vector [accel] with srp-acceleration vector [a_srp,a_srp_wo]
                  accel = accel - a_srp + a_srp_wo        ! "virtual" accel for full sunlight (f = 1)

                  ! get new differences
                  do i=1,3

                    ! first column is acceleration itself
                    newdiff(i,1) = accel(i)

                    do m=2,this%k+1
                      ! following columns are the back differences of the complete accelerations. The dont need to be modified,
                      ! because their real position is still in the sun + the virtual acceleration at the leading node is
                      ! matched to include SRP (as the states before)
                      newdiff(i,m) = newdiff(i,m-1) - this%diff(i,m-1)

                    enddo!m=2,k+1

                  end do!i=1,3

                  ! accelerations and differences modified successfully!

                end if!if betaSRP = 1/-1

                ! reset errors just in case
                terrs = 0.0d0
                terrd = 0.0d0

                do i = 1,3

                  ! corrector
                  pos(i) = pos(i) + stepsize2 * ( g(2,this%k+1) + ratio*gp(2,this%k+1) ) * newdiff(i,this%k+1)
                  vel(i) = vel(i) + this%stepsize * g(1,this%k+1) * newdiff(i,this%k+1)

                  ! get error estimate
                  !wts(i) = abs(this%lastvel(i)) * this%releps + this%abseps / VB (10-APR-2015): Should be 'abs(vel)' here !!?
                  wts(i) = abs(vel(i)) * this%releps + this%abseps
                  !wtd(i) = abs(this%lastpos(i)) * this%releps + this%abseps / VB (10-APR-2015): Should be 'abs(pos)' here !!?
                  wtd(i) = abs(pos(i)) * this%releps + this%abseps
                  terrd  = terrd + (newdiff(i,this%k+1) / wtd(i))**2

                end do

                ! Error estimation
                terrs = sqrt(terrs)
                terrd = sqrt(terrd)
                this%errd  = abs( stepsize2 * ( g(2,this%k+1) - g(2,this%k) + ratio * ( gp(2,this%k+1) - gp(2,this%k) ) ) ) * terrd
                this%errs  = abs( this%stepsize * ( g(1,this%k+1) - g(1,this%k) ) ) * terrs

              else if (border == 1) then ! Two boundaries crossed

                ! Integration node in sunlight or umbra
                pos     = this%lastpos                        ! resetting position to last node
                vel     = this%lastvel                        ! resetting velocity to last node
                this%lastpos = pos2back                       ! resetting past position to pos2back nodes


                ! 1) Predict pos and vel for n+1
                forall(i = 1:3)
                  pos(i) = (1.d0 + ratio) * this%lastpos(i) - ratio*pos2back(i) + stepsize2 * sum2
                  vel(i) = this%lastvel(i) + this%stepsize * sum1
                end forall

                currtime = this%inttime + this%stepsize

                ! 2) Apply method of modified back differences
                call radiation_model%getSrpAcceleration(    satellite_model,    &
                                                            solarsystem_model,  &
                                                            reduction,          &
                                                            pos,                &
                                                            vel,                &
                                                            this%start_epoch%mjd + currtime/sec_per_day, &
                                                            a_srp)

                newsrpdiff(1:3,1) = a_srp

                ! SRP acceleration differences are "added to" OR "subbracted from" the overall acceleration difference,
                ! because the integrator must not be confronted with a rapid change in acceleration

                ! Method of modified back differences applied to the last k positions
                call derivatives_model%deriv(                       &
                                                gravity_model,      &
                                                atmosphere_model,   &
                                                manoeuvres_model,   &
                                                radiation_model,    &
                                                satellite_model,    &
                                                solarsystem_model,  &
                                                thirdbody_model,    &
                                                tides_model,        &
                                                reduction,          &
                                                pos,                &
                                                vel,                &
                                                currtime,           &
                                                accel)                          ! overall acceleration vector
                if(hasFailed()) return

                newdiff(:,1) = accel(:)

                if (betaSRP == 1) then ! umbra -> sunlight (current f = 1)
                  ! all past k acceleration-differences-vectors will be modified by adding
                  ! the srp-acc-differences vector to match the SRP-level
                  ! on the leading node!
                  do m = 2,this%k+1
                    newdiff(1:3,m) = newdiff(1:3,m) + newsrpdiff(1:3,m)
                  end do

                else ! sunlight -> umbra (current f = 0)
                  ! all past k acceleration-differences-vectors will be modified by
                  ! subtracting the srp-acc-differences vector to match the SRP-level
                  ! on the leading node!
                  do m = 2,this%k+1
                    newdiff(1:3,m) = newdiff(1:3,m) - newsrpdiff(1:3,m)
                  end do

                end if!betaSRP = 1/-1

                !resetting erros just in case
                terrs = 0.0d0
                terrd = 0.0d0

                do i = 1,3

                  ! corrector
                  pos(i) = pos(i) + stepsize2 * ( g(2,this%k+1) + ratio*gp(2,this%k+1) ) * newdiff(i,this%k+1)
                  vel(i) = vel(i) + this%stepsize * g(1,this%k+1) * newdiff(i,this%k+1)

                  ! get error estimate
                  !wts(i) = abs(this%lastvel(i)) * this%releps + this%abseps/ VB (10-APR-2015): Should be 'abs(pos)' here !!?
                  wts(i) = abs(vel(i)) * this%releps + this%abseps
                  !wtd(i) = abs(this%lastpos(i)) * this%releps + this%abseps/ VB (10-APR-2015): Should be 'abs(pos)' here !!?
                  wtd(i) = abs(pos(i)) * this%releps + this%abseps

                  terrs  = terrs + (newdiff(i,this%k+1) / wts(i))**2
                  terrd  = terrd + (newdiff(i,this%k+1) / wtd(i))**2

                end do

                ! Error estimation
                terrs = sqrt(terrs)
                terrd = sqrt(terrd)
                this%errd  = abs( stepsize2 * ( g(2,this%k+1) - g(2,this%k) + ratio * ( gp(2,this%k+1) - gp(2,this%k) ) ) ) * terrd
                this%errs  = abs( this%stepsize * ( g(1,this%k+1) - g(1,this%k) ) ) * terrs

              end if!border=0/1

              if(border == -1) then
                reset = 1 ! step failed but NOT due to SRP discontinuity...
              else
                reset = 0 ! do NOT reset integrator and remain order
                this%fail  = 0 ! reset fail counter
              end if
            else ! no boundary...
              reset = 1
            end if
          else ! no shadow correction
            reset = 1
          end if
        endif

        !** failed, but fail is less than 3....
        if (this%fail < 3) then

          pos     = this%lastpos             ! resetting position
          vel     = this%lastvel             ! resetting velocity
          this%lastpos = pos2back            ! resetting past position

          this%diff(1:3,1:this%k) = diffsave(1:3,1:this%k)          ! resetting acc differences
          if(this%correctSrp .and. this%flag_srp) then
            srpdiff(1:3,1:this%k) = srpdiffsave(1:3,1:this%k)  ! resetting srp acc differences
          end if

          currtime = savecurrtime        ! resetting currtime, because hasnt changed
          delt     = 0.d0                ! do not advance in time_mjd_sec

        end if

        if(isControlled()) then
          call checkOut(csubid)
        end if
        ! exit routine
        return

      endif

      ! step succeded, reset fail count
      this%fail = 0

      if (this%k == kmax) then

        ! no longer in variable-order mode
        ! calculate step for next time
        sigma = 1.0d0

        do i = 2,10
          sigma = (i-1) * alpha(i-1) * sigma
        enddo

        erks = abs( this%stepsize * gammastars(this%k) * sigma ) * terrs
        erkd = abs( stepsize2 * gammastard(this%k) * sigma ) * terrd

        ! compute this%r for single and double integration
        rs = (this%eps/(2.d0*erks)) ** (1.d0/dble(this%k+1))
        rd = (this%eps/(2.d0*erkd)) ** (1.d0/dble(this%k+2))

        ! use smallest r
        this%r = min(rs,rd)

        ! bound r
        if (this%r > 2.d0)  this%r = 2.d0
        if (this%r < 0.5d0) this%r = 0.5d0

        ! set last used value of k
        this%lk = this%k

      else

        ! variable-order start-up mode
        ! Evaluate again, increment k and double step
        ! evaluate force model
        call derivatives_model%deriv(                       &
                                        gravity_model,      &
                                        atmosphere_model,   &
                                        manoeuvres_model,   &
                                        radiation_model,    &
                                        satellite_model,    &
                                        solarsystem_model,  &
                                        thirdbody_model,    &
                                        tides_model,        &
                                        reduction,          &
                                        pos,                &
                                        vel,                &
                                        currtime,           &
                                        accel)
        if(hasFailed()) return

        ! get new differences
        do i=1,3
          newdiff(i,1) = accel(i)
          do m=2,this%k+1
            newdiff(i,m) = newdiff(i,m-1) - this%diff(i,m-1)
          enddo
        enddo

        ! set last used value of k, in case an interpolation
        ! is needed now
        this%lk = this%k
        this%k = this%k+1
        this%r = 2.0d0

      endif

      ! Get ready for next step: set integration time and state to
      ! current values, and copy latest differences from newdiff
      ! to diff.

      this%inttime = currtime

      this%diff(1:3,1:this%lk+1) = newdiff(1:3,1:this%lk+1)
      this%intpos = pos
      this%intvel = vel

    endif

    !=======================================================================================
    !
    ! End of integration step, check to see if INTERPOLATION is needed
    !
    !------------------------------------------------------------------------

    if ( (this%inttime > rqtime .and. this%intdir == 1) .or. (this%inttime < rqtime .and. this%intdir == -1) ) then

      if(logInfo == 'UNK') logInfo = 'INT'

      ! integrated too far, interpolate to request time
      ! set interpolation step
      hI     = rqtime - this%inttime
      ratio  = hI/ this%steps(this%lk)
      invrat = 1.0d0 / ratio

      ! compute gamma values
      gamma1(1) = hI / this%psi(1)

      do i=2,this%lk
        gamma1(i) = (hI + this%psi(i-1)) / this%psi(i)
      enddo

      gammap(1) = -1.0d0
      gammap(2) = 0.0d0

      do i=3,this%lk
        gammap(i) = this%psin(i-2) / this%psi(i)
      enddo

      ! calculate coefficients
      do i=1,this%lk+1
        do q=1,this%lk+3-i

          if (i == 1) then
            gint(q,i) = 1.0d0 / q
            gpint(q,i) = 1.0d0 / q * (-invrat)**q
          else
            gint(q,i) = gamma1(i-1) * gint(q,i-1) - hI/this%psi(i-1) * gint(q+1,i-1)
            gpint(q,i) = gammap(i-1) * gpint(q,i-1) - hI/this%psi(i-1) * gpint(q+1,i-1)
          endif

        enddo
      enddo

      ! compute interpolation values
      do i=1,3

        sum1 = 0.0d0
        sum2 = 0.0d0

        do m=1,this%lk+1
          sum1 = sum1 + gint(1,m) * this%diff(i,m)
          sum2 = sum2 + (gint(2,m) + ratio*gpint(2,m)) * this%diff(i,m)
        enddo

        vel(i) = this%intvel(i) + hI * sum1
        pos(i) = (1.0d0 + ratio) * this%intpos(i) - ratio * this%lastpos(i) + hI**2 * sum2

      enddo

      ! set current time to the request time
      currtime = rqtime

    !else if ((rqtime - this%inttime) < epsilon(1.0)) then
       ! Edited by CHK because the integrator somehow got caught in a loop where
       !  this%inttime was equal to rqtime but currtime was not moved forward for some reason ?!
    !  currtime = rqtime
    endif ! end of interpolation

    !====================================================================================
    !
    ! 4) Finish
    !
    !-------------------------------------------------------------------------------------

    ! set the change in time
    delt = currtime - savecurrtime

    !** log stepsize, number of backpoints, local error
    if (this%flog) then
      if(logInfo=='UNK') then
        write(*,*) this%inttime, rqtime, this%intdir
        stop
      end if
      call mjd2gd(currtime/86400.d0,year,month,day,fracday)
      call dayFraction2HMS(fracday,hour,minute,second)
      write(this%intlog,100) year, month, day, hour, minute, second, currtime/86400.d0, &
                        this%stepsize, this%k, this%errd, this%errs, this%countCallsVarstormcow, logInfo
    end if

    !** done
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

100 format(4x,i4.4,2("-",i2.2),1x,2(i2.2,":"),i2.2,7x,f17.10,11x, f12.5,20x,i1,6x,2(14x,e13.6e2),11x,i8,3x,a3)

  end subroutine varstormcow



!========================================================================
!
!> @anchor      fixedgaussjackson
!!
!> @brief       Gauss-Jackson integrator (fixed step)
!> @brief      Andre Horstmann
!!
!> @date        <ul>
!!                <li> 06.2015 (initial design)
!!                <li> 24.01.2015 (changed the way the counter works - now via module variable, cf. varstormcow)</li>
!!              </ul>
!!
!> @detail      The subroutine fixedgaussjackson performs one step of the fixed-step
!!              Gauss-Jackson algorithm. If the integration step takes the
!!              current time beyond the request time, the routine will
!!              interpolate to the request time. The routine changes the value
!!              of current time, position, and velocity during the call. The
!!              routine sets the value of delt, which is the amount the
!!              current time changed.
!!              Call with reset = 1 to start the integrator. The routine will
!!              normally return with reset = 0. If the routine returns
!!              reset = 1, the method must restart because a step has failed
!!              after 3 tries.
!!              The algorithms used in this code are explained in:
!!
!------------------------------------------------------------------------
  subroutine fixedgaussjackson(     &
                          rqtime,   &                                           ! <--  DBL    requested time
                          currtime, &                                           ! <--> DBL    current time
                          pos,      &                                           ! <--> DBL()  position vector
                          vel,      &                                           ! <--> DBL()  velocity vector
                          reset,    &                                           ! <--> INT    1 = restarting
                                                                                !             0 = normal return
                          delt      &                                           ! -->  DBL    change in current time
                        )
    !** interface
    !------------------------------------------------------------------------
    real(dp),                 intent(in)    :: rqtime   ! requested time
    real(dp),                 intent(inout) :: currtime ! current time
    real(dp), dimension(3),   intent(inout) :: pos      ! position vector
    real(dp), dimension(3),   intent(inout) :: vel      ! velocity vector
    integer,                  intent(inout) :: reset    ! 1 if restarting
    real(dp),                 intent(out)   :: delt     ! change in currtime
    !------------------------------------------------------------------------

    ! during call
    ! note: the units of rqtime, currtime, pos, vel, and accel must
    ! all be compatible

!   integer, parameter :: kmax   = 8 ! maximum value of k

!   !------------------------------------------------------------------------
!   ! local variables

    character(len=*), parameter :: csubid = "fixedgaussjackson"
!   character(len=3) :: logInfo
!   real(dp), dimension(3)             :: a_srp         ! current SRP acceleration
!   real(dp), dimension(3)             :: a_srp_wo      ! SRP acceleration without shadow! (f=1)
!   real(dp)                           :: stepsize2     ! step size squared
!   real(dp), dimension(3)             :: accel         ! acceleration vector
!   real(dp), dimension(3)             :: alpha_e       ! angle satellite -> earth
!   real(dp), dimension(3)             :: alpha_s       ! angle satellite -> sun
!   real(dp), dimension(3,kmax+1), save             :: diff           ! modified divided differences
!   real(dp), dimension(3,kmax)                     :: diffsave       ! diff from last step
!   real(dp), save                     :: r             ! amount to change the step by on next step
!   real(dp), dimension(kmax), save    :: steps         ! step size history
!   real(dp), dimension(kmax), save    :: psi           ! sums of steps
!   real(dp), dimension(kmax)          :: alpha         ! ratio of current step to psi
!   real(dp), dimension(0:kmax-2)      :: psinm1        ! sums of steps,
!                                                 ! starting 2 steps back
!   real(dp), dimension(kmax-1), save  :: psin          ! sums of steps,
!                                                 ! starting 1 step back
!   real(dp), dimension(kmax+2,kmax+1) :: g             ! integration coefficients
!   real(dp), dimension(kmax+2,kmax+1) :: gp            ! integration coefficients
!   real(dp), dimension(kmax)          :: beta          ! ratio of phi’s
!   real(dp), dimension(3)             :: pos2back      ! position from 2 steps ago
!   real(dp), dimension(3), save                        :: lastpos       ! position from last step
!   real(dp), dimension(3), save              :: lastvel       ! velocity from last step
!   real(dp), dimension(3,kmax+1)             :: newdiff       ! differences after evaluation
!   real(dp)                           :: savecurrtime  ! last value of current time
!   real(dp), save                     :: errd          ! error estimate on position
!   real(dp), save                     :: errs          ! error estimate on velocity
!   real(dp)                           :: erkd          ! position error estimate
!   real(dp)                           :: erks          ! velocity error estimate
!   real(dp)                           :: terrd         ! temp value for errd
!   real(dp)                           :: terrs         ! temp value for errs
!   real(dp)                           :: rs            ! r value for single integration
!   real(dp)                           :: rd            ! r value for double integration
!   real(dp)                           :: sigma         ! parameter used in error estimate
!   real(dp)                           :: savestep      ! saved step size of
!                                                 ! furthest backpoint
!   real(dp)                           :: ratio         ! ratio of latest step sizes
!   real(dp)                           :: invrat        ! inverse of ratio
!   real(dp)                           :: sum1          ! integration sum for velocity
!   real(dp)                           :: sum2          ! integration sum for position
!   real(dp), dimension(kmax)          :: gammastars = (/-0.5d0,                 -0.08333333333333333d0,      &
!                                                        -0.04166666666666667d0, -0.02638888888888891d0,      &
!                                                        -0.01874999999999999d0, -0.01426917989417986d0,      &
!                                                        -0.01136739417989419d0, -0.009356536596119958d0/)
!                                                       ! difference of single fixed-step coefficients
!   real(dp), dimension(kmax)          :: gammastard = (/-1.d0, 0.08333333333333333d0, 0.0d0, &
!                                                        -0.004166666666666666d0, -0.004166666666666666d0,    &
!                                                        -0.003654100529100521d0, -0.003141534391534404d0,    &
!                                                        -0.002708608906525564d0/)
!                                                       ! difference of double fixed-step coefficients
!   real(dp), save                     :: eps           ! tolerance
!   real(dp), save                     :: releps        ! relative error / eps
!   real(dp), save                     :: abseps        ! absolute error / eps
!   real(dp)                           :: macheps       ! machine epsilon
!   real(dp), dimension(3)             :: wts           ! error weighting for velocity
!   real(dp), dimension(3)             :: wtd           ! error weighting for position
!   real(dp)                           :: h1s           ! initial step for single integration
!   real(dp)                           :: h1d           ! initial step for double integration
!   real(dp), dimension(3), save             :: intpos        ! last position integrated to
!   real(dp), dimension(3), save             :: intvel        ! last velocity integrated to
!   real(dp), save                           :: inttime       ! last time integrated to
!                                                       ! inttime is the same as currtime
!                                                       ! unless the last step had an
!                                                       ! interpolation
!   real(dp)       :: fracday
!   real(dp), dimension(3)             :: fshadow       ! shadow factors [0:1] for three subsequent points
!   real(dp), dimension(kmax+2,kmax+1) :: gint          ! interp. coefficients
!   real(dp), dimension(kmax+2,kmax+1) :: gpint         ! interp. coefficients
!   real(dp), dimension(kmax)          :: gamma1        ! gamma associated with gint
!   real(dp), dimension(kmax)          :: gammap        ! gamma associated with gpint
!   real(dp)                           :: hI            ! interpolation step
!   real(dp), dimension(3,kmax+1)      :: srpdiff       ! SRP acceleration difference
!   real(dp), dimension(3,kmax+1)      :: newsrpdiff    ! SRP acceleration difference new
!   real(dp), dimension(3,3)           :: posl2         ! matrix containing three position vectors: pos, lastpos, pos2back
!   real(dp), dimension(3,kmax+1)      :: srpdiffsave   ! SRP acceleration saved
!   real(dp), save                     :: stepsize      ! step size
!   real(dp), dimension(2,3)           :: T             ! T-function for SRP correction
!   real(dp), dimension(3)             :: theta         ! angle theta

!   integer                          :: betaSRP       ! telling whether object moves from shadow to sun (=1) or from sun to shadow (=-1)
!   integer                          :: border        ! telling how many boundaries where crossed since last step (0=one boundary, 1=two boundaries, -1=no crossings)
!   integer, save                    :: intdir        ! integration direction
!                                                     !  1 = forward
!                                                     ! -1 = backward
!   integer                          :: q             ! coefficient index
!   integer                          :: i,m           ! loop control
!   integer, save                    :: k             ! number of backpoints
!   integer, save                    :: lk            ! k in last integration step for interpolator
!   integer, save                    :: fail          ! number of consecutive failures
!   integer                          :: num
!   integer                          :: year, month, day, hour, minute, second
!
!   logical                          :: again         ! true if repeating initial step
!   logical                          :: isBoundary    ! true if shadow boundary crossing happened
!   logical, save                    :: flag_srp = .false.  ! becomes true as soon as SRP perturbation has been activated and varstormcow is initialized

    !------------------------------------------------------------------------

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

   ! write(*,*) "entered 'fixedgaussjackson' subroutine'"

    !------------------------------------------------------------------------
    !
    !  1. Initialising (8th order)
    !
    !------------------------------------------------------------------------

    ! 1.1 Evaluate initial acceleration and propagate one step. Evaluate again and define polynomial

    delt = 0.0d0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkOut(csubid)
    end if

  end subroutine fixedgaussjackson

!========================================================================
!
!> @anchor      srpCorrection
!!
!> @brief       Shadow boundary correction algorithm for failed integration step
!> @brief      Andre Horstmann, Vitali Braun
!!
!> @date        <ul>
!!                <li> 10.2013 (initial design)  </li>
!!                <li> 11.2013 (optimization)    </li>
!!                <li> 14.03.2014 (optimization) </li>
!!              </ul>
!!
!> @detail      The Lundberg algorithm in this subroutine performs a necessary
!!              correction to a state vector which is exposed to
!!              integration errors due to the abrupt change in acceleration caused by
!!              the solar radiation pressure while crossing shadow boundaries.
!!              After not beeing able to integrate the equation
!!              of motion with the third bisection of the stepsize, the step is taken
!!              nonetheless, maintaining the present stepsize, but with the application
!!              of a state correction at the leading node, which results in a
!!              mitigation of integration errors. Additionally, the required
!!              computation time is reduced, because the integrator does not need to
!!              re-initialise.
!!
!-------------------------------------------------------------------------------------
  subroutine srpCorrection(                     &
                            this,               &
                            radiation_model,    &
                            satellite_model,    &
                            solarsystem_model,  &
                            reduction,          &
                            pos,                &
                            vel,                &
                            currtime,           &
                            stepsize,           &
                            T,                  & ! steps, T, fshadow, &
                            theta,              &
                            alpha_e,            &
                            alpha_s,            &
                            a_srp_wo,           &
                            beta,               &
                            border              &
                            )
  !dpos, dvel, steps_,stepsize_,beta_,a_srp_wo_,border_)

    ! -- initialization -----------------------------------------------------
    class(Numint_class),intent(inout)       :: this
    type(Radiation_class),intent(inout)     :: radiation_model
    type(Satellite_class),intent(inout)     :: satellite_model
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)     :: reduction
    real(dp),dimension(3,3),  intent(inout) :: pos                              ! current position
    real(dp),                 intent(in)    :: currtime                         ! current time
    real(dp),dimension(3),    intent(inout) :: vel                              ! current velocity
    !real(dp),dimension(8),    intent(in)    :: steps                           ! stepsize history
    real(dp),                 intent(in)    :: stepsize                         ! current stepsize
    real(dp), dimension(2,3), intent(in)    :: T                                ! T-function
    !real(dp), dimension(3),   intent(in)    :: fshadow                         ! shadow function
    real(dp), dimension(3),   intent(in)    :: theta                            ! Angle Theta
    real(dp), dimension(3),   intent(in)    :: alpha_e                          ! Angle Alpha earth
    real(dp), dimension(3),   intent(in)    :: alpha_s                          ! Angle Alpha
    real(dp), dimension(3),   intent(out)   :: a_srp_wo                         ! acc w/o shadow
    integer,                  intent(out)   :: beta                             ! deciding factor
    integer,                  intent(out)   :: border                           ! current node in penumbra = 0, in sun =-1, or umbra = 1
    ! -----------------------------------------------------------------------
    ! local variables

    character(len=*), parameter :: csubid = "srpCorrection"

    real(dp), dimension(3)   :: dpos           ! position correction
    real(dp), dimension(3)   :: dvel           ! velocity correction

    real(dp)                 :: lp_T1          ! angle differences
    real(dp)                 :: lp2_T1         ! angle differences
    real(dp)                 :: t_c1           ! first crossing time
    real(dp)                 :: t_c2           ! second crossing time
    real(dp)                 :: delta_t1       ! time correction 1
    real(dp)                 :: delta_t2       ! time correction 2

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !===================================================
    !
    ! Compute supposed crossing times
    !
    !-------------------------------------

    ! a) sunlight/penumbra boundary (factor -1)
    !write (*,*) T

    lp_T1  = T(1,1) - T(1,2)
    lp2_T1 = T(1,1) - 2.d0*T(1,2) + T(1,3)

    !write (*,*) lp_T1, lp2_T1

    ! determine correction direction:
    if (lp_T1 > 0.d0) then
      beta = 1        ! sat going: shadow -> sun
    else
      beta = -1       ! sat going: sun -> shadow
    end if

    ! time_1 difference calculation
    delta_t1 = stepsize/(2.d0*lp2_T1)*(2.d0*lp_T1 + lp2_T1 - dble(beta)*((2.d0*lp_T1 + lp2_T1)**2.d0 - 8.d0*lp2_T1*T(1,1))**0.5d0)
    t_c1     = currtime - delta_t1    ! crossing time 1: sunlight / penumbra


    ! b) penumbra/umbra boundary (factor +1)

    lp_T1  = T(2,1) - T(2,2)
    lp2_T1 = T(2,1) - 2.d0*T(2,2) + T(2,3)

    ! determine correction direction
    if (lp_T1 > 0.d0) then
      beta = 1         ! sat going: shadow -> sun
    else
      beta = -1        ! sat going: sun -> shadow
    end if

    ! time_2 difference calculation
    delta_t2 = stepsize/(2.d0*lp2_T1)*(2.d0*lp_T1 + lp2_T1 - beta*((2.d0*lp_T1 + lp2_T1)**2.d0 - 8.d0*lp2_T1*T(2,1))**0.5d0)
    t_c2     = currtime - delta_t2    ! crossing time 2: penumbra / umbra

    !========================================================================
    !
    ! Check crossing times
    !
    !------------------------------------------------------------------------

    if (ieee_is_nan(t_c1) .or. ieee_is_nan(t_c2)) then ! only one crossing time found

      ! check angles:
      if ((theta(1) < alpha_e(1) + alpha_s(1)) .and. &
          (theta(1) > alpha_e(1) - alpha_s(1))) then      ! penumbra

        border = 0

      else

        call setNeptuneError(E_SRP_CORRECTION, FATAL)
        return

      end if

    else if(ieee_is_nan(t_c1) .and. ieee_is_nan(t_c2)) then ! no crossing time found - thus integrator fail is not due to SRP...

      border = -1

    else ! two crossing times found - correct state vector

      ! getting maximum SRP-acceleration (without shadow -> f = 1): a_srp_wo
      a_srp_wo = radiation_model%getSrpWithoutShadow(   satellite_model,    &
                                                        solarsystem_model,  &
                                                        reduction,          &
                                                        pos(1,:),           &
                                                        vel,                &
                                                        this%start_epoch%mjd + currtime/sec_per_day)

      ! calculating correction
      dpos = -beta*(0.5d0*(currtime - t_c2)**2.d0 + ((t_c2 - t_c1)**2)/6.d0)*a_srp_wo
      dvel = -beta*((currtime - t_c2) + 0.5d0*(t_c2 - t_c1))*a_srp_wo

      ! applying calculated correction
      pos(1,:) = pos(1,:) + dpos
      vel      = vel      + dvel

      ! marker for sunlight or umbra
      border = 1

    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine srpCorrection

end module numint
