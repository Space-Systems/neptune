!===============================================================================================
!
!> @anchor      derivatives
!!
!> @brief       Force model evaluation for routines for state vector and covariance propagation
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull
!!
!> @date        <ul>
!!                <li>VB:  11.11.2013 (architecture change to module)</li>
!!                <li>VB:  10.06.2014 (added SRP variational equations for covariance matrix) </li>
!!                <li>VB:  06.10.2015 (added support for additional third bodies: all planets of the solar system) </li>
!!                <li>CHK: 13.11.2015 (updated to use libslam) </li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>VB:  14.07.2017 (added function to set start epoch)</li>
!!                <li>CHK: 02.01.2018 (updated to use atmosphere class) </li>
!!                <li>CHK: 03.01.2018 (created derivatives class) </li>
!!              </ul>
!!
!> @details     This module contains the functions and parameters required for the evaluation of
!!              the force function, which provides the total acceleration in the inertial frame (GCRF)
!!              for the numerical integration of the state vector (deriv), as well as the derivative
!!              of the state transition matrix (deriv_cov).
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      derivatives
!!
!!------------------------------------------------------------------------------------------------
module derivatives

  use slam_astro,             only: getEarthGravity
  use slam_astro_conversions, only: getGeodeticLatLon
  use slam_io,                only: C_ON, C_OFF
  use slam_time,              only: time_t
  use atmosphere,             only: Atmosphere_class
  use neptune_error_handling, only: E_MIN_ALTITUDE, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, &
                                    E_OUT_OF_BOUNDS, checkIn, checkOut
  use gravity,                only: Gravity_class
  use slam_types,             only: dp
  use neptuneParameters,      only: C_SOLID_TIDES, C_OCEAN_TIDES, C_SRP, C_ATMOSPHERE, C_MANEUVERS, C_ALBEDO, C_GEOPOTENTIAL, &
                                    C_MERCURY, C_VENUS, C_MARS, C_JUPITER, C_SATURN, C_URANUS, C_NEPTUNE, C_MOON, C_SUN, C_WIND
  use maneuvers,              only: Manoeuvres_class
  use slam_math,              only: mag
  use radiation,              only: Radiation_class
  use satellite,              only: Satellite_class
  use slam_reduction_class,   only: Reduction_type
  use solarsystem,            only: Solarsystem_class, ID_SUN, ID_MOON, ID_JUPITER, ID_VENUS, ID_MARS, ID_MERCURY, ID_SATURN, ID_URANUS, ID_NEPTUNE
  use tides,                  only: Tides_class, SOLID_TIDES, OCEAN_TIDES
  use thirdbody,              only: ThirdBody_class

  implicit none

  private

    !** parameters
    !-----------------------------------------------------------
    integer, parameter :: LEN_PERTURBATION_NAME = 32  !< length of the name of a perturbation
    integer, parameter :: NUM_PERTURBATIONS     = 17  !< number of perturbations

    ! NOTE: the integer values below are referenced as array indices. Make sure they start with 1 and there is no gap in between.
    ! NOTE: if changes occur to the values below, please also check in derivatives.f90, whether the pertSwitch array setup is still valid.
    integer, parameter, public :: PERT_GRAVITY     = 1   !< geopotential perturbations
    integer, parameter, public :: PERT_DRAG        = 2   !< drag perturbations
    integer, parameter, public :: PERT_SUN         = 3   !< sun gravity perturbations
    integer, parameter, public :: PERT_MOON        = 4   !< moon gravity perturbations
    integer, parameter, public :: PERT_SRP         = 5   !< solar radiation pressure
    integer, parameter, public :: PERT_ALBEDO      = 6   !< albedo
    integer, parameter, public :: PERT_SOLID_TIDES = 7   !< solid tides
    integer, parameter, public :: PERT_OCEAN_TIDES = 8   !< ocean tides
    integer, parameter, public :: PERT_WIND        = 9   !< horizontal wind
    integer, parameter, public :: PERT_MANEUVERS   = 10  !< orbital maneuvers
    integer, parameter, public :: PERT_MERCURY     = 11  !< mercury gravity perturbations
    integer, parameter, public :: PERT_VENUS       = 12  !< venus gravity perturbations
    integer, parameter, public :: PERT_MARS        = 13  !< mars gravity perturbations
    integer, parameter, public :: PERT_JUPITER     = 14  !< jupiter gravity perturbations
    integer, parameter, public :: PERT_SATURN      = 15  !< saturn gravity perturbations
    integer, parameter, public :: PERT_URANUS      = 16  !< uranus gravity perturbations
    integer, parameter, public :: PERT_NEPTUNE     = 17  !< neptune gravity perturbations


    type, public :: perturbation_t
        integer                              :: id                  !< perturbation ID
        logical                              :: switch              !< if .true., perturbation is considered
        character(len=LEN_PERTURBATION_NAME) :: name                !< perturbation name
    end type

    ! The atmospheric model class (type defintion)
    type, public :: Derivatives_class

      !** flags
      !-----------------------------------------------------------
      logical :: initialisedPerturbations  !< initialisation status of this module
      logical :: is_start_epoch_set        !< to check if start epoch has been set

      type(perturbation_t), dimension(:), allocatable :: neptunePerturbations   !< containing all information about the perturbations

      real(dp) :: last_epoch        !< last epoch, e.g. for re-entries / MJD
      real(dp) :: last_altitude     !< last altitude, e.g. for re-entries / km
      real(dp) :: last_latitude    !< last latitude, e.g. for re-entries / rad
      real(dp) :: last_longitude     !< last longitude, e.g. for re-entries / rad

      type(time_t) :: start_epoch   !< epoch of initial state vector at propagation start

    contains

      !** public methods
      !----------------------------------------------------------
      procedure :: deriv
      procedure :: deriv_cov

      !** setter
      procedure :: setPertSwitch

      !** getter
      procedure :: getLastPositionECEF
      procedure :: getPertSwitch
      procedure :: get_neptune_perturbation
      procedure :: get_neptune_perturbation_number
      procedure :: set_start_epoch

      procedure :: initPerturbations
      procedure :: destroy

    end type Derivatives_class

    ! Constructor
    interface Derivatives_class
        module procedure constructor
    end interface Derivatives_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !>  @author     Christopher Kebschull
    !>  @date       <ul>
    !!                  <li>ChK: 03.01.2018 (initial implementation)</li>
    !!              </ul>
    !>  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Derivatives_class) function constructor()
        constructor%initialisedPerturbations = .false.                          !< initialisation status of this module
        constructor%is_start_epoch_set = .false.                                !< to check if start epoch has been set

        constructor%last_epoch     = -1.d0                                      !< last epoch, e.g. for re-entries / MJD
        constructor%last_altitude  = -1.d0                                      !< last altitude, e.g. for re-entries / km
        constructor%last_latitude  = -1.d0                                      !< last latitude, e.g. for re-entries / rad
        constructor%last_longitude = -1.d0                                      !< last longitude, e.g. for re-entries / rad

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
        class(Derivatives_class)    :: this

        if (allocated(this%neptunePerturbations)) deallocate(this%neptunePerturbations)

    end subroutine destroy

  !=====================================================================
  !
  !> @brief     Set the start epoch which is added to the offset during propagation
  !!
  !> @author    Vitali Braun
  !!
  !> @param[in] epoch   Start epoch
  !!
  !> @date      <ul>
  !!              <li>2017-07-14 (initial design) </li>
  !!            </ul>
  !!
  !> @anchor    set_start_epoch
  !!
  !-----------------------------------------------------------------
    subroutine set_start_epoch(this,epoch)

        implicit none

        class(Derivatives_class)    :: this
        type(time_t), intent(in)    :: epoch

        this%start_epoch = epoch
        this%is_start_epoch_set = .true.

        return
    end subroutine set_start_epoch

  !=====================================================================
  !
  !> @brief     Return the number of perturbations modelled by NEPTUNE
  !!
  !> @author    Vitali Braun
  !!
  !> @return    Number of functions
  !!
  !> @date      <ul>
  !!              <li>13.03.2017 (initial design) </li>
  !!            </ul>
  !!
  !! @anchor    get_neptune_perturbation_number
  !!
  !-----------------------------------------------------------------
    pure integer function get_neptune_perturbation_number(this) result(res)

        implicit none

        class(Derivatives_class),intent(in) :: this
        res = NUM_PERTURBATIONS

        return
    end function get_neptune_perturbation_number

  !=====================================================================
  !
  !> @brief     Return the perturbation for a given ID
  !!
  !> @author    Vitali Braun
  !!
  !> @param[in] gravity_model   the gravity model
  !> @param[in] id              ID of perturbation
  !!
  !> @return    the associated perturbation as perturbation_t
  !!
  !> @date      <ul>
  !!              <li>13.03.2017 (initial design) </li>
  !!            </ul>
  !!
  !> @anchor   get_neptune_perturbation
  !!
  !------------------------------------
    type(perturbation_t) function get_neptune_perturbation(this,gravity_model,id) result(res)

         implicit none

        class(Derivatives_class)            :: this
        type(Gravity_class),intent(inout)   :: gravity_model
        integer, intent(in)                 :: id

        character(len=*), parameter :: csubid = "get_neptune_perturbation"

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        ! check initialisation status
        if(.not. this%initialisedPerturbations) then
            call this%initPerturbations(gravity_model)
        end if

        if(id > 0 .and. id <= NUM_PERTURBATIONS) then
            res = this%neptunePerturbations(id)
            if(isControlled()) then
                call checkOut(csubid)
            end if
        else
            call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
        end if
        return

    end function get_neptune_perturbation

  !==========================================================
  !
  !> @brief     Initialise perturbations module
  !!
  !> @author    Vitali Braun
  !!
  !> @param[in] gravity_model   the gravity model
  !!
  !> @date      <ul>
  !!              <li>13.03.2017 (initial design) </li>
  !!            </ul>
  !!
  !------------------------------------
    subroutine initPerturbations(this,gravity_model)

        implicit none

        class(Derivatives_class)            :: this
        type(Gravity_class),intent(inout)   :: gravity_model
        character(len=*), parameter         :: csubid = "initPerturbations"

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(csubid)
        end if

        ! allocate perturbations array
        ! NOTE: you have to change the parameter NUM_PERTURBATIONS, if perturbations are added or removed!
        if(allocated(this%neptunePerturbations)) deallocate(this%neptunePerturbations)
        allocate(this%neptunePerturbations(NUM_PERTURBATIONS))
        this%neptunePerturbations(1:NUM_PERTURBATIONS)%switch = .false.

        this%neptunePerturbations(PERT_GRAVITY)%name   = C_GEOPOTENTIAL
        this%neptunePerturbations(PERT_GRAVITY)%id     = PERT_GRAVITY
        this%neptunePerturbations(PERT_GRAVITY)%switch = .true.
        call gravity_model%setGeoDegree(6)    ! the default value for the geopotential

        this%neptunePerturbations(PERT_DRAG)%name   = C_ATMOSPHERE
        this%neptunePerturbations(PERT_DRAG)%id     = PERT_DRAG
        this%neptunePerturbations(PERT_DRAG)%switch = .true.

        this%neptunePerturbations(PERT_SUN)%name   = C_SUN
        this%neptunePerturbations(PERT_SUN)%id     = PERT_SUN
        this%neptunePerturbations(PERT_SUN)%switch = .true.

        this%neptunePerturbations(PERT_MOON)%name   = C_MOON
        this%neptunePerturbations(PERT_MOON)%id     = PERT_MOON
        this%neptunePerturbations(PERT_MOON)%switch = .true.

        this%neptunePerturbations(PERT_SRP)%name   = C_SRP
        this%neptunePerturbations(PERT_SRP)%id     = PERT_SRP
        this%neptunePerturbations(PERT_SRP)%switch = .true.

        this%neptunePerturbations(PERT_ALBEDO)%name   = C_ALBEDO
        this%neptunePerturbations(PERT_ALBEDO)%id     = PERT_ALBEDO
        this%neptunePerturbations(PERT_ALBEDO)%switch = .true.

        this%neptunePerturbations(PERT_SOLID_TIDES)%name   = C_SOLID_TIDES
        this%neptunePerturbations(PERT_SOLID_TIDES)%id     = PERT_SOLID_TIDES
        this%neptunePerturbations(PERT_SOLID_TIDES)%switch = .true.

        this%neptunePerturbations(PERT_OCEAN_TIDES)%name   = C_OCEAN_TIDES
        this%neptunePerturbations(PERT_OCEAN_TIDES)%id     = PERT_OCEAN_TIDES
        this%neptunePerturbations(PERT_OCEAN_TIDES)%switch = .false.

        this%neptunePerturbations(PERT_WIND)%name   = C_WIND
        this%neptunePerturbations(PERT_WIND)%id     = PERT_WIND
        this%neptunePerturbations(PERT_WIND)%switch = .true.

        this%neptunePerturbations(PERT_MANEUVERS)%name   = C_MANEUVERS
        this%neptunePerturbations(PERT_MANEUVERS)%id     = PERT_MANEUVERS
        this%neptunePerturbations(PERT_MANEUVERS)%switch = .false.

        this%neptunePerturbations(PERT_MERCURY)%name   = C_MERCURY
        this%neptunePerturbations(PERT_MERCURY)%id     = PERT_MERCURY
        this%neptunePerturbations(PERT_MERCURY)%switch = .false.

        this%neptunePerturbations(PERT_VENUS)%name   = C_VENUS
        this%neptunePerturbations(PERT_VENUS)%id     = PERT_VENUS
        this%neptunePerturbations(PERT_VENUS)%switch = .false.

        this%neptunePerturbations(PERT_MARS)%name   = C_MARS
        this%neptunePerturbations(PERT_MARS)%id     = PERT_MARS
        this%neptunePerturbations(PERT_MARS)%switch = .false.

        this%neptunePerturbations(PERT_JUPITER)%name   = C_JUPITER
        this%neptunePerturbations(PERT_JUPITER)%id     = PERT_JUPITER
        this%neptunePerturbations(PERT_JUPITER)%switch = .false.

        this%neptunePerturbations(PERT_SATURN)%name   = C_SATURN
        this%neptunePerturbations(PERT_SATURN)%id     = PERT_SATURN
        this%neptunePerturbations(PERT_SATURN)%switch = .false.

        this%neptunePerturbations(PERT_URANUS)%name   = C_URANUS
        this%neptunePerturbations(PERT_URANUS)%id     = PERT_URANUS
        this%neptunePerturbations(PERT_URANUS)%switch = .false.

        this%neptunePerturbations(PERT_NEPTUNE)%name   = C_NEPTUNE
        this%neptunePerturbations(PERT_NEPTUNE)%id     = PERT_NEPTUNE
        this%neptunePerturbations(PERT_NEPTUNE)%switch = .false.

        this%initialisedPerturbations = .true.

        if(isControlled()) then
            call checkOut(csubid)
        end if
        return

   end subroutine initPerturbations

  !==========================================================
  !
  !> @brief     Returns last position in ECEF (e.g. for re-entries)
  !!
  !! @author    Vitali Braun
  !!
  !> @return    4 element vector containing epoch, altitude, latitude and longitude
  !!
  !! @date      <ul>
  !!              <li> 16.10.2014 (initial design) </li>
  !!            </ul>
  !!
  !! @detail    A four-component real-valued vector is returned containing the
  !!            last epoch (MJD), last altitude (km), last latitude (deg) and last longitude (deg)
  !!
  !! @anchor    getLastPositionECEF
  !!
  !------------------------------------
 function getLastPositionECEF(this)

    class(Derivatives_class)    :: this
    real(dp), dimension(4)      :: getLastPositionECEF

    getLastPositionECEF(1) = this%last_epoch
    getLastPositionECEF(2) = this%last_altitude
    getLastPositionECEF(3) = this%last_latitude
    getLastPositionECEF(4) = this%last_longitude

    return

  end function getLastPositionECEF

  !==========================================================
  !
  !> @brief     Returns perturbation switch status
  !!
  !> @author    Vitali Braun
  !!
  !> @return    .true. when a perturbation is enabled
  !!
  !> @date      <ul>
  !!              <li> 17.11.2013 (initial design) </li>
  !!              <li> 13.03.2017 (adapted to new perturbations structure - and changed return value to boolean) </li>
  !!            </ul>
  !!
  !> @anchor    getPertSwitch
  !!
  !------------------------------------
  logical function getPertSwitch(this,iswitch)

    class(Derivatives_class)    :: this
    integer, intent(in)         :: iswitch

    if(iswitch > 0 .and. iswitch <= NUM_PERTURBATIONS) then
      getPertSwitch = this%neptunePerturbations(iswitch)%switch
    else
      getPertSwitch = .false.
    end if
    return

  end function getPertSwitch

  !==========================================================
  !
  !> @brief     Sets perturbation switch status
  !!
  !> @author    Vitali Braun
  !!
  !> @param[in] gravity_model   the gravity model
  !> @param[in] iswitch         Perturbation ID
  !> @param[in] val             Switch value
  !!
  !> @date      <ul>
  !!              <li> 17.11.2013 (initial design) </li>
  !!              <li> 13.03.2017 (adapted to new structure - and changed to boolean) </li>
  !!            </ul>
  !!
  !! @anchor    setPertSwitch
  !!
  !------------------------------------
  subroutine setPertSwitch(this,gravity_model,iswitch,val)

    class(Derivatives_class)            :: this
    type(Gravity_class),intent(inout)   :: gravity_model
    integer, intent(in)                 :: iswitch
    logical, intent(in)                 :: val

    character(len=*), parameter :: csubid = "setPertSwitch"

    if(isControlled()) then
        if(hasToReturn()) return
        call checkIn(csubid)
    end if

    ! check initialisation status
    if(.not. this%initialisedPerturbations) then
        call this%initPerturbations(gravity_model)
    end if

    if(iswitch < 1 .or. iswitch > NUM_PERTURBATIONS) then
        call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
        return
    else
        this%neptunePerturbations(iswitch)%switch = val
    end if

    !** done!
    if(isControlled()) then
        call checkOut(csubid)
    end if
    return

  end subroutine setPertSwitch

  !===============================================================================================
  !
  !> @anchor      deriv
  !!
  !> @brief       Force model evaluation for state vector propagation
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
  !> @param[in] r_in                radius vector
  !> @param[in] v_in                velocity vector
  !> @param[in] time_offset         time offset to start epoch
  !> @param[out] acc                resulting acceleration
  !!
  !! @date        <ul>
  !!                <li> 01.10.2012 (initial design)               </li>
  !!                <li> 04.06.2013 (code optimization)            </li>
  !!              </ul>
  !!
  !! @details     This routine provides the evaluation of the force function, which gives the total
  !!              acceleration in the inertial frame (GCRF).
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine deriv( this,               &
                    gravity_model,      &                                       ! <-> TYPE gravity model
                    atmosphere_model,   &                                       ! <-> TYPE atmosphere model
                    manoeuvres_model,   &                                       ! <-> TYPE manoeuvres model
                    radiation_model,    &                                       ! <-> TYPE radiation model
                    satellite_model,    &                                       ! <-> TYPE satellite model
                    solarsystem_model,  &                                       ! <-> TYPE solarsystem model
                    thirdbody_model,    &                                       ! <-> TYPE thirdbody model
                    tides_model,        &                                       ! <-> TYPE tides model
                    reduction,          &                                       ! <-> TYPE reduction
                    r_in,               &                                       ! <-- DBL() radius vector / km
                    v_in,               &                                       ! <-- DBL() velocity vector / km/s
                    time_offset,        &                                       ! <-- DBL   time offset to start epoch / s
                    acc                 &                                       ! --> DBL() acceleration in GCRF / km/s**2
                  )

    !** interface
    !--------------------------------------------
    class(Derivatives_class)                    :: this
    type(Gravity_class),intent(inout)           :: gravity_model
    type(Atmosphere_class),intent(inout)        :: atmosphere_model
    type(Manoeuvres_class),intent(inout)        :: manoeuvres_model
    type(Radiation_class),intent(inout)         :: radiation_model
    type(Satellite_class),intent(inout)         :: satellite_model
    type(Solarsystem_class),intent(inout)       :: solarsystem_model
    type(ThirdBody_class),intent(inout)         :: thirdbody_model
    type(Tides_class),intent(inout)             :: tides_model
    type(Reduction_type),intent(inout)          :: reduction
    real(dp), dimension(3), intent(in)          :: r_in
    real(dp), dimension(3), intent(in)          :: v_in
    real(dp),               intent(in)          :: time_offset

    real(dp), dimension(3), intent(out)         :: acc
    !--------------------------------------------

    character(len=*), parameter :: csubid = "deriv"
    integer :: ierr ! error flag
    integer :: i
    integer, parameter :: NUMBER_THIRD_BODIES = 9
    integer, dimension(NUMBER_THIRD_BODIES), parameter :: bodiesIndex  = (/ID_SUN, ID_MOON, ID_JUPITER, ID_VENUS, ID_MARS, &
                                                                           ID_MERCURY, ID_SATURN, ID_URANUS, ID_NEPTUNE/)
    integer, dimension(NUMBER_THIRD_BODIES), parameter :: bodiesSwitch = (/PERT_SUN, PERT_MOON, PERT_JUPITER, PERT_VENUS,  &
                                                                           PERT_MARS, PERT_MERCURY, PERT_SATURN, PERT_URANUS, PERT_NEPTUNE/)

    real(dp), dimension(3) :: r, v
    real(dp), dimension(3) :: acc_atmosphere ! accelerations due to atmosphere perturbations in inertial frame
    real(dp), dimension(3) :: acc_gravity    ! accelerations due to gravity potential perturbations in inertial frame
    real(dp), dimension(NUMBER_THIRD_BODIES,3) :: acc_3b         ! accelerations due to third bodies gravity perturbations in inertial frame
    real(dp), dimension(3) :: acc_3b_temp    ! accelerations due to one third body in inertial frame
    real(dp), dimension(3) :: acc_srp        ! accelerations due to solar radiation pressure in inertial frame
    real(dp), dimension(3) :: acc_alb        ! accelerations due to solar radiation pressure in inertial frame
    real(dp), dimension(3) :: acc_std        ! accelerations due to solid Earth tides in inertial frame
    real(dp), dimension(3) :: acc_otd        ! accelerations due to ocean tides in inertial frame
    real(dp), dimension(3) :: acc_man        ! accelerations due to orbital maneuvers in inertial frame

    real(dp)               :: altitude       ! altitude for given radius vector in ITRF (km)
    real(dp)               :: dtemp1, dtemp2 ! temporary doubles
    real(dp), dimension(3) :: r_ITRF         ! radius vector in ITRF
    real(dp)               :: time_mjd       ! MJD in days
    real(dp), dimension(3) :: v_ITRF         ! velocity vector in ITRF

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    ! convert passed time offset to actual epoch
    ! NOTE: It is important that the start_epoch has been set before - otherwise, this routine will not work.
    !       NEPTUNE takes care of this, while calling this routine from somewhere else should be done with care!
    time_mjd = this%start_epoch%mjd + time_offset/86400.d0
    !write(*,*) "TIME MJD = "
    !write(*,*) time_mjd, start_epoch%mjd
    !read(*,*)
    r = r_in
    v = v_in

    !** convert from GCRF to ITRF
    !-----------------------------------------------------
    call reduction%inertial2earthFixed(r, v, time_mjd, r_itrf, v_itrf)
    if(hasFailed()) return
!    write(*,*) "r_in = ", r_in, mag(r_in)
!    write(*,*) "v_in = ", v_in
!
!    write(*,*)
!    write(*,*) "r_itrf = ", r_itrf
!    write(*,*) "v_itrf = ", v_itrf

    call gravity_model%getGravityAcceleration(&
                                 reduction,   &
                                 r_itrf,      &   ! <-- DBL(3) radius vector in ITRF / km
                                 v_itrf,      &   ! <-- DBL(3) velocity vector in ITRF / km/s
                                 time_mjd,    &   ! <-- DBL    time in MJD
                                 acc_gravity  &   ! --> DBL(3) acceleration vector in GCRF / km/s**2
                               )
    if(hasFailed()) return

    if(this%neptunePerturbations(PERT_DRAG)%switch) then
      call atmosphere_model%getAtmosphereAcceleration(      &
                                       gravity_model,       &  ! <-> TYPE gravity model
                                       satellite_model,     &  ! <-> TYPE satellite model
                                       solarsystem_model,   &  ! <-> TYPE solarsystem model
                                       reduction,           &
                                       r,                   &  ! <-- DBL(3) radius vector in GCRF
                                       v,                   &  ! <-- DBL(3) velocity vector in GCRF
                                       r_itrf,              &  ! <-- DBL(3) radius vector in ITRF
                                       v_itrf,              & ! <-- DBL(3) velocity vector in ITRF
                                       time_mjd,            &  ! <-- DBL    current time (MJD)
                                       acc_atmosphere       &  ! --> DBL(3) acceleration vector in inertial frame
                                    )

      if(hasFailed()) then  ! e.g. due to re-entry
        ! save information on epoch, altitude, longitude and latitude
        call getGeodeticLatLon(r_itrf, this%last_altitude, this%last_latitude, this%last_longitude)
        this%last_epoch     = time_mjd
        acc_atmosphere = 0.d0
        return
      end if

    else

      !** check altitude for lower boundary
      call getGeodeticLatLon(r_itrf, altitude, dtemp1, dtemp2)

      if(altitude < atmosphere_model%getMinAltitude()) then
        call setNeptuneError(E_MIN_ALTITUDE, FATAL)
        ! save information on epoch, altitude, longitude and latitude
        call getGeodeticLatLon(r_itrf, this%last_altitude, this%last_latitude, this%last_longitude)
        this%last_epoch    = time_mjd
        return
      else
        !** set current altitude
        ierr = gravity_model%setCurrentAltitude(altitude)
        !** initialize accelerations
        acc_atmosphere = 0.d0
      end if

    end if

    ! all (currently only: Sun and Moon) third body perturbations are handled in the following loop
    do i = 1, 2 ! NUMBER_THIRD_BODIES
      if(this%neptunePerturbations(bodiesSwitch(i))%switch) then
        call thirdbody_model%getThirdBodyAcceleration(      &
                                       solarsystem_model,   &
                                       r,                   &    ! <-- DBL(3) radius vector in GCRF
                                       time_mjd,            &    ! <-- DBL    current time (MJD)
                                       bodiesIndex(i),      &    ! <-- CHR()  body identifier
                                       acc_3b_temp          &    ! --> DBL(3) acceleration vector in inertial frame
                                     )
        acc_3b(i,:) = acc_3b_temp
      else
        acc_3b(i,:) = 0.d0
      end if
    end do

    if(this%neptunePerturbations(PERT_SRP)%switch) then
      call radiation_model%getSrpAcceleration(      &
                                satellite_model,    &  ! <-> TYPE satellite model
                                solarsystem_model,  &  ! <-> TYPE solarsystem model
                                reduction,          &  ! <-> TYPE reduction
                                r,                  &  ! <-- DBL(3) GCRF radius of spacecraft
                                v,                  &  ! <-- DBL(3) GCRF velocity of spacecraft
                                time_mjd,           &  ! <-- DBL    current MJD
                                acc_srp             &  ! --> DBL(3) acceleration in inertial frame
                             )
    else
      acc_srp = 0.d0
    end if

    if(this%neptunePerturbations(PERT_ALBEDO)%switch) then
      call radiation_model%getAlbedoAcceleration(       &
                                    satellite_model,    &  ! <-> TYPE satellite model
                                    solarsystem_model,  &  ! <-> solarsystem model
                                    reduction,          &  ! <-> reduction
                                    r,                  &  ! <-- DBL(3) GCRF radius of spacecraft
                                    v,                  &  ! <-- DBL(3) GCRF velocity of spacecraft
                                    time_mjd,           &  ! <-- DBL    current MJD
                                    acc_alb             &  ! --> DBL(3) acceleration in inertial frame
                                )
    else
      acc_alb = 0.d0
    end if

    if(this%neptunePerturbations(PERT_SOLID_TIDES)%switch) then
      call tides_model%getTidesAcceleration(            &
                                  solarsystem_model,    &  ! <-> TYPE   solarsystem model
                                  reduction,            &  ! <-> TYPE   reduction
                                  r,                    &  ! <-- DBL(3) radius vector in GCRF frame
                                  r_itrf,               &  ! <-- DBL(3) radius vector in ITRF frame
                                  v_itrf,               &  ! <-- DBL(3) velocity vector in ITRF frame
                                  time_mjd,             &  ! <-- DBL    current MJD
                                  SOLID_TIDES,          &  ! <-- INT    tide type (1=solid, 2=ocean)
                                  acc_std               &  ! --> DBL(3) acceleration vector in inertial frame
                               )
    else
      acc_std = 0.d0
    end if

    if(this%neptunePerturbations(PERT_OCEAN_TIDES)%switch) then
      call tides_model%getTidesAcceleration(            &
                                  solarsystem_model,    &  ! <-> TYPE   solarsystem model
                                  reduction,            &  ! <-> TYPE   reduction
                                  r,                    &  ! <-- DBL(3) radius vector in GCRF frame
                                  r_itrf,               &  ! <-- DBL(3) radius vector in ITRF frame
                                  v_itrf,               &  ! <-- DBL(3) velocity vector in ITRF frame
                                  time_mjd,             &  ! <-- DBL    current MJD
                                  OCEAN_TIDES,          &  ! <-- INT    tide type (1=solid, 2=ocean)
                                  acc_otd               &  ! --> DBL(3) acceleration vector in inertial frame
                               )
    else
      acc_otd = 0.d0
    end if

    if(this%neptunePerturbations(PERT_MANEUVERS)%switch) then
      call manoeuvres_model%get_maneuver_acceleration(&
                                     reduction,     &
                                     r,             &  ! <-- DBL(3) radius vector in GCRF / km
                                     v,             &  ! <-- DBL(3) velocity vector in GCRF / km/s
                                     time_mjd,      &  ! <-- DBL    current MJD
                                     acc_man        &  ! --> DBL(3) acceleration vector in inertial frame
                                  )
    else
      acc_man = 0.d0
    end if

!if((time_mjd - 55276.75001d0) < eps15) then
! write(*,*) "+++ Date: ", time_mjd, " +++"
! write(*,*) "accgrav = ", acc_gravity, mag(acc_gravity)
! write(*,*) "atmo    = ", acc_atmosphere
! write(*,*) "sun     = ", acc_3b(1,:)
! write(*,*) "moon    = ", acc_3b(2,:)
! write(*,*) "srp     = ", acc_srp
! write(*,*) "solid   = ", acc_std
! write(*,*) "ocean   = ", acc_otd
! write(*,*) "albedo  = ", acc_alb


!end if
    if(.not. hasFailed()) then
      acc = acc_gravity + acc_atmosphere + acc_srp + acc_alb + acc_std + acc_otd + acc_man
      !** add also third body contributions
      do i = 1, 2 !NUMBER_THIRD_BODIES
        acc = acc + acc_3b(i,:)
        !write(*,*) "acc  = ", acc_3b(i,:), i
      end do
    else
      acc = 0.d0
      return
    end if

    !write(*,*) "acc  = ", acc
    !read(*,*)
    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine deriv

  !===============================================================================================
  !
  !> @anchor      deriv_cov
  !!
  !! @brief       State error transition matrix derivative evaluation for given state vector
  !! @author      Vitali Braun
  !! @author      Christopher Kebschull
  !!
  !> @param[in] gravity_model       the gravity model
  !> @param[in] atmosphere_model    the atmosphere model
  !> @param[in] manoeuvres_model    the manoeuvres model
  !> @param[in] radiation_model     the gravity model
  !> @param[in] satellite_model     the satellite model
  !> @param[in] solarsystem_model   the solar system model
  !> @param[in] thirdbody_model     the third body model
  !> @param[in] reduction           the reduction model
  !> @param[in] r_eci               radius vector
  !> @param[in] v_eci               velocity vector
  !> @param[in] set_matrix          state error transition matrix
  !> @param[in] time_mjd            Epoch in MJD
  !> @param[out] dset_matrix        derivatives of the state error transition matrix
  !!
  !! @date        <ul>
  !!                <li>VB:  04.07.2013 (initial design)</li>
  !!                <li>VB:  10.06.2014 (added variational equations for SRP) </li>
  !!                <li>CHK: 02.12.2018 (update to use atmosphere class) </li>
  !!              </ul>
  !!
  !! @details     This subroutine provides the derivative of the state error transition matrix
  !!              for the numerical integration routine.
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine deriv_cov( this,               &
                        gravity_model,      &                                   ! <-> TYPE gravity model
                        atmosphere_model,   &                                   ! <-> TYPE atmosphere model
                        radiation_model,    &                                   ! <-> TYPE radiation_model
                        satellite_model,    &                                   ! <-> TYPE satellite model
                        solarsystem_model,  &                                   ! <-> TYPE solarsystem model
                        thirdbody_model,    &                                   ! <-> TYPE thirdbody model
                        reduction,          &                                   ! <-> TYPE reduction
                        r_eci,              &                                   ! <-- DBL() radius vector [km]
                        v_eci,              &                                   ! <-- DBL() velocity vector [km/s]
                        setMatrix,          &                                   ! <-- DBL() state error transition matrix
                        time_mjd,           &                                   ! <-- DBL() current MJD
                        dsetMatrix          &                                   ! --> DBL() derivative of state error transition matrix
                      )

    !** interface
    !--------------------------------------------------------------------
    class(Derivatives_class)                    :: this
    type(Gravity_class),intent(inout)           :: gravity_model
    type(Atmosphere_class),intent(inout)        :: atmosphere_model
    type(Radiation_class),intent(inout)         :: radiation_model
    type(Satellite_class),intent(inout)         :: satellite_model
    type(Solarsystem_class),intent(inout)       :: solarsystem_model
    type(ThirdBody_class),intent(inout)         :: thirdbody_model
    type(Reduction_type),intent(inout)         :: reduction
    real(dp), dimension(3),   intent(in)        :: r_eci
    real(dp), dimension(3),   intent(in)        :: v_eci
    real(dp),                 intent(in)        :: time_mjd
    real(dp), dimension(:,:), intent(in)        :: setMatrix

    real(dp), dimension(size(setMatrix,1),size(setMatrix,1)), intent(out) :: dsetMatrix
    !--------------------------------------------------------------------

    character(len=*), parameter :: csubid = "deriv_cov"

    real(dp)                                                    :: mu           ! Earth's gravity constant
    real(dp), dimension(size(setMatrix,1),size(setMatrix,1))    :: pdm          ! partial derivative matrix
    real(dp), dimension(3)                                      :: r            ! radius vector in body-fixed frame
    real(dp), dimension(3)                                      :: v            ! velocity vector in body-fixed frame
    real(dp)                                                    :: r3           ! r**3
    real(dp)                                                    :: r5           ! r**5
    real(dp), dimension(3,3)                                    :: tempMat

    integer                                                     :: i            ! loop counter

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** convert ECI state to ECEF, as acceleration are evaluated in body-fixed frame
    call reduction%inertial2earthFixed(r_eci, v_eci, time_mjd, r, v)

    mu = getEarthGravity()

    r3 = mag(r)**3.d0
    r5 = mag(r)**5.d0

    !===========================================================
    !
    ! Compute elements of partial derivative matrix (for 2-body contribution)
    !
    ! + - - - + - - - +
    ! |       |       |
    ! |   0   |   E   |   E = identity matrix
    ! |       |       |
    ! + - - - + - - - +
    ! |       |       |
    ! |   G   |   0   |   G = partial derivatives (geopotential)
    ! |       |       |
    ! + - - - + - - - +
    !
    !----------------------------------------------------
    pdm   = 0.d0
    forall(i=1:3) pdm(i,i+3) = 1.d0

    pdm(4,1) = 3.d0*mu*r(1)*r(1)/r5 - mu/r3
    pdm(5,1) = 3.d0*mu*r(1)*r(2)/r5
    pdm(6,1) = 3.d0*mu*r(1)*r(3)/r5

    pdm(4,2) = pdm(5,1)
    pdm(5,2) = 3.d0*mu*r(2)*r(2)/r5 - mu/r3
    pdm(6,2) = 3.d0*mu*r(2)*r(3)/r5

    pdm(4,3) = pdm(6,1)
    pdm(5,3) = pdm(6,2)
    pdm(6,3) = 3.d0*mu*r(3)*r(3)/r5 - mu/r3


    !============================================================
    !
    ! Add contribution due to geopotential
    !
    !-----------------------------------------
    if(gravity_model%getGeoCovDegree() >= 2) then

      !write(*,*) "pdm wo j2"
      !do i=4,6
      !  write(*,*) (pdm(i,k), k=1,3)
      !end do

      pdm(4:6,1:3) = pdm(4:6,1:3) + gravity_model%getGravityCovariance(r, v)

      !write(*,*) "pdm w j2"
      !do i=4,6
      !  write(*,*) (pdm(i,k), k=1,3)
      !end do
      !read(*,*)

    end if

    !** convert to inertial frame, as other contributions are computed according to the inertial frame:
    tempMat = reduction%getRotationMatrixRC2IT(time_mjd)

    pdm(4:6,1:3) = matmul(matmul(transpose(tempMat), pdm(4:6,1:3)), tempMat)


    !============================================================
    !
    ! Add contribution due to atmospheric drag
    !
    !-----------------------------------------
    if(atmosphere_model%getDragCovFlag()) then

!     write(*,*) "Drag variational eqs...."
!
!     write(*,*) "pdm wo drag"
!     do i=4,6
!       write(*,'(6(e22.15e2))') (pdm(i,k), k=1,6)
!     end do

      pdm(4:6,1:6) = pdm(4:6,1:6) + atmosphere_model%getDragCovariance( gravity_model,      &
                                                                        satellite_model,    &
                                                                        solarsystem_model,  &
                                                                        reduction,          &
                                                                        r_eci,              &
                                                                        v_eci,              &
                                                                        time_mjd)

!    write(*,*) "pdm with drag"
!    do i=4,6
!      write(*,'(6(e22.15e2))') (pdm(i,k), k=1,6)
!    end do
!    read(*,*)

    end if

    if(hasFailed()) return

    !============================================================
    !
    ! Add contribution due to Sun and Moon
    !
    !-----------------------------------------
    if(thirdbody_model%getThirdBodyCovFlag(ID_SUN)) then

     !write(*,*) "pdm wo sun"
     !do i=4,6
     !  write(*,*) (pdm(i,k), k=1,3)
     !end do

      !** for Sun
      pdm(4:6,1:3) = pdm(4:6,1:3) + thirdbody_model%getThirdBodyCovariance(solarsystem_model, r_eci, time_mjd, ID_SUN)

     !write(*,*) "pdm w sun"
     !do i=4,6
     !  write(*,*) (pdm(i,k), k=1,3)
     !end do
     !read(*,*)

    end if

    if(hasFailed()) return ! either 'thirdbody_model%getThirdBodyCovFlag' or 'thirdbody_model%getThirdBodyCovariance'

    if(thirdbody_model%getThirdBodyCovFlag(ID_MOON)) then

      !write(*,*) "pdm wo moon"
      !do i=4,6
      !  write(*,*) (pdm(i,k), k=1,3)
      !end do

      !** for Moon
      pdm(4:6,1:3) = pdm(4:6,1:3) + thirdbody_model%getThirdBodyCovariance(solarsystem_model, r_eci, time_mjd, ID_MOON)

      !write(*,*) "pdm w moon"
      !do i=4,6
      !  write(*,*) (pdm(i,k), k=1,3)
      !end do
      !read(*,*)

    end if

    if(hasFailed()) return ! either 'thirdbody_model%getThirdBodyCovFlag' or 'thirdbody_model%getThirdBodyCovariance'

    !============================================================
    !
    ! Add contribution due to solar radiation pressure
    !
    !-----------------------------------------
    if(radiation_model%getSrpCovFlag()) then

      !write(*,*) "SRP variational eqs...."

!      write(*,*) "pdm wo SRP"
!      do i=4,6
!        write(*,'(6(e22.15e2))') (pdm(i,k), k=1,6)
!      end do

      pdm(4:6,1:3) = pdm(4:6,1:3) + radiation_model%getSrpCovariance(   satellite_model,    &
                                                                        solarsystem_model,  &
                                                                        reduction,          &
                                                                        r_eci,              &
                                                                        v_eci,              &
                                                                        time_mjd)

!      write(*,*) "pdm with SRP"
!      do i=4,6
!        write(*,'(6(e22.15e2))') (pdm(i,k), k=1,6)
!      end do
!      read(*,*)

    end if

    if(hasFailed()) return
    !write(50,'(9(e14.7e2,x))') pdm(4:6,1:3)

    !============================================================
    !
    ! Compute derivative of partial derivative matrix
    !
    !---------------------------------------------------
    dsetMatrix = matmul(pdm, setMatrix)

!    write(*,*) "dset with drag"
!    do i=1,6
!      write(*,'(6(e22.15e2))') (pdm(i,k), k=1,6)
!    end do
!    read(*,*)

!    write(51,'(36(e14.7e2,x))') dsetMatrix


    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine deriv_cov

end module derivatives
