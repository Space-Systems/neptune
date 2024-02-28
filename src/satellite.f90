!==================================================================================================
!!
!> @brief       Satellite macro model representation
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  13.03.2013 (initial design)</li>
!!                <li>VB:  04.06.2013 (code optimization) </li>
!!                <li>VB:  26.01.2013 (added oriented surfaces)</li>
!!                <li>VB:  31.03.2015 (added ESH feature)</li>
!!                <li>CHK: 13.11.2015 (updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 04.01.2018 (created Satellite_class)
!!              </ul>
!!
!> @details     This module contains type definitions and routines required for the representation
!!              of the satellite's oriented surfaces. It always holds the current orientation of
!!              all surfaces and performs an update as soon as one of the update functions is called,
!!              as long as the surfaces are rotating.
!!
!!              \par Includes
!!
!!              <ul>
!!                <li> neptune_error_handling.f90   </li>
!!                <li> slam_io.f90                  </li>
!!                <li> slam_math.f90                </li>
!!                <li> slam_reduction.f90           </li>
!!                <li> solarsystem.f90              </li>
!!                <li> slam_types.f90               </li>
!!              </ul>
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      satellite
!!
!!------------------------------------------------------------------------------------------------
module satellite

  use neptune_error_handling, only: E_NEGATIVE_MASS, E_MAX_SURFACES, E_INVALID_SATELLITE, &
                                    E_UNDEFINED_SATELLITE, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, WARNING, FATAL, hasFailed, &
                                    E_ID_MATCH, E_ANGLE_OUT_OF_RANGE, checkIn, checkOut
  use slam_io,                only: cdelimit, SEQUENTIAL, IN_FORMATTED, LOG_AND_STDOUT, openFile, message, nxtbuf
  use slam_math,              only: pi, halfPi, mag, deg2rad, eps9, cross, undefined, redang
  use slam_reduction_class,   only: Reduction_type
  use solarsystem,            only: Solarsystem_class, ID_SUN, getBodyPosition
  use slam_types,             only: dp

  implicit none

  private

  !** Flags
  integer, parameter, public :: MODE_CANNON_BALL = 0    ! Cannon ball mode (spherical)
  integer, parameter, public :: MODE_ORIENTED    = 1    ! Oriented surface mode

  integer, parameter, public :: ORIENT_SPHERE = 0  ! Sphere flag (no orientation required)
  integer, parameter, public :: ORIENT_EARTH  = 1  ! Earth oriented surface flag
  integer, parameter, public :: ORIENT_SUN    = 2  ! Sun oriented surface flag
  integer, parameter, public :: ORIENT_SPACE  = 3  ! Inertial oriented surface flag

  integer, parameter, public :: ID_ATMOSPHERE = 1 !> ID which shall be used if hasChangedSatellite() is called for atmosphere perturbations
  integer, parameter, public :: ID_SRP        = 2 !> ID which shall be used if hasChangedSatellite() is called for SRP perturbations
  integer, parameter, public :: ID_ALBEDO     = 3 !> ID which shall be used if hasChangedSatellite() is called for Albedo perturbations


  !** parameters
  integer, parameter, public :: MAX_SURFACES = 30   ! maximum number of surfaces which can be defined

  !** surfaces
  !-------------------------------------------------------------------------------
  type surface_t

    integer                :: id            ! surface identifier (0 = not specified)
    integer                :: orientation   ! surface orientation
    real(dp)               :: area          ! surface area in km**2
    real(dp), dimension(2) :: normal_angle  ! normal unit vector of the area described by two angles in radians
                                            !  1 = Right ascension, Azimuth or Sun incidence angle
                                            !  2 = Declination or Elevation
    real(dp), dimension(3) :: normal_GCRF   ! normal unit vector in GCRF
    real(dp)               :: reflSpec      ! spectral reflectivity coefficient
    real(dp)               :: reflDiff      ! diffuse reflectivity coefficient
    real(dp)               :: esh           ! equivalent solar hours of surface / h

  end type
  !-------------------------------------------------------------------------------

  !** macro-model of a satellite
  !-------------------------------------------------------------------------------
  type, public :: Satellite_class

    logical                                      :: eshFlag                     ! if true, the Equivalent Solar Hours are computed for the individual surfaces
    logical                                      :: cannon_ball_mode            ! becomes true if any of the specified surfaces has the Sphere flag set (ORIENT_SPHERE = 0)

    !** Input file
    character(len=512)                           :: csurfile                    ! surfaces definition file

    integer                                      :: nsurfaces                   ! number of different surfaces for a satellite
    logical                                      :: multiOriented               ! flag indicating that there are Earth- AND Sun-oriented surfaces
    real(dp)                                     :: cdrag                       ! drag coefficient of the satellite
    real(dp)                                     :: mass                        ! object mass
    real(dp)                                     :: mjd                         ! epoch for which data is valid (negative as long as no epoch is available)
    type(surface_t),dimension(MAX_SURFACES)      :: surface                     ! surfaces
    logical                                      :: isSet                       ! initialization flag

    logical, dimension(3)                        :: hasChanged                  ! indicator as soon as satellite properties change, e.g. mass, drag coefficient, etc.
                                                                                ! the array supports a flag for atmosphere (=1), SRP (=2) and albedo (=3) perturbations as all of them will
                                                                                ! access the data supported by the satellite type
    logical                                      :: check_coefficien_boundaries ! Whether or not it should be checked if c_r an c_d are in realistic ranges

    logical                                      :: first_ESH                   ! for initialization
    logical                                      :: first_orientation           ! for initialization
    real(dp)                                     :: startMJD                    ! first MJD required to sum up time
    real(dp), dimension(30)                      :: lastMJD                     ! MJD in last call for any surface - necessary to sum only in forward direction

    contains

        !** getter
        procedure :: getESH
        procedure :: getESHFlag
        procedure :: getObjectBallisticCoefficient
        procedure :: getObjectCrossSection
        procedure :: getObjectDragCoefficient
        procedure :: getObjectMass
        procedure :: getOrientationMode
        procedure :: getSurfaceDefinitionFileName
        procedure :: getSurfaceNormal
        procedure :: getSurfaceNumber
        procedure :: getSurfaceArea
        procedure :: getSurfaceSpecRefl
        procedure :: getSurfaceDiffRefl
        procedure :: hasChangedSatellite

        !** setter
        procedure :: resetObject
        procedure :: setObjectMass
        procedure :: setSurfaceDefinitionFileName
        procedure :: setBoundaryCheck

        procedure :: initObject
        procedure :: readSurfaceDefinitionFile
        procedure :: updateESH

        procedure, private  :: updateOrientation

    end type Satellite_class

    ! Constructor
    interface Satellite_class
        module procedure constructor
    end interface Satellite_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !!  @author     Christopher Kebschull
    !!  @date       <ul>
    !!                  <li>ChK: 24.12.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Satellite_class) function constructor()
        constructor%eshFlag             = .false.                               ! if true, the Equivalent Solar Hours are computed for the individual surfaces
        constructor%cannon_ball_mode    = .false.                               ! becomes true if any of the specified surfaces has the Sphere flag set (ORIENT_SPHERE = 0)
        constructor%csurfile            = 'neptune.srf'                         ! surfaces definition file
        constructor%multiOriented       = .false.                               ! flag indicating that there are Earth- AND Sun-oriented surfaces
        constructor%mjd                 = -1.d0                                 ! epoch for which data is valid (negative as long as no epoch is available)
        constructor%isSet               = .false.
        constructor%nsurfaces           = 1
        constructor%first_esh           = .true.
        constructor%first_orientation   = .true.

        constructor%surface(:)%id   = 0
        constructor%surface(:)%esh  = 0.0
        constructor%check_coefficien_boundaries = .true.

    end function constructor

!=============================================================================
!
!> @anchor      initObject
!!
!> @brief       Initialise the object
!> @author      Vitali Braun
!!
!> @param[in]   obj   object of type satellit_t, which initializes the module variable 'object'
!!
!> @date        <ul>
!!                <li> 07.11.2013 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  subroutine initObject(this)

    class(Satellite_class)      :: this

    character(len=*), parameter :: csubid = 'initObject'
    character(len=255)          :: cmess

    integer                     :: i    ! loop counter

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%isSet) then   ! initialization already done -> re-initialize!

      call this%resetObject()

    end if

    !=============================================
    !
    ! Check values for validity
    !
    !--------------------------------------

    !** mass and drag coefficient
    if(this%mass <= 0.d0) then
      cmess = "Satellite's mass negative or zero."
      call setNeptuneError(E_NEGATIVE_MASS, FATAL, (/cmess/))
      return

    else if(this%cdrag <= 0.d0) then

      cmess = "Satellite's drag coefficient is negative or zero."
      call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
      return

    end if

    !** for all surfaces, check that values are >= 0 and sum of Crs + Crd <= 1, also angles can have limits, e.g. 90° for declination/elevation
    do i=1,this%nsurfaces

      if(this%surface(i)%area < 0.d0) then

        write(cmess,'(a,i3)') "Area is negative for surface #", this%surface(i)%id
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return

      else if(this%surface(i)%reflSpec + this%surface(i)%reflDiff > 1.d0 .and. this%check_coefficien_boundaries) then

        cmess = "Sum of specular and diffuse reflectivity > 1."
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return

      else if(this%surface(i)%orientation /= ORIENT_SPHERE .and. &
              this%surface(i)%orientation /= ORIENT_EARTH  .and. &
              this%surface(i)%orientation /= ORIENT_SUN    .and. &
              this%surface(i)%orientation /= ORIENT_SPACE)  then

        write(cmess,'(a,i2,a)') "Unsupported orientation mode: '", this%surface(i)%orientation, "'."
        call setNeptuneError(E_INVALID_SATELLITE, FATAL, (/cmess/))
        return

      else if(abs(this%surface(i)%normal_angle(2)) > halfPi) then

        write(cmess,'(a,i2)') " Elevation can not be higher than +-90 deg. Set to '+90 deg'."
        call setNeptuneError(E_ANGLE_OUT_OF_RANGE, WARNING, (/cmess/))

        this%surface(i)%normal_angle(2) = halfPi

      end if


    end do

    !** correct area to km**2
    this%surface%area = this%surface%area*1.d-6 ! in km**2

    if(any(this%surface%orientation == ORIENT_SPHERE)) then

      this%cannon_ball_mode = .true.   ! a spherical shape is assumed.

    else if(any(this%surface%orientation == ORIENT_EARTH) .and. any(this%surface%orientation == ORIENT_SUN)) then

      this%multiOriented = .true.   ! e.g. if Earth oriented satellite with solar arrays is simulated

    end if

    this%hasChanged = .true.   ! new data available for all accessing modules
    this%isSet      = .true.   ! initialized object is available

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine initObject

!=============================================================================
!
!> @anchor      getESH
!!
!> @brief       Returns the ESH for a surface
!> @author      Vitali Braun
!!
!> @param[in]   id_in    Surface ID
!!
!> @return      Equivalent solar hours / h
!!
!> @date        <ul>
!!                <li> 31.03.2015 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function getESH(this,id_in)

    class(Satellite_class)  :: this
    integer, intent(in)     :: id_in

    getESH = this%surface(max(min(this%nsurfaces,id_in),0))%esh

    return

  end function getESH


!=============================================================================
!
!> @anchor      getESHFlag
!!
!> @brief       Returns the ESH flag
!> @author      Vitali Braun
!!
!> @return      .true. when the ESH is calculated
!!
!> @date        <ul>
!!                <li> 31.03.2015 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  logical function getESHFlag(this)

    class(Satellite_class)  :: this
    getESHFlag = this%eshFlag
    return

  end function getESHFlag

!===================================================================================
!
!> @anchor      updateESH
!!
!> @brief       Update the current counters of the ESH for the individual surfaces
!> @author      Vitali Braun
!!
!> @param[in]   timeMJD   Current datetime / MJD
!> @param[in]   cosinc    Solar incidence angle / deg
!> @param[in]   id        Surface ID
!!
!> @date        <ul>
!!                <li> 31.03.2015 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  subroutine updateESH(this, timeMJD, id, cosinc)

    class(Satellite_class)  :: this
    real(dp), intent(in)    :: timeMJD
    real(dp), intent(in)    :: cosinc
    integer,  intent(in)    :: id

    character(len=*), parameter :: csubid = 'updateESH'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%first_esh) then
      this%startMJD    = timeMJD
      this%lastMJD(:)  = timeMJD
      this%first_esh   = .false.
    end if

    !** do nothing and leave if MJD is less than last MJD
    if(timeMJD <= this%lastMJD(id)) then
      if(isControlled()) call checkOut(csubid)
      return
    end if

    if(this%isSet) then
      this%surface(id)%esh = this%surface(id)%esh + max(cosinc,0.d0)*(timeMJD - this%lastMJD(id))*24.d0 ! in hours
      this%lastMJD(id) = timeMJD
    end if

    if(isControlled()) call checkOut(csubid)
    return

  end subroutine updateESH

!=============================================================================
!
!> @anchor      getObjectMass
!!
!> @brief       Returns the mass of the object
!> @author      Vitali Braun
!!
!> @return      mass / kg
!!
!> @date        <ul>
!!                <li> 07.11.2013 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function getObjectMass(this)

    class(Satellite_class)      :: this
    character(len=*), parameter :: csubid = 'getObjectMass'

    getObjectMass = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    getObjectMass = this%mass

    if(isControlled()) then
      call checkOut(csubid)
    end if


  end function getObjectMass

!=============================================================================
!
!> @anchor      getObjectDragCoefficient
!!
!> @brief       Returns the drag coefficient of the object
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 27.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function getObjectDragCoefficient(this)

    class(Satellite_class)      :: this
    character(len=*), parameter :: csubid = 'getObjectDragCoefficient'

    getObjectDragCoefficient = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    getObjectDragCoefficient = this%cdrag

    if(isControlled()) then
      call checkOut(csubid)
    end if

  end function getObjectDragCoefficient

!=============================================================================
!
!> @anchor      getSurfaceNumber
!!
!> @brief       Returns the number of surfaces
!> @author      Vitali Braun
!!
!> @return      number of surfaces
!!
!> @date        <ul>
!!                <li> 27.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  integer function getSurfaceNumber(this)

    class(Satellite_class)      :: this
    character(len=*), parameter :: csubid = 'getSurfaceNumber'

    getSurfaceNumber = -1

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    getSurfaceNumber = this%nsurfaces

    if(isControlled()) then
      call checkOut(csubid)
    end if


  end function getSurfaceNumber
!=============================================================================
!
!> @anchor      getSurfaceNormal
!!
!> @brief       Returns the normal vector for a requested surface
!> @author      Vitali Braun
!!
!> @param[in]   id          ID of surface
!> @param[in]   time_mjd    MJD for vector orientation
!> @param[in]   r           radius vector in GCRF
!> @param[in]   v           velocity vector in GCRF
!!
!> @return      normal vector
!!
!> @date        <ul>
!!                <li> 28.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  function getSurfaceNormal(this, solarsystem_model, reduction, id, time_mjd, r, v)

    real(dp), dimension(3)                  :: getSurfaceNormal
    class(Satellite_class)                  :: this
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)     :: reduction
    real(dp),               intent(in)      :: time_mjd
    real(dp), dimension(3), intent(in)      :: r
    real(dp), dimension(3), intent(in)      :: v
    integer,                intent(in)      :: id

    character(len=*), parameter :: csubid = 'getSurfaceNormal'
    integer                     :: i
    logical                     :: found

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    found = .false.

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    do i=1,this%nsurfaces

      if(id == this%surface(i)%id) then

        !** perform an update
        call this%updateOrientation(solarsystem_model,reduction,time_mjd, r, v)

        getSurfaceNormal(:) = this%surface(i)%normal_GCRF(:)
        found = .true.
        exit

      end if

    end do

    if(.not. found) then
      call setNeptuneError(E_ID_MATCH, FATAL)
      getSurfaceNormal = 0.d0
      return
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if


  end function getSurfaceNormal


!=============================================================================
!
!> @anchor      getSurfaceSpecRefl
!!
!> @brief       Returns the specular reflectivity for requested surface.
!> @author      Vitali Braun
!!
!> @param[in]   id    Surface id
!!
!> @return      specular reflectivity
!!
!> @date        <ul>
!!                <li> 28.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function getSurfaceSpecRefl(this,id)

    class(Satellite_class)  :: this
    integer, intent(in)     :: id

    character(len=*), parameter :: csubid = 'getSurfaceSpecRefl'
    integer :: i    ! counter
    logical :: found

    getSurfaceSpecRefl = undefined
    found = .false.

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    do i=1,this%nsurfaces
      if(id == this%surface(i)%id) then
        getSurfaceSpecRefl = this%surface(i)%reflSpec
        found = .true.
        exit
      end if

    end do

    if(.not. found) then
      call setNeptuneError(E_ID_MATCH, FATAL)
      getSurfaceSpecRefl = -1.d0
      return
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSurfaceSpecRefl

!=============================================================================
!
!> @anchor      getSurfaceDiffRefl
!!
!> @brief       Returns the diffuse reflectivity for requested surface.
!> @author      Vitali Braun
!!
!> @param[in]   id    Surface id
!!
!> @return      diffuse reflectivity
!!
!> @date        <ul>
!!                <li> 28.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function getSurfaceDiffRefl(this,id)

    class(Satellite_class)  :: this
    integer, intent(in)     :: id

    character(len=*), parameter :: csubid = 'getSurfaceDiffRefl'
    integer :: i    ! counter
    logical :: found = .false.

    getSurfaceDiffRefl = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    do i=1,this%nsurfaces

      if(id == this%surface(i)%id) then
        getSurfaceDiffRefl = this%surface(i)%reflDiff
        found = .true.
        exit
      end if

    end do

    if(.not. found) then
      call setNeptuneError(E_ID_MATCH, FATAL)
      getSurfaceDiffRefl = -1.d0
      return
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSurfaceDiffRefl

!=============================================================================
!
!> @anchor      getSurfaceArea
!!
!> @brief       Returns the area for requested surface.
!> @author      Vitali Braun
!!
!> @param[in]   id    Surface id
!!
!> @return      area / m**2
!!
!> @date        <ul>
!!                <li> 28.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function getSurfaceArea(this, id)

    class(Satellite_class)  :: this
    integer, intent(in)     :: id

    character(len=*), parameter :: csubid = 'getSurfaceArea'
    integer :: i    ! counter
    logical :: found = .false.

    getSurfaceArea = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then
      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return
    end if

    do i=1,this%nsurfaces

      if(id == this%surface(i)%id) then
        getSurfaceArea = this%surface(i)%area
        found = .true.
        exit
      end if

    end do

    if(.not. found) then
      call setNeptuneError(E_ID_MATCH, FATAL)
      getSurfaceArea = -1.d0
      return
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSurfaceArea

!==================================================================================================
!
!> @anchor      getObjectBallisticCoefficient
!!
!> @brief       Returns the ballistic coefficient of the object for the given epoch in velocity direction
!> @author      Vitali Braun
!!
!> @param[in]   solarsystem_model The solar system model
!> @param[in]   reduction         reduction handler
!> @param[in]   time_mjd          Modified Julian Day of epoch
!> @param[in]   r                 Radius vector Earth -> Satellite in GCRF
!> @param[in]   v                 Velocity vector, actually being the vector the cross section is computed for
!!
!> @return      ballistic coefficient / kg/m**2
!!
!> @date        <ul>
!!                <li> 02.06.2014 (initial design)</li>
!!              </ul>
!!
!> @details     The ballistic coefficient (in kg/m**2) is computed via the cross-section as seen
!!              from the velocity direction without taking into account the shadowing effects
!!              which some surfaces may experience. Also the current mass and the drag parameter are
!!              taken into account. Note that the ballistic coefficient only takes into account
!!              the mass and the cross-section (not Cd or Cr).
!!
!---------------------------------------------------------------------------------------------------
  real(dp) function getObjectBallisticCoefficient(this, solarsystem_model, reduction, time_mjd, r, v) result(bc)

    class(Satellite_class)                  :: this
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)      :: reduction
    real(dp),               intent(in)      :: time_mjd
    real(dp), dimension(3), intent(in)      :: r
    real(dp), dimension(3), intent(in)      :: v

    character(len=*), parameter :: csubid = 'getObjectBallisticCoefficient'

    bc = 0.d0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then

      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return

    else

      bc = this%getObjectMass()/(this%getObjectCrossSection(solarsystem_model, reduction, time_mjd, r, v)*1.d6) ! 1.d6 to correct for 'm**2'
      if(hasFailed()) return

    end if

    !** done
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getObjectBallisticCoefficient

!==================================================================================================
!
!> @anchor      getObjectCrossSection
!!
!> @brief       Returns the cross section of the object for the given epoch in velocity direction
!> @author      Vitali Braun
!!
!> @param[in]   solarsystem_model The solar system model
!> @param[in]   reduction         reduction handler
!> @param[in]   time_mjd          Modified Julian Day of epoch
!> @param[in]   r                 Radius vector Earth -> Satellite in GCRF
!> @param[in]   v                 Velocity vector, actually being the vector the cross section is computed for
!!
!> @return      cross section / km**2
!!
!> @date        <ul>
!!                <li> 27.01.2014 (initial design)</li>
!!              </ul>
!!
!> @details     The cross-section as seen from the velocity direction is computed without taking
!!              into account the shadowing effects which some surfaces may experience.
!!
!---------------------------------------------------------------------------------------------------
  real(dp) function getObjectCrossSection(this,solarsystem_model,reduction,time_mjd, r, v)

    class(Satellite_class)                  :: this
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)     :: reduction
    real(dp),               intent(in)      :: time_mjd
    real(dp), dimension(3), intent(in)      :: r
    real(dp), dimension(3), intent(in)      :: v

    character(len=*), parameter :: csubid = 'getObjectCrossSection'
    integer  :: i       ! loop counter
    real(dp) :: temp    ! temporary

    getObjectCrossSection = 0.d0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%isSet) then

      call setNeptuneError(E_UNDEFINED_SATELLITE, FATAL)
      return

    else if(this%surface(1)%orientation == ORIENT_SPHERE) then    ! return cross-section for 'sphere' mode

      getObjectCrossSection = this%surface(1)%area
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return

    end if

    !** perform update of satellite's orientation
    call this%updateOrientation(solarsystem_model, reduction, time_mjd, r, v)


    !write(*,*) "nsurfaces = ", this%nsurfaces
    !write(*,*) "surface area = ", this%surface(1)%area, this%surface(2)%area
    !write(*,*) "normal 1: ", this%surface(1)%normal_GCRF
    !write(*,*) "normal 2: ", this%surface(2)%normal_GCRF

    !read(*,*)
    !** compute cross-section as sum of projected surfaces with cos < 90 deg.
    do i=1,this%nsurfaces

      temp = dot_product(this%surface(i)%normal_GCRF, v)

      if(temp > 0.d0) then

        getObjectCrossSection = getObjectCrossSection + temp/mag(v)*this%surface(i)%area

      end if

    end do

    !** done
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getObjectCrossSection


!========================================================================================================
!
!> @anchor      hasChangedSatellite
!!
!> @brief       Returns TRUE if there was any change to the satellite setup between different calls
!> @author      Vitali Braun
!!
!> @param[in]   id    Subroutine ID
!!
!> @return       .true. when there has been changes to the satellite
!!
!> @date        <ul>
!!                <li> 27.01.2014 (initial design)</li>
!!              </ul>
!!
!> @details     This routine returns a flag, which tells, whether there have been changes to parameters
!!              in the satellite configuration, which are specific to the calling routine. E.g. for
!!              atmospheric perturbations, the drag coefficient and the mass are important, while
!!              for SRP, it is the specular and diffuse reflectivity, as well as the mass.
!---------------------------------------------------------------------------------------------------------
  logical function hasChangedSatellite(this, id)

    class(Satellite_class)  :: this
    integer, intent(in)     :: id

    select case(id)
      case(ID_ATMOSPHERE)
        hasChangedSatellite = this%hasChanged(1)
        this%hasChanged(1)      = .false.
      case(ID_SRP)
        hasChangedSatellite = this%hasChanged(2)
        this%hasChanged(2)      = .false.
      case(ID_ALBEDO)
        hasChangedSatellite = this%hasChanged(3)
        this%hasChanged(3)      = .false.
      case default  ! for an unknown id, always TRUE will be returned, so that a re-acquiring of data is forced in case of some error which might have led here...
        hasChangedSatellite     = .true.

    end select

  end function hasChangedSatellite


!=============================================================================
!
!> @anchor      readSurfaceDefinitionFile
!!
!> @brief       Reads the surfaces definition file ('neptune.srf' as default)
!> @author      Vitali Braun
!> @author      Christopher Kebschull
!!
!> @param[in]   cpath   Path where the surfaces definition file is located.
!!
!> @date        <ul>
!!                <li>VB:  26.01.2014 (initial design)</li>
!!                <li>CHK: 04.01.2017 (simplified routine by applying object-oriented approach)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  subroutine readSurfaceDefinitionFile(this,cpath)

    class(Satellite_class)          :: this
    character(len=*), intent(in)    :: cpath

    character(len=255)                                          :: cbuf         ! character buffer
    character(len=(len_trim(cpath) + len(this%csurfile) + 1))   :: cfile
    character(len=*), parameter                                 :: csubid = 'readSurfaceDefinitionFile'

    integer                                     :: ich                          ! data file channel
    integer                                     :: itemp                        ! temporary
    type(surface_t),dimension(MAX_SURFACES)     :: tmp_surface                  ! temporary surfaces

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** open file
    cfile = trim(adjustl(cpath))//cdelimit//trim(adjustl(this%csurfile))
    ich = openFile(cfile, SEQUENTIAL, IN_FORMATTED)

    call message(' - Reading surfaces definition file...', LOG_AND_STDOUT)

    !=========================================
    !
    ! Read file
    !
    !-------------------------------

    !** Satellite mass
    call nxtbuf('#', 0, ich, cbuf)
    read(cbuf,*) this%mass

    !** Drag coefficient
    call nxtbuf('#', 0, ich, cbuf)
    read(cbuf,*) this%cdrag

    !** Equivalent Solar Hours feature
    call nxtbuf('#', 0, ich, cbuf)
    read(cbuf,*) itemp

    if(itemp /= 0) then
      this%eshFlag = .true.
    end if

    !** Read individual surfaces and count them

    this%nsurfaces = 0

    do

      call nxtbuf('#', 0, ich, cbuf)
      if(len_trim(cbuf) == 0) exit  ! end of file

      this%nsurfaces = this%nsurfaces + 1
      read(cbuf,*) tmp_surface(this%nsurfaces)%area,               &
                   tmp_surface(this%nsurfaces)%orientation,        &
                   tmp_surface(this%nsurfaces)%normal_angle(1),    &
                   tmp_surface(this%nsurfaces)%normal_angle(2),    &
                   tmp_surface(this%nsurfaces)%reflSpec,           &
                   tmp_surface(this%nsurfaces)%reflDiff

      this%surface(this%nsurfaces)%id = this%nsurfaces

    end do

    this%surface(1:this%nsurfaces) = tmp_surface(1:this%nsurfaces)
    this%surface(1:this%nsurfaces)%normal_angle(1) = tmp_surface(1:this%nsurfaces)%normal_angle(1)*deg2rad  ! convert to radians
    this%surface(1:this%nsurfaces)%normal_angle(2) = tmp_surface(1:this%nsurfaces)%normal_angle(2)*deg2rad  ! convert to radians

    !** now init object
    call this%initObject()

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine readSurfaceDefinitionFile


!=============================================================================
!
!> @anchor      setObjectMass
!!
!> @brief       Sets the mass of the object
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 08.11.2013 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  subroutine setObjectMass(this,mass)

    class(Satellite_class)  :: this
    real(dp), intent(in)    :: mass

    character(len=*), parameter :: csubid = 'setObjectMass'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(mass < 0.d0) then
      call setNeptuneError(E_NEGATIVE_MASS, FATAL)
      return
    else
      this%mass = mass
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine setObjectMass


  !=============================================================================
  !
  !> @anchor      setBoundaryCheck
  !!
  !> @brief       Sets whether or not to check coefficient boundaries
  !> @author      Daniel Lück
  !!
  !> @date        <ul>
  !!                <li> 28.02.2024 (initial design)</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------
    subroutine setBoundaryCheck(this,boundaryCheck)
  
      class(Satellite_class)  :: this
      logical, intent(in)     :: boundaryCheck
  
      character(len=*), parameter :: csubid = 'setBoundaryCheck'
  
      if(isControlled()) then
        if(hasToReturn()) return
        call checkIn(csubid)
      end if
  
      
      this%check_coefficien_boundaries = boundaryCheck
      
  
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
  
    end subroutine setBoundaryCheck


!=============================================================================
!
!> @anchor      resetObject
!!
!> @brief       Resets the surface ids and sets the init flag to .false.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 07.11.2013 (initial design)</li>
!!                <li> 04.01.2018 (surfaces are re-allocated)</li>
!!                <li> 09.03.2018 (surfaces are no longer dynamically allocated)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  subroutine resetObject(this)

    class(Satellite_class)  :: this

    this%surface(:)%id=0
    this%surface(:)%esh = 0.0

    this%isSet = .false.
    return

  end subroutine resetObject

!=============================================================================
!
!> @anchor      updateOrientation
!!
!> @brief       Updates the orientation of the satellite for the given epoch
!> @author      Vitali Braun
!!
!> @param[in]   time_mjd    MJD for which update has to be performed
!> @param[in]   r           Radius vector Earth -> satellite (GCRF)
!> @param[in]   v           Velocity vector (GCRF)
!!
!> @date        <ul>
!!                <li> 28.01.2014 (initial design)</li>
!!                <li> 02.02.2015 (fixed bug for multi-oriented surfaces in minimum incidence evaluation) </li>
!!              </ul>
!!
!> @details     This routine performs the update of all surface normals for
!!              the defined satellite configuration depending on the orientation
!!              mode. It is private and can thus only be called from within this
!!              module. Note that it does not have any means of error handling,
!!              as this routine assumes that the configuration of the satellite
!!              is available as soon as it is called.
!!
!-----------------------------------------------------------------------------
  subroutine updateOrientation(this, solarsystem_model, reduction, time_mjd, r, v)

    class(Satellite_class)                  :: this
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)     :: reduction
    real(dp), intent(in)                    :: time_mjd
    real(dp), dimension(3), intent(in)      :: r
    real(dp), dimension(3), intent(in)      :: v

    integer :: i    ! loop counter

    real(dp)               :: de_norm   ! declination of the normal vector in radian
    real(dp)               :: de_sun    ! declination of the sun in radian
    real(dp), dimension(3) :: h         ! orbit normal vector
    real(dp)               :: habs      ! magnitude of h
    real(dp)               :: incidence ! sun incidence angle for surface
    real(dp)               :: minIncidence ! minimum sun incidence angle for surface (due to orbit geometry)
    real(dp), dimension(3) :: n_uvw     ! normal vector in UVW coordinates
    real(dp)               :: ra_norm   ! right ascension of the normal vector in radian
    real(dp)               :: ra_sun    ! right ascension of the sun in radian
    real(dp), dimension(3) :: rsun      ! Sun vector in GCRF
    real(dp)               :: rsun_abs  ! magnitude of rsun
    real(dp)               :: temp      ! temporary

    !** check whether update is required
    !--------------------------------------------------------------------
    if(abs(this%mjd - time_mjd) < eps9) then    ! no update required

      return

    else if(this%mjd < 0.d0) then   ! called for first time, so that time_mjd actually is the initial epoch

      this%first_orientation      = .true.

    end if

    !** perform update
    !-------------------------------------------------------------------

    this%mjd = time_mjd

    do i=1,this%nsurfaces

      select case(this%surface(i)%orientation)

        case(ORIENT_EARTH)

          !** coordinates in UVW system
          n_uvw(1) = sin(this%surface(i)%normal_angle(2))                                         ! u = sin(El)
          n_uvw(2) = cos(this%surface(i)%normal_angle(2))*cos(this%surface(i)%normal_angle(1))  ! v = cos(El)*cos(Az)
          n_uvw(3) = cos(this%surface(i)%normal_angle(2))*sin(this%surface(i)%normal_angle(1))  ! w = cos(El)*sin(Az)

          !** convert to inertial frame
          call reduction%uvw2eci(r, v, n_uvw, this%surface(i)%normal_GCRF)

        case(ORIENT_SUN)

          !** get Sun orientation
          rsun     = solarsystem_model%getBodyPosition(time_mjd, ID_SUN)
          rsun_abs = mag(rsun)

          !** compute Sun's right ascension and declination
          ra_sun = atan2(rsun(2), rsun(1))
          de_sun = asin(rsun(3)/rsun_abs)

          !** now add surface orientation - first right ascension:
          ra_norm = ra_sun + this%surface(i)%normal_angle(1)

          !** correct to interval [0:360] deg
          ra_norm = redang(ra_norm, 2, 1, .false.)

          de_norm = de_sun + this%surface(i)%normal_angle(2)

          !** for declination, it is important to note, that a correction is required, as soon as
          !   the difference in right ascension is in the interval [90:270] deg.
          if((ra_norm > halfPi) .and. (ra_norm - halfPi < pi)) then
            de_norm = -1.d0*de_norm
          end if

          !** get minimum incidence angle if Earth-oriented satellite with Sun-oriented surfaces
          if(this%multiOriented) then

            !** check whether the values are above minimum incidence angle
            incidence = acos(cos(ra_norm)*cos(de_norm))

            h    = cross(r,v)
            habs = mag(h)

            temp = acos(dot_product(rsun, h)/rsun_abs/habs)

            if(temp > halfPi) then
              temp = pi - temp
            end if

            minIncidence = halfPi - temp

            if(minIncidence > incidence) then ! correct Delta right ascension or Delta elevation

              if(minIncidence > this%surface(i)%normal_angle(1)) then   ! add some declination

                temp    = acos(cos(minIncidence)/cos(this%surface(i)%normal_angle(1)))
                de_norm = de_norm + temp

              else  ! add some right ascension

                temp    = acos(cos(minIncidence)/cos(this%surface(i)%normal_angle(2)))
                ra_norm = ra_norm + temp

              end if

            end if

          end if

          !** convert to GCRF
          this%surface(i)%normal_GCRF(1) = cos(ra_norm)*cos(de_norm)
          this%surface(i)%normal_GCRF(2) = sin(ra_norm)*cos(de_norm)
          this%surface(i)%normal_GCRF(3) = sin(de_norm)

        case(ORIENT_SPACE)

          if(this%first_orientation) then

            this%surface(i)%normal_GCRF(1) = cos(this%surface(i)%normal_angle(1))*cos(this%surface(i)%normal_angle(2))
            this%surface(i)%normal_GCRF(2) = sin(this%surface(i)%normal_angle(1))*cos(this%surface(i)%normal_angle(2))
            this%surface(i)%normal_GCRF(3) = sin(this%surface(i)%normal_angle(2))

          end if

        case default
          cycle

      end select

    end do


  end subroutine updateOrientation

!=============================================================================
!
!> @anchor      getOrientationMode
!!
!> @brief       Get the orientation mode which is used in the simulation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 03.03.2014: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  integer function getOrientationMode(this)

    class(Satellite_class)  :: this
    if(this%cannon_ball_mode) then
      getOrientationMode = MODE_CANNON_BALL
    else
      getOrientationMode = MODE_ORIENTED
    end if

    return

  end function getOrientationMode

!=============================================================================
!
!> @anchor      getSurfaceDefinitionFileName
!!
!> @brief       Get the name of the surfaces definition file
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 26.01.2014: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  character(len=255) function getSurfaceDefinitionFileName(this)

    class(Satellite_class)  :: this
    getSurfaceDefinitionFileName = this%csurfile
    return

  end function getSurfaceDefinitionFileName


!=============================================================================
!
!> @anchor      setSurfaceDefinitionFileName
!!
!> @brief       Sets the name of the surfaces definition file
!> @author      Vitali Braun
!!
!> @param[in]   val     Requested name of the SD file
!!
!> @date        <ul>
!!                <li> 26.01.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------

  subroutine setSurfaceDefinitionFileName(this,val)

    class(Satellite_class)          :: this
    character(len=*), intent(in)    :: val

    this%csurfile = val(1:min(len(val),len(this%csurfile)))
    return

  end subroutine setSurfaceDefinitionFileName


end module satellite


