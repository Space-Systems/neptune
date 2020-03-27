!>-----------------------------------------------------------------------------------------------
!!
!> @brief       Radiation pressure modeling
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  12.03.2013 (initial design)</li>
!!                <li>VB:  04.06.2013 (code optimization) </li>
!!                <li>VB:  28.01.2014 (added oriented surfaces) </li>
!!                <li>VB:  10.06.2014 (added variational equations due to SRP) </li>
!!                <li>VB:  31.03.2014 (added ESH feature, which only works in combination with SRP perturbations)</li>
!!                <li>CHK: 13.11.2015 (updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 03.01.2018 (created Radiation_class)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for radiation
!!              pressure perturbations modeling, which are caused by the Sun and the Earth,
!!              specifically in the context of numerical integration of a satellite trajectory.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      radiation
!!
!! @todo        Computation of Cr from surface properties (Doornbos PhD, p.66)
!! @todo        PRIORITY: check results for HAMR cases
!!
!!------------------------------------------------------------------------------------------------
module radiation

  use slam_types,             only: dp
  use slam_astro,             only: getEarthRadius, getSunRadius, getAstronomicalUnit
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, checkIn, checkOut
  use slam_math,              only: mag, pi, halfPi, twopi, eps9, identity_matrix, outerproduct, undefined
  use slam_reduction_class,   only: Reduction_type
  use satellite,              only: Satellite_class, MAX_SURFACES, ID_ALBEDO, ID_SRP, MODE_CANNON_BALL
  use solarsystem,            only: Solarsystem_class, ID_SUN
  use slam_strings,           only: toUppercase
  use slam_time,              only: sec_per_day

  implicit none

  private

    !** parameters
    !---------------------------------------------------------------------------------------------------------------------
    character(len=*), parameter :: C_CONICAL = "CONICAL"

    integer, parameter :: CONICAL_SHADOW = 1                                    ! conical shadow model
    integer, parameter :: NO_SHADOW      = 0                                    ! no shadow

    integer, parameter :: nrings = 2                                            ! number of rings used for the Earth plate discretisation in the Albedo model

    real(dp), parameter :: length_of_year = 365.2425d0
    real(dp), dimension(21), parameter :: aphelions = (/                        & ! array containing the aphelion dates (MJD) for the years 2000 - 2020
                                            51729.d0,         52094.583d0,      &
                                            52461.16666667d0, 52824.25d0,       &
                                            53191.45833333d0, 53556.20833333d0, &
                                            53919.95833333d0, 54288.d0,         &
                                            54651.33333333d0, 55016.08333333d0, &
                                            55383.5d0,        55746.625d0,      &
                                            56113.16666667d0, 56478.625d0,      &
                                            56842.d0,         57209.83333333d0, &
                                            57573.66666667d0, 57937.83333333d0, &
                                            58305.70833333d0, 58668.91666666d0, &
                                            59034.5d0/)

    type, public :: Radiation_class

        !** variables
        !-------------------------------------------------------------------------------
        integer :: ringElements                                                 ! total number of ring elements based on the number of rings (will be set in the initialization)
        integer :: shadowMode                                                   ! conical shadow model as default...

        logical :: AlbInitialized                                               ! for initialization check
        logical :: covSrp                                                       ! considering SRP variational equations in covariance propagation
        logical :: SrpInitialized                                               ! for initialization check

        logical :: first_srp                                                    ! first call ops


        real(dp), dimension(MAX_SURFACES) :: areaSRP                            ! array containing the area of each surface in km**2
        real(dp), dimension(MAX_SURFACES) :: crsSRP                             ! specular reflectivity coefficient
        real(dp), dimension(MAX_SURFACES) :: crdSRP                             ! diffuse reflectivity coefficient


        ! Formerly saved variables
        real(dp), private                           :: last_srp_mjd
        real(dp), dimension(3), private             :: last_srp_r
        real(dp), dimension(3), private             :: last_srp_acc
        real(dp), private                           :: saved_mass               ! satellite mass
        integer, private                            :: apidx                    ! index of aphelion array, which is saved to accelerate array index search
        real(dp), private                           :: last_press               ! saving value for last call
        real(dp), private                           :: last_mjd                 ! saving last MJD

        logical, private                            :: first_albedo
        real(dp), dimension(MAX_SURFACES), private  :: area_albedo              ! array containing the area of each surface in km**2
        real(dp), dimension(MAX_SURFACES), private  :: crs_albedo               ! specular reflectivity coefficient
        real(dp), dimension(MAX_SURFACES), private  :: crd_albedo               ! diffuse reflectivity coefficient

        real(dp), private                           :: mass_albedo              ! object mass (kg)
        real(dp), private                           :: mass_srp                 ! object mass (kg)




    contains
        !** get
        procedure :: getAlbedoAcceleration
        procedure :: getSrpAcceleration
        procedure :: getSrpCovariance
        procedure :: getSrpCovFlag
        procedure :: getSrpWithoutShadow

        !** set
        procedure :: setAlbedoInitFlag
        procedure :: setShadowModel
        procedure :: setSrpCovFlag
        procedure :: setSrpInitFlag

        !** others
        procedure :: checkBoundaryCrossing
        procedure :: initAlbedo

        procedure, private :: getShadow
        procedure, private :: getSolarPressure

    end type Radiation_class

  ! Constructor
    interface Radiation_class
        module procedure constructor
    end interface Radiation_class

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
    type(Radiation_class) function constructor()
        constructor%shadowMode      = CONICAL_SHADOW
        constructor%AlbInitialized  = .false.                                   ! for initialization check
        constructor%covSrp          = .false.                                   ! considering SRP variational equations in covariance propagation
        constructor%SrpInitialized  = .false.                                   ! for initialization check

        constructor%first_srp       = .true.                                    ! first call ops

        constructor%last_srp_mjd    = -1.d0
        constructor%last_srp_r      =  0.d0
        constructor%last_srp_acc    =  0.d0
        constructor%apidx           = 1
        constructor%last_press      = -1.d0
        constructor%last_mjd        = -1.d0
        constructor%first_albedo = .true.
    end function constructor

!!------------------------------------------------------------------------------------------------
!> @anchor      getSrpCovFlag
!!
!> @brief       Gets the appropriate flag covSrp
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 10.06.2014 (initial design)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  logical function getSrpCovFlag(this)

    class(Radiation_class)  :: this

    getSrpCovFlag = this%covSrp
    return

  end function getSrpCovFlag

!!------------------------------------------------------------------------------------------------
!> @anchor      setSrpCovFlag
!!
!> @brief       Sets the appropriate flag covSrp
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial design)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  subroutine setSrpCovFlag(this,val)

    class(Radiation_class)  :: this
    logical, intent(in)     :: val

    this%covSrp = val
    return

  end subroutine setSrpCovFlag

!------------------------------------------------------------------------------------------------
!
!> @anchor      getSrpCovariance
!!
!> @brief       Gets the contributions to the variational equations due to solar radiation pressure
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 10.06.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   r_gcrf            radius vector of satellite in GCRF
!> @param[in]   v_gcrf            velocity vector of satellite in GCRF
!> @param[in]   time_mjd          MJD
!!
!!------------------------------------------------------------------------------------------------
  function getSrpCovariance(                    &
                            this,               &
                            satellite_model,    &
                            solarsystem_model,  &
                            reduction,          &
                            r_gcrf,             &
                            v_gcrf,             &
                            time_mjd            &
                            ) result(cov)

    implicit none

    !** interface
    !---------------------------------------------------
    class(Radiation_class)              :: this
    type(Satellite_class)               :: satellite_model
    type(Solarsystem_class)             :: solarsystem_model
    type(Reduction_type)               :: reduction
    real(dp), dimension(3), intent(in)  :: r_gcrf
    real(dp), dimension(3), intent(in)  :: v_gcrf
    real(dp),               intent(in)  :: time_mjd

    real(dp), dimension(3,3) :: cov     ! result in 3x3 matrix
    !---------------------------------------------------

    character(len=*), parameter :: csubid =  "getSrpCovariance"

    real(dp)                 :: cosinc        ! cosine of incidence angle / rad
    real(dp)                 :: fshadow       ! shadow factor
    real(dp), dimension(3,3) :: mati          ! identity matrix
    real(dp), dimension(3)   :: r_sun         ! GCRF radius of Sun
    real(dp), dimension(3)   :: ros           ! radius vector from object to Sun
    real(dp)                 :: ros_abs       ! magnitude of vector ros
    real(dp)                 :: ros_abs2      ! ros_abs*ros_abs

    real(dp)                 :: rsun_abs      ! magnitude of Sun radius vector
    real(dp)                 :: rsun_abs2     ! rsun_abs*rsun_abs
    real(dp)                 :: solarPressure ! normalized solar flux at 1 AU
    real(dp)                 :: temp          ! temporary variable to sum all surfaces
    real(dp), dimension(3)   :: vecn          ! surface normal vector

    integer :: i
    integer :: surf   ! number of surfaces

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !=======================================================
    !
    ! Get SRP coefficient and mass
    !
    !-------------------------------------------------------
    surf = satellite_model%getSurfaceNumber()
    if(hasFailed()) return

    if(satellite_model%hasChangedSatellite(ID_SRP) .or. this%first_srp)  then    ! only required, if configuration has changed or called for first time

      this%saved_mass = satellite_model%getObjectMass()
      if(hasFailed()) return

      do i=1,surf

        this%crsSRP(i)  = satellite_model%getSurfaceSpecRefl(i)
        if(hasFailed()) return

        this%crdSRP(i)  = satellite_model%getSurfaceDiffRefl(i)
        if(hasFailed()) return

        this%areaSRP(i) = satellite_model%getSurfaceArea(i)
        if(hasFailed()) return

      end do

      this%first_srp = .false.
      this%crdSRP(:) = this%crdSRP(:)/3.d0  ! only this factor will be required!

    end if

    !** initialize identity matrix
    call identity_matrix(mati)

    !=====================================================
    !
    ! Get Sun position
    !
    !----------------------------

    r_sun = solarsystem_model%getBodyPosition(time_mjd, ID_SUN)

    !=====================================================
    !
    ! Compute shadow factor
    !
    !----------------------------

    if(this%shadowMode == NO_SHADOW) then
      fshadow = 1.d0
    else
      fshadow = this%getShadow(     &
                            r_sun,  &  ! <-- DBL(3)  GCRF radius of sun
                            r_gcrf  &  ! <-- DBL(3)  GCRF radius of object
                         )
      if(hasFailed()) return
    end if

    !=====================================================
    !
    ! Get solar irradiance
    !
    !------------------------------------

    solarPressure = this%getSolarPressure(time_mjd)
    if(hasFailed()) return

    !=====================================================
    !
    ! Compute SRP acceleration
    !
    !----------------------------

    !** vector object -> sun
    ros       = r_sun - r_gcrf
    ros_abs   = mag(ros)
    ros_abs2  = ros_abs*ros_abs
    rsun_abs  = mag(r_sun)
    rsun_abs2 = rsun_abs*rsun_abs

    !======================================================
    !
    ! Build result matrix
    !
    !-----------------------------------------------------
    if(fshadow > eps9) then

      if(satellite_model%getOrientationMode() == MODE_CANNON_BALL) then

        temp = this%areaSRP(1)*(1.d0 + this%crsSRP(1))

      else

        temp = 0.d0

        do i=1,surf

          !** get normal vector in GCRF
          vecn = satellite_model%getSurfaceNormal(solarsystem_model, reduction, i, time_mjd, r_gcrf, v_gcrf)
          if(hasFailed()) return

          cosinc = dot_product(vecn,ros)

          temp = temp + this%areaSRP(i)*cosinc*(1.d0 + this%crsSRP(i))

        end do

      end if

      cov(:,:) = -fshadow*temp*solarPressure*rsun_abs2/(this%saved_mass*ros_abs2*ros_abs)*(mati(:,:) - 3.d0*outerproduct(ros, ros)/ros_abs2)

    else

      cov(:,:) = 0.d0

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSrpCovariance


!=================================================================================================
!
!> @anchor      setAlbedoInitFlag
!!
!> @brief       Set initialization flag for albedo perturbations to .false.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 05.06.2013 (initial design) </li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  subroutine setAlbedoInitFlag(this)

    class(Radiation_class)  :: this

    this%AlbInitialized = .false.

  end subroutine setAlbedoInitFlag


!=================================================================================================
!
!> @anchor      setSrpInitFlag
!!
!> @brief       Set initialization flag for SRP perturbations to .false.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 05.06.2013 (initial design) </li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  subroutine setSrpInitFlag(this)

    class(Radiation_class)  :: this

    this%SrpInitialized = .false.

  end subroutine setSrpInitFlag

!==================================================================================================
!
!> @anchor      initAlbedo
!!
!> @brief       Initialization of Earth radiation pressure routine
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 24.01.2014 (initial design)</li>
!!              </ul>
!!
!> @details     This routine initializes the Earth radiation pressure (Albedo) routine.
!!              Global parameters are set and Earth's surface data is read from an input file.
!!
!!------------------------------------------------------------------------------------------------

  subroutine initAlbedo(this)

    class(Radiation_class)      :: this
    character(len=*), parameter :: csubid = "initAlbedo"

    integer :: i    ! loop counter

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%AlbInitialized) then
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
    end if

    !** compute number of ring elements based on the number of rings
    this%ringElements = 1

    do i=1,nrings

      this%ringElements = this%ringElements + 6*i

    end do

    !** initialization complete
    this%AlbInitialized = .true.

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine initAlbedo

!==================================================================================================
!
!> @anchor      initSRP
!!
!> @brief       Initialization of solar radiation pressure routine
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.03.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!> @details     This routine initializes the solar radiation pressure (SRP) routine.
!!              Global parameters, e.g. the solar constant, are set and surface data
!!              is read from an input file into the specific surfaces array.
!!
!!------------------------------------------------------------------------------------------------

! subroutine initSRP()
!
!   character(len=*), parameter :: csubid = "initSRP"

!   if(isControlled()) then
!     if(hasToReturn()) return
!     call checkIn(csubid)
!   end if

!   !** TBD

!   this%SrpInitialized = .true.

!   !** done!
!   if(isControlled()) then
!     call checkOut(csubid)
!   end if

!   return

! end subroutine initSRP

!======================================================================================
!
!> @anchor      setShadowModel
!!
!> @brief       Set shadow model to be used in propagation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 14.03.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   mode_name     Selected mode name
!!
!> @details     This function sets the shadow model to be used. It accepts only documented
!!              input. For unknown entries, no shadow will be used.
!!
!!------------------------------------------------------------------------------------------------
  subroutine setShadowModel(this,mode_name)

    class(Radiation_class)          :: this
    character(len=*), intent(in)    :: mode_name

    if(toUppercase(mode_name) == C_CONICAL) then
      this%shadowMode = CONICAL_SHADOW
    else
      this%shadowMode = NO_SHADOW
    end if

  end subroutine setShadowModel

!======================================================================================
!
!> @anchor      getShadow
!!
!> @brief       Returns shadow factor indicating shadowing conditions Earth/Sun
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.03.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!> @param[in]   r_sun      Radius vector of the Sun in GCRF
!> @param[in]   r_sat      Radius vector of the Satellite in GCRF
!!
!> @details     This function computes a shadowing factor f between 0 and 1, indicating whether
!!              the object is in umbra (f = 1), penumbra (0 < f < 1) or in full sunlight (f = 0).
!!              An algorithm from Montenbruck (2000) is used.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Compute vector sun<>satellite                                    </li>
!!                <li> Compute apparent radii of Sun and Earth as seen by the satellite </li>
!!                <li> Occultation check and factor computation </li>
!!                <li> Finish.           </li>
!!              </ol>
!!
!! @todo        Check Montenbruck's shadow function and add Earth's flattening
!!
!!------------------------------------------------------------------------------------------------

  real(dp) function getShadow(          &
                              this,     &
                              r_sun,    &  ! <-- DBL(3)  radius of the sun in GCRF
                              r_sat     &  ! <-- DBL(3)  radius of the sat in GCRF
                             )

    !** interface
    !----------------------------------------------
    class(Radiation_class)              :: this
    real(dp), dimension(3), intent(in)  :: r_sun
    real(dp), dimension(3), intent(in)  :: r_sat
    !----------------------------------------------

    character(len=*), parameter :: csubid = "getShadow"

    real(dp)               :: a         ! apparent radius of Sun as seen by satellite
    real(dp)               :: area_occ  ! occulted area of the Sun
    real(dp)               :: b         ! apparent radius of Earth as seen by satellite
    real(dp)               :: c         ! apparent separation of the Sun and Earth center
    real(dp)               :: res_abs   ! magnitude of vector r_sun
    real(dp), dimension(3) :: rss       ! vector satellite => sun
    real(dp)               :: rss_abs   ! magnitude of rss
    real(dp)               :: x1        ! auxiliary
    real(dp)               :: y1        ! auxiliary

    getShadow = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** vector satellite => sun
    rss = r_sun - r_sat

    !** distance satellite-sun [km]
    rss_abs = mag(rss)

    !** distance satellite-earth [m]
    res_abs = mag(r_sat)

    !** apparent radius of sun as seen from satellite
    a = asin(getSunRadius()/rss_abs)

    !** apparent radius of earth as seen from satellite
    b = asin(getEarthRadius()/res_abs)

    !** apparent separation of the sun and earth center
    c = acos(-dot_product(r_sat,rss)/(res_abs*rss_abs))

    !** check if occultation takes place
    if(c >= (a+b)) then

      getShadow = 1.d0

    else if(c < abs(a-b)) then

      getShadow = 0.d0

    else

      x1 = (c*c + a*a - b*b)/(2*c)
      y1 = sqrt(a*a - x1*x1)

      !** occulted area of the sun
      area_occ = a*a*acos(x1/a) + b*b*acos((c-x1)/b) - c*y1

      !** remaining fraction of sunlight
      getShadow = 1.d0 - area_occ/(pi*a*a)

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getShadow

!==============================================================================
!
!> @anchor      getSrpWithoutShadow
!!
!> @brief       Returns solar radiation pressure acceleration as if there was no shadow
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 14.03.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   r           GCRF radius vector / km
!> @param[in]   v           GCRF velocity vector / km/s
!> @param[in]   time_mjd    Current time / MJD
!!
!> @return      acceleration vector
!!
!> @details     This function returns the current solar radiation pressure acceleration
!!              as if there was no shadow.
!!
!!------------------------------------------------------------------------------------------------
  function getSrpWithoutShadow(this, &
                                satellite_model,    &
                                solarsystem_model,  &
                                reduction,          &
                                r,                  &
                                v,                  &
                                time_mjd)           &
                                result(awos)

    real(dp), dimension(3)             :: awos                                  ! acceleration without shadow
    class(Radiation_class)             :: this
    type(Satellite_class)              :: satellite_model
    type(Solarsystem_class)            :: solarsystem_model
    type(Reduction_type)              :: reduction
    real(dp), dimension(3), intent(in) :: r
    real(dp), dimension(3), intent(in) :: v
    real(dp),               intent(in) :: time_mjd

    character(len=*), parameter :: csubid = 'getSrpWithoutShadow'
    integer :: currentShadowMode

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** store current shadow configuration...
    currentShadowMode = this%shadowMode

    !** change to no shadow...
    this%shadowMode = NO_SHADOW

    !** get SRP acceleration without shadow
    call this%getSrpAcceleration(satellite_model, solarsystem_model, reduction, r, v, time_mjd, awos)
    if(hasFailed()) return

    !** reset shadow mode
    this%shadowMode = currentShadowMode

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSrpWithoutShadow

!==============================================================================
!
!> @anchor      getSolarPressure
!!
!> @brief       Returns solar pressure in kg/km/s**2
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 13.03.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!              </ul>
!!
!> @param[in]   time_mjd   Current time (MJD)
!!
!> @details     This function computes the solar irradiance as a function of
!!              the Earth's position along its orbit. Then, solar pressure is computed
!!              by dividing it by the speed of light. The formula as given by
!!              Wertz, "Spacecraft Attitude Determination and Control", 1978, is used.
!!              The aphelion dates of the Earth's orbit were taken from the USNO website:
!!              http://www.usno.navy.mil/USNO/astronomical-applications/data-services/earth-seasons
!!              and are available for the years 2000 through 2020. If the passed MJD is
!!              earlier than 2000 or beyond 2020 then the aphelion is assumed to be repeated every 365.2524 days.
!!
!! @todo        Test function 'getSolarPressure' for time-varying solar 'constant'
!!------------------------------------------------------------------------------------------------

  real(dp) function getSolarPressure(this, time_mjd)     ! <-- DBL  current time MJD

    !** interface
    !----------------------------------------------
    class(Radiation_class)  :: this
    real(dp), intent(in)    :: time_mjd
    !----------------------------------------------

    character(len=*), parameter :: csubid = "getSolarPressure"

    integer :: nyears            ! number of years between time_mjd and first/last aphelion date

    real(dp)            :: ap_end               ! coming aphelion date
    real(dp)            :: ap_start             ! last aphelion date
    real(dp)            :: cosd                 ! cosine of the distance from the aphelion of Earth's orbit.
    real(dp), parameter :: fac = 4.5298d-3      ! factor resulting from dividing 1358 W/m**2 by the speed of light and converting it to kg/s**2/km

    getSolarPressure = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** if already available for given epoch, then return that value immediately...
    if(time_mjd == this%last_mjd) then

      getSolarPressure = this%last_press

    else

      !** find index in 'aphelions' array
      if(time_mjd < aphelions(1)) then

        nyears   = int((aphelions(1) - time_mjd)/length_of_year) + 1
        ap_start = aphelions(1) - nyears*length_of_year
        ap_end   = aphelions(1) - (nyears - 1)*length_of_year

      else if(time_mjd > aphelions(size(aphelions))) then

        nyears   = int((time_mjd - aphelions(size(aphelions)))/length_of_year) + 1
        ap_start = aphelions(size(aphelions)) + (nyears - 1)*length_of_year
        ap_end   = aphelions(size(aphelions)) + nyears*length_of_year

      else

        this%apidx = 2

        do

          if(aphelions(this%apidx-1) >= time_mjd .and. aphelions(this%apidx) < time_mjd) then
            this%apidx = this%apidx - 1
            exit
          end if

          this%apidx = this%apidx + 1
          if(this%apidx > size(aphelions)) then
            this%apidx = 1 ! set to '1' for each epoch not within defined range
            exit
          end if

        end do

        ap_start = aphelions(this%apidx)
        ap_end   = aphelions(this%apidx + 1)

      end if

      cosd = cos(twopi*(time_mjd - ap_start)/(ap_end - ap_start))

      getSolarPressure = fac/(1.004d0 + 3.34d-2*cosd)

      this%last_press = getSolarPressure   ! save for multiple calls
      this%last_mjd   = time_mjd

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSolarPressure

!================================================================================
!
!> @anchor      getAlbedoAcceleration
!!
!> @brief       Return acceleration due to Earth radiation pressure
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 30.01.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   r_gcrf      Radius vector of spacecraft in GCRF
!> @param[in]   v_gcrf      Velocity vector of spacecraft in GCRF
!> @param[in]   time_mjd    Current MJD
!!
!> @param[out]  acc_srp     Acceleration vector in inertial frame
!!
!> @details     This routine computes the acceleration due to Earth radiation pressure.
!!              The model is based on Knocke(1989). The position of the sun is computed
!!              using JPL ephemerides.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Finish.                  </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------

  subroutine getAlbedoAcceleration(   this,             &
                                      satellite_model,  &  ! <-> TYPE   Satellite model
                                      solarsystem_model,&  ! <-> TYPE   Solarsystem model
                                      reduction,        &  ! <-> TYPE   Reduction
                                      r_gcrf,           &  ! <-- DBL(3) GCRF radius of spacecraft
                                      v_gcrf,           &  ! <-- DBL(3) GCRF velocity of spacecraft
                                      time_mjd,         &  ! <-- DBL    current MJD
                                      acc_alb           &  ! --> DBL(3) acceleration in inertial frame
                                  )

    !** interface
    !-------------------------------------------------
    class(Radiation_class)              :: this
    type(Satellite_class)               :: satellite_model
    type(Solarsystem_class)             :: solarsystem_model
    type(Reduction_type)               :: reduction
    real(dp), dimension(3), intent(in)  :: r_gcrf
    real(dp), dimension(3), intent(in)  :: v_gcrf
    real(dp),               intent(in)  :: time_mjd

    real(dp), dimension(3), intent(out) :: acc_alb
    !-------------------------------------------------

    character(len=*), parameter :: csubid = "getAlbedoAcceleration"

    integer :: i,j                                                              ! loop counter
    integer :: k                                                                ! iteration counter
    integer :: nsegments                                                        ! number of segments in a ring
    integer :: satsurf                                                          ! holding the number of satellite surfaces

    real(dp)                                :: albedo                           ! albedo of an Earth plate
    real(dp)                                :: albFlux                          ! quantity, computed via Es/c/m, where Es is the sun irradiance, c is the speed of light and m is the satellite's mass
    real(dp), dimension(nrings+1)           :: beta                             ! array containing the angles for the Earth rings as seen from the Earth
    real(dp), dimension(nrings+1)           :: betac                            ! array containing the centers for the Earth rings as seen from the Earth
    real(dp)                                :: cosinc                           ! cosine of incidence angle between sun and surface normal
    real(dp)                                :: cosjd                            ! factor, which takes into account the annual albedo variation
    real(dp)                                :: costheta                         ! cosine of incidence angle of the sun wrt. Earth plate
    real(dp)                                :: coszeta                          ! cosine of zeta(1)
    real(dp)                                :: emiss                            ! emissivity of an Earth plate
    real(dp)                                :: emissFlux                        ! quantity, computed via Mb/c/m, where Mb is the exitance and equal to Es/4 (black body), c is the speed of light and m is the satellite's mass
    real(dp), dimension(nrings)             :: gama                             ! array containing auxiliary angles
    real(dp)                                :: lambda                           ! azimuth of Earth plate relative to transverse component (positive counterclockwise)
    real(dp)                                :: lat                              ! sine of latitude
    real(dp)                                :: plateArea                        ! attenuated weighted area of each Earth plate
    real(dp)                                :: rabs                             ! absolute value for r_gcrf
    real(dp)                                :: re                               ! Earth's radius
    real(dp), dimension(3)                  :: rseg                             ! radius vector of segment in GCRF
    real(dp), dimension(3)                  :: rseg_uvw                         ! radius vector of segment in UVW frame
    real(dp), dimension(3)                  :: rss                              ! vector from segment to satellite surface
    real(dp), dimension(3)                  :: rsun                             ! sun vector
    real(dp)                                :: segmWidth                        ! angular width of segments in a ring
    real(dp)                                :: solarPressure                    ! radiation pressure (kg/km/s**2)
    real(dp), dimension(3)                  :: tempvec                          ! auxiliary vector
    real(dp)                                :: tmp, tmp2                        ! auxiliary
    real(dp), dimension(3)                  :: unit_r                           ! unit vector of r_gcrf
    real(dp), dimension(3)                  :: unit_rseg                        ! unit vector of rseg in GCRF
    real(dp), dimension(3)                  :: unit_rseg_uvw                    ! unit vector of rseg in UVW frame
    real(dp), dimension(3)                  :: unit_rss                         ! unit vector of rss
    real(dp), dimension(3)                  :: unit_rsun                        ! unit vector of rsun
    real(dp), dimension(MAX_SURFACES,3)     :: vecn                             ! normal vector for each satellite surface
    real(dp), dimension(nrings)             :: zeta                             ! array containing the angles for the Earth rings as seen from the satellite
    real(dp)                                :: zetaM                            ! maximum surface visibility angle (rad)

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

!    write(*,*) "computing albedo...", time_mjd

    rabs    = mag(r_gcrf)
    unit_r  = r_gcrf/rabs
    re      = getEarthRadius()
    satsurf = satellite_model%getSurfaceNumber()                                ! number of satellite surfaces
    if(hasFailed()) return

    cosjd  = cos(twopi/365.25d0*(time_mjd - 44960.d0))                          ! for start epoch at Dec 22, 1981 (Knocke, 1989)

    !===========================================
    !
    ! Get Sun's position
    !
    !-----------------------------------------
    rsun      = solarsystem_model%getBodyPosition(time_mjd, ID_SUN)
    unit_rsun = rsun/mag(rsun)

    !===========================================
    !
    ! Get radiation pressure
    !
    !-----------------------------------------
    solarPressure = this%getSolarPressure(time_mjd)
    if(hasFailed()) return

    !====================================================
    !
    ! Retrieve satellite properties (the mass is
    ! only required once and as soon as it changes...)
    !
    !--------------------------------------------------
    if(satellite_model%hasChangedSatellite(ID_ALBEDO) .or. this%first_albedo)  then

      this%mass_albedo = satellite_model%getObjectMass()
      if(hasFailed()) return

      do i=1,satsurf

        this%crs_albedo(i)  = satellite_model%getSurfaceSpecRefl(i)
        if(hasFailed()) return

        this%crd_albedo(i)  = satellite_model%getSurfaceDiffRefl(i)
        if(hasFailed()) return

        this%area_albedo(i) = satellite_model%getSurfaceArea(i)
        if(hasFailed()) return

      end do

      this%first_albedo = .false.
      this%crd_albedo(:) = this%crd_albedo(:)/3.d0  ! only this factor will be required!

    end if

    !** compute combined quantities
    albFlux   = solarPressure/this%mass_albedo
    emissFlux = 0.25d0*albFlux      ! black body assumption

    !** maximum Earth visibility angle zetaM (from satellite)
    !write(*,*) "ZETAM"
    if(re > rabs) then
      zetaM = 1.d0
    else
      zetaM = asin(re/rabs)
    end if

    !** maximum central angle (angle between satellite radius and Earth radius pointing to tangent point)
    beta(nrings+1) = halfPi - zetaM

    !=====================================================
    !
    ! Compute central cap
    !
    !---------------------------------------------
    coszeta = (this%ringElements - 1 + cos(zetaM))/this%ringElements

    zeta(1) = acos(coszeta)
    !write(*,*) "gama1"
    gama(1) = asin(rabs*sin(zeta(1))/re)
    !write(*,*) "ok"
    beta(1) = gama(1) - zeta(1)


    !====================================================
    !
    ! Compute zeta and beta angles for other rings
    !
    !---------------------------------------------
    k = 0

    do i=2,nrings

      k       = k + 6*(i - 1)
      tmp     = k*coszeta - k + 1
      zeta(i) = acos(tmp)
      !write(*,*) "gama", i
      gama(i) = asin(rabs*sin(zeta(i))/re)
      !write(*,*) "ok"
      beta(i) = gama(i) - zeta(i)

    end do

    !write(*,*) "gama = ", gama
    !write(*,*) "zeta = ", zeta

    !====================================================
    !
    ! Compute ring centers
    !
    !----------------------------------------------
    do i=2,nrings+1

      betac(i) = 0.5d0*(beta(i) + beta(i-1))

    end do

    !===============================================================
    !
    ! Compute central cap contribution
    !
    !-------------------------------------------------

    !** weighted plate area
    plateArea = 2.d0*(1.d0 - coszeta)

    !** angle between satellite and sun as seen from Earth's center
    costheta = dot_product(unit_r,unit_rsun)

    if(costheta < 0.d0) costheta = 0.d0    ! central cap is in shadow as soon as sun is perpendicular to satellite's radius

    !** sine of latitude of central cap
    lat = unit_r(3)

    tmp  = cosjd*lat
    tmp2 = (3.d0*lat*lat - 1.d0)

    !** albedo and emissivity
    albedo = 0.34d0 + 0.10d0*tmp + 0.145d0*tmp2
    emiss  = 0.68d0 - 0.07d0*tmp - 0.180d0*tmp2

    !** now compute central cap contribution...
    if(satellite_model%getOrientationMode() == MODE_CANNON_BALL) then

      acc_alb(1:3) = this%area_albedo(1)*((1.d0 + this%crs_albedo(1))*unit_r(1:3) + 2.d0*this%crd_albedo(1)*unit_r(1:3))

    else

      !** accumulate contributions due to central cap
      acc_alb = 0.d0

      do i=1,satsurf

        vecn(i,:) = satellite_model%getSurfaceNormal(solarsystem_model, reduction, i, time_mjd, r_gcrf, v_gcrf)
!       write(*,*) "normal vector = ", vecn(i,:)
!       write(*,*) "unit_r = ", unit_r
        if(hasFailed()) return

        cosinc  = dot_product(vecn(i,:),-unit_r)

        if(cosinc < 0.d0) cosinc = 0.d0 ! surface not seeing central cap

!       write(*,*) "cosinc = ", cosinc
!       write(*,*) "area   = ", this%area_albedo(i)

        acc_alb(:) = acc_alb(:) + cosinc*this%area_albedo(i)*(2.d0*(this%crd_albedo(i) + this%crs_albedo(i)*cosinc)*vecn(i,:) + (1.d0 - this%crs_albedo(i))*unit_r(:))

      end do

    end if

    acc_alb(:) = acc_alb(:)*(albedo*albFlux*costheta + emiss*emissFlux)*plateArea

!    write(*,*) "albflux = ", albFlux, albedo
!    write(*,*) "emissFlux = ", emissFlux, emiss
!    write(*,*) "plateArea = ", plateArea
!     write(*,*) "acc_alb = ", acc_alb
     !read(*,*)

    !==============================================================
    !
    ! Compute ring contributions
    !
    !--------------------------------------------------

    !** do for each ring
    do i=2,nrings+1

      nsegments = 6*(i-1)
      segmWidth = twopi/nsegments

      !** do for each ring segment
      do j=1,6

        lambda = (j-1)*segmWidth

        !** compute radius vector for segment in UVW system
        unit_rseg_uvw(1) = cos(betac(i))
        unit_rseg_uvw(2) = sin(betac(i))*cos(lambda)
        unit_rseg_uvw(3) = sin(betac(i))*sin(lambda)

        !write(*,*) "betac = ", betac
        !write(*,*) "lambda = ", lambda
        rseg_uvw = unit_rseg_uvw*re

        !** convert to inertial frame
        call reduction%uvw2eci(r_gcrf, v_gcrf, rseg_uvw, rseg)

        unit_rseg = rseg/mag(rseg)

        !** compute direction from segment to satellite
        rss      = r_gcrf - rseg
        unit_rss = rss/mag(rss)

        !** sine of latitude
        lat = unit_rseg(3)

        tmp  = cosjd*lat
        tmp2 = (3.d0*lat*lat - 1.d0)

        !** albedo and emissivity
        albedo = 0.34d0 + 0.10d0*tmp + 0.145d0*tmp2
        emiss  = 0.68d0 - 0.07d0*tmp - 0.180d0*tmp2

        !** angle between satellite and sun as seen from Earth's center
        costheta = dot_product(unit_rseg,unit_rsun)

        if(costheta < 0.d0) costheta = 0.d0    ! segment is in shadow as soon as sun is perpendicular to radius vector


        !** contribution of segment to each satellite surface
        if(satellite_model%getOrientationMode() == MODE_CANNON_BALL) then

          tempvec(:) = this%area_albedo(1)*((1.d0 + this%crs_albedo(1))*unit_rss(1:3) + 2.d0*this%crd_albedo(1)*unit_rss(1:3))

        else

          tempvec = 0.d0

          do k = 1, satsurf

            cosinc     = dot_product(vecn(i,:),-unit_rss)
            if(cosinc < 0.d0) cosinc = 0.d0   ! surface not seeing element

            tempvec(:) = tempvec(:) + cosinc*this%area_albedo(i)*(2.d0*(this%crd_albedo(i) + this%crs_albedo(i)*cosinc)*vecn(i,:) + (1.d0 - this%crs_albedo(i))*unit_rss(:))

          end do

        end if

        acc_alb(:) = acc_alb(:) + tempvec(:)*(albedo*albFlux*costheta + emiss*emissFlux)*plateArea

      end do

    end do

    !write(*,*) "acc_alb = ", acc_alb, albedo, albFlux, costheta, emiss, emissFlux, plateArea, this%area_albedo(1), this%crd_albedo(1), this%crs_albedo(1)
!    write(*,*) "done.", acc_alb
    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine getAlbedoAcceleration

!========================================================================
!
!> @anchor      checkBoundaryCrossing
!!
!> @brief       Check whether there have been a shadow boundary crossing within the last steps
!> @author      Andre Horstmann, Vitali Braun
!!
!> @date        <ul>
!!                <li> 11.2013 (initial design)             </li>
!!                <li> 14.03.2014 (migrated to radiation module and did some optimization)
!!              </ul>
!!
!! @detail      This subroutine evaluates several geometric parameters
!!              to calculate the function T(t) that is defined in the
!!              SRP correction algorithm. It then checks whether there was a boundary
!!              crossing or not for the passed position vectors.
!!
!------------------------------------------------------------------------
  subroutine checkBoundaryCrossing(this, solarsystem_model, pos, time_mjd, steps, bound, T, fshadow, &
                                   theta, alpha_e, alpha_s)

    ! interface
    !----------------------------------------------------------------------------
    class(Radiation_class)                :: this
    type(Solarsystem_class),intent(inout) :: solarsystem_model
    real(dp), dimension(3,3), intent(in)  :: pos         ! current (=1) and 2 prior positions (2,3)
    real(dp),                 intent(in)  :: time_mjd    ! current time / MJD
    real(dp), dimension(:),   intent(in)  :: steps       ! step history
    logical,                  intent(out) :: bound       ! flag indicating that boundary has been crossed
    real(dp), dimension(2,3), intent(out) :: T           ! T-function
    real(dp), dimension(3),   intent(out) :: fshadow     ! shadow function
    real(dp), dimension(3),   intent(out) :: theta       ! Angle Theta
    real(dp), dimension(3),   intent(out) :: alpha_e     ! Angle Alpha earth
    real(dp), dimension(3),   intent(out) :: alpha_s     ! Angle Alpha
   !-----------------------------------------------------------------------------

    real(dp), dimension(3,3) :: r_e                             ! satellite -> earth
    real(dp)                 :: Rearth                          ! radius earth
    real(dp)                 :: Rsun                            ! radius sun
    real(dp), dimension(3,3) :: r_s                             ! satellite -> sun
    real(dp), dimension(3,3) :: r_sun                           ! earth -> sun

    integer :: i    ! loop counter

    Rearth = getEarthRadius()                                   ! constant
    Rsun   = getSunRadius()                                     ! constant

    r_e(:,:)   = -pos(:,:)                                ! satellite -> earth
    r_sun(:,1) = solarsystem_model%getBodyPosition(time_mjd,                                   ID_SUN)   ! earth -> sun
    r_sun(:,2) = solarsystem_model%getBodyPosition(time_mjd - steps(9)/sec_per_day,            ID_SUN)
    r_sun(:,3) = solarsystem_model%getBodyPosition(time_mjd - (steps(9)+steps(8))/sec_per_day, ID_SUN)

    r_s(:,:)   = r_sun(:,:) - pos(:,:) ! satellite -> sun

    forall(i = 1:3)
      theta(i)   = acos(dot_product(r_e(:,i),r_s(:,i))/(mag(r_e(:,i))*mag(r_s(:,i))))
      alpha_e(i) = asin(Rearth/mag(r_e(:,i)))
      alpha_s(i) = asin(Rsun/mag(r_s(:,i)))
    end forall

    do i=1,3
      fshadow(i) = this%getShadow(r_sun(:,i), pos(:,i))
      if(hasFailed()) return
    end do

    ! calculate values of T

    T(1,:) = theta(:) - alpha_e(:) - alpha_s(:)
    T(2,:) = theta(:) - alpha_e(:) + alpha_s(:)

    if ((theta(1) > alpha_e(1) + alpha_s(1)) .and. .not.(theta(3) > alpha_e(3) + alpha_s(3))) then
      bound =.true.
    elseif ((alpha_e(1) + alpha_s(1) > theta(1)) .and. &
            (theta(1) > alpha_e(1) - alpha_s(1)) .and. &
            .not. ((alpha_e(3) + alpha_s(3) > theta(3)) .and.  &
                   (theta(3) > alpha_e(3) - alpha_s(3)))) then
      bound = .true.
    elseif ((theta(1) < alpha_e(1) - alpha_s(1)) .and. .not.(theta(3) < alpha_e(3) - alpha_s(3))) then
      bound = .true.
    else
      bound = .false.
    end if

  end subroutine checkBoundaryCrossing


!================================================================================
!
!> @anchor      getSrpAcceleration
!!
!> @brief       Return acceleration due to solar radiation pressure
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.03.2013 (initial design)</li>
!!                <li> 04.06.2013 (code optimization) </li>
!!                <li> 28.01.2013 (added velocity as input parameter, in order to get oriented surfaces working.)</li>
!!                <li> 10.06.2014 (correction factor to account for distances <> 1 AU, as previously, the solar
!!                                 constant was computed for the distance 1 AU.)</li>
!!                <li> 13.02.2015 (fixed bug: correction factor (see above) had to be inverted...)</li>
!!                <li> 30.03.2015 (added new feature: ESH computation for oriented surfaces)</li>
!!              </ul>
!!
!> @param[in]   r_gcrf      Radius vector of spacecraft in GCRF
!> @param[in]   v_gcrf      Velocity vector of spacecraft in GCRF
!> @param[in]   time_mjd    Current MJD
!!
!> @param[out]  acc_srp     Acceleration vector in inertial frame
!!
!> @details     This routine computes the acceleration due to solar radiation pressure.
!!              A conic shadow model is assumed. The position of the sun is computed
!!              using JPL ephemerides.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Get Sun position         </li>
!!                <li> Compute shadow factor    </li>
!!                <li> Compute SRP acceleration </li>
!!                <li> Finish.                  </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------

  subroutine getSrpAcceleration(  this,             &
                                  satellite_model,  &  ! <-> TYPE Satellite model
                                  solarsystem_model,&  ! <-> TYPE Solarsystem model
                                  reduction,        &  ! <-> TYPE Reduction
                                  r_gcrf,           &  ! <-- DBL(3) GCRF radius of spacecraft
                                  v_gcrf,           &  ! <-- DBL(3) GCRF velocity of spacecraft
                                  time_mjd,         &  ! <-- DBL    current MJD
                                  acc_srp           &  ! --> DBL(3) acceleration in inertial frame
                               )

    !** interface
    !-------------------------------------------------
    class(Radiation_class)              :: this
    type(Satellite_class)               :: satellite_model
    type(Solarsystem_class)             :: solarsystem_model
    type(Reduction_type)               :: reduction
    real(dp), dimension(3), intent(in)  :: r_gcrf
    real(dp), dimension(3), intent(in)  :: v_gcrf
    real(dp),               intent(in)  :: time_mjd

    real(dp), dimension(3), intent(out) :: acc_srp
    !-------------------------------------------------

    character(len=*), parameter :: csubid = "getSrpAcceleration"

    integer  :: i                             ! loop counter
    integer  :: surf                          ! surface number

    real(dp) :: cosinc   ! cosine of incidence angle between sun and surface normal
    real(dp) :: fshadow                       ! shadow factor:
                                              ! 1 = no shadow
                                              ! 0 = umbra
                                              ! 0 < fshadow < 1 = penumbra
    real(dp), dimension(3) :: r_sun           ! Sun position vector in GCRF
    real(dp), dimension(3) :: ros             ! radius vector object -> sun
    real(dp)               :: ros_abs         ! magnitude of ros
    real(dp)               :: solarPressure   ! solar pressure in N/m**2
    real(dp)               :: srCorrect       ! correction factor for solar constant, the latter being determined for 1 AU
    real(dp), dimension(3) :: vecn            ! normal vector of surface

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check if already available (in those cases where called twice...)
    if((abs(time_mjd - this%last_srp_mjd) < eps9) .and. (mag(r_gcrf - this%last_srp_r) < eps9)) then

      acc_srp = this%last_srp_acc

      if(isControlled()) then
        call checkOut(csubid)
      end if
      return

    end if

    !=====================================================
    !
    ! Get Sun position
    !
    !----------------------------

    r_sun = solarsystem_model%getBodyPosition(time_mjd, ID_SUN)

    !=====================================================
    !
    ! Compute shadow factor
    !
    !----------------------------

    if(this%shadowMode == NO_SHADOW) then
      fshadow = 1.d0
    else
      fshadow = this%getShadow(     &
                            r_sun,  &  ! <-- DBL(3)  GCRF radius of sun
                            r_gcrf  &  ! <-- DBL(3)  GCRF radius of object
                         )
      if(hasFailed()) return
    end if

    !=====================================================
    !
    ! Get solar irradiance
    !
    !------------------------------------

    !solarPressure = getSolarPressure(time_mjd)  this function is not required anymore, as a correction factor is used
    !                                            for the 1 ua - based solar constant

    solarPressure = 4.5391d-3 ! constant for 1 ua in kg/s**2/km

    if(hasFailed()) return

    !=====================================================
    !
    ! Compute SRP acceleration
    !
    !----------------------------

    surf = satellite_model%getSurfaceNumber()
    if(hasFailed()) return

    !** vector object -> sun
    ros     = r_sun - r_gcrf
    ros_abs = mag(ros)
    ros     = ros/ros_abs   ! unit vector

    !** compute correction for solar radiation constant
    srCorrect = (getAstronomicalUnit()/ros_abs)**2.d0

    if(satellite_model%hasChangedSatellite(ID_SRP) .or. this%first_srp)  then    ! only required, if configuration has changed or called for first time

      this%mass_srp = satellite_model%getObjectMass()
      if(hasFailed()) return

      do i=1,surf

        this%crsSRP(i)  = satellite_model%getSurfaceSpecRefl(i)
        if(hasFailed()) return

        this%crdSRP(i)  = satellite_model%getSurfaceDiffRefl(i)
        if(hasFailed()) return

        this%areaSRP(i) = satellite_model%getSurfaceArea(i)
        if(hasFailed()) return

      end do

      this%first_srp = .false.
      this%crdSRP(:) = this%crdSRP(:)/3.d0  ! only this factor will be required!

    end if

    if(fshadow > eps9) then

      if(satellite_model%getOrientationMode() == MODE_CANNON_BALL) then

        acc_srp(1:3) = this%areaSRP(1)*((1.d0 + this%crsSRP(1))*ros(1:3) + 2.d0*this%crdSRP(1)*ros(1:3))

      else

        acc_srp = 0.d0

        do i=1,surf

          !** get normal vector in GCRF
          vecn   = satellite_model%getSurfaceNormal(solarsystem_model, reduction, i, time_mjd, r_gcrf, v_gcrf)

          if(hasFailed()) return
          cosinc = dot_product(vecn,ros)

          if(cosinc <= 0.d0) cycle ! skip surface, as it is in shadow

          !** sum up equivalent solar hours
          if(satellite_model%getESHFlag()) call satellite_model%updateESH(time_mjd, i, cosinc)

          acc_srp(1:3) = acc_srp(1:3) + this%areaSRP(i)*cosinc*(2.d0*(this%crdSRP(i) + this%crsSRP(i)*cosinc)*vecn(1:3) + (1.d0 - this%crsSRP(i))*ros(1:3))

        end do

      end if

      acc_srp(1:3) = -fshadow*acc_srp(1:3)*solarPressure/this%mass_srp*srCorrect

    else

      acc_srp = 0.d0

    end if

   !write(*,*) "acc_srp = "
   !write(*,*) fshadow, "fshadow"
   !write(*,*) solarPressure, "press"
   !write(*,*) mass, "this%mass_srp"
   !write(*,*) srCorrect, "srCorrect"
   !write(*,*) acc_srp, "acc_srp"

    !** store data in module variable
    this%last_srp_acc = acc_srp
    this%last_srp_mjd = time_mjd
    this%last_srp_r   = r_gcrf

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine getSrpAcceleration

end module radiation
