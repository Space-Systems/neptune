!=================================================================================================
!!
!> @brief       Solid and ocean tides modeling
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!> @author      Andrea Turchi (ATU)
!!
!> @date        <ul>
!!                <li>VB:  16.07.2013 (initial design)    </li>
!!                <li>CHK: 13.11.2015 (updated to use with libslam)    </li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 04.01.2018 (Using Solarsystem_class)</li>
!!                <li>CHK: 05.01.2017 (created Tides_class)
!!                <li>ATU: 01.09.2022 (added fes2004 model)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for Earth's
!!              solid and ocean tides modeling, specifically in the context of numerical
!!              integration of a satellite trajectory.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      tides
!!
!!------------------------------------------------------------------------------------------------
module tides

  use slam_astro,             only: getEarthRadius, getEarthGravity, getEarthMass, getGMST
  use slam_astro_conversions, only: getGeodeticLatLon, getRadiusLatLon
  use neptune_error_handling, only: E_SOLAR_SYSTEM_INIT, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, checkIn, checkOut
  use slam_math,              only: mag, pi, eps9, rad2deg, factorial, deg2rad, as2deg
  use slam_reduction_class,   only: Reduction_type
  use solarsystem,            only: Solarsystem_class, ID_SUN, ID_MOON
  use slam_types,             only: dp
  use slam_io,                only: openFile, closeFile, SEQUENTIAL, IN_FORMATTED, cdelimit, message, LOG_AND_STDOUT
  use slam_time,              only: jd2000, jd245, julian_century

  implicit none

  private

    integer, parameter, public          :: SOLID_TIDES = 1
    integer, parameter, public          :: OCEAN_TIDES = 2

    type, public :: Tides_class

        logical, public :: tidesInitialized
        real(dp), dimension(1:18,2:6,0:6), public :: dC_p      ! 3-D array of the prograde C coefficients of all tides
        real(dp), dimension(1:18,2:6,0:6), public :: dS_p      ! 3-D array of the prograde S coefficients of all tides
        real(dp), dimension(1:18,2:6,0:6), public :: dC_m      ! 3-D array of the retrograde C coefficients of all tides
        real(dp), dimension(1:18,2:6,0:6), public :: dS_m      ! 3-D array of the retrograde C coefficients of all tides
        real(dp), dimension(1:18), public         :: doodson   ! vector of the Doodson numbers of all tides
        integer, public                           :: max_l     ! maximum degree of harmonic coefficients

        character(len=255) :: dataPath                         ! Path where input data files are located
        character(len=255) :: fesDataFile                      ! Tides constituents data file
        
    contains

        procedure :: initTides
        procedure :: getTidesAcceleration
        procedure :: setTidesInitFlag
        procedure :: get_Delaunay_arg
        procedure :: get_Doodson_arg
        procedure :: init_FES2004
        procedure :: get_FES2004_corrections

    end type Tides_class

    ! Constructor
    interface Tides_class
        module procedure constructor
    end interface Tides_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !!  @author     Christopher Kebschull
    !!  @date       <ul>
    !!                  <li>ChK: 05.01.2018 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Tides_class) function constructor()

        constructor%dataPath         = "data"                 ! Path where input data files are located
        constructor%fesDataFile      = "fes2004_Cnm-Snm.dat"  ! Tides constituents data file

        ! other variables
        constructor%max_l            = 6                      ! maximum degree of harmonic coefficients

        constructor%tidesInitialized = .false.

    end function constructor

  !==============================================================================
  !
  !> @anchor      initTides
  !!
  !> @brief       Initialization of tides module
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 16.07.2013 (initial design)    </li>
  !!              </ul>
  !!
  !!-----------------------------------------------------------------------------
  subroutine initTides( this, &
                        cpath)

    class(Tides_class)            :: this
    character(len=*), intent(in)  :: cpath

    character(len=*), parameter   :: csubid = "initTides"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%tidesInitialized) then
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
    end if

    call message(' - Initializing tides model...', LOG_AND_STDOUT)

    !** set data path in module scope
    this%dataPath = trim(adjustl(cpath(1:min(len(this%dataPath),len(cpath)))))

    call init_FES2004(this)
    !update the initialization status
    this%tidesInitialized = .true.

    !** in geopotential module, activate function to provide computation parameters
    !call setTideSupport(.true.)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine initTides

  !============================================================================
  !
  !> @anchor      getTidesAcceleration
  !!
  !> @brief       Providing accelerations due to ocean and solid Earth tides
  !> @author      Vitali Braun
  !!
  !> @param[in]   solarsystem_model   the solar system model
  !> @param[in]   reduction           reduction handler
  !> @param[in]   r_gcrf              position vector (gcrf)
  !> @param[in]   r_itrf              position vector (itrf)
  !> @param[in]   v_itrf              velocity vector (itrf)
  !> @param[in]   time_mjd            current MJD
  !> @param[in]   tidetype            1 = solid; 2 = ocean
  !> @param[out]  accel               acceleration vector in inertial frame
  !!
  !> @date        <ul>
  !!                <li> 16.07.2013 (initial design) </li>
  !!                <li> 11.02.2015 (changed to a stable method for Legendre functions recursion)</li>
  !!                <li> 12.02.2015 (corrected ocean tides and added solid earth pole tides)</li>
  !!              </ul>
  !!
  !> @details     This routine computes the acceleration due to Earth's ocean
  !!              and solid tides in an inertial frame (ECI).
  !!
  !!              \par Overview
  !!
  !!              <ol>
  !!                <li> Compute Legendre polynomials </li>
  !!                <li> Determine partial derivatives of disturbing potential </li>
  !!                <li> Compute acceleration in GCRF </li>
  !!                <li> Finish. </li>
  !!              </ol>
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine getTidesAcceleration(                          &
                                        this,               &
                                        solarsystem_model,  &
                                        reduction,          &
                                        r_gcrf,             &  ! <-- DBL(3) radius vector in GCRF frame
                                        r_itrf,             &  ! <-- DBL(3) radius vector in ITRF frame
                                        v_itrf,             &  ! <-- DBL(3) velocity vector in ITRF frame
                                        time_mjd,           &  ! <-- DBL    current MJD
                                        tidetype,           &  ! <-- INT    tide type (1=solid, 2=ocean)
                                        accel               &  ! --> DBL(3) acceleration vector in inertial frame
                                      )


    !** interface
    !----------------------------------------------
    class(Tides_class)                      :: this
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)      :: reduction
    real(dp), dimension(3), intent(in)      :: r_itrf, r_gcrf
    real(dp), dimension(3), intent(in)      :: v_itrf
    real(dp),               intent(in)      :: time_mjd
    integer,                intent(in)      :: tidetype


    real(dp), dimension(3), intent(out) :: accel
    !----------------------------------------------

    character(len=*), parameter :: csubid = "getTidesAcceleration" ! subroutine id

    integer, parameter     :: SUN  = 1
    integer, parameter     :: MOON = 2
    integer                :: dm              ! coefficient in transformation between normalized and unnormalized quantities
    integer                :: i,l,m
    integer                :: lmax            ! maximum order of potential

    real(dp), dimension(2)         :: body_lat
    real(dp), dimension(2)         :: body_lon
    real(dp), dimension(0:6)       :: costerm, sinterm                            ! cos(ml) and sin(ml)
    real(dp), dimension(2:6,0:6)   :: dC, dS                                      ! deviations in coefficients C and S due to tides
    real(dp), dimension(2:3,0:3)   :: k = reshape((/0.29525d0, 0.093d0, &         ! nominal love numbers for degree and order
                                                    0.29470d0, 0.093d0, &
                                                    0.29801d0, 0.093d0, &
                                                    0.d0,      0.094d0/), (/2,4/))
    real(dp), dimension(0:2)       :: kp =(/-0.00087d0, -0.00079d0, -0.00057d0/)  ! nominal love numbers for correction of fourth degree term (k+nm)
    real(dp), dimension(2,0:6,0:6) :: lp                                          ! legendre polynomials for Sun and Moon
    real(dp), dimension(0:7,0:7)   :: lpsat                                       ! legendre polynomials for satellite
    real(dp), dimension(2)         :: mu                                          ! GM parameter of Sun and Moon
    real(dp), dimension(2)         :: pom, pomAvg                                 ! polar motion variables xp,yp (current and running average) / rad
    real(dp), dimension(2)         :: rabs_body
    real(dp), dimension(3,2)       :: rBodyGCRF
    real(dp), dimension(3,2)       :: rBodyITRF

    real(dp) :: const                   ! constant factor
    real(dp) :: dudlambda               ! dU/d(lambda)
    real(dp) :: dudphi                  ! dU/d(phi_gc)
    real(dp) :: dudr                    ! dU/d(rabs)
    real(dp) :: fac                     ! factorial term for the conversion to unnormalized coefficients
    real(dp) :: insig1, insig2, insig3  ! cumulated sums
    real(dp) :: lambda                  ! longitude
    real(dp) :: m1, m2                  ! auxiliaries to account for pole tide
    real(dp) :: muEarth                 ! Earth's gravity constant
    real(dp) :: oorabs                  ! 1/rabs
    real(dp) :: oorabs2                 ! 1/rabs2
    real(dp) :: oorabs3                 ! 1/rabs3
    real(dp) :: oosqrt_r1r2             ! 1/sqrt(r1r2)
    real(dp) :: phi_gc                  ! latitude geocentric
    real(dp) :: r1r2                    ! radius(1)**2 + radius(2)**2
    real(dp) :: rabs                    ! magnitude of radius vector
    real(dp) :: rabs2                   ! squared radius magnitude
    real(dp) :: rabs3                   ! cubed radius magnitude
    real(dp) :: rekm                    ! Earth's radius in km
    real(dp) :: rrfac                   ! temporary
    real(dp) :: sqrt_r1r2               ! sqrt(r1r2)
    real(dp) :: tanphi                  ! tan(phi_gc)
    real(dp) :: temp, templ, temp2      ! temporary
    real(dp) :: temp_t1                 ! temporary (used to support tides)
    real(dp) :: temp_t2                 ! temporary (used to support tides)
    real(dp) :: temp_t3                 ! temporary (used to support tides)
    real(dp) :: temp_t4                 ! temporary (used to support tides)
    real(dp) :: temp_t5                 ! temporary (used to support tides)

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !=========================================================================
    !
    ! Get position of Sun and Moon in ITRF (latitude and longitude)
    !
    !-------------------------------------------------------------------------
    if(.not. solarsystem_model%getSolarSystemInitFlag()) then
      call setNeptuneError(E_SOLAR_SYSTEM_INIT, FATAL)
      return
    end if

    !** set maximum order
    select case(tidetype)
      case(SOLID_TIDES)
        lmax = 4
      case(OCEAN_TIDES)
        lmax = 6
    end select

    !** position in GCRF
    rBodyGCRF(:,SUN)  = solarsystem_model%getBodyPosition(time_mjd, ID_SUN)
    if(hasFailed()) return
    rBodyGCRF(:,MOON) = solarsystem_model%getBodyPosition(time_mjd, ID_MOON)
    if(hasFailed()) return

    rabs_body(SUN)  = mag(rBodyGCRF(:,SUN))
    rabs_body(MOON) = mag(rBodyGCRF(:,MOON))

    !** convert to ITRF
    call reduction%inertial2earthFixed(rBodyGCRF(:,SUN),  time_mjd, rBodyITRF(:,SUN))
    if(hasFailed()) return
    call reduction%inertial2earthFixed(rBodyGCRF(:,MOON), time_mjd, rBodyITRF(:,MOON))
    if(hasFailed()) return

    !** compute longitude and latitude
    call getGeodeticLatLon(rBodyITRF(:,SUN),  temp, body_lat(SUN),  body_lon(SUN))
    if(hasFailed()) return
    call getGeodeticLatLon(rBodyITRF(:,MOON), temp, body_lat(MOON), body_lon(MOON))
    if(hasFailed()) return


    !=========================================================================
    !
    ! Get GM parameters for Sun and Moon
    !
    !-------------------------------------------
    mu(SUN)  = solarsystem_model%getBodyGM(ID_SUN)
    if(hasFailed()) return
    mu(MOON) = solarsystem_model%getBodyGM(ID_MOON)
    if(hasFailed()) return


    !=========================================================================
    !
    ! Compute Legendre polynomials
    !
    !------------------------------------------
    do i = SUN, MOON

      lp(i,0,0) = 1.d0
      lp(i,0,1) = 0.d0
      lp(i,1,0) = sin(body_lat(i))
      lp(i,1,1) = cos(body_lat(i))

      !** determine legendre polynomials recursively
      do m = 0, lmax
        do l = max(2,m), lmax
          if(l == m) then
            lp(i,m,m) = (2*m-1)*lp(i,1,1)*lp(i,m-1,m-1)
          else if(l == m + 1) then
            lp(i,l,m) = (2*m+1)*lp(i,1,0)*lp(i,l-1,m)
          else
            lp(i,l,m) = ((2*l-1)*lp(i,1,0)*lp(i,l-1,m) - (l+m-1)*lp(i,l-2,m))/(l-m)
          end if
        end do
      end do
    end do ! Bodies

    !=========================================================================
    !
    ! Compute coefficients deviations for solid Earth tides
    !
    !----------------------------------------------------------------------

    dC = 0.d0
    dS = 0.d0

    muEarth = getEarthGravity()
    rekm    = getEarthRadius()

    if(tidetype == SOLID_TIDES) then

      do l = 2,lmax-1
        do m = 0,l
          if(m == 0) then
            dm = 1
          else
            dm = 2
          end if

          ! consider transformation for unnormalized coefficients
          fac   = sqrt(real(factorial(l-m)/factorial(l+m))) 
          templ = k(l,m)*sqrt(real(dm))*fac/(muEarth*sqrt(2.d0*l+1.d0))

          do i = SUN, MOON
            temp    = templ*mu(i)*(rekm/rabs_body(i))**(l+1.d0)*lp(i,l,m)
            dC(l,m) = dC(l,m) + temp*cos(m*body_lon(i))
            dS(l,m) = dS(l,m) + temp*sin(m*body_lon(i))
          end do
        end do
      end do

      !** correct changes in degree 4 coefficients due to degree 2 solid tides
      !----------------------------------------------------------------------------
      do m = 0,2
        if(m == 0) then
          dm = 1
        else
          dm = 2
        end if

        templ = kp(m)*sqrt(real(dm*factorial(2-m)/factorial(2+m)))/(muEarth*sqrt(5.d0))

        do i = SUN, MOON
          temp    = templ*mu(i)*(rekm/rabs_body(i))**3*lp(i,2,m)
          dC(4,m) = dC(4,m) + temp*cos(m*body_lon(i))
          dS(4,m) = dS(4,m) + temp*sin(m*body_lon(i))
        end do
      end do

      !** pole tide correction to 2,1 terms (only if EOP selected and already initialized)
      !--------------------------------------------------------------------------------------
      if(reduction%getEopInitFlag()) then

        pom    = reduction%getPolarMotion(time_mjd)
        pomAvg = reduction%getPolarMotionAvg(time_mjd)

        m1 = pom(1) - pomAvg(1)
        m2 = pomAvg(2) - pom(2)

        dC(2,1) = dC(2,1) - 1.721d-9*(m1 - 0.0115*m2)
        dS(2,1) = dS(2,1) - 1.721d-9*(m2 + 0.0115*m1)

      end if

    !=========================================================================
    !
    ! Compute coefficients deviations for ocean tides
    !
    !----------------------------------------------------------------------
    else if(tidetype == OCEAN_TIDES) then

      ! get all ocean tides corrections on harmonic coefficients (including the pole tide ones) 
      call get_FES2004_corrections(this, time_mjd, lmax, reduction, dC, dS)

    end if

    !===============================================================================================
    !
    ! Compute required quantities required for the accelerations (or take from geoopotential...)
    !
    !---------------------------------------------------------------------------------------------
    insig1 = 0.d0
    insig2 = 0.d0
    insig3 = 0.d0

    !** get radius, geocentric longitude and latitude
    call getRadiusLatLon(r_itrf, v_itrf, rabs, phi_gc, lambda)

    !** orbital radius
    rabs2      = rabs*rabs
    rabs3      = rabs2*rabs
    oorabs     = 1.d0/rabs
    oorabs2    = oorabs*oorabs
    oorabs3    = oorabs2*oorabs

    lpsat(0,0) = 1.d0
    lpsat(0,1) = 0.d0
    lpsat(1,0) = sin(phi_gc)
    lpsat(1,1) = cos(phi_gc)

    !** determine legendre polynomials recursively
    do m = 0, lmax

      do l = max(2,m), lmax

        if(l == m) then
          lpsat(m,m) = (2*m-1)*lpsat(1,1)*lpsat(m-1,m-1)
        else if(l == m + 1) then
          lpsat(l,m) = (2*m+1)*lpsat(1,0)*lpsat(l-1,m)
        else
          lpsat(l,m) = ((2*l-1)*lpsat(1,0)*lpsat(l-1,m) - (l+m-1)*lpsat(l-2,m))/(l-m)
        end if

      end do

    end do

    ! determine partial derivatives of the disturbing potential

    costerm(0) = 1.d0
    costerm(1) = cos(lambda)

    sinterm(0) = 0.d0
    sinterm(1) = sin(lambda)

    tanphi = tan(phi_gc)

    insig1 = 0.d0
    insig2 = 0.d0
    insig3 = 0.d0

    do l = 2,lmax

      !** recursive computation of sin(ml) and cos(ml)
      !--------------------------------------------------
      costerm(l) = 2.d0*costerm(1)*costerm(l-1)-costerm(l-2)
      sinterm(l) = 2.d0*costerm(1)*sinterm(l-1)-sinterm(l-2)

      ! determine pre-factor for expression inside the sigmas
      rrfac        = (rekm*oorabs)**l
      lpsat(l,l+1) = 0.d0

      do m = 0,l !MIN(l,maxord)

        temp_t1 = rrfac*(l+1)*lpsat(l,m)
        temp_t2 = rrfac*(lpsat(l,m+1) - (m * tanphi * lpsat(l,m)))
        temp_t3 = rrfac* m * lpsat(l,m)

        ! radial partial derivative of the potential
        insig1 = insig1 + temp_t1 * (dC(l,m)*costerm(m) + dS(l,m)*sinterm(m))

        ! phi (latitudal) derivative of the potential
        insig2 = insig2 + temp_t2 * (dC(l,m)*costerm(m) + dS(l,m)*sinterm(m))

        ! lambda (longitudal) derivative of the potential
        insig3 = insig3 + temp_t3 * (dS(l,m)*costerm(m) - dC(l,m)*sinterm(m))

      end do

    end do

    temp_t1 = -muEarth*oorabs3
    temp_t2 =  muEarth*oorabs

    write(*,*) "temp_t = ", temp_t1, temp_t2, muEarth, oorabs3, oorabs

    ! compute the radial partial derivative of the potential
    dudr      = temp_t1 * insig1

    ! compute the latitutal partial derivative of the potential
    dudphi    =  temp_t2  * insig2

    ! compute the longitudal partial derivative of the potential
    dudlambda =  temp_t2  * insig3

    write(*,*) "du = ", dudr, dudphi, dudlambda

    write(*,*) "r_gcrf = ", r_gcrf

    ! pre-compute terms which are used for the acceleration components
    r1r2        = r_gcrf(1)*r_gcrf(1) + r_gcrf(2)*r_gcrf(2)
    temp_t3     = 1.d0/r1r2
    sqrt_r1r2   = sqrt(r1r2)
    oosqrt_r1r2 = 1.d0/sqrt_r1r2
    temp_t4     = r_gcrf(3)*oorabs2*oosqrt_r1r2
    temp_t5     = oorabs2*sqrt_r1r2
    temp2       = dudlambda*temp_t3
    temp        = dudr - temp_t4*dudphi

    write(*,*) "temp = ", temp, temp2, temp_t5, dudr, dudphi

    !==========================================================================
    !
    ! Finally, compute the non-spherical perturbative accelerations in the GCRF
    !
    !----------------------------------------------------------

    ! i-direction [km/s^2]
    accel(1) = temp * r_gcrf(1) - temp2 * r_gcrf(2)

    ! j-direction [km/s^2]
    accel(2) = temp * r_gcrf(2) + temp2 * r_gcrf(1)

    ! k-direction [km/s^2]
    accel(3) = dudr * r_gcrf(3) + temp_t5 * dudphi

    write(*,*) "accel = ", accel, temp, temp2, temp_t5, dudr
    read(*,*)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

  end subroutine getTidesAcceleration

!=========================================================================
!
!> @anchor      setTidesInitFlag
!!
!> @brief       Set initialization flag to .false.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 16.07.2013 (initial design) </li>
!!              </ul>
!!
!!-----------------------------------------------------------------------
  subroutine setTidesInitFlag(this)

    class(Tides_class)  :: this
    this%tidesInitialized = .false.

  end subroutine setTidesInitFlag

!=========================================================================
!
!> @anchor      getDelaunay_arg
!!
!> @brief       Compute the Delaunay arguments (IERS chapter 5)
!> @author      Andrea Turchi (ATU)
!!
!> @date        <ul>
!!                <li> 17.08.22 (initial design) </li>
!!              </ul>
!!
!!-----------------------------------------------------------------------
  subroutine get_Delaunay_arg(this, time_mjd, F_vect)

    implicit none
    class(Tides_class)                  :: this
    real(dp), intent(in)                :: time_mjd  ! time in modified julian date
    real(dp)                            :: time_jd   ! time in julian date
    real(dp)                            :: jd_cent   ! time in julian centuries
    real(dp)                            :: l         ! mean anomaly of the Moon
    real(dp)                            :: l_prime   ! mean anomaly of the Sun
    real(dp)                            :: F         ! mean longitude of the Moon minus Omega
    real(dp)                            :: D         ! mean elongation of the Moon from the Sun
    real(dp)                            :: Omega     ! mean longitude of the ascending node of the Moon
    real(dp), dimension(5), intent(out) :: F_vect    ! Delaunay arguments vector

    !get input time in Julian centuries
    time_jd = time_mjd + jd245
    jd_cent = (time_jd - jd2000)/julian_century
    
    !compute Delaunay arguments in degrees 
    l = 134.96340251d0 + (1717915923.217800d0 * jd_cent + 31.879200d0 * jd_cent**2.d0 &
        + 0.05163500d0 * jd_cent**3.d0 - 0.0002447000d0 * jd_cent**4.d0) * as2deg
    l_prime = 357.52910918d0 + (129596581.048100d0 * jd_cent - 0.553200d0 * jd_cent**2.d0 &
              + 0.00013600d0 * jd_cent**3.d0 - 0.0000114900d0 * jd_cent**4.d0) * as2deg
    F = 93.27209062d0 + (1739527262.847800d0 * jd_cent - 12.751200d0 * jd_cent**2.d0 &
        - 0.00103700d0 * jd_cent**3.d0 + 0.0000041700d0 * jd_cent**4.d0) * as2deg
    D = 297.85019547d0 + (1602961601.209000d0 * jd_cent - 6.370600d0 * jd_cent**2.d0 &
        + 0.00659300d0 * jd_cent**3.d0 - 0.0000316900d0 * jd_cent**4.d0) * as2deg
    Omega = 125.04455501d0 + (- 6962890.543100d0 * jd_cent + 7.472200d0 * jd_cent**2.d0 &
            + 0.00770200d0 * jd_cent**3.d0 - 0.0000593900d0 * jd_cent**4.d0) * as2deg

    ! insert Delaunay arguments inside a vector
    F_vect = (/l, l_prime, F, D, Omega/)                           

  end subroutine get_Delaunay_arg

!=========================================================================
!
!> @anchor      getDoodson_arg
!!
!> @brief       Compute the Doodson arguments (it has been used the same approach as shown in Orekit at the following link:
!               https://github.com/CS-SI/Orekit/blob/develop/src/main/java/org/orekit/forces/gravity/potential/OceanTidesWave.java)
!> @author      Andrea Turchi (ATU)
!!
!> @date        <ul>
!!                <li> 23.08.22 (initial design) </li>
!!              </ul>
!!
!!-----------------------------------------------------------------------
  subroutine get_Doodson_arg(this, doodson, ctheta_g, cl, cl_prime, cF, cD, cOmega)

    implicit none
    class(Tides_class)   :: this
    real(dp), intent(in) :: doodson                                 ! Doodson number of the ocean tide
    integer              :: cPs, cNPrime, cP, cH, cS, cTau          ! coefficients extracted from the Doodson number
    integer, intent(out) :: ctheta_g, cl, cl_prime, cF, cD, cOmega  ! Doodson arguments

    ! extract coefficients from the tide's Doodson number
    cPs     = mod(doodson, 10.d0) - 5
    cNPrime = mod((doodson / 10.d0), 10.d0) - 5
    cP      = mod((doodson / 100.d0), 10.d0) - 5
    cH      = mod((doodson / 1000.d0), 10.d0) - 5
    cS      = mod((doodson / 10000.d0), 10.d0) - 5
    cTau    = mod((doodson / 100000.d0), 10.d0)

    ! compute Doodson arguments
    ctheta_g = cTau
    cl       = -cP
    cl_prime = -cPs
    cF       = -cTau + cS + cH + cP + cPs
    cD       = -cH - cPs
    cOmega   = -cTau + cS + cH + cP - cNPrime + cPs

  end subroutine get_Doodson_arg

!=========================================================================
!
!> @anchor      initFES2004
!!
!> @brief       Initialize the FES2004 ocean tides model by collecting
!!              harmonic coefficients data from an external file
!!
!> @author      Andrea Turchi (ATU)
!!
!> @date        <ul>
!!                <li> 25.08.22 (initial design) </li>
!!              </ul>
!!
!!------------------------------------------------------------------------
  subroutine init_FES2004(this)

    implicit none
    class(Tides_class)                             :: this
    integer                                        :: ich, ios      ! variables for data loading          
    integer                                        :: temp_l        ! temporary degree value
    integer                                        :: temp_m        ! temporary order value
    integer                                        :: ind           ! index
    character(len=255)                             :: cbuf          ! variable for data loading
    character(len=255)                             :: Darw          ! Darwin symbol of the ocean tide
    real(dp)                                       :: temp_dCp      ! temporary prograde C coefficient
    real(dp)                                       :: temp_dSp      ! temporary prograde S coefficient
    real(dp)                                       :: temp_dCm      ! temporary retrograde C coefficient
    real(dp)                                       :: temp_dSm      ! temporary retrograde S coefficient
    real(dp)                                       :: temp_Doodson  ! temporary Doodson number of the ocean tide

    !read data from external file
    ich = openFile(trim(adjustl(this%dataPath))//cdelimit//trim(adjustl(this%fesDataFile)), SEQUENTIAL, IN_FORMATTED)
    do ind = 1, 59462
      read(ich, '(a)', iostat=ios) cbuf
      read(cbuf, *) temp_Doodson, Darw, temp_l, temp_m, temp_dCp, temp_dSp, temp_dCm, temp_dSm
      !collect data for the desired l values
      if (temp_l >= 2 .and. temp_l <= this%max_l) then
        !collect data for each tide constituent
        if (temp_Doodson >= 55.564d0 .and. temp_Doodson <= 55.566d0) then
          this%dC_p(1, temp_l, temp_m)  = temp_dCp
          this%dS_p(1, temp_l, temp_m)  = temp_dSp
          this%dC_m(1, temp_l, temp_m)  = temp_dCm
          this%dS_m(1, temp_l, temp_m)  = temp_dSm
          this%doodson(1)               = temp_Doodson
        else if (temp_Doodson >= 55.574d0 .and. temp_Doodson <= 55.576d0) then
          this%dC_p(2, temp_l, temp_m)  = temp_dCp
          this%dS_p(2, temp_l, temp_m)  = temp_dSp
          this%dC_m(2, temp_l, temp_m)  = temp_dCm
          this%dS_m(2, temp_l, temp_m)  = temp_dSm
          this%doodson(2)               = temp_Doodson
        else if (temp_Doodson >= 56.553d0 .and. temp_Doodson <= 56.555d0) then
          this%dC_p(3, temp_l, temp_m)  = temp_dCp
          this%dS_p(3, temp_l, temp_m)  = temp_dSp
          this%dC_m(3, temp_l, temp_m)  = temp_dCm
          this%dS_m(3, temp_l, temp_m)  = temp_dSm
          this%doodson(3)               = temp_Doodson
        else if (temp_Doodson >= 57.554d0 .and. temp_Doodson <= 57.556d0) then
          this%dC_p(4, temp_l, temp_m)  = temp_dCp
          this%dS_p(4, temp_l, temp_m)  = temp_dSp
          this%dC_m(4, temp_l, temp_m)  = temp_dCm
          this%dS_m(4, temp_l, temp_m)  = temp_dSm
          this%doodson(4)               = temp_Doodson
        else if (temp_Doodson >= 65.454d0 .and. temp_Doodson <= 65.456d0) then
          this%dC_p(5, temp_l, temp_m)  = temp_dCp
          this%dS_p(5, temp_l, temp_m)  = temp_dSp
          this%dS_p(5, temp_l, temp_m)  = temp_dSp
          this%dS_p(5, temp_l, temp_m)  = temp_dSp
          this%dC_m(5, temp_l, temp_m)  = temp_dCm
          this%dS_m(5, temp_l, temp_m)  = temp_dSm
          this%doodson(5)               = temp_Doodson
        else if (temp_Doodson >= 75.554d0 .and. temp_Doodson <= 75.556d0) then
          this%dC_p(6, temp_l, temp_m)  = temp_dCp
          this%dS_p(6, temp_l, temp_m)  = temp_dSp
          this%dC_m(6, temp_l, temp_m)  = temp_dCm
          this%dS_m(6, temp_l, temp_m)  = temp_dSm
          this%doodson(6)               = temp_Doodson
        else if (temp_Doodson >= 85.454d0 .and. temp_Doodson <= 85.456d0) then
          this%dC_p(7, temp_l, temp_m)  = temp_dCp
          this%dS_p(7, temp_l, temp_m)  = temp_dSp
          this%dC_m(7, temp_l, temp_m)  = temp_dCm
          this%dS_m(7, temp_l, temp_m)  = temp_dSm
          this%doodson(7)               = temp_Doodson
        else if (temp_Doodson >= 93.554d0 .and. temp_Doodson <= 93.556d0) then
          this%dC_p(8, temp_l, temp_m)  = temp_dCp
          this%dS_p(8, temp_l, temp_m)  = temp_dSp
          this%dC_m(8, temp_l, temp_m)  = temp_dCm
          this%dS_m(8, temp_l, temp_m)  = temp_dSm
          this%doodson(8)               = temp_Doodson
        else if (temp_Doodson >= 135.654d0 .and. temp_Doodson <= 135.656d0) then
          this%dC_p(9, temp_l, temp_m)  = temp_dCp
          this%dS_p(9, temp_l, temp_m)  = temp_dSp
          this%dC_m(9, temp_l, temp_m)  = temp_dCm
          this%dS_m(9, temp_l, temp_m)  = temp_dSm
          this%doodson(9)               = temp_Doodson
        else if (temp_Doodson >= 145.554d0 .and. temp_Doodson <= 145.556d0) then
          this%dC_p(10, temp_l, temp_m) = temp_dCp
          this%dS_p(10, temp_l, temp_m) = temp_dSp
          this%dC_m(10, temp_l, temp_m) = temp_dCm
          this%dS_m(10, temp_l, temp_m) = temp_dSm
          this%doodson(10)              = temp_Doodson
        else if (temp_Doodson >= 163.554d0 .and. temp_Doodson <= 163.556d0) then
          this%dC_p(11, temp_l, temp_m) = temp_dCp
          this%dS_p(11, temp_l, temp_m) = temp_dSp
          this%dC_m(11, temp_l, temp_m) = temp_dCm
          this%dS_m(11, temp_l, temp_m) = temp_dSm
          this%doodson(11)              = temp_Doodson
        else if (temp_Doodson >= 165.554d0 .and. temp_Doodson <= 165.556d0) then
          this%dC_p(12, temp_l, temp_m) = temp_dCp
          this%dS_p(12, temp_l, temp_m) = temp_dSp
          this%dC_m(12, temp_l, temp_m) = temp_dCm
          this%dS_m(12, temp_l, temp_m) = temp_dSm
          this%doodson(12)              = temp_Doodson
        else if (temp_Doodson >= 235.754d0 .and. temp_Doodson <= 235.756d0) then
          this%dC_p(13, temp_l, temp_m) = temp_dCp
          this%dS_p(13, temp_l, temp_m) = temp_dSp
          this%dC_m(13, temp_l, temp_m) = temp_dCm
          this%dS_m(13, temp_l, temp_m) = temp_dSm
          this%doodson(13)              = temp_Doodson
        else if (temp_Doodson >= 245.654d0 .and. temp_Doodson <= 245.656d0) then
          this%dC_p(14, temp_l, temp_m) = temp_dCp
          this%dS_p(14, temp_l, temp_m) = temp_dSp
          this%dC_m(14, temp_l, temp_m) = temp_dCm
          this%dS_m(14, temp_l, temp_m) = temp_dSm
          this%doodson(14)              = temp_Doodson
        else if (temp_Doodson >= 255.554d0 .and. temp_Doodson <= 255.556d0) then
          this%dC_p(15, temp_l, temp_m) = temp_dCp
          this%dS_p(15, temp_l, temp_m) = temp_dSp
          this%dC_m(15, temp_l, temp_m) = temp_dCm
          this%dS_m(15, temp_l, temp_m) = temp_dSm
          this%doodson(15)              = temp_Doodson
        else if (temp_Doodson >= 273.554d0 .and. temp_Doodson <= 273.556d0) then
          this%dC_p(16, temp_l, temp_m) = temp_dCp
          this%dS_p(16, temp_l, temp_m) = temp_dSp
          this%dC_m(16, temp_l, temp_m) = temp_dCm
          this%dS_m(16, temp_l, temp_m) = temp_dSm
          this%doodson(16)              = temp_Doodson
        else if (temp_Doodson >= 275.554d0 .and. temp_Doodson <= 275.556d0) then
          this%dC_p(17, temp_l, temp_m) = temp_dCp
          this%dS_p(17, temp_l, temp_m) = temp_dSp
          this%dC_m(17, temp_l, temp_m) = temp_dCm
          this%dS_m(17, temp_l, temp_m) = temp_dSm
          this%doodson(17)              = temp_Doodson
        else if (temp_Doodson >= 455.554d0 .and. temp_Doodson <= 455.556d0) then
          this%dC_p(18, temp_l, temp_m) = temp_dCp
          this%dS_p(18, temp_l, temp_m) = temp_dSp
          this%dC_m(18, temp_l, temp_m) = temp_dCm
          this%dS_m(18, temp_l, temp_m) = temp_dSm
          this%doodson(18)              = temp_Doodson
        end if
      end if
    end do
    ich = closeFile(ich)
    !assign the right order of magnitude to harmonic coefficients
    this%dC_p = this%dC_p * 1.d-11
    this%dC_m = this%dC_m * 1.d-11
    this%dS_p = this%dS_p * 1.d-11
    this%dS_m = this%dS_m * 1.d-11

  end subroutine init_FES2004

!=========================================================================
!
!> @anchor      getFES2004_corrections
!!
!> @brief       Compute harmonic coefficients corrections produced by the
!!              FES2004 ocean tides model
!!
!> @author      Andrea Turchi (ATU)
!!
!> @date        <ul>
!!                <li> 25.08.22 (initial design) </li>
!!              </ul>
!!
!!------------------------------------------------------------------------
  subroutine get_FES2004_corrections(this, time_mjd, lmax, reduction, dC, dS)

    implicit none
    class(Tides_class)                          :: this
    real(dp), intent(in)                        :: time_mjd                               ! time in modified julian date
    integer, intent(in)                         :: lmax                                   ! maximum degree of harmonic coefficients
    type(Reduction_type), intent(inout)         :: reduction
    real(dp), dimension(5)                      :: F_vect                                 ! Delaunay arguments vector
    real(dp)                                    :: theta_g                                ! Greenwich mean sidereal time in radians
    real(dp)                                    :: theta_f                                ! argument of the tide constituent
    real(dp)                                    :: m1, m2                                 ! auxiliaries to account for pole tide
    real(dp)                                    :: fac                                    ! unnormalization factor
    integer                                     :: i                                      ! index
    integer                                     :: l                                      ! coefficients degree
    integer                                     :: m                                      ! coefficients order
    integer                                     :: dm                                     ! coefficient in transformation between normalized and unnormalized quantities
    integer                                     :: ctheta_g, cl, cl_prime, cF, cD, cOmega ! Doodson arguments
    real(dp), dimension(2)                      :: pom, pomAvg                            ! polar motion variables xp,yp (current and running average) / rad
    real(dp), dimension(2:6,0:6), intent(inout) :: dC                                     ! matrix with corrections on the C harmonic coefficients
    real(dp), dimension(2:6,0:6), intent(inout) :: dS                                     ! matrix with corrections on the S harmonic coefficients

    !get Delaunay arguments in radians for the current time
    call get_Delaunay_arg(this, time_mjd, F_vect)
    F_vect = F_vect*deg2rad

    !get Greenwich Mean Sidereal Time
    theta_g = getGMST(time_mjd)

    do l = 2, lmax
      do m = 0, l
        if (m == 0) then
          dm = 1
        else
          dm = 2
        end if 
        !compute the corresponding unnormalization factor
        fac = sqrt(factorial(l - m)*dm*(2.d0*l + 1.d0)/factorial(l + m))
        do i = 1, 18
          !compute the argument for each tide constituent            
          call get_Doodson_arg(this, this%doodson(i), ctheta_g, cl, cl_prime, cF, cD, cOmega)
          theta_f = ctheta_g * theta_g + cl * F_vect(1) + cl_prime * F_vect(2) + cF * F_vect(3) + &
                    cD * F_vect(4) + cOmega * F_vect(5)

          !get the produced gravity field corrections
          dC(l, m) = dC(l, m) + fac*((this%dC_p(i, l, m) + this%dC_m(i, l, m))*cos(theta_f) + &
                     (this%dS_p(i, l, m) + this%dS_m(i, l, m))*sin(theta_f))
          if (m == 0) then
            dS(l, m) = 0
          else
            dS(l, m) = dS(l, m) + fac*((this%dS_p(i, l, m) - this%dS_m(i, l, m))*cos(theta_f) - &
                       (this%dC_p(i, l, m) - this%dC_m(i, l, m))*sin(theta_f))
          end if
        end do
      end do
    end do

    !Ocean pole tide correction
    if(reduction%getEopInitFlag()) then

      pom    = reduction%getPolarMotion(time_mjd)
      pomAvg = reduction%getPolarMotionAvg(time_mjd)

      m1 = pom(1) - pomAvg(1)
      m2 = pomAvg(2) - pom(2)

      dC(2,1) = dC(2,1) - 2.1778d-10*(m1 - 0.01724*m2)*sqrt(5.d0/3.d0)
      dS(2,1) = dS(2,1) - 1.7232d-10*(m2 - 0.03365*m1)*sqrt(5.d0/3.d0)

    end if

  end subroutine get_FES2004_corrections    

end module tides
