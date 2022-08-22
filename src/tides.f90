!=================================================================================================
!!
!> @brief       Solid and ocean tides modeling
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  16.07.2013 (initial design)    </li>
!!                <li>CHK: 13.11.2015 (updated to use with libslam)    </li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 04.01.2018 (Using Solarsystem_class)</li>
!!                <li>CHK: 05.01.2017 (created Tides_class)
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
!> @todo        Find bug in Solid tides (comparing to 2-body trajectory -> discontinuity?!)
!> @todo        Ocean tide model implementation incorrect? No reference and Vallado seems to be wrong. Consider implementing a REAL IERS model
!!------------------------------------------------------------------------------------------------
module tides

  use slam_astro,             only: getEarthRadius, getEarthGravity, getEarthMass
  use slam_astro_conversions, only: getGeodeticLatLon, getRadiusLatLon
  use neptune_error_handling, only: E_SOLAR_SYSTEM_INIT, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, checkIn, checkOut
  use slam_math,              only: mag, pi, eps9, rad2deg, factorial
  use slam_reduction_class,   only: Reduction_type
  use solarsystem,            only: Solarsystem_class, ID_SUN, ID_MOON
  use slam_types,             only: dp

  implicit none

  private

    integer, parameter, public :: SOLID_TIDES = 1
    integer, parameter, public :: OCEAN_TIDES = 2

    type, public :: Tides_class

        logical, public :: tidesInitialized

    contains

        procedure :: initTides
        procedure :: getTidesAcceleration
        procedure :: setTidesInitFlag

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
  subroutine initTides(this)

    class(Tides_class)  :: this

    character(len=*), parameter  :: csubid = "initTides"

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
    type(Reduction_type),intent(inout)     :: reduction
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
    real(dp), dimension(0:6)       :: costerm, sinterm    ! cos(ml) and sin(ml)
    real(dp), dimension(2:6,0:6)   :: dC, dS      ! deviations in coefficients C and S due to tides
    real(dp), parameter            :: densWater = 1025.d9             ! water density in kg/km**3
    real(dp), dimension(2:3,0:3)   :: k = reshape((/0.29525d0, 0.093d0, &     ! nominal love numbers for degree and order
                                                    0.29470d0, 0.093d0, &
                                                    0.29801d0, 0.093d0, &
                                                    0.d0,      0.094d0/), (/2,4/))
    real(dp), dimension(2:6)       :: kld = (/-0.3075d0, -0.1950d0, & ! load deformation coefficients
                                              -0.1320d0, -0.1032d0, &
                                              -0.0892d0/)
    real(dp), dimension(0:2)       :: kp =(/-0.00087d0, -0.00079d0, -0.00057d0/)  ! nominal love numbers for correction of fourth degree term (k+nm)
    real(dp), dimension(2,0:6,0:6) :: lp          ! legendre polynomials for Sun and Moon
    real(dp), dimension(0:7,0:7)   :: lpsat       ! legendre polynomials for satellite
    real(dp), dimension(2)         :: mu          ! GM parameter of Sun and Moon
    real(dp), dimension(2)         :: pom, pomAvg ! polar motion variables xp,yp (current and running average) / rad
    real(dp), dimension(2)         :: rabs_body
    real(dp), dimension(3,2)       :: rBodyGCRF
    real(dp), dimension(3,2)       :: rBodyITRF

    real(dp) :: const         ! constant factor
    real(dp) :: dudlambda     ! dU/d(lambda)
    real(dp) :: dudphi        ! dU/d(phi_gc)
    real(dp) :: dudr          ! dU/d(rabs)
    real(dp) :: fac           ! factorial term for the conversion to unnormalized coefficients
    real(dp) :: insig1, insig2, insig3  ! cumulated sums
    real(dp) :: lambda              ! longitude
    real(dp) :: m1, m2              ! auxiliaries to account for pole tide
    real(dp) :: muEarth             ! Earth's gravity constant
    real(dp) :: oorabs              ! 1/rabs
    real(dp) :: oorabs2             ! 1/rabs2
    real(dp) :: oorabs3             ! 1/rabs3
    real(dp) :: oosqrt_r1r2         ! 1/sqrt(r1r2)
    real(dp) :: phi_gc              ! latitude geocentric
    real(dp) :: r1r2                ! radius(1)**2 + radius(2)**2
    real(dp) :: rabs                ! magnitude of radius vector
    real(dp) :: rabs2               ! squared radius magnitude
    real(dp) :: rabs3               ! cubed radius magnitude
    real(dp) :: rekm                ! Earth's radius in km
    real(dp) :: rrfac               ! temporary
    real(dp) :: sqrt_r1r2           ! sqrt(r1r2)
    real(dp) :: tanphi              ! tan(phi_gc)
    real(dp) :: temp, templ, temp2  ! temporary
    real(dp) :: temp_t1             ! temporary (used to support tides)
    real(dp) :: temp_t2             ! temporary (used to support tides)
    real(dp) :: temp_t3             ! temporary (used to support tides)
    real(dp) :: temp_t4             ! temporary (used to support tides)
    real(dp) :: temp_t5             ! temporary (used to support tides)


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

          fac   = sqrt(real(factorial(l-m)/factorial(l+m))) ! consider transformation for unnormalized coefficients
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

        !templ = kp(m)*dm*sqrt(1.8d0*(4-m)*(3-m)/(4.d0+m)/(3.d0+m))*factorial(2-m)/factorial(2+m)/muEarth
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

        !dC(2,1) = dC(2,1) - 1.721d-9*(m1 - 0.0115*m2)
        dC(2,1) = dC(2,1) - 1.333d-9*(m1 - 0.0115*m2)*sqrt(5.d0/3.d0)
        !dS(2,1) = dS(2,1) - 1.721d-9*(m2 + 0.0115*m1)
        dS(2,1) = dS(2,1) - 1.333d-9*(m2 + 0.0115*m1)*sqrt(5.d0/3.d0)

      end if

    !=========================================================================
    !
    ! Compute coefficients deviations for ocean tides
    !
    !----------------------------------------------------------------------
    else if(tidetype == OCEAN_TIDES) then

      const = densWater/muEarth*4.d0*pi*rekm**2.d0/getEarthMass()

      do l = 2,lmax

        templ = const/(2.d0*l + 1.d0)*(1.d0 + kld(l))

        do m = 0,l
          do i = SUN, MOON

            temp = templ*mu(i)*(rekm/rabs_body(i))**(l+1.d0)*lp(i,l,m)

            dC(l,m) = dC(l,m) + temp*cos(m*body_lon(i))
            dS(l,m) = dS(l,m) + temp*sin(m*body_lon(i))

          end do
        end do
      end do

      !** DEBUG
!     write(*,*) "---", time_mjd, "----"
!     do l = 2, lmax
!       do m=0,l
!         write(*,*) "l, m, dC, dS = ", l, m, dC(l,m), dS(l,m)
!         if(l==6 .and. m==5) then
!           write(*,*) "lp = ", lp(1,l,m), lp(2,l,m)
!         end if
!       end do
!     end do
!     write(*,*) "--------------"

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
    rabs2        = rabs*rabs
    rabs3        = rabs2*rabs
    oorabs       = 1.d0/rabs
    oorabs2      = oorabs*oorabs
    oorabs3      = oorabs2*oorabs

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

    ! compute the radial partial derivative of the potential
    dudr      = temp_t1 * insig1

    ! compute the latitutal partial derivative of the potential
    dudphi    =  temp_t2  * insig2

    ! compute the longitudal partial derivative of the potential
    dudlambda =  temp_t2  * insig3

    ! pre-compute terms which are used for the acceleration components
    r1r2        = r_gcrf(1)*r_gcrf(1) + r_gcrf(2)*r_gcrf(2)
    temp_t3     = 1.d0/r1r2
    sqrt_r1r2   = sqrt(r1r2)
    oosqrt_r1r2 = 1.d0/sqrt_r1r2
    temp_t4     = r_gcrf(3)*oorabs2*oosqrt_r1r2
    temp_t5     = oorabs2*sqrt_r1r2
    temp2       = dudlambda*temp_t3
    temp        = dudr - temp_t4*dudphi

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

    !write(*,*) "accel = ", accel
    !read(*,*)

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

end module tides
