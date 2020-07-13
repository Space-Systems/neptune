!=================================================================================================
!!!
!> @brief       Third body perturbations
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  28.02.2013 (initial design)</li>
!!                <li>VB:  03.06.2013 (code optimization) </li>
!!                <li>VB:  17.07.2013 (extended routine 'getThirdBodyPosition' to include moon) </li>
!!                <li>VB:  27.01.2014 (moved functions like 'initThirdBody' or 'getThirdBodyPosition' to newly created module 'solarsystem')</li>
!!                <li>CHK: 16.11.2015 (updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 05.01.2017 (created ThirdBody_class)
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for the
!!              estimation of perturbations due to third body gravitational attraction.
!!
!>  @todo Add perturbations due to Jupiter and Venus (or even all solar system bodies?!), as these introduce errors in the 10 to 1 m regime (comparable to tides or Albedo)
!!
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      thirdBody
!!
!!------------------------------------------------------------------------------------------------
module thirdBody

  use slam_astro,             only: getEarthGravity
  use neptune_error_handling, only: E_WRONG_THIRD_BODY, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, checkIn, checkOut
  use slam_math,              only: mag, identity_matrix, outerproduct
  use solarsystem,            only: Solarsystem_class, ID_SUN, ID_MOON
  use slam_types,             only: dp

  implicit none

  private

    type, public :: ThirdBody_class

        !** variables
        !----------------------------------------------------
        integer                 :: last_body_id                                 ! the body for which the getThirdBodyCovariance routine has been called in the previous call
        logical, dimension(2)   :: cov3b                                        ! switches for Sun (ID_SUN) and Moon (ID_MOON), which can be considered for the covariance matrix computation individually
        real(dp)                :: last_time_mjd                                ! the MJD for which the getThirdBodyCovariance routine has been called in the previous call

    contains

        !** getter
        procedure :: getThirdBodyAcceleration
        procedure :: getThirdBodyCovariance
        procedure :: getThirdBodyCovFlag

        !** setter
        procedure :: setThirdBodyCovFlag

    end type ThirdBody_class

    ! Constructor
    interface ThirdBody_class
        module procedure constructor
    end interface ThirdBody_class

contains

    ! ====================================================================
    !!
    !> @brief      Constructor that should be used to initialize variables.
    !!
    !> @author     Christopher Kebschull
    !> @date       <ul>
    !!                  <li>ChK: 24.12.2017 (initial implementation)</li>
    !!              </ul>
    !> @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(ThirdBody_class) function constructor()
        constructor%cov3b = .false.   ! switches for Sun (ID_SUN) and Moon (ID_MOON), which can be considered for the covariance matrix computation individually
        constructor%last_time_mjd = -1.0
    end function constructor

! =======================================================================================================
!
!> @anchor      getThirdBodyAcceleration
!!
!> @brief       Computes acceleration due to gravity perturbations caused by third bodies
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 28.02.2013 (initial design)</li>
!!                <li> 27.01.2014 (changed input variable to be a third body's ID instead of a character)</li>
!!              </ul>
!!
!> @param[in]   solarsystem_model the solar system model
!> @param[in]   r_gcrf            radius vector in GCRF
!> @param[in]   time_mjd          current time MJD
!> @param[in]   body_id           third body identifier (e.g. SUN = 1 or MOON = 2)
!> @param[out]  acc_thirdbody     acceleration vector in inertial frame
!!
!! @details     This routine computes the acceleration in the GCRF due to third body
!!              gravity perturbations. It is based on SPICE routine SPKPOS to compute the
!!              position in the GCRF frame (called J2000 in SPICE) for the given time. However
!!              the data is interpolated using chebyshev polynomials.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> ... </li>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  subroutine getThirdBodyAcceleration(                      &
                                       this,                &
                                       solarsystem_model,   &
                                       r_gcrf,              &  ! <-- DBL(3) radius vector in GCRF
                                       time_mjd,            &  ! <-- DBL    current time MJD
                                       body_id,             &  ! <-- INT    third body identifier (e.g. 'SUN' = 1, or 'MOON' = 2)
                                       acc_thirdbody        &  ! --> DBL(3) acceleration vector in inertial frame
                                     )


    !** interface
    !---------------------------------------------------
    class(ThirdBody_class)                  :: this
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    real(dp), dimension(3), intent(in)      :: r_gcrf
    real(dp),               intent(in)      :: time_mjd
    integer,                intent(in)      :: body_id

    real(dp), dimension(3), intent(out)     :: acc_thirdbody
    !---------------------------------------------------

    character(len=*), parameter :: csubid =  "getThirdBodyAcceleration"

    real(dp)  :: mu3      ! mass GM of third body
    real(dp)  :: OOr33    ! one over mag(r_gcrf)
    real(dp)  :: OOrobj33 ! one over robj**3
    real(dp), dimension(3) :: pos3    ! position vector of third body in GCRF frame
    real(dp)  :: robj3    ! magnitude of position vector satellite -> third body
    real(dp), dimension(3) :: vec_obj3  ! position vector satellite -> third body

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !=====================================================================
    !
    ! Get gravitational constant for passed body
    !
    !------------------------------------------------

    mu3 = solarsystem_model%getBodyGM(body_id)
    if(hasFailed()) return

    !=======================================================================
    !
    ! get third body position from SPICE
    !
    !-----------------------------------------------------

    pos3 = solarsystem_model%getBodyPosition(time_mjd, body_id)
    if(hasFailed()) return

    !=======================================================================
    !
    ! compute accelerations
    !
    !------------------------------------------------------

    !** compute vector from satellite to third body
    vec_obj3 = pos3 - r_gcrf

    !** magnitude of vec_obj3
    robj3 = mag(vec_obj3)

    OOrobj33 = 1.d0/robj3**3.d0
    OOr33    = 1.d0/mag(pos3)**3.d0

    !** finally compute acceleration
    acc_thirdbody(1:3) = mu3*(vec_obj3(1:3)*OOrobj33 - pos3(1:3)*OOr33)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine getThirdBodyAcceleration

!!------------------------------------------------------------------------------------------------
!> @anchor      setThirdBodyCovFlag
!!
!> @brief       Sets the appropriate flags to consider contributions of a third body in the covariance matrix propagation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 03.02.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   body_id           third body identifier (e.g. SUN = 1 or MOON = 2)
!> @param[in]   val               Value the cov3d flag is set to for the passed body
!!
!!------------------------------------------------------------------------------------------------
  subroutine setThirdBodyCovFlag(                &
                                  this,          &
                                  body_id,       &  ! <-- INT    third body identifier (e.g. 'SUN' = 1, or 'MOON' = 2)
                                  val            &  ! <-- LOG    switch value
                                )

    !** interface
    !---------------------------------------------------
    class(ThirdBody_class)  :: this
    integer, intent(in)     :: body_id
    logical, intent(in)     :: val
    !---------------------------------------------------

    character(len=*), parameter :: csubid =  "setThirdBodyCovFlag"
    character(len=255) :: cmess

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(body_id > size(this%cov3b ) .or. body_id < 0) then

      write(cmess,'(a,i3)') "id: ", body_id
      call setNeptuneError(E_WRONG_THIRD_BODY, FATAL, (/cmess/))
      return

    else

      this%cov3b (body_id) = val

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine setThirdBodyCovFlag

!!------------------------------------------------------------------------------------------------
!> @anchor      getThirdBodyCovFlag
!!
!> @brief       Gets the appropriate flag cov3b
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 03.02.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   body_id           third body identifier (e.g. SUN = 1 or MOON = 2)
!!
!!------------------------------------------------------------------------------------------------
  logical function getThirdBodyCovFlag(                &
                                        this,          &
                                        body_id        &  ! <-- INT    third body identifier (e.g. 'SUN' = 1, or 'MOON' = 2)
                                      )

    !** interface
    !---------------------------------------------------
    class(ThirdBody_class)  :: this
    integer, intent(in)     :: body_id
    !---------------------------------------------------

    character(len=*), parameter :: csubid =  "getThirdBodyCovFlag"
    character(len=255) :: cmess

    getThirdBodyCovFlag = .false.

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(body_id > size(this%cov3b ) .or. body_id < 0) then

      write(cmess,'(a,i3)') "id: ", body_id
      call setNeptuneError(E_WRONG_THIRD_BODY, FATAL, (/cmess/))
      return

    else

      getThirdBodyCovFlag = this%cov3b (body_id)

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getThirdBodyCovFlag

!!------------------------------------------------------------------------------------------------
!> @anchor      getThirdBodyCovariance
!!
!> @brief       Gets the contributions to the variational equations due to third bodies (Sun and Moon)
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 03.02.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   r_gcrf            radius vector of satellite in GCRF
!> @param[in]   time_mjd          MJD
!> @param[in]   body_id           third body identifier (e.g. SUN = 1 or MOON = 2)
!!
!!------------------------------------------------------------------------------------------------
  function getThirdBodyCovariance(  this,               &
                                    solarsystem_model,  &
                                    r_gcrf,             &
                                    time_mjd,           &
                                    body_id             &
                                    ) result(cov)

    !** interface
    !---------------------------------------------------
    class(ThirdBody_class)                  :: this
    type(Solarsystem_class), intent(inout)  :: solarsystem_model
    real(dp), dimension(3), intent(in)      :: r_gcrf
    real(dp),               intent(in)      :: time_mjd
    integer,                intent(in)      :: body_id

    real(dp), dimension(3,3) :: cov     ! result in 3x3 matrix
    !---------------------------------------------------

    character(len=*), parameter :: csubid =  "getThirdBodyCovariance"
    character(len=255) :: cmess

    real(dp)               :: mu        ! gravitational constant of the Earth
    real(dp)               :: mu3       ! gravitational constant of third body
    real(dp)               :: rabs      ! magnitude of r_gcrf
    real(dp), dimension(3) :: rbody     ! position of third body in GCRF
    real(dp), dimension(3) :: rdiff     ! difference vector: rbody - rgcrf
    real(dp)               :: temp      ! auxiliary

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(body_id > size(this%cov3b ) .or. body_id < 0) then

      write(cmess,'(a,i3)') "id: ", body_id
      call setNeptuneError(E_WRONG_THIRD_BODY, FATAL, (/cmess/))
      return

    end if

    !** get constants
    mu   = getEarthGravity()
    mu3  = solarsystem_model%getBodyGM(body_id)
    if(hasFailed()) return

    rabs = mag(r_gcrf)

    if(time_mjd == this%last_time_mjd  .and. body_id /= this%last_body_id) then    ! the contribution due to the central body has only to be applied once, if Sun and Moon both are considered.
      temp = 0.d0
    else
      temp = mu/rabs**3
    end if

    !** get position of third body
    rbody = solarsystem_model%getBodyPosition(time_mjd, body_id)
    if(hasFailed()) return

    rdiff = rbody - r_gcrf

    !** initialize result matrix (identity matrix)
    call identity_matrix(cov)

    !================================================
    !
    ! Start computation
    !
    !-----------------------------------------

    ! 1) part with identity matrix (r**3)
    !cov = -(temp + mu3/(mag(rdiff))**3.d0)*cov
    cov = -(mu3/(mag(rdiff))**3.d0)*cov

    ! 2) part with r**5
    !cov = cov + 3.d0*(temp/rabs**2.d0*outerproduct(r_gcrf, r_gcrf) + mu3*outerproduct(rdiff, rdiff)/(mag(rdiff))**5.d0)
    cov = cov - 3.d0*mu3*outerproduct(rdiff, rdiff)/(mag(rdiff))**5.d0

    !write(*,*) "mu3 = ", mu3
    !write(*,*) "rdiff = ", rdiff

    this%last_time_mjd  = time_mjd
    this%last_body_id  = body_id

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getThirdBodyCovariance


end module thirdBody
