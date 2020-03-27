!> @anchor      astroConvert
!!
!> @brief       Provides simple conversion functions for astronomical applications
!!
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  06.06.2014 (initial design)</li>
!!                <li>CHK: 13.11.2015 (moved to neptune from the lib and updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 04.01.2018 (Using Solarsystem_class)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions related to the
!!              simple conversion of different quantities as listed below:
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      astroConvert
!!
!!------------------------------------------------------------------------------------------------
module astroConvert

  use slam_math,              only: deg2rad, twopi, mag, rad2deg, undefined
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, checkIn, checkOut
  use solarsystem,            only: Solarsystem_class, ID_SUN
  use slam_types,             only: dp

  implicit none

  private

  public :: raan2solarTime
  public :: solarTime2Raan

contains
  
  !=========================================================================
  !
  !>  @brief  Converts an angle given as solar time to RAAN
  !!
  !>  @author Vitali Braun
  !!
  !>  @param[inout] solarsystem_model  The solar system model
  !>  @param[in]    mjd   Modified Julian Day
  !>  @param[in]    lst   Local solar time
  !!
  !>  @return RAAN
  !!
  !>  @date   <ul>
  !!            <li>06.06.2014 (initial design)</li>
  !!          </ul>
  !!
  !!------------------------------------------------------------
  real(dp) function solarTime2Raan(solarsystem_model, mjd, lst) result(raan)

    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    real(dp), intent(in)                    :: mjd
    real(dp), intent(in)                    :: lst


    character(len=*), parameter :: csubid = 'solarTime2Raan'

    real(dp), dimension(3) :: rsun
   
    raan = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    ! get Sun's position
    rsun = solarsystem_model%getBodyPosition(mjd, ID_SUN)

    if(hasFailed()) return
    
    ! compute RAAN
    raan = atan2(rsun(2), rsun(1)) + 15.d0*(lst - 12.d0)*deg2rad

    ! convert to angle between 0...360 deg
    if(raan < 0.d0) then
      raan = raan + twopi
    end if

    ! done
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function solarTime2Raan
  
  !=========================================================================
  !
  !>  @brief  Converts a given RAAN to solar time in hours
  !!
  !!  @author Vitali Braun
  !!
  !>  @param[inout] solarsystem_model  The solar system model
  !>  @param[in]    mjd   Modified Julian Day
  !>  @param[in]    raan  Right ascension of ascending node / rad
  !!
  !>  @return RAAN
  !!
  !!  @date   <ul>
  !!            <li>06.06.2014 (initial design)</li>
  !!          </ul>
  !!
  !!------------------------------------------------------------
  real(dp) function raan2solarTime(solarsystem_model, mjd, raan) result(lst)

    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    real(dp), intent(in)                    :: mjd
    real(dp), intent(in)                    :: raan

    character(len=*), parameter :: csubid = 'raan2solarTime'

    real(dp), dimension(3) :: rsun
   
    lst = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    ! get Sun's position
    rsun = solarsystem_model%getBodyPosition(mjd, ID_SUN)

    if(hasFailed()) return
   
    ! compute RAAN
    lst = (raan - atan2(rsun(2), rsun(1)))*rad2deg/15.d0 + 12.d0

    ! convert to hours between 0...24 h
    if(lst < 0.d0) then
      lst = lst + 24.d0
    else if(lst >= 24.d0) then
      lst = lst - 24.d0
    end if

    ! done
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function raan2solarTime
 
end module astroConvert
