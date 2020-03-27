!>-----------------------------------------------------------------------------------------------
!> @brief       Averaging techniques to obtain mean orbital elements
!!
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB: 30.07.2013 (initial design) </li>
!!                <li>VB: 16.01.2014 (now using 'units' module) </li>
!!                <li>CHK: 13.11.2015 (updated to use libslam, added some more documentation) </li>
!!              </ul>
!!
!> @details     This module contains all routines and parameters for the application of
!!              averaging techniques and the providing of mean elements.
!!
!!              \par Includes
!!
!!              <ul>
!!                <li> neptune_error_handling.f90 </li>
!!                <li> slam_math.f90 </li>
!!                <li> slam_astro_conversions.f90 </li>
!!                <li> slam_orbit_types.f90 </li>
!!                <li> slam_types.f90 </li>
!!                <li> slam_units.f90 </li>
!!              </ul>
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      averaging
!!
!!------------------------------------------------------------------------------------------------
module averaging

  use slam_astro_conversions,  only: getOrbitType
  use slam_error_handling,     only: FATAL, E_KEPLER_UNITS, isControlled, hasToReturn, checkIn, checkOut
  use neptune_error_handling,  only: setNeptuneError, E_MEAN_ELS_MISSING
  use slam_math,               only: deg2rad
  use slam_orbit_types,        only: kepler_t, assignment(=)
  use slam_types,              only: dp
  use slam_units,              only: UNIT_DEG, UNIT_RAD, UNIT_KM, UNIT_M

  implicit none

  private

  logical :: isAvailable = .false.    ! flag for the availability of a current mean state (kepler elements)

  type(kepler_t) :: meanState

! Averaging is not a class yet, but meanState and isAvailable should be private per thread
!$omp threadprivate(isAvailable,meanState)

  !===================================
  !
  ! public methods
  !
  !------------------
  public :: getMeanElements
  public :: setMeanElements

contains

!>-----------------------------------------------------------------------------------------------
!> @brief       Flag for the availability of mean elements as global variable of this module
!!
!> @author      Vitali Braun
!!
!> @return      .true. when mean elements are available
!!
!> @date        <ul>
!!                <li> 30.07.2013 (initial design) </li>
!!              </ul>
!!
!> @anchor      isAvailableMeanState
!!
!!------------------------------------------------------------------------------------------------
  logical function isAvailableMeanState()

    isAvailableMeanState = isAvailable

  end function isAvailableMeanState

!>-----------------------------------------------------------------------------------------------
!> @brief       Returns the current mean elements
!!
!> @author      Vitali Braun
!!
!> @return      mean keplerian elements as kepler_t
!!
!! @date        <ul>
!!                <li> 30.07.2013 (initial design) </li>
!!              </ul>
!!
!! @anchor      getMeanElements
!!
!!------------------------------------------------------------------------------------------------
  type(kepler_t) function getMeanElements()
    
    character(len=*), parameter :: csubid = "getMeanElements"
    real(dp) :: zero = 0.0_dp 
   
    getMeanElements = zero

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(isAvailable) then
      getMeanElements = meanState

    else
      call setNeptuneError(E_MEAN_ELS_MISSING, FATAL)
      getMeanElements = zero
      return

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getMeanElements

!>-----------------------------------------------------------------------------------------------
!> @brief       Set mean elements
!!
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 30.07.2013 (initial design) </li>
!!              </ul>
!!
!> @param[in]   kep     Kepler elements
!!
!> @anchor      setMeanElements
!!
!!------------------------------------------------------------------------------------------------
  subroutine setMeanElements(kep)

    type(kepler_t), intent(in) :: kep

    character(len=*), parameter :: csubid = 'setMeanElements'
    integer :: orbtype
    
    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    meanState = kep

    !** check for RAD/DEG units
    select case(meanState%angles_unit)

      case(UNIT_DEG) ! convert to radian
        meanState%inc     = meanState%inc*deg2rad       ! inclination
        meanState%raan    = meanState%raan*deg2rad      ! right ascension of ascending node
        meanState%aop     = meanState%aop*deg2rad       ! argument of pericenter
        meanState%man     = meanState%man*deg2rad       ! mean anomaly
        meanState%ecan    = meanState%ecan*deg2rad      ! eccentric anomaly
        meanState%tran    = meanState%tran*deg2rad      ! true anomaly
        meanState%arglat  = meanState%arglat*deg2rad    ! argument of true latitude
        meanState%lonper  = meanState%lonper*deg2rad    ! longitude of perigee
        meanState%truelon = meanState%truelon*deg2rad   ! true longitude

        meanState%angles_unit = UNIT_RAD

      case(UNIT_RAD)
        ! nothing happens

      case default
        call setNeptuneError(E_KEPLER_UNITS, FATAL)
        return
 
    end select

    !** check for SMA unit
    select case(meanState%sma_unit)

      case(UNIT_KM)
        ! nothing happens..
      
      case(UNIT_M)
        meanState%sma = meanState%sma/1.d3
        meanState%sma_unit = UNIT_KM

      case default
        call setNeptuneError(E_KEPLER_UNITS, FATAL)
        return
 
    end select

    !** determine orbit type and compute missing values
    call getOrbitType(meanState,orbtype)
    
    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    isAvailable = .true.
    return

  end subroutine setMeanElements

end module averaging


