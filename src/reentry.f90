!==================================================================================================
!!
!> @brief       Enabling IADC re-entry campaigns
!> @author      Vitali Braun
!> @author      Christopher Kebschull
!!
!> @date        <ul>
!!                <li>VB: 16.10.2014 (initial design)</li>
!!                <li>VB: 20.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 02.01.2018 (regards NEPTUNE class)
!!                <li>CHK: 04.01.2018 (created Reentry_class)
!!              </ul>
!!
!> @details     This module contains type definitions, parameters and routines required for
!!              performing re-entry analysis, designed specifically to support IADC re-entry
!!              campaigns using NEPTUNE.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      reentry
!!
!!------------------------------------------------------------------------------------------------
module reentry

  use slam_math,    only: rad2deg
  use neptuneClass, only: Neptune_class
  use slam_types,   only: dp
  use slam_time,    only: time_t, mjd2gd, date2string, getDateTimeNowUtc

  implicit none

  private

    character(len=*), parameter :: coriginator = 'Wiedemann'                    !> originator of re-entry analysis
    integer,parameter :: LEN_REENTRY_MESSAGE = 1000                             !> length of IADC REDB reentry message


    type, public :: Reentry_class

        logical     :: reEntryInitialised                                       !> initialisation flag
        real(dp)    :: uncertaintyInterval                                      !> 20% uncertainty

    contains

        procedure :: initReEntry
        procedure :: getReEntryMessage
        procedure :: isInitialisedReEntry

    end type Reentry_class

    ! Constructor
    interface Reentry_class
        module procedure constructor
    end interface Reentry_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !!  @author     Christopher Kebschull
    !!  @date       <ul>
    !!                  <li>ChK: 04.01.2018 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Reentry_class) function constructor()
        constructor%reEntryInitialised = .false.                                !> initialisation flag
        constructor%uncertaintyInterval = 0.2d0                                 !> 20% uncertainty

    end function constructor

!=============================================================================
!
!> @anchor      initReEntry
!!
!> @brief       Initialises this module
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 16.10.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  subroutine initReEntry(this)

    class(Reentry_class)    :: this

    !** set initialisation flag
    this%reEntryInitialised = .true.

  end subroutine initReEntry

!=============================================================================
!
!> @anchor      getReEntryMessage
!!
!> @brief       Provides a re-entry report message, which can directly be used for the IADC REDB
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 16.10.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  character(len=LEN_REENTRY_MESSAGE) function getReEntryMessage(this,neptune) result(mess)

    class(Reentry_class)                :: this
    type(Neptune_class),intent(inout)   :: neptune

    character(len=160)      :: cdata   ! data string containing reentry information
    integer                 :: ierr
    real(dp), dimension(4)  :: temp_vec  ! containing re-entry epoch, altitude, latitude and longitude
   
    type(time_t) :: minCOIW_epoch ! start of COIW 
    type(time_t) :: maxCOIW_epoch ! end of COIW
    type(time_t) :: orig_epoch    ! originator epoch
    type(time_t) :: orbit_epoch   ! orbit epoch
    type(time_t) :: reentry_epoch ! re-entry epoch

    ! receive current time
    !--------------------------------------------
    orig_epoch = getDateTimeNowUtc()

    ! get orbit epoch
    !---------------------------------------------
    call neptune%getStartEpoch(orbit_epoch, ierr)

    ! get COIW epoch
    !---------------------------------------------
    temp_vec = neptune%derivatives_model%getLastPositionECEF()

    reentry_epoch%mjd = temp_vec(1)
    call mjd2gd(reentry_epoch)

    minCOIW_epoch%mjd = reentry_epoch%mjd - this%uncertaintyInterval*(reentry_epoch%mjd-orbit_epoch%mjd)
    maxCOIW_epoch%mjd = reentry_epoch%mjd + this%uncertaintyInterval*(reentry_epoch%mjd-orbit_epoch%mjd)

    call mjd2gd(minCOIW_epoch)
    call mjd2gd(maxCOIW_epoch)

100 format(a,';',5(a,';'),f6.2,';',f6.2,';',f6.2)
    write(cdata,100) coriginator, date2string(orig_epoch), date2string(orbit_epoch), &
                     date2string(minCOIW_epoch), date2string(reentry_epoch), date2string(maxCOIW_epoch), &
                     temp_vec(4)*rad2deg, temp_vec(3)*rad2deg, neptune%atmosphere_model%getMinAltitude()


    mess = new_line('a')//repeat('=',160)//new_line('a')//new_line('a')// &
           '  IADC Re-entry campaign - REDB format string:'//new_line('a')//new_line('a')// &
           '  Originator;Origination Epoch;Orbit Epoch;COIW start;COIW;COIW end;longitude/deg;latitude/deg;altitude/km'//new_line('a')// &
           repeat('-',160)//new_line('a')// &
           '  '//cdata//new_line('a')//new_line('a')//repeat('=',160)
           
  end function getReEntryMessage

!=============================================================================
!
!> @anchor      isInitialisedReEntry
!!
!> @brief       Gives the init flag for this module
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 16.10.2014 (initial design)</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  logical function isInitialisedReEntry(this)

    class(Reentry_class)    :: this
    isInitialisedReEntry = this%reEntryInitialised

  end function isInitialisedReEntry

end module


