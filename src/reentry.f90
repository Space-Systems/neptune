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
  use slam_time,    only: time_t, mjd2gd

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

    character(len=140)      :: cdata   ! data string containing reentry information
    integer, dimension(8)   :: date ! date string
    integer                 :: ierr
    real(dp), dimension(4)  :: temp_vec  ! containing re-entry epoch, altitude, latitude and longitude
   
    type(time_t) :: minCOIW_epoch ! start of COIW 
    type(time_t) :: maxCOIW_epoch ! end of COIW
    type(time_t) :: orig_epoch    ! originator epoch
    type(time_t) :: orbit_epoch   ! orbit epoch
    type(time_t) :: reentry_epoch ! re-entry epoch

    ! receive current time
    !--------------------------------------------
    call date_and_time(values=date)
    orig_epoch%year   = date(1)
    orig_epoch%month  = date(2)
    orig_epoch%day    = date(3)
    orig_epoch%hour   = date(5)
    orig_epoch%minute = date(6)
    orig_epoch%second = date(7)

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

100 format(5(x,i4.4,'/',i2.2,'/',i2.2,x,i2.2,':',i2.2,':',i2.2),x,f7.2,x,f6.2,x,i3,x,i4.4,'/',i2.2,'/',i2.2,x,i2.2,':',i2.2,':',i2.2)

    write(cdata,100) orig_epoch%year, orig_epoch%month, orig_epoch%day, orig_epoch%hour, orig_epoch%minute, int(orig_epoch%second), &
                     orbit_epoch%year, orbit_epoch%month, orbit_epoch%day, orbit_epoch%hour, orbit_epoch%minute, int(orbit_epoch%second), &
                     minCOIW_epoch%year, minCOIW_epoch%month, minCOIW_epoch%day, minCOIW_epoch%hour, minCOIW_epoch%minute, int(minCOIW_epoch%second), &
                     reentry_epoch%year, reentry_epoch%month, reentry_epoch%day, reentry_epoch%hour, reentry_epoch%minute, int(reentry_epoch%second), &
                     maxCOIW_epoch%year, maxCOIW_epoch%month, maxCOIW_epoch%day, maxCOIW_epoch%hour, maxCOIW_epoch%minute, int(maxCOIW_epoch%second), &
                     temp_vec(4)*rad2deg, temp_vec(3)*rad2deg, int(neptune%atmosphere_model%getMinAltitude()), &
                     reentry_epoch%year, reentry_epoch%month, reentry_epoch%day, reentry_epoch%hour, reentry_epoch%minute, int(reentry_epoch%second)

    mess = new_line('a')//repeat('=',160)//new_line('a')//new_line('a')// &
           '  IADC Re-entry campaign - REDB format string:'//new_line('a')//new_line('a')// &
           '  Originator;Origination Epoch;Orbit Epoch;COIW start;COIW;COIW end;longitude/deg;latitude/deg;altitude/km'//new_line('a')// &
           repeat('-',160)//new_line('a')// &
           '  '//coriginator//cdata//new_line('a')//new_line('a')//repeat('=',160)
           
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


