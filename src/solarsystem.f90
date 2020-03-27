!==================================================================================================
!!
!> @brief       Provides physical properties and ephemerides of bodies in the solar system
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  27.01.2014 (initial design)</li>
!!                <li>VB:  06.10.2015 (supporting all planets of the solar system)</li>
!!                <li>CHK: 13.11.2015 (updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 04.01.2017 (created Solarsystem_class)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions related to the
!!              bodies in our solar system.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      solarsystem
!!
!!------------------------------------------------------------------------------------------------
module solarsystem

  use neptune_error_handling, only: E_SOLAR_SYSTEM_INIT, E_WRONG_THIRD_BODY, E_THIRD_BODY_SERIES, E_CHEB_COEFF_MISSING, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, WARNING, FATAL, E_SPECIAL, checkIn, checkOut
  use slam_math,              only: undefined
  use slam_io,                only: LOG_AND_STDOUT, message
  use slam_time,              only: time_t, jd245, jd2000, mjd2gd, delta_AT
  use slam_types,             only: dp

  implicit none

    !** IDs (NOTE: the ID's are used as array indices, so there always should be an ID '1', which is
    !              the lowest an, while all the other ID's shall be subsequent natural numbers.)
    !-------------------------------------------------------
    integer, parameter :: ID_SUN     = 1
    integer, parameter :: ID_MOON    = 2
    integer, parameter :: ID_JUPITER = 3
    integer, parameter :: ID_VENUS   = 4
    integer, parameter :: ID_MARS    = 5
    integer, parameter :: ID_MERCURY = 6
    integer, parameter :: ID_SATURN  = 7
    integer, parameter :: ID_URANUS  = 8
    integer, parameter :: ID_NEPTUNE = 9

    integer, parameter :: n_supported_bodies = 9

    !** ephemerides model parameters
    !--------------------------------------------------
    integer, parameter :: DE_405 = 1    ! this one lasts up to 2050
    integer, parameter :: DE_421 = 2    ! this one lasts up to 2053
    integer, parameter :: DE_430 = 3    ! this one lasts from 1850 up to 2150
    character(len=*), parameter :: C_DE_405 = 'DE-405'
    character(len=*), parameter :: C_DE_421 = 'DE-421'
    character(len=*), parameter :: C_DE_430 = 'DE-430'
    character(len=*), dimension(3), parameter :: ephemModel     = (/C_DE_405, C_DE_421, C_DE_430/)
    character(len=*), dimension(3), parameter :: ephemModelFile = (/'de405.bsp', 'de421.bsp', 'de430.bsp'/)

    character(len=*), parameter :: C_FILE_LEAP  ='naif0012.tls'

    character(len=*), dimension(n_supported_bodies), parameter :: body_name = (/'SUN    ', 'MOON   ', 'JUPITER', 'VENUS  ', 'MARS   ', &
                                                                              'MERCURY', 'SATURN ', 'URANUS ', 'NEPTUNE'/)

    integer,parameter                                   :: cheby_degree = 20    ! Degree of the polynomials
    real(dp),dimension(:,:,:,:),allocatable             :: chebyshev_polynomials! Chebyshev polynomials for solar system bodies
    real(dp),dimension(:,:,:),allocatable               :: chebyshev_epochs     ! Chebyshev epochs for solar system bodies

    type, public :: Solarsystem_class

        integer                                     :: selectedModel
        logical, dimension(n_supported_bodies)      :: isAvailableState         ! indicating that position available for the supported body for given epoch
        real(dp)                                    :: last_time_mjd            ! saving the MJD of last call of a time conversion computation
                                                                                ! in order to prevent the same time conversions for the same point in time

        real(dp), dimension(n_supported_bodies)     :: mu                       ! mass parameter GM for all supported bodies (km**3/s**2)
        real(dp), dimension(n_supported_bodies,3)   :: state_gcrf               ! position of the body in the GCRF

        logical                                     :: solarSystemInitialized

    contains

        !** getter
        procedure           :: getBodyGM
        procedure           :: getBodyPosition
        procedure           :: getPlanetaryEphemFileName
        procedure           :: getLeapSecondFileName
        procedure           :: getSolarSystemInitFlag

        !** setter
        procedure           :: setSolarSystemInitFlag

        !** other
        procedure,public    :: initSolarSystem
        procedure,private   :: extractBodyInformation
        procedure           :: destroy

    end type Solarsystem_class

    ! Constructor
    interface Solarsystem_class
        module procedure constructor
    end interface Solarsystem_class

contains

    ! ====================================================================
    !!
    !> @brief      Constructor that should be used to initialize variables.
    !!
    !> @author     Christopher Kebschull
    !> @date       <ul>
    !!                  <li>ChK: 04.01.2018 (initial implementation)</li>
    !!              </ul>
    !> @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Solarsystem_class) function constructor()
        constructor%isAvailableState = .false.                                  ! indicating that position available for the supported body for given epoch
        constructor%solarSystemInitialized = .false.

    end function constructor

    !===========================================================================
    !!
    !> @anchor     destroy
    !!
    !> @brief      Destroys all that needs destruction
    !> @author     Christopher Kebschull
    !!
    !!
    !> @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Solarsystem_class)    :: this

        if (allocated(chebyshev_polynomials)) deallocate(chebyshev_polynomials)
        if (allocated(chebyshev_epochs)) deallocate(chebyshev_epochs)

    end subroutine destroy

! =======================================================================================================
!
!> @anchor      getLeapSecondFileName
!!
!> @brief       Returns the SPICE Kernel file leap seconds are read from
!> @author      Vitali Braun
!!
!> @return      file name as string
!!
!> @date        <ul>
!!                <li>VB: 24.10.2016 (initial design)</li>
!!              </ul>
!!
! -----------------------------------------------------------------------------------------
  character(len=len(C_FILE_LEAP)) function getLeapSecondFileName(this) result(fname)

    class(Solarsystem_class)    :: this

    fname = C_FILE_LEAP
    return
  end function


! =======================================================================================================
!
!> @anchor      getPlanetaryEphemFileName
!!
!> @brief       Returns the file planetary ephemerides are read from
!> @author      Vitali Braun
!!
!> @return      file name as string
!!
!> @date        <ul>
!!                <li>VB: 21.10.2016 (initial design)</li>
!!              </ul>
!
! -----------------------------------------------------------------------------------------
  character(len=len(ephemModelFile(1))) function getPlanetaryEphemFileName(this) result(fname)

    class(Solarsystem_class)    :: this

    if(.not. this%solarSystemInitialized) then
        fname = ''
    else
        fname = ephemModelFile(this%selectedModel)
    end if
    return
  end function

!===========================================================================
!
!> @anchor      getSolarSystemInitFlag
!!
!> @brief       Get initialization flag
!> @author      Vitali Braun
!!
!> @return      .true. when the init has happend
!!
!> @date        <ul>
!!                <li>VB: 27.01.2014 (initial design) </li>
!!              </ul>
!!
!!--------------------------------------------------------------------------
  logical function getSolarSystemInitFlag(this)

    class(Solarsystem_class)    :: this

    getSolarSystemInitFlag = this%solarSystemInitialized
    return
  end function getSolarSystemInitFlag


!===========================================================================
!
!> @anchor      setSolarSystemInitFlag
!!
!> @brief       Set initialization flag to given value
!> @author      Vitali Braun
!!
!! @param[in]   flag    Value which the init flag is set to
!!
!> @date        <ul>
!!                <li>VB: 27.01.2014 (initial design) </li>
!!              </ul>
!!
!!--------------------------------------------------------------------------
  subroutine setSolarSystemInitFlag(this,flag)

    class(Solarsystem_class)    :: this
    logical, intent(in)         :: flag

    this%solarSystemInitialized = flag
    return
  end subroutine setSolarSystemInitFlag


!===========================================================================
!
!> @anchor      initSolarSystem
!!
!> @brief       Initialization of solar system bodies
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li>VB: 28.02.2013 (initial design)</li>
!!                <li>VB: 03.06.2013 (code optimization) </li>
!!                <li>VB: 27.01.2014 (moved to 'solarsystem' module and renamed)</li>
!!                <li>VB: 06.10.2015 (all planets of the solar system are now initialized)</li>
!!              </ul>
!!
!! @param[in]   cpath         Path to read ephemeris data from
!! @param[in]   ephem_series  String identifier which DE-series to use, e.g. 'DE-421'
!! @param[in]   date_end      End date, identifying maximum epoch for ephemeris initialization
!!
!> @details     This routine initializes the solar system module. Global parameters are set and
!!              the SPICE kernel is loaded.
!!
!!------------------------------------------------------------------------------------------------
  subroutine initSolarSystem(               &
                              this,         &
                              cpath,        &  ! <-- CHR()   path to read data from
                              ephem_series, &  ! <-- CHR()   string identifier which DE-series to use, e.g 'DE-421'
                              date_end      &  ! <-- TYP     end date, identifying maximum epoch for ephemeris initialization
                            )
    use slam_io,    only: cdelimit, SEQUENTIAL, IN_FORMATTED, openFile, closeFile, nxtbuf
    !** interface
    !-------------------------------------------
    class(Solarsystem_class)                    :: this
    character(len=*),               intent(in)  :: cpath
    character(len=len(ephemModel)), intent(in)  :: ephem_series
    type(time_t),                   intent(in)  :: date_end
    !------------------------------------------
    integer                             :: ibody,isegment
    integer                             :: num_bodies,num_segments,final_segments
    character(len=1000)                 :: file_name
    integer                             :: file_unit
    logical                             :: lexists
    integer                             :: l
    real(dp)                            :: poly_start,poly_end,oldpoly_start
    real(dp),dimension(0:cheby_degree)  :: polys
    character(len=2000)                 :: cbuf                                 ! Line buffer
    character(len=9)                    :: cfile                                ! ephemeris data file name
    character(len=255)                  :: cmess                                ! message string
    character(len=255)                  :: de_file                              ! file name to load SPICE kernel (JPL ephemerides)
    character(len=255)                  :: leap_file                            ! file name to load SPICE kernel (leap seconds)
    character(len=*),parameter          :: csubid = "initSolarSystem"           ! subroutine ID

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%solarSystemInitialized) then   ! return if already initialized
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
    end if

    call message(' - Reading planetary ephemerides data...', LOG_AND_STDOUT)

    if(.not. any(ephemModel == ephem_series)) then
      call setNeptuneError(E_THIRD_BODY_SERIES, FATAL, (/ephem_series/))
      return
    end if

    !==========================================================================
    select case (ephem_series)
      case (C_DE_421)
        this%selectedModel = DE_421
        ! check max. epoch for DE-421, which is 2053-10-09 or JD = 2471184.5
        if(date_end%jd  >= 2471184.5d0) then
            cmess = "Third body ephemerides not available in JPL series '"//C_DE_421//"' for given epoch. Switching to '"//C_DE_405//"'."
            call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
            this%selectedModel = DE_405
        end if

      case (C_DE_430)
          !This one goes until 2150, maybe we need to do a check here as well
          this%selectedModel = DE_430
      case default
        ! DE_405 last only until  2050 JAN 01 00:01:51.623.
        ! maybe we need to check for that as well
        if(date_end%jd  >= 2469807.577083) then
          cmess = "Third body ephemerides not available in JPL series '"//C_DE_405//"' for given epoch. Crashing."
            call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
        end if

      end select


    !==========================================================================
    !
    ! set gravitational constants
    !
    !--------------------------------------------
    if(this%selectedModel == DE_421) then

      ! data from Folkner W.M., Williams J.G., Boggs D.H.,
      ! The Planetary and Lunar Ephemeris DE 421, Memorandum IOM 343R-08-003,
      ! Jet Propulsion Laboratory, 31 March 2008
      this%mu(ID_SUN)     = 132712440040.944d0
      this%mu(ID_MOON)    = 4902.800076d0
      this%mu(ID_JUPITER) = 126712764.8d0
      this%mu(ID_VENUS)   = 324858.592d0
      this%mu(ID_MARS)    = 42828.375214d0
      this%mu(ID_MERCURY) = 22032.09d0
      this%mu(ID_SATURN)  = 37940585.2d0
      this%mu(ID_URANUS)  = 5794548.6d0
      this%mu(ID_NEPTUNE) = 6836535.d0

      cfile = ephemModelFile(DE_421)

    else if(this%selectedModel == DE_405) then

      ! data from Standish E.M., JPL Planetary and Lunar Ephemerides DE405/LE405,
      ! Memorandum IOM 312.F-98-048, Jet Propulsion Laboratory, 26 August 1998
      this%mu(ID_SUN)     = 132712440017.987d0
      this%mu(ID_MOON)    = 4902.801d0
      this%mu(ID_JUPITER) = 126712767.863d0
      this%mu(ID_VENUS)   = 324858.599d0
      this%mu(ID_MARS)    = 42828.314d0
      this%mu(ID_MERCURY) = 22032.08d0
      this%mu(ID_SATURN)  = 37940626.063d0
      this%mu(ID_URANUS)  = 5794549.007d0
      this%mu(ID_NEPTUNE) = 6836534.064d0

      cfile = ephemModelFile(DE_405)

    else if (this%selectedModel == DE_430) then

      ! data from https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf
      ! William M. Folkner, James G. Williams, Dale H. Boggs, Ryan S. Park, and Petr Kuchynka
      ! IPN Progress Report 42-196 February 15, 2014
      ! Table 8, page 49
      this%mu(ID_SUN)     = 132712440041.939400d0
      this%mu(ID_MOON)    = 4902.800066d0
      this%mu(ID_JUPITER) = 126712764.800000d0
      this%mu(ID_VENUS)   = 324858.592000d0
      this%mu(ID_MARS)    = 42828.375214d0
      this%mu(ID_MERCURY) = 22031.780000d0
      this%mu(ID_SATURN)  = 37940585.200000d0
      this%mu(ID_URANUS)  = 5794548.600000d0
      this%mu(ID_NEPTUNE) = 6836527.100580d0

      cfile = ephemModelFile(DE_430)

    end if

    !==========================================================================
    !
    ! load SPICE kernel
    !
    !---------------------------------------------

    !** build data file string
    de_file   = trim(adjustl(cpath))//cdelimit//cfile
    leap_file = trim(adjustl(cpath))//cdelimit//C_FILE_LEAP

    ! Switching from SPICE to interpolatino of body positions
    !** load kernel
    !call furnsh(leap_file)
    !call furnsh(de_file)

    num_bodies = 0
    final_segments = 0
    do ibody=1, n_supported_bodies
        file_name = trim(adjustl(cpath))//cdelimit//ephem_series//"_Earth_J2000_"//trim(body_name(ibody))//".dat"
        inquire (file = file_name , exist = lexists)
        ! Determine whether file is available
        if (lexists) then
            num_bodies = num_bodies + 1
            ! Determine the number of available segments
            file_unit = openFile(file_name, SEQUENTIAL, IN_FORMATTED)
            num_segments = 0
            oldpoly_start = 0.0
            do ! Read station data
                call nxtbuf('#', 0, file_unit, cbuf)
                ! Read until the end of the file
                if (.not. cbuf .eq. '') then
                    read(cbuf,*) l, poly_start, poly_end, polys(0:cheby_degree)
                    if (oldpoly_start /= poly_start) then
                        num_segments = num_segments + 1
                    end if
                    oldpoly_start = poly_start
                else
                    exit
                end if
            end do
            ! Choose the largest segment count for the array
            if (num_segments > final_segments) then
                final_segments = num_segments
            end if
            file_unit = closeFile(file_unit)
        end if
    end do

    ! Allocate the chebyshev polynomial coefficient and epoch arrays
    if (final_segments > 0 .and. num_bodies > 0) then
        if (allocated(chebyshev_polynomials)) deallocate(chebyshev_polynomials)
        if (allocated(chebyshev_epochs)) deallocate(chebyshev_epochs)
        allocate(chebyshev_polynomials(1:n_supported_bodies,1:final_segments,0:cheby_degree,1:3))
        allocate(chebyshev_epochs(1:n_supported_bodies,1:final_segments,1:2))
        chebyshev_polynomials(:,:,:,:) = 0.d0
        chebyshev_epochs(:,:,:) = 0.d0
    else
        call setNeptuneError(E_CHEB_COEFF_MISSING, FATAL)
        return
    end if

    ! Now lets read the polynomial coefficients
    do ibody=1, n_supported_bodies
        file_name = trim(adjustl(cpath))//cdelimit//ephem_series//"_Earth_J2000_"//trim(body_name(ibody))//".dat"
        inquire (file = file_name , exist = lexists)
        ! Determine whether file is available
        if (lexists) then
            ! Determine the number of available segments
            file_unit = openFile(file_name, SEQUENTIAL, IN_FORMATTED)
            isegment = 1
            oldpoly_start = 0
            l = 1
            do ! Read station data
                call nxtbuf('#', 0, file_unit, cbuf)
                ! Read until the end of the file
                if (.not. cbuf .eq. '') then
                    read(cbuf,*) l, chebyshev_epochs(ibody,isegment,1), chebyshev_epochs(ibody,isegment,2), chebyshev_polynomials(ibody,isegment,0:cheby_degree,l)
                    l = l + 1
                    if (l > 3) then
                        l = 1
                        isegment = isegment + 1
                    end if
                else
                    exit
                end if
            end do
            file_unit = closeFile(file_unit)
        end if
    end do

    this%last_time_mjd          = -1.d0
    this%solarSystemInitialized = .true.

    ! This is just a temporary operation to extract the position of the solarsystem
    ! bodies.
    !call this%extractBodyInformation(this%selectedModel)
    !stop

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine initSolarSystem



!===========================================================================
!
!> @anchor      extractBodyInformation
!!
!> @brief       Extract solar system body information from kernels
!> @author      Christopher Kebschull
!!
!> @param[in]   selectedModel  indicates which kernels to load
!!
!>  @date       <ul>
!!                <li>23.02.2018 (initial implementation)</li>
!!              </ul>
!!
!> @details     In order to get rid of SPICE we extract the body information directly relative to Earth
!!  and write them into data files per solar system body
!! retrieving from DE_430 from Jan. 1st, 1961 to Oct 09, 2053
!! JD_S = 2437300.5d0
!! JD_E = 2471183.5d0
!! Using 32 days for Chebyshev interpolation (from https://ipnpr.jpl.nasa.gov/progress_report/42-196/196C.pdf)
!!
!!------------------------------------------------------------------------------------------------
  subroutine extractBodyInformation(this, selectedModel)
    use slam_io,                only: SEQUENTIAL, OUT_FORMATTED_OVERWRITE, openFile, closeFile
    use slam_interpolation,     only: get_chebyshev_polynomials, chebyshev_interpolation
    use slam_strings,           only: toString

    class(Solarsystem_class)    :: this
    integer                     :: selectedModel
    character(len=6)            :: ephFile
    integer                     :: segment_interval                             ! Interval of the chebyshev polynomials coefficients is days
    integer                     :: ibody,iday,isegment,max_data
    real(dp)                    :: rday
    character(len=1000)         :: file_name                                    !< file name
    integer                     :: file_unit
    real(dp)                    :: start_jd = 2437300.5d0
    real(dp)                    :: end_jd   = 2471183.5d0
    !real(dp)                    :: start_mjd,end_mjd,step_mjd
    real(dp)                    :: time_mjd
    real(dp), dimension(3)      :: pos3
    character(len=80)           :: csep                                         ! separator
    character(len=8)            :: cdate                                        ! system time
    character(len=10)           :: ctime                                        ! system time
    integer                     :: steps
    real(dp),dimension(:,:),allocatable     :: interpol_state
    real(dp),dimension(:),allocatable       :: interpol_epoch
    real(dp), dimension(:,:),allocatable    :: interpol_ephemeris
    integer                                 :: l
    real(dp),allocatable,dimension(:,:)     :: cheby_polynomial                 ! Chebyshev polynomials

    select case (selectedModel)
        case (DE_405)
            ephFile = C_DE_405
        case (DE_421)
            ephFile = C_DE_421
        case (DE_430)
            ephFile = C_DE_430
    end select

    csep = '#'//repeat('-',79)

    ! Interval for polynomial fit is 32 days
    segment_interval = 32

    do ibody=1, n_supported_bodies
        ! There is no data on Jupiter prior to the year 2001, we need to account for that
        if (ibody == ID_JUPITER) then
            start_jd = 2451910.5d0  ! Jan 1st, 2001
        else
            start_jd = 2437300.5d0  ! Jan 1st, 1961
        end if
        file_name = ephFile//"_Earth_J2000_"//trim(body_name(ibody))//".dat"
        file_unit = openFile(file_name, SEQUENTIAL, OUT_FORMATTED_OVERWRITE)
        ! Write a header here
        !** Get current date and time
        call date_and_time(cdate,ctime)
        write(file_unit,'(A)')  "#"
        write(file_unit,'(A)')  "#     _/     _/  _/_/_/_/  _/_/_/   _/_/_/_/_/  _/   _/  _/      _/  _/_/_/_/"
        write(file_unit,'(A)')  "#    _/_/   _/  _/        _/    _/     _/      _/   _/  _/_/    _/  _/"
        write(file_unit,'(A)')  "#   _/ _/  _/  _/_/_/    _/_/_/       _/      _/   _/  _/  _/  _/  _/_/_/"
        write(file_unit,'(A)')  "#  _/   _/_/  _/        _/           _/      _/   _/  _/    _/_/  _/"
        write(file_unit,'(A)')  "# _/     _/  _/_/_/_/  _/           _/        _/_/   _/      _/  _/_/_/_/"
        write(file_unit,'(A)')  "#"
        write(file_unit,'(A)') trim(csep)
        write(file_unit,'(A)') "# This file was generated by the program NEPTUNE on"
        write(file_unit,'(A)') "# "//cdate(1:4)//"-"//cdate(5:6)//"-"//cdate(7:8)//"T"//ctime(1:2)//":"//ctime(3:4)//":"//ctime(5:6)
        write(file_unit,'(A)') trim(csep)
        write(file_unit,'(A)') "# Extraction of "//trim(body_name(ibody))//" position realtive to Earth"
        write(file_unit,'(A)') "#   from JPL Planetary and Lunar Ephemerides "//ephFile//"."
        write(file_unit,'(A)') "# Data presented in the GCRF frame using chebyshev polynomials using degree "//toString(cheby_degree)//"."
        write(file_unit,'(A)') trim(csep)

        do isegment = 0, int(end_jd - start_jd), segment_interval

            ! Adjust the segment interval size toward the end... 
            !   there might by less than 32 days left in the last segment
            if (start_jd + isegment + segment_interval > end_jd) then
                segment_interval = segment_interval - int(start_jd + isegment + segment_interval - end_jd)
            end if

            steps = 10
            max_data = segment_interval * steps

            allocate(interpol_state(1:max_data,1:3))
            allocate(interpol_epoch(1:max_data))
            allocate(interpol_ephemeris(1:max_data,1:2))
            allocate(cheby_polynomial(0:cheby_degree,1:3))

            interpol_state(:,:) = 0.d0
            interpol_epoch(:)   = 0.d0

            do iday=0, max_data -1
                ! Go through the time
                rday = dble(iday)/dble(max_data-1)*(segment_interval)
                time_mjd = start_jd + isegment + rday - 2400000.5d0
                pos3 = this%getBodyPosition(time_mjd,ibody)
                interpol_state(iday+1,1:3) = pos3(1:3)
                interpol_epoch(iday+1) = time_mjd
                !write(file_unit,'(F12.5,1X,F18.6,1X,F18.6,1X,F18.6)') time_mjd, pos3
            end do
            interpol_ephemeris(:,1) = interpol_epoch(:)
            do l = 1,3
                interpol_ephemeris(:,2) = interpol_state(:,l)
                call chebyshev_interpolation(cheby_degree,5,max_data,interpol_ephemeris(:,1:2),cheby_polynomial(:,l))
                write(file_unit,'(I1,2(1X,F12.5),31(1X,F18.6))') l, interpol_epoch(1), interpol_epoch(max_data), cheby_polynomial(:,l)
            end do

!            ! Checking that it actually works
!            start_mjd = start_jd + isegment -  2400000.5d0
!            end_mjd = start_mjd + segment_interval
!            step_mjd = start_mjd + (segment_interval)/dble(2)*dble(0)
!            pos3 = interpolate(step_mjd, start_mjd, end_mjd, cheby_degree, cheby_polynomial)
!            write (*,'(A24,F10.4,3(1X,F18.6))') "Interpolated (0) state error: ", step_mjd, pos3-this%getBodyPosition(step_mjd,ibody)
!            step_mjd = start_mjd + (segment_interval)/dble(2)*dble(1)
!            pos3 = interpolate(step_mjd, start_mjd, end_mjd, cheby_degree, cheby_polynomial)
!            write (*,'(A24,F10.4,3(1X,F18.6))') "Interpolated (1) state error: ", step_mjd, pos3-this%getBodyPosition(step_mjd,ibody)
!            step_mjd = start_mjd + (segment_interval)/dble(2)*dble(2)
!            pos3 = interpolate(step_mjd, start_mjd, end_mjd, cheby_degree, cheby_polynomial)
!            write (*,'(A24,F10.4,3(1X,F18.6))') "Interpolated (2) state error: ", step_mjd, pos3-this%getBodyPosition(step_mjd,ibody)

            deallocate(interpol_state)
            deallocate(interpol_epoch)
            deallocate(interpol_ephemeris)
            deallocate(cheby_polynomial)

        end do
        file_unit = closeFile(file_unit)

    end do

  end subroutine extractBodyInformation

!===========================================================================
!
!> @anchor      getBodyGM
!!
!> @brief       Returns the mass parameter GM for a solar system body
!> @author      Vitali Braun
!!
!! @param[in]   id    Body ID for which the GM is requested
!!
!> @return      GM
!!
!> @date        <ul>
!!                <li>VB: 17.07.2013 (initial design) </li>
!!                <li>VB: 27.01.2014 (moved to 'solarsystem' module and renamed after merging of getSunGM and getMoonGM)
!!              </ul>
!!
!!--------------------------------------------------------------------------
  real(dp) function getBodyGM(this, id)

    class(Solarsystem_class)    :: this
    integer, intent(in)         :: id

    character(len=255) :: cmess
    character(len=*), parameter :: csubid = 'getBodyGM'

    getBodyGM = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(id > size(this%mu) .or. id < 0) then
      write(cmess,'(a,i3)') "id: ", id
      call setNeptuneError(E_WRONG_THIRD_BODY, FATAL, (/cmess/))
      return
    end if

    getBodyGM = this%mu(id)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getBodyGM

!!------------------------------------------------------------------------------------------------
!> @anchor      getBodyPosition
!!
!> @brief       Computes the position of a solar system body in GCRF
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li>VB: 12.03.2013 (initial design 'getSunPosition')   </li>
!!                <li>VB: 03.06.2013 (code optimization)                 </li>
!!                <li>VB: 17.07.2013 (extended to also include the moon) </li>
!!                <li>VB: 23.01.2014 (fixed bug, which did not allow to get Moon's position when working without lunar perturbations.</li>
!!                <li>VB: 27.01.2014 (moved to 'solarsystem' module and renamed) </li>
!!                <li>VB: 15.05.2016 (added hasFailed checks for date conversions) </li>
!!              </ul>
!!
!> @param[in]   time_mjd    current time MJD
!> @param[in]   body_id     ID of requested body
!!
!> @return      body position / km
!!
!> @details     This routine computes the position of (currently) the sun or the moon. It uses the SPICE routine
!!              SPKPOS to compute the position in the GCRF frame (called J2000 in SPICE)
!!              for the given time. If this routine is called by an external routine
!!              and the Sun's or the Moon's position is already available for the given epoch,
!!              then it is returned without a new computation
!!
!!------------------------------------------------------------------------------------------------
  function getBodyPosition(this, time_mjd, body_id)

    !** interface
    !-------------------------------------------
    class(Solarsystem_class)    :: this
    real(dp), intent(in)        :: time_mjd
    integer,  intent(in)        :: body_id

    real(dp), dimension(3)      :: getBodyPosition
    !-------------------------------------------

    character(len=*), parameter :: csubid = "getBodyPosition"
    character(len=255)          :: cmess

    integer :: day          ! day
    integer :: jerr         ! error flag
    integer :: mon          ! month
    integer :: year         ! year

    real(dp)  :: dat        ! delta atomic time
    real(dp)  :: dfrac      ! day fraction
    !real(dp)  :: dtemp     ! temporary
    real(dp)  :: TT_j2000   ! TT in seconds since J2000.0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(.not. this%solarSystemInitialized) then
      call setNeptuneError(E_SOLAR_SYSTEM_INIT, FATAL)
      return
    end if

    if(body_id < 0 .or. body_id > n_supported_bodies) then
      write(cmess,'(a,i2,a)') "Unknown celestial body with id '", body_id, "'."
      call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
      return
    end if

    if(time_mjd < 0 .or. time_mjd > 115178) then ! Data available till epoch 2174 MAR 23 04:46:14.022.
      write(cmess,'(a,f10.1,a)') "Insufficient data from spice kernel for mjd '", time_mjd, "'."
      call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
      return
    end if

    if(time_mjd == this%last_time_mjd .and. this%isAvailableState(body_id)) then
      getBodyPosition = this%state_gcrf(body_id,:)
    else
      this%isAvailableState = .false.

      !** get terrestrial time TT
      !-----------------------------------------------
      call mjd2gd(time_mjd, year, mon, day, dfrac)
      if(hasFailed()) return

      call delta_AT (year, mon, day, dfrac, dat, jerr)
      if(hasFailed()) return

      !** compute TT as seconds from J2000.0:
      TT_j2000 = (time_mjd + jd245 - jd2000)*86400.d0 + dat + 32.184d0
      !------------------------------------------------------------------
      !** compute geometric state of body relative to earth
      !call spkpos(trim(body_name(body_id)), TT_j2000, 'J2000', 'NONE', 'EARTH', this%state_gcrf(body_id,:), dtemp)
      !write (*,*) "Delta: ", this%state_gcrf(body_id,:) - interpolateBodyPosition(time_mjd,body_id)
      this%state_gcrf(body_id,:) = interpolateBodyPosition(time_mjd,body_id)
      this%isAvailableState(body_id) = .true.

      getBodyPosition = this%state_gcrf(body_id,:)

    end if

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getBodyPosition

    !=========================================================================
    !!
    !> @anchor     interpolateBodyPosition
    !!
    !! @brief      Function use chebyshev polynomials to interpolate an orbit
    !!
    !! @author     Christopher Kebschull
    !!
    !! @date       <ul>
    !!                <li>23.02.2018 (initial implementation)</li>
    !!              </ul>
    !!
    !> @param[in]  stepEpoch - interpolation epoch [mjd]
    !!
    !> @param[in]  body_id - identifies the solar system body
    !!
    !> @result state - a state vector [km;km/s]
    !!
    !-------------------------------------------------------------------------
    function interpolateBodyPosition(stepEpoch,body_id)
        use slam_strings,   only: toString

        real(dp),intent(in)     :: stepEpoch
        integer,intent(in)      :: body_id
        real(dp),dimension(3)   :: interpolateBodyPosition

        integer                 :: num_segments,segment_index
        real(dp)                :: first_poly_epoch,last_poly_epoch
        character(len=255)      :: cmess                                        ! message string
        real(dp),dimension(:,:),allocatable  :: coeff

        num_segments = size(chebyshev_epochs,2)
        first_poly_epoch = chebyshev_epochs(body_id,1,1)
        last_poly_epoch = chebyshev_epochs(body_id,num_segments,2)
        if (stepEpoch >= first_poly_epoch .and. stepEpoch <= last_poly_epoch) then
            ! determine segment in the array, where the coefficients are stored for the requested time
            ! Polynomial fit span is 32 days
            segment_index = int((stepEpoch - first_poly_epoch)/32) + 1
            allocate(coeff(0:cheby_degree,1:3))
            coeff(0:cheby_degree,1:3) = chebyshev_polynomials(body_id,segment_index,0:cheby_degree,1:3)
            interpolateBodyPosition = interpolate(                                                      &
                                                    stepEpoch,                                          &
                                                    chebyshev_epochs(body_id,segment_index,1),          &
                                                    chebyshev_epochs(body_id,segment_index,2),          &
                                                    cheby_degree,                                       &
                                                    coeff                                               &
                                                    )
            deallocate(coeff)
        else
            cmess = "Third body ephemerides for "//body_name(body_id)//" cannot be interpolated for given epoch "//toString(stepEpoch)//"."
            call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
        end if

    end function

    !=========================================================================
    !!
    !>  @anchor     interpolate
    !!
    !!  @brief      Function use chebyshev polynomials to interpolate an orbit
    !!
    !!  @author     Christopher Kebschull
    !!
    !!  @date       <ul>
    !!                <li>23.02.2018 (initial implementation)</li>
    !!              </ul>
    !!
    !!  @param[in]  stepEpoch - interpolation epoch [mjd]
    !!
    !!  @param[in]  startEpoch - start epoch of the interpolation frame [mjd]
    !!
    !!  @param[in]  endepoch -  end epoch of the interpolation frame [mjd]
    !!
    !!  @param[in]  cheby_degree -  Degree of the polynomial [-]
    !!
    !!  @param[in]  cheby_coefficients -  Coefficients of the polynomials [-]
    !!
    !!  @result state - a state vector [km]
    !!
    !-------------------------------------------------------------------------
    function interpolate(stepEpoch, startEpoch, endEpoch, cheby_degree, cheby_coefficients) result(state)
        use slam_types,             only: dp
        use slam_interpolation,     only: get_chebyshev_polynomials

        implicit none

        ! function parameters
        real(dp),intent(in)                             :: stepEpoch
        real(dp),intent(in)                             :: startEpoch
        real(dp),intent(in)                             :: endEpoch
        integer,intent(in)                              :: cheby_degree
        real(dp),intent(in),dimension(:,:),allocatable  :: cheby_coefficients

        ! working variables
        real(dp),dimension(3)                           :: state
        real(dp)                                        :: xx
        integer                                         :: l,p
        real(dp),dimension(:),allocatable               :: tcheby

        ! Retrieve interpolation details
        allocate(tcheby(0:cheby_degree))
        ! Loop over state vector
        xx = (2.d0 * stepEpoch - (startEpoch + endEpoch)) &
                / (endEpoch - startEpoch)                                       ! transform to interval [-1, 1]
        do l=1,3
            call get_chebyshev_polynomials(cheby_degree+1, xx, tcheby)
            state(l) = 0.d0
            do p = 0, cheby_degree
                state(l) = state(l) + cheby_coefficients(p,l)*tcheby(p)
            end do
        end do

        deallocate(tcheby)

    end function interpolate

end module solarsystem
