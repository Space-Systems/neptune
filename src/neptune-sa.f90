!-----------------------------------------------------------------------------------------------
!
!> @brief NEPTUNE main routine (standalone)
!!
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  18.12.2012 (Geopotential + EOP)       </li>
!!                <li>VB:  24.01.2013 (Atmosphere)               </li>
!!                <li>VB:  28.02.2013 (Third bodies)             </li>
!!                <li>VB:  14.03.2013 (Solar radiation pressure) </li>
!!                <li>VB:  20.04.2013 (Interface/API redesign)   </li>
!!                <li>VB:  04.06.2013 (Code optimization)        </li>
!!                <li>VB:  10.11.2014 (Added timing)             </li>
!!                <li>VB:  31.01.2015 (added logo printing)      </li>
!!                <li>CHK: 16.11.2015 (updated to use libslam)   </li>
!!                <li>VB:  14.06.2016 (removed autoconfiguration functionality)   </li>
!!              </ul>
!!
!> @details     This is the main driver for the NEPTUNE (<b>N</b>PI <b>E</b>phemeris
!!              <b>P</b>ropagation <b>T</b>ool with <b>Un</b>certainty <b>E</b>stimation) routines.
!!              It reads the input file and performs the
!!              initialization of NEPTUNE as well as the propagation. The NEPTUNE library
!!              is thus used as a standalone propagation tool using this driver.
!!
!> @see         rdinp()
!> @see         libneptune::init_neptune()
!> @see         libneptune::propagate()
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      neptune_sa
!!
!!------------------------------------------------------------------------------------------------
program neptune_sa

  use slam_astro
  use derivatives,              only: PERT_SRP
  use neptune_error_handling,   only: getNeptuneErrorMessage, E_MIN_ALTITUDE, getLatestError, hasFailed, getLogfileChannel, & 
                                      setLogfileName, getLogfileName, initErrorHandler, ERRORS, setLogVerbosity, setCliVerbosity, setLogfileChannel, printTrace
  use slam_io,                  only: LOG_AND_STDOUT, openFile, STREAM, SEQUENTIAL, SWITCHED_OFF, IN_UNFORMATTED, closeFile, message
  use libneptune,               only: init_neptune, propagate
  use neptuneClass,             only: Neptune_class
  use slam_orbit_types,         only: state_t, covariance_t
  use reentry,                  only: Reentry_class
  use satellite,                only: Satellite_class
  use slam_time,                only: time_t, date2string
  use slam_timer,               only: startTimer, getElapsedTime, resetTimer
  use slam_types,               only: dp, sp
  use parsetle,                 only: readTLEfile, tle_t

  implicit none

    type(Neptune_class)           :: neptune                                    ! NEPTUNE instance
    type(Reentry_class)           :: reentry_model                              ! Reentry model

    character(len=255)            :: cmess                                      ! message string
    character(len=*), parameter   :: csubid = "neptune-sa"                      !  program id
    character(len=4)              :: ctemp4
    character(len=12)             :: ctemp
    character(len=200)            :: f_input                                    ! input file name
    character(len=200)            :: f_logo                                     ! logo file name

    integer                       :: ierr                                       ! error flag
    integer                       :: itemp                                      ! temporary integer

    real(dp)                      :: dtemp
    real(sp)                      :: totalTime
    real(sp), dimension(2)        :: timerArray

    type(state_t)                 :: state, state_out                           ! state vector
    type(covariance_t)            :: covar, covar_out                           ! covariance matrix
    type(time_t), dimension(2)    :: epoch                                      ! begin and end date of propagation
                                                                                !   1 = begin date
                                                                                !   2 = end date


  !========================================================================
  !
  ! Create NEPTUNE instance
  !
  !-----------------------------------------
  neptune = Neptune_class()
  reentry_model = Reentry_class()

  !========================================================================
  !
  ! Command Line Options
  !
  !----------------------------------------------------
  call parseCommandLine(neptune, reentry_model, f_input, f_logo)

  !========================================================================
  !
  ! Open log file
  !
  !----------------------------------------------------
  ierr = setLogfileChannel(openFile(getLogfileName(), SEQUENTIAL, -1))

  !========================================================================
  !
  ! Initialize error handling
  !
  !-----------------------------------------
  call initErrorHandler(control = 'YES', errAction = 'RETURN', traceback = 'YES')



  !========================================================================
  !
  ! Welcome message
  !
  !----------------------------------------------------
  call neptune%version_model%printNeptuneLogo(f_logo)

  call message('       '//neptune%version_model%getNeptuneAbbreviation(), LOG_AND_STDOUT)
  call message('                       ('//trim(neptune%version_model%to_string(neptune%version_model%get_neptune_version()))//' :: '//trim(neptune%version_model%getVersionDate(neptune%version_model%to_string(neptune%version_model%get_neptune_version())))//')', LOG_AND_STDOUT)
  call message(' _____________________________ IRAS -- ESA __________________________________', LOG_AND_STDOUT)
  call message (" ", LOG_AND_STDOUT)

  !========================================================================
  !
  ! Read input file (which actually sets NEPTUNE initialization parameters...)
  !
  !----------------------------------------------------
  call rdinp( neptune,    &      ! <-> CLASS  NEPTUNE class
              f_input,    &      ! <--  CHR() input file name
              state,      &      ! -->  TYP   initial cartesian state vector
              covar,      &      ! -->  TYP   initial covariance matrix
              epoch       &      ! -->  TYP   initial epoch
            )

  if(hasFailed()) then
    write(*,'(a,i3,a)') " Error No. ", getLatestError(), " occurred. Initialization stopped."
    write(*,*)          repeat("-",66)
    call getNeptuneErrorMessage(getLatestError(), cmess)
    write(*,'(a)')          " Error Message: "//trim(cmess)
    write(*,'(a)')          " +++ Program terminated! +++"
    stop
  end if

  !========================================================================
  !
  ! Propagation initialization
  !
  !----------------------------------------------------
  call message(' - Initializing propagation...', LOG_AND_STDOUT)
  ierr = init_neptune(neptune,state)

  if(hasFailed()) then
    write(ctemp4,'(i4)') getLatestError()
    write(*,'(a)') "  Error No. "//ctemp4//" occurred. NEPTUNE initialization stopped."
    write(*,'(2x,a)') repeat("-",66)
    call getNeptuneErrorMessage(getLatestError(), cmess)
    write(*,'(a)') "  Error Message: "//trim(cmess)
    write(*,'(a)') "  +++ Program terminated! +++"
    stop
  end if

  !** dump all the inputs
  call neptune%dump_input("another_dump.out")

  !========================================================================
  !
  ! Propagation
  !
  !----------------------------------------------------
  call message(' - Starting propagation...', LOG_AND_STDOUT)

  call propagate(          &
                neptune,   &   ! <-> CLASS  NEPTUNE class
                state,     &   ! <-- TYP    input state vector (GCRF)
                covar,     &   ! <-- TYP    input covariance matrix (GCRF)
                epoch,     &   ! <-- TYP    propagation start and end epoch
                               !            epoch(1) : start epoch
                               !            epoch(2) : end epoch
                state_out, &   ! --> TYP    output state vector
                covar_out, &   ! --> TYP    output covariance matrix
                .true.     &   ! <-- LOG    reset flag
              )

  if(hasFailed()) then
    if(reentry_model%isInitialisedReEntry() .and. getLatestError() == E_MIN_ALTITUDE) then ! re-entry mode, object re-entered
      write(*,*) reentry_model%getReEntryMessage(neptune)
    else
      write(ctemp4,'(i4)') getLatestError()
      write(*,'(a)') "  Error No. "//ctemp4//" occurred. Propagation stopped."
      write(*,'(2x,a)') repeat("-",66)
      call getNeptuneErrorMessage(getLatestError(), cmess)
      write(*,'(a,a)') "  Error Message: ", trim(cmess)
      call printTrace()
      write(*,'(a)')   "  +++ Program terminated! +++"
      write(*,'(2x,a)') repeat("-",66)
    end if
    stop
  end if

  call message(' NEPTUNE finished nominally.', LOG_AND_STDOUT)
  call message(' - Elapsed time for propagation: ', LOG_AND_STDOUT)
  call dtime(timerArray, totalTime) ! There seems to be an issue on ifort with using totalTime ...
  write(ctemp,'(f12.4)') (timerArray(1) + timerArray(2)) ! ... so we just use both of them
  call message('    * total time: '//trim(adjustl(ctemp))//' s', LOG_AND_STDOUT)
  write(ctemp,'(f12.4)') timerArray(1)
  call message('    * user time : '//trim(adjustl(ctemp))//' s', LOG_AND_STDOUT)
  write(ctemp,'(f12.4)') timerArray(2)
  call message('    * system time : '//trim(adjustl(ctemp))//' s', LOG_AND_STDOUT)
  if(neptune%satellite_model%getESHFlag()) then
    if(neptune%derivatives_model%getPertSwitch(PERT_SRP)) then
      write(*,'(A)') '  * Note that Equivalent Solar Hours (ESH) are only estimated for SRP perturbation being activated.'
    else
      write(*,'(A)') '  * Equivalent Solar Hours (ESH) for individual surfaces: '
      write(*,'(A)') '        _ID____ESH____'

      do itemp = 1, neptune%satellite_model%getSurfaceNumber()
        dtemp = neptune%satellite_model%getESH(itemp)
        write(*,'(a,i2,x,f10.2)')  '        ', itemp, dtemp
      end do
    end if
  end if

  !call destroy_neptune(neptune)

  ierr = setLogfileChannel(closeFile(getLogfileChannel()))

contains

!==============================================================================
!
!> @brief   Routine for the parsing of the optional command line parameters
!!
!> @author  Vitali Braun
!> @date    <ul>
!!            <li>VB:  04.02.2014 (initial design)</li>
!!            <li>VB:  14.06.2016 (removed autoconfiguration options)</li>
!!            <li>VB:  13.10.2016 (added new options to export and import input to and from file,
!!                                 and changed from iargc/getarg to 2003 syntax.)</li>
!!          </ul>
!!
!> @param[inout]  neptune         Neptune class instance
!> @param[inout]  reentry_model   The re-entry model
!> @param[inout]  f_input         Name of NEPTUNE input file
!> @param[inout]  f_logo          Name of NEPTUNE logo file
!!
!-------------------------------------------------------------------------------
  subroutine parseCommandLine(neptune, reentry_model, f_input, f_logo)
    use slam_time,           only: gd2mjd
    use slam_types,          only:
    !use neptuneInput,        only: read_input_from_dump, dump_input
    use slam_error_handling, only: isControlled, hasToReturn, hasFailed, checkIn, checkOut

    implicit none

    type(Neptune_class)             :: neptune                                  ! NEPTUNE instance
    type(Reentry_class)             :: reentry_model                            ! Reentry model
    character(len=*), intent(inout) :: f_input
    character(len=*), intent(inout) :: f_logo

    character(len=*), parameter :: csubid = 'parseCommandLine'
    character(len=200) :: cin       ! command line string containing one argument
    character(len=600) :: cmd       ! command line string containing the entire command
    character(len=20)  :: ctemp     ! temporary
    character(len=200) :: f_dump_input  ! dump file to read input from
    character(len=200) :: f_dump_output ! dump file to read input from

    integer :: cliStatus        ! command line parsing status flag
    integer :: i                ! loop counter
    integer :: ierr
    integer :: ios
    integer :: flag_skip        ! telling how many of the following command line options are skipped

    logical :: actParsNotGiven
    logical :: flag_dump_input  ! read neptune input dump
    logical :: flag_dump_output ! dump neptune input
    logical :: flag_input_file  ! neptune input file
    logical :: flag_log_file    ! neptune log file
    logical :: flag_logo_file   ! neptune logo file
    logical :: flag_verbosity   ! verbosity
    logical, dimension(8) :: parSet

    if(isControlled()) then
        if(hasToReturn()) return
        call checkIn(csubid)
    end if

    !** reset flags and set some initial parameters
    !--------------------------------
    flag_dump_input  = .false.
    flag_dump_output = .false.
    flag_input_file  = .false.
    flag_logo_file   = .false.
    flag_log_file    = .false.
    flag_skip        = 0
    flag_verbosity   = .false.

    actParsNotGiven = .true.
    parSet(:)       = .false.

    ierr = setLogVerbosity(3)
    !--------------------------------

    ! check first, if the program was launched with the -h (help) flag - if so, display CLI help screen and stop
    call get_command(command=cmd)
    if(index(cmd, ' -h ') /= 0 .or. index(cmd, ' --help ') /= 0) then
        call get_command_argument(0, cin)   ! get program name
        write(*,*) ' === NEPTUNE v'//trim(neptune%version_model%to_string(neptune%version_model%get_neptune_version()))//' -- '//trim(neptune%version_model%getVersionDate(neptune%version_model%to_string(neptune%version_model%get_neptune_version())))//' === '
        write(*,*) ' Usage: '//trim(cin)//' [OPTION]...'
        write(*,*) '   --export-dump FILE'
        write(*,*) '        Will dump all input parameters to a specified file (see also --import-dump)'
        write(*,*) '   -h   --help'
        write(*,*) '        Display this help'
        write(*,*) '   --import-dump FILE'
        write(*,*) '        Import input parameters that have been previously dumped (see --export-dump)'
        write(*,*) '   --log LOGFILE'
        write(*,*) '        The path for the logfile'
        write(*,*) '   --logo FILE'
        write(*,*) '        Provide the path for the NEPTUNE logo (ASCII), which is displayed in the CL'
        write(*,*) '        Example: --logo=data/neptune_logo.dat'
        write(*,*) '   --re-entry'
        write(*,*) '        For re-entry analysis, NEPTUNE will extract impact window information and paste it in the format used for IADC campaigns.'
        write(*,*) '   -v N '
        write(*,*) '        Set command line verbosity:'
        write(*,*) '            0 = error messages only'
        write(*,*) '            1 = errors and warnings'
        write(*,*) '            2 = errors, warnings and remarks'
        write(*,*) '            3 = all messages (incl. file open/close, debug)'
        write(*,*) '   --version'
        write(*,*) '        Display the NEPTUNE version'
        stop
    end if

    do i=1,command_argument_count()

        call get_command_argument(i,cin,status=cliStatus)
        if(cliStatus == -1) then
            write(ctemp,'(i2)') cliStatus
            write(*,*) 'The command line argument no. '//ctemp//' was truncated to: '
            write(*,*) '   "'//cin//'"'
            write(*,*) ' ->> Please fix in command line parsing (by providing more characters) or see if you can shorten it.'
            cycle
        else if(cliStatus /= 0) then
            write(*,*) 'There was an error with the command line argument no. '//ctemp//' ('//cin//'). Ignoring.'
            cycle
        end if

        if(flag_skip > 0) then
            flag_skip = flag_skip - 1
            cycle
        end if

        if(index(cin,'--version ') == 1) then ! '--version' writes current version of NEPTUNE
            write(*,*) "NEPTUNE v"//trim(neptune%version_model%to_string(neptune%version_model%get_neptune_version()))//" -- "//trim(neptune%version_model%getVersionDate(neptune%version_model%to_string(neptune%version_model%get_neptune_version())))
            stop

        else if(index(cin,'--export-dump ') == 1) then
            if(i+1 <= command_argument_count()) then
                call get_command_argument(i+1,f_dump_output)
                flag_skip        = 1 ! skip next option, as it is already processed
                flag_dump_output = .true.
            end if

        else if(index(cin,'--import-dump ') == 1) then
            if(i+1 <= command_argument_count()) then
                call get_command_argument(i+1,f_dump_input)
                flag_skip       = 1 ! skip next option, as it is already processed
                flag_dump_input = .true.
            end if

        else if(index(cin,'-f ') == 1) then    ! '-f neptune.inp' providing input file name to be used
            if(i+1 <= command_argument_count()) then
                call get_command_argument(i+1,f_input)
                flag_skip       = 1 ! skip next option, as it is already processed ('neptune.inp')
                flag_input_file = .true.
            end if

        else if(index(cin,'--logo ') == 1) then  ! '--logo=<LOGO_FILE>' providing logo file name to be used
            if(i+1 <= command_argument_count()) then
                call get_command_argument(i+1,f_logo)
                flag_skip      = 1
                flag_logo_file = .true.
            end if

        else if(index(cin,'-v ') == 1) then
            if(i+1 <= command_argument_count()) then
                call get_command_argument(i+1,ctemp)
                read(ctemp,*,iostat=ios) itemp
                if(ios/=0) then
                    write(*,*) "Unknown command line option: '-v "//trim(adjustl(ctemp))//"'."
                    write(*,*) "Program terminated."
                    stop
                else
                    ierr = setCliVerbosity(itemp)
                    flag_verbosity = .true.
                    flag_skip      = 1
                end if
            end if

        else if(index(cin,'--log ') == 1) then

            ierr = setLogfileName(cin(index(cin,"log")+4:len_trim(cin)))
            if(ierr /= 0) then
                write(*,*) "Logfile could not be initialized..."
                write(*,*) "Program terminated."
                stop
            end if
            flag_log_file = .true.

        else if(index(cin,'--re-entry ') == 1) then  ! IADC re-entry campaign mode providing COIW information after re-entry of object
            call reentry_model%initReEntry()

        else
            write(*,*) "Unknown command line option: '"//trim(adjustl(cin))//"'."
            write(*,*) "Program terminated."
            stop
        end if

    end do

    if(flag_dump_input) then
        call neptune%read_input_from_dump(f_dump_input)
        if(hasFailed()) return
    end if

    if(flag_dump_output) then
        call neptune%dump_input(f_dump_output)
    end if

    if(.not. flag_input_file) then
        f_input = "input/neptune.inp"
    end if

    if(.not. flag_logo_file) then
        f_logo = "data/neptune_logo_small.dat"
    end if

    if(.not. flag_verbosity) then
        ierr = setCliVerbosity(ERRORS)
    end if

    if(.not. flag_log_file) then
        ierr = setLogfileName("neptune.log")
        if(ierr /= 0) then
            write(*,*) "Logfile could not be initialized..."
            write(*,*) "Program terminated."
            stop
        end if
    end if

    !** done!
    if(isControlled()) then
        call checkOut(csubid)
    end if
    return

  end subroutine parseCommandLine

  end program neptune_sa
