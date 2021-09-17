!===============================================================================================
!
!> @anchor      version
!!
!> @brief       Providing the versioning for NEPTUNE
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li>VB:  04.02.2014 (initial implementation)</li>
!!                <li>VB:  10.04.2015 (changed versioning for beta version to scheme 0.9.x with x being incremented)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!              </ul>
!!
!! @details     This module provides the functions and parameters required for the version
!!              handling in NEPTUNE. By calling the function 'getNeptuneVersion' one obtains
!!              the current version number. The function 'getVersionDate(version)' provides the date
!!              associated with the given version.
!!              A 'semantic versioning' is used within NEPTUNE, which means, that each version
!!              consists of a three-digit number: major.minor.patch
!!              The patch number is incremented for changes and bug-fixes, which do not
!!              change the API of NEPTUNE. The minor version is incremented for releases which add
!!              new, but backward-compatible, API features, the major version is incremented for
!!              API changes, which are not backward-compatible (http://semver.org/spec/v2.0.0.html).
!!
!!              Version 1.1.2-release (17-SEP-2021)
!!              <ul>
!!                <li> Feature: Support of ITRF input state</li>
!!                <li> Feature: Support updated JB2008 solmag files</li>
!!                <li> Fix: integrator infinite loop when overstepping requested epoch</li>
!!              </ul>
!!              Version 1.1.1-release (15-APR-2021)
!!              <ul>
!!                <li> Feature: Support for Flang (AMD) compiler</li>
!!                <li> Feature: OpenMP works with gfortran-8</li>
!!              </ul>
!!              Version 1.1.0-release (07-MAR-2021)
!!              <ul>
!!                <li> Feature: JB2008 atmopsheric model</li>
!!                <li> Fix: Covariance propagation shows proper behaviour over longer timespans</li>
!!              </ul>
!!              Version 1.0.0-release (26-MAR-2020)
!!              <ul>
!!                <li> Note: Cleaned code</li>
!!                <li> Note: Add more documentation</li>
!!              </ul>
!!              Version 0.9.18-beta (30-APR-2018)
!!              <ul>
!!                <li> Feature: Object oriented NEPTUNE core</li>
!!                <li> Feature: Contains back-ported NRLMSISE-00 from a ported C-version ;-)</li>
!!                <li> Feature: No longer depends on SPICE for SUN and MOON ephemerides</li>
!!                <li> Feature: New data files with Chebyshev coefficients for SUN and MOON ephemerides</li>
!!                <li> Feature: Reads all data files entirely - no need for re-reading any files post initialization</li>
!!                <li> Feature: Contains new OpenMP parallelized test-executable openmp-test-sa.f90</li>
!!                <li> Note: Needs new libslam with object oriented Reduction_class</li>
!!                <li> Note: Removed PFUnit for now as there were isssues downloading the source several times</li>
!!                <li> Note: Deactivate trace-back when using OpenMP as libslam cannot handle multiple calls to the error handling yet</li>
!!              </ul>
!!
!!              Version 0.9.17-beta (27-NOV-2017)
!!              <ul>
!!                <li> Feature: DE-430 kernel support for gravitatal perturbation modeling</li>
!!                <li> Feature: Ability to use from the API</li>
!!                <li> Feature: Tests based on PFunit, which can be run before releasing </li>
!!                <li> Fixed: Several issues stemming from the switch to a 0.0s start point in the integrator</li>
!!                <li> Fixed: Benchmark output after running neptune-sa updated with more information and support for ifort</li>
!!              </ul>
!!
!!              Version 0.9.16-beta (28-APR-2017)
!!              <ul>
!!                <li> Fixed: checking for optional parameter in neptune error handling</li>
!!              </ul>
!!
!!              Version 0.9.15-beta (12-APR-2017)
!!              <ul>
!!                <li> Feature: input dump and export and import</li>
!!                <li> Feature: progress file writing</li>
!!                <li> Fixed: error handling for error codes delegated to SLAM</li>
!!              </ul>
!!
!!              Version 0.9.14-beta (14-JUN-2016)
!!              <ul>
!!                <li> Feature: removed autoconfiguration, as it is not used any longer</li>
!!                <li> Feature: backward propagation is now working (in principle - thorough testing still required...)</li>
!!              </ul>
!!
!!              Version 0.9.13-beta (26-JAN-2016)
!!              <ul>
!!                <li> Fixed: recovery after reset of integrator - not jumping back to last output state, if not required.</li>
!!                <li> Feature: removed some sections which passed non-contiguous array sections (performance issue)</li>
!!              </ul>
!!
!!              Version 0.9.12-beta (24-JAN-2016)
!!              <ul>
!!                <li> Fixed: string length check in neptuneInput for file names</li>
!!                <li> Fixed: deactivated CSSI space weather file read, due to reported bug, for the time being.</li>
!!                <li> Fixed: for fap_day read, if Kp's end, daily Ap value will be used (constant for all 3hr coefficients)</li>
!!                <li> Fixed: the routine getDensityMSIS2000 now uses correct values for latitude and longitude</li>
!!              </ul>
!!
!!              Version 0.9.11-beta (24-JAN-2016)
!!              <ul>
!!                <li> neptune() now checks the reference frames for both, state and covariance</li>
!!                <li> Fixed: for numint log, the header wrote an uninitialised string for the covariance matrix. Now there is a NOTE, if this happens.</li>
!!                <li> Fixed: issue in varstormcow, when requested time was exactly equal to integration time. Although a very unlikely case, this happened in simulations</li>
!!                     and NEPTUNE was caught in an infinite loop.</li>
!!                <li> Fixed: counter of number of calls to the numerical integration routines are now handled by module variables. This makes sure, that the count of the
!!                     integration steps starts at 1 again, as soon as NEPTUNE is called with reset=.true.</li>
!!              </ul>
!!              Version 0.9.10-beta (20-JAN-2016)
!!              <ul>
!!                <li> Fixed: uninitialised value in epoch check for start and end epoch</li>
!!                <li> Fixed: in numint.f90 the integration logfile is now only closed if it was open before</li>
!!                <li> Upgrade: geopotential can now be run with individual harmonics switched on or off</li>
!!              </ul>
!!
!!              Version 0.9.9-beta (11-JAN-2016)
!!              <ul>
!!                <li> Fixed: start and epoch reset now leads also to re-initialisation of EOP and SOLMAG arrays</li>
!!                <li> Fixed: leap second for June 2015</li>
!!                <li> Fixed: integration logfile is now closed after neptune call</li>
!!              </ul>
!!
!!              Version 0.9.8-beta (28-AUG-2015)
!!              <ul>
!!                <li> Corrected timer for runtime output: now user and total time are provided.</li>
!!              </ul>
!!
!!              Version 0.9.7-beta (23-AUG-2015)
!!              <ul>
!!                <li> Introducing an exponential atmosphere model - as well as the option to select the atmosphere model.</li>
!!              </ul>
!!
!!              Version 0.9.6-beta (08-JUN-2015)
!!              <ul>
!!                <li> Introducing the framework for a new integrator (the latter TBD).</li>
!!              </ul>
!!
!!              Version 0.9.5-beta (17-MAY-2015)
!!              <ul>
!!                <li> Fixed in Albedo acceleration routine: handling cases, where r < R_E (close to re-entry).</li>
!!              </ul>
!!
!!              Version 0.9.4-beta (04-MAY-2015)
!!              <ul>
!!                <li> Numerical integration logfile now provides two additional flags, indicating whether integrator call resulted in an integration or just an interpolatin.</li>
!!              </ul>
!!
!!              Version 0.9.3-beta (01-MAY-2015)
!!              <ul>
!!                <li> Fixed some bugs in covariance propagation related to reference frames.</li>
!!              </ul>
!!
!!              Version 0.9.2-beta (20-APR-2015)
!!              <ul>
!!                <li> Added covariance matrix store in neptuneGet module.</li>
!!              </ul>
!!
!!              Version 0.9.1-beta (10-APR-2015)
!!              <ul>
!!                <li> Fixed bug in varstormcow routine, where the weights wts/wtd in the stepsize control were (incorrectly) computed with lastvel/lastpos..</li>
!!              </ul>
!!
!!              Version 1.0.0-beta (06/2014)
!!              <ul>
!!                <li> Feature-complete.
!!                <li> Added Lundberg algorithm for SRP corrections in numerical integration.</li>
!!                <li> Added autoconfiguration mode.</li>
!!                <li> Fixed various bugs.</li>
!!              </ul>
!!
!!              Version 1.0.0-alpha (02/2014)
!!              <ul>
!!                <li> First version of NEPTUNE with basic propagation functionalities.</li>
!!              </ul>
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!!------------------------------------------------------------------------------------------------
module version

  use slam_io, only: openFile, closeFile, LOG_AND_STDOUT, SEQUENTIAL, IN_FORMATTED, message
  use slam_error_handling, only: hasToReturn, isControlled, E_PARSE, FATAL, checkIn, checkOut
  use neptune_error_handling, only: setNeptuneError

  implicit none

  private

    integer, parameter :: LEN_VERSION_STRING = 15 !< max length of version string

    !** abbreviation NEPTUNE
    character(len=*), parameter :: neptuneAbbrev = 'NPI Ephemeris Propagation Tool with Uncertainty Extrapolation'


    type, public :: version_t
        integer :: major  = 0
        integer :: minor  = 0
        integer :: patch  = 0
        character(len=20) :: suffix = ''
    contains
        procedure :: version_gt
        generic :: operator(>) => version_gt
    end type

    type, public :: Version_class

        !** current version of NEPTUNE
        type(version_t)                     :: current_version
        integer                             :: logoLineCounter                   ! number of lines in NEPTUNE logo data file
        logical                             :: loadedLogo                       ! flag for having loaded NEPTUNE logo
        character(len=255), dimension(100)  :: neptuneLogo ! array containing NEPTUNE logo

    contains
        !** getter
        procedure :: getNeptuneAbbreviation
        procedure :: get_neptune_version
        procedure :: getVersionDate
        procedure :: printNeptuneLogo
        procedure :: parse_version_string
        procedure :: to_string

    end type Version_class

    ! Constructor
    interface Version_class
        module procedure constructor
    end interface Version_class

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
    type(Version_class) function constructor()
        constructor%current_version = version_t(1,1,2,'release')
        constructor%logoLineCounter = 0                                         ! number of lines in NEPTUNE logo data file
        constructor%loadedLogo      = .false.
    end function constructor

!===============================================================================================
!
!> @anchor      version_gt
!!
!> @brief       Compare two version for 'v1 > v2'
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 14.04.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    pure function version_gt(v1,v2) result(b)
        class(version_t), intent(in)     :: v1,v2

        logical :: b

        b = .false.
        if(v1%major > v2%major) then
            b = .true.
            return
        else if(v1%major == v2%major) then
            if(v1%minor > v2%minor) then
                b = .true.
            else if(v2%minor == v2%minor) then
                if(v1%patch > v2%patch) then
                    b = .true.
                end if
            end if
        end if
        return
    end function

!===============================================================================================
!
!> @anchor      getNeptuneAbbreviation
!!
!> @brief       Gets the meaning of the abbreviation 'NEPTUNE'
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 03.11.2014 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------

  character(len=len(neptuneAbbrev)) function getNeptuneAbbreviation(this)
    class(Version_class)    :: this
    getNeptuneAbbreviation = neptuneAbbrev
    return

  end function getNeptuneAbbreviation

!===============================================================================================
!
!> @anchor      get_neptune_version
!!
!> @brief       Gets the current version of NEPTUNE
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial implementation)</li>
!!                <li> 14.04.2017 (changed to derived type for version)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
    pure type(version_t) function get_neptune_version(this)
        class(Version_class),intent(in)    :: this
        get_neptune_version = this%current_version
        return
    end function

!===============================================================================================
!
!> @anchor      to_string
!!
!> @brief       Gets the current version of NEPTUNE as a string
!> @author      Vitali Braun
!!
!! @param[in]   version     Version to convert to string
!!
!> @date        <ul>
!!                <li> 14.04.2017 (initial implementation)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------

  character(len=LEN_VERSION_STRING) function to_string(this,version) result(v)

    class(Version_class)        :: this
    type(version_t), intent(in) :: version
    integer :: itemp
    character(len=4)  :: c_major, c_minor, c_patch
    character(len=30) :: format_string

    if(version%major == 0) then
        itemp = 0
    else
        itemp = int(log10(version%major*1.0)) + 1
    end if
    write(c_major,'(i4)') itemp

    if(version%minor == 0) then
        itemp = 0
    else
        itemp = int(log10(version%minor*1.0)) + 1
    end if
    write(c_minor,'(i4)') itemp

    if(version%patch == 0) then
        itemp = 0
    else
        itemp = int(log10(version%patch*1.0)) + 1
    end if
    write(c_patch,'(i4)') itemp

    format_string = '(i'//trim(adjustl(c_major))//',".",i'//trim(adjustl(c_minor))//',".",i'//trim(adjustl(c_patch))//')'
    write(v,trim(format_string)) version%major, version%minor, version%patch
    if(len(version%suffix) > 0) then
        v = trim(v)//'-'//trim(version%suffix)
    end if
    return

  end function to_string

!===============================================================================================
!
!> @anchor      printNeptuneLogo
!!
!> @brief       Print NEPTUNE logo on screen and to logfile
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 31.01.2015 (initial implementation)</li>
!!                <li> 01.02.2015 (added optional prefix for each logo line)</li>
!!              </ul>
!!
!! @param[in]   f_logo    Name of logo file
!! @param[in]   prefix    (optional) prefix for each line, e.g. '#' for a commented line (max. 2 characters are used)
!!
!!------------------------------------------------------------------------------------------------
  subroutine printNeptuneLogo(this, f_logo, prefix)

    class(Version_class)                    :: this
    character(len=*),           intent(in)  :: f_logo
    character(len=*), optional, intent(in)  :: prefix

    integer            :: i         ! counter
    integer            :: ich_logo  ! logo file channel
    integer            :: ios       ! I/O status
    character(len=255) :: cbuf      ! buffer
    character(len=2)   :: c_pre     ! filtered prefix to be used in logo printing

    if(.not. this%loadedLogo) then ! read logo from data file

      ich_logo = openFile(f_logo, SEQUENTIAL, IN_FORMATTED)

      do

        read(ich_logo,'(a)', iostat=ios) cbuf

        if(ios /= 0) exit ! end of file

        this%logoLineCounter = this%logoLineCounter + 1
        if(this%logoLineCounter > size(this%neptuneLogo)) then
          write(*,*) "Something wrong with NEPTUNE logo read... Check number of lines in data file."
          stop
        end if

        read(cbuf,'(a)') this%neptuneLogo(this%logoLineCounter)

      end do

      ich_logo   = closeFile(ich_logo)
      this%loadedLogo = .true.

    end if

    !** check for prefix
    if(present(prefix)) then
      c_pre = prefix(1:min(2,len(prefix)))
    else
      c_pre = ''
    end if

    call message(trim(c_pre), LOG_AND_STDOUT)

    do i = 1, this%logoLineCounter
      call message(trim(c_pre)//trim(this%neptuneLogo(i)), LOG_AND_STDOUT)
    end do
!    call message('     _/     _/  _/_/_/_/  _/_/_/   _/_/_/_/_/  _/   _/  _/      _/  _/_/_/_/', LOG_AND_STDOUT)
!    call message('    _/_/   _/  _/        _/    _/     _/      _/   _/  _/_/    _/  _/',        LOG_AND_STDOUT)
!    call message('   _/ _/  _/  _/_/_/    _/_/_/       _/      _/   _/  _/  _/  _/  _/_/_/',     LOG_AND_STDOUT)
!    call message('  _/   _/_/  _/        _/           _/      _/   _/  _/    _/_/  _/',          LOG_AND_STDOUT)
!    call message(' _/     _/  _/_/_/_/  _/           _/        _/_/   _/      _/  _/_/_/_/',     LOG_AND_STDOUT)
    call message(trim(c_pre), LOG_AND_STDOUT)
    return

  end subroutine printNeptuneLogo

!===============================================================================================
!
!> @anchor      getVersionDate
!!
!> @brief       Gets the date for a given version of NEPTUNE
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial implementation)</li>
!!                <li> 10.04.2014 (changed version date to scheme DD-MON-YYYY)</li>
!!              </ul>
!!
!! @param[in]   version     The version for which the date is requested
!!
!!------------------------------------------------------------------------------------------------


  character(len=11) function getVersionDate(this,version)

    class(Version_class)            :: this
    character(len=*), intent(in)    :: version

    select case(version)
      case('1.1.2-release')
        getVersionDate = '17-Sep-2021'
      case('1.1.1-release')
        getVersionDate = '15-Apr-2021'
      case('1.1.0-release')
        getVersionDate = '07-Mar-2021'
      case('1.0.0-release')
        getVersionDate = '26-Mar-2020'
      case('0.9.18-beta')
        getVersionDate = '30-Apr-2018'
      case('0.9.17-beta')
        getVersionDate = '27-Nov-2017'
      case('0.9.16-beta')
        getVersionDate = '28-Apr-2017'
      case('0.9.15-beta')
        getVersionDate = '12-Apr-2017'
      case('0.9.14-beta')
        getVersionDate = '14-JUN-2016'
      case('0.9.13-beta')
        getVersionDate = '26-JAN-2016'
      case('0.9.12-beta')
        getVersionDate = '24-JAN-2016'
      case('0.9.11-beta')
        getVersionDate = '24-JAN-2016'
      case('0.9.10-beta')
        getVersionDate = '20-JAN-2016'
      case('0.9.9-beta')
        getVersionDate = '11-JAN-2016'
      case('0.9.8-beta')
        getVersionDate = '28-AUG-2015'
      case('0.9.7-beta')
        getVersionDate = '23-AUG-2015'
      case('0.9.6-beta')
        getVersionDate = '08-JUN-2015'
      case('0.9.5-beta')
        getVersionDate = '17-MAY-2015'
      case('0.9.4-beta')
        getVersionDate = '04-MAY-2015'
      case('0.9.3-beta')
        getVersionDate = '01-MAY-2015'
      case('0.9.2-beta')
        getVersionDate = '20-APR-2015'
      case('0.9.1-beta')
        getVersionDate = '10-APR-2015'
      case('1.0.0-alpha')
        getVersionDate = '01-FEB-2014'
      case('1.0.0-beta')
        getVersionDate = '01-JUN-2014'
      case default
        getVersionDate = '00-XXX-0000'
    end select

    return

  end function getVersionDate

!===============================================================================================
!
!> @anchor      parse_version_string
!!
!> @brief       Parsing the version from a string
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.04.2017 (initial implementation)</li>
!!              </ul>
!!
!! @param[in]   version_string  A string in format MAJOR.MINOR.PATCH[-SUFFIX] (suffix may be, e.g. alpha or beta)
!!
!!------------------------------------------------------------------------------------------------
    type(version_t) function parse_version_string(this,version_string) result(v)

        class(Version_class)         :: this
        character(len=*), intent(in) :: version_string

        character(len=*), parameter :: routine_id = 'parse_version_string'

        integer :: dash, dot_one, dot_two   ! positions in string of first and second dot and the optional dash
        integer :: ios

        character(len=*), parameter :: format_error = 'Version parsing format error. '// &
                                                      'Required format MAJOR.MINOR.PATH[-SUFFIX]'

        if(isControlled()) then
            if(hasToReturn()) return
            call checkIn(routine_id)
        end if

        ! check if there is a suffix
        dash    = index(version_string,'-')
        dot_one = index(version_string,'.')
        if(dot_one == 0) then
            call setNeptuneError(E_PARSE, FATAL, (/format_error//'1'/))
            return
        end if
        dot_two = index(version_string(dot_one+1:),'.') + dot_one
        if(dot_two == 0) then
            call setNeptuneError(E_PARSE, FATAL, (/format_error//'2'/))
            return
        end if

        read(version_string(:dot_one-1),*,iostat=ios) v%major
        if(ios /= 0) then
            call setNeptuneError(E_PARSE, FATAL, (/format_error//'3'/))
            return
        end if
        read(version_string(dot_one+1:dot_two-1),*,iostat=ios) v%minor
        if(ios /= 0) then
            call setNeptuneError(E_PARSE, FATAL, (/format_error//'4'/))
            return
        end if
        if(dash > 0) then
            read(version_string(dot_two+1:dash-1),*, iostat=ios) v%patch
            v%suffix = version_string(dash+1:min(len(v%suffix),len(version_string)))
        else
            read(version_string(dot_two+1:),*, iostat=ios) v%patch
            v%suffix = ''
        end if

        if(isControlled()) then
            call checkOut(routine_id)
        end if
        return

    end function


end module version
