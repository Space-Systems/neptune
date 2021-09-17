!> @anchor      atmosphere
!!
!> @brief       Earth atmosphere modeling
!!
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!> @author      Daniel Lück (DLU)
!!
!> @date        <ul>
!!                <li>VB:  17.01.2013 (initial design)</li>
!!                <li>VB:  30.09.2013 (added HWM07)   </li>
!!                <li>VB:  03.11.2013 (added ESA solmag data file handling) </li>
!!                <li>VB:  23.08.2015 (added US Standard Atmosphere) </li>
!!                <li>CHK: 13.11.2015 (updated to use libslam) </li>
!!                <li>VB:  24.01.2016 (CSSI space weather cannot be used currently, bug reported)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 26.12.2017 (transformed to atmopshere class)
!!                <li>CHK: 05.01.2018 (Moved variables with save statements to Atmosphere class)
!!                <li>CHK: 16.02.2018 (Moving to object-oriented nrlmsise-00 implementation)
!!                <li>DLU: 12.01.2021 (Added JB2008 atmosphere model)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for Earth's
!!              atmosphere perturbations modeling, specifically in the context of numerical
!!              integration of a satellite trajectory. Also a horizontal wind model (HWM07) is
!!              implemented.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      atmosphere
!!
!!------------------------------------------------------------------------------------------------
module atmosphere

    use slam_types,             only: dp, sp
    use slam_astro,             only: EGM96, getEarthRotation
    use slam_astro_conversions, only: getGeodeticLatLon, getRadiusLatLon
    use neptune_error_handling, only: setNeptuneError, E_SGA_INPUT, E_SGA_UNSUPPORTED_TYPE, E_SGA_TIME_TAG, &
                                    E_SGA_UNSUPPORTED_VERSION, &
                                    E_SGA_NO_VERSION, E_SGA_MISSING, E_MIN_ALTITUDE, E_ATMOSPHERE_INIT, &
                                    getNeptuneErrorMessage
    use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, REMARK, &
                                    E_UNKNOWN_PARAMETER, E_SPECIAL, checkIn, checkOut
    !use gravity,                only: setCurrentAltitude
    use gravity,                only: Gravity_class
    use slam_io,                only: openFile, closeFile, SEQUENTIAL, IN_FORMATTED, cdelimit, LOG_AND_STDOUT, message
    use slam_math,              only: eps9, eps15, mag, rad2deg, deg2rad, outerproduct, identity_matrix, infinite, undefined, pi
    use slam_interpolation,     only: lagrange_interpolation
    use slam_reduction_class,   only: Reduction_type
    use satellite,              only: Satellite_class, ID_ATMOSPHERE
    use solarsystem,            only: Solarsystem_class
    use slam_statistics,        only: mean
    use slam_time,              only: getMonthNumber, time_t, mjd2yyddd, mjd2daySeconds, getLocalSolarTime, date2string, &
                                    gd2mjd, mjd2gd, getDateTimeNow
    use nrlmsise00Class,        only: Nrlmsise00_class, nrlmsise_flags
    use hwm07Class,             only: Hwm07_class
    use JB2008module
    implicit none

    private

    integer,public,parameter :: FILE_CSSI   = 1
    integer,public,parameter :: FILE_FAPDAY = 2
    integer,public,parameter :: FILE_JB08   = 3
    integer,public,parameter :: DAILY       = 1
    integer,public,parameter :: MONTHLY     = 2

    integer,public,parameter :: n_supported_sga_file_types = 3                  ! Number of different supported data file types for solar and geomagnetic activity

    integer,public,parameter :: JB08        = 3
    integer,public,parameter :: MSIS2000    = 2
    integer,public,parameter :: EXPONENTIAL = 1
    integer,public,parameter :: n_supported_models      = 3                     ! Number of supported atmosphere models
    integer,public,parameter :: n_supported_wind_models = 1                     ! Number of supported horizontal wind models

    ! Model names
    character(len=*), dimension(n_supported_models),      parameter :: modelName     = (/"Exponential","NRLMSISE-00","JB2008     "/)
    character(len=*), dimension(n_supported_wind_models), parameter :: windModelName = (/"HWM07"/)
    character(len=*),public,parameter,dimension(n_supported_sga_file_types) :: fileIdent    = (/"CssiS", "d/mm/", " F10,"/)       ! characteristic identifiers of sga data file
    character(len=*),public,parameter,dimension(n_supported_sga_file_types) :: suppSgaFiles = (/"CSSI (v1.2)", "ESA daily  ", "JB Indices "/)    ! supported solmag data files

    ! Solar and geomagnetic activity (SGA) type for Nrmmsise
    type sga_t

        real(dp) :: mjd                                                         ! Modified julian day of SGA
        integer, dimension(8) :: ap                                             ! 3hr planetary amplitude Ap values
        real(dp) :: f81ctr                                                      ! Observed F10.7, 81-day centered
        real(dp) :: f107                                                        ! Observed F10.7

        ! JB08 specific indices
        integer, dimension(24) :: dtc                                           ! 1hr temperature change from dtcfile
        real(dp) :: f10b_jb                                                     ! F10.7 index, 81-day centered
        real(dp) :: f10_jb                                                      ! F10.7 index
        real(dp) :: s10b_jb                                                     ! EUV Index, 81-day centered
        real(dp) :: s10_jb                                                      ! EUV index
        real(dp) :: m10b_jb                                                     ! MG2 index, 81-day centered
        real(dp) :: m10_jb                                                      ! MG2 index
        real(dp) :: y10b_jb                                                     ! x-ray and lya index, 81-day centered
        real(dp) :: y10_jb                                                      ! x-ray and lya index

    end type sga_t

    ! Solar and geomagnetic activity (SGA) data array for Nrmmsise
    type(sga_t), dimension(:), allocatable :: sga_data                          ! Containing solar and geomagnetic activity data (SGA)

    ! The atmospheric model class (type defintion)
    type, public :: Atmosphere_class

        !=======================================================================
        !
        ! Models
        !
        !-----------------------------------------------------------------------
        type(Nrlmsise00_class)     :: nrlmsise00_model                          ! NRLMSISE-00 model
        type(Hwm07_class)          :: hwm07_model                               ! HWM07 model

        real(dp) :: p_altitude_min                                              ! Lower boundary for altitude wrt. reference ellipsoid (km), 10 km being the default

        real(dp) :: cdom                                                        ! Ratio: cdrag over mass
        real(dp) :: cdrag                                                       ! Drag coefficient
        real(dp) :: dmass                                                       ! Mass of the satellite

        logical :: first                                                        ! If drag routine(s) called first time, some values have to be acquired once for the satellite configuration
        !--------------------------------------------------------------------------

        integer            :: nmodel                                            ! Atmosphere model to be used (default: NRLMSISE-00)
        integer            :: nwmodel                                           ! Wind model to be used

        !==============================================================================================================
        !
        ! Data files
        !
        !--------------------------------------------------------------------------------------------------------------
        character(len=255) :: dataPath                                          ! Path where input data files are located
        character(len=255) :: sgaDataFile                                       ! Solar and geomagnetic activity data file
        character(len=255) :: sgaMonthlyDataFile                                ! Monthly solar and geomagnetic activity data file
        character(len=255) :: hwmDataFile                                       ! Name of data file for Horizontal Wind Model (HWM07 as default)
        character(len=255) :: distWindDataFile                                  ! Name of data file for Disturbance Winds
        character(len=255) :: qdGridDataFile                                    ! Name of pre-computed interpolation grid for quasi-dipole coordinates

        character(len=3)  :: solmagFileVersion
        type(time_t)      :: solmagFileDateDaily
        type(time_t)      :: solmagFileDateMonthly

        !---------------------------------------------------------------------------------------------------------------

        ! other variables/parameters
        integer  :: ap_predict                                                  ! Geomagnetic activity value to use for long-term forecasts
        real(dp) :: sol_predict                                                 ! Solar activity value to use for long-term forecasts

        logical :: atmosphereInitialized                                        ! For initialization check
        logical :: monthlyDataRequired                                          ! Will be true, as soon as NEPTUNE starts reading monthly data
        logical :: considerWind                                                 ! Considering horizontal wind contributions if true
        logical :: covDrag                                                      ! Considering drag variational equations in covariance propagation

        !** saved variables for more efficiency
        real(dp), dimension(3) :: last_r
        real(dp)               :: last_time_mjd
        real(dp)               :: last_rho
        real(dp)               :: rho                                           ! atmospheric density in kg/km**3
        type(time_t)           :: date_start                                    ! First date solmag data is available for
        type(time_t)           :: date_end                                      ! Last date solmag data is available for

    contains
        procedure :: initAtmosphere
        procedure :: readCssiDataFile
        procedure :: read_jb2008_files
        procedure :: get_jb2008_date_interval
        procedure :: readFapDayDataFile
        procedure :: readFapMonDataFile
        procedure :: parseCSSIDataLine
        procedure :: kp2ap
        procedure :: destroy

        !** setter
        procedure :: setMinAltitude
        procedure :: setAtmosphereInitFlag
        procedure :: setAtmosphereModel
        procedure :: setApForecast
        procedure :: setSolarFluxForecast
        procedure :: setHorizontalWind
        procedure :: setHorizontalWindFileName
        procedure :: setDisturbanceWindFileName
        procedure :: setDragCovFlag
        procedure :: setQDGridFileName
        procedure :: setSolmagFileName

        !** getter
        procedure :: getAtmosphereModelName
        procedure :: getAtmosphereInitFlag
        procedure :: getAtmosphereAcceleration
        procedure :: getAtmosphericDensity
        procedure :: getApForecast
        procedure :: getDateForSolarActivity
        procedure :: getDragCovariance
        procedure :: getDragCovFlag
        procedure :: getMinAltitude
        procedure :: getSolarFluxForecast
        procedure :: getWindModelName
        procedure :: getHorizontalWindFileName
        procedure :: getDisturbanceWindFileName
        procedure :: getQDGridFileName
        procedure :: getRelativeVelocity
        procedure :: getSolarActivityForDate
        procedure :: getSolmagFileName
        procedure :: getSolmagFileVersion
        procedure :: getSgaDataFileType
        procedure :: getDensityExponential
        procedure :: getDensityMSIS2000
        procedure :: getDensityJB2008

    end type

    ! Constructor
    interface Atmosphere_class
        module procedure constructor
    end interface Atmosphere_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !>  @author     Christopher Kebschull
    !>  @date       <ul>
    !!                  <li>ChK: 24.12.2017 (initial implementation)</li>
    !!              </ul>
    !>  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Atmosphere_class) function constructor()

        constructor%nrlmsise00_model = Nrlmsise00_class()                       ! Instatiate NRLMSISE-00
        constructor%hwm07_model = Hwm07_class()

        constructor%p_altitude_min  = 10.d0                                     ! Lower boundary for altitude wrt. reference ellipsoid (km), 10 km being the default
        constructor%first           = .true.                                    ! If drag routine(s) called first time, some values have to be acquired once for the satellite configuration

        constructor%nmodel          = 2                                         ! Atmosphere model to be used (default: NRLMSISE-00)
        constructor%nwmodel         = 1                                         ! Wind model to be used

        constructor%dataPath           = "data"                                 ! Path where input data files are located
        constructor%sgaDataFile        = "fap_day.dat"                          ! Solar and geomagnetic activity data file
        constructor%sgaMonthlyDataFile = "fap_mon.dat"                          ! Monthly solar and geomagnetic activity data file
        constructor%hwmDataFile        = "hwm071308e.dat"                       ! Name of data file for Horizontal Wind Model (HWM07 as default)
        constructor%distWindDataFile   = "dwm07b_104i.dat"                      ! Name of data file for Disturbance Winds
        constructor%qdGridDataFile     = "apexgrid.dat"                         ! Name of pre-computed interpolation grid for quasi-dipole coordinates

        constructor%solmagFileVersion   = '---'

        ! other variables/parameters
        constructor%ap_predict = 10                                             ! Geomagnetic activity value to use for long-term forecasts
        constructor%sol_predict = 80.d0                                         ! Solar activity value to use for long-term forecasts

        constructor%atmosphereInitialized = .false.                             ! For initialization check
        constructor%monthlyDataRequired   = .false.                             ! Will be true, as soon as NEPTUNE starts reading monthly data
        constructor%considerWind          = .false.                             ! Considering horizontal wind contributions if true
        constructor%covDrag               = .false.                             ! Considering drag variational equations in covariance propagation
    end function constructor

    !===========================================================================
    !!
    !>  @anchor     destroy
    !!
    !>  @brief      Destroys all that needs destruction
    !>  @author     Christopher Kebschull
    !!
    !>  @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Atmosphere_class)    :: this

        call this%hwm07_model%destroy()
        if (allocated(sga_data)) deallocate(sga_data)
    end subroutine destroy

!==========================================================================
!
!> @anchor      setAtmosphereInitFlag
!!
!> @brief       Set initialization flag to .false.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 05.06.2013 (initial design) </li>
!!              </ul>
!!
!!-------------------------------------------------------------------------
  subroutine setAtmosphereInitFlag(this)
    class(Atmosphere_class) :: this
    this%atmosphereInitialized = .false.
    return
  end subroutine setAtmosphereInitFlag

!==========================================================================
!
!> @anchor      getAtmosphereInitFlag
!!
!> @brief       Get initialization flag
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.03.2014 (initial design) </li>
!!              </ul>
!!
!!-------------------------------------------------------------------------
  logical function getAtmosphereInitFlag(this)
    class(Atmosphere_class) :: this
    getAtmosphereInitFlag = this%atmosphereInitialized

  end function getAtmosphereInitFlag

!===========================================================================
!
!> @anchor      setHorizontalWind
!!
!> @brief       Set wind flag to .true.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 25.10.2013 (initial design) </li>
!!              </ul>
!!
!!---------------------------------------------------------------------------
  subroutine setHorizontalWind(this,val)
    class(Atmosphere_class) :: this
    logical, intent(in) :: val
    this%considerWind = val

  end subroutine setHorizontalWind


!=============================================================================
!
!> @anchor      initAtmosphere
!!
!> @brief       Initialization of atmosphere module
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li> 17.01.2013 (VB: initial design)</li>
!!                <li> 03.06.2013 (VB: code optimization)</li>
!!                <li> 23.08.2015 (VB: added exponential atmosphere handling)</li>
!!                <li> 28.02.2018 (CHK: Always reading the full file)</li>
!!              </ul>
!!
!> @param[in]   cpath       Path to read space solar and geomagnetic data from
!!
!> @details     This routine initializes the atmosphere module. Global parameters are set and
!!              solar and geomagnetic activity data is read from a current data file into the specific arrays.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Allocate solmag data array          </li>
!!                <li> Check which data file shall be read </li>
!!                <li> Call appropriate routine
!!                <li> Read data </li>
!!                <li> Perform Lagrange interpolation to get daily data,
!!                     if only monthly data is available </li>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------

  subroutine initAtmosphere(  this, &
                              cpath)

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    character(len=*), intent(in)        :: cpath
    !------------------------------------------

    character(len=*), parameter :: csubid = "initAtmosphere"
    integer                     :: ich                                          ! input channel
    integer                     :: idays                                        ! number of days for which data is required
    integer                     :: sgaDataType                                  ! data file type for solmag input read (e.g. 1=CSSI, 2=FAP_DAY)
    real(dp)                    :: mjd_end                                      ! MJD of last SGA data entry
    real(dp)                    :: mjd_start                                    ! MJD of first SGA data entry
    character(len=255)          :: cmess                                        ! message string


    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%atmosphereInitialized .or. this%nmodel == EXPONENTIAL) then         ! return if already initialized or exponential atmosphere is used
      this%atmosphereInitialized = .true.                                            ! (for the latter, no init is required)
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
    end if

    !** set data path in module scope
    this%dataPath = trim(adjustl(cpath(1:min(len(this%dataPath),len(cpath)))))

    !** open solar and geomagnetic activity data file
    !---------------------------------------------------------------------------------------------------
    ich = openFile(trim(adjustl(this%dataPath))//cdelimit//trim(adjustl(this%sgaDataFile)), SEQUENTIAL, IN_FORMATTED)

    call message(' - Reading solar and geomagnetic activity data...', LOG_AND_STDOUT)

    !** find out type of data file
    !-----------------------------------------------------------------
    sgaDataType = this%getSgaDataFileType(ich)
    if(hasFailed()) return

    if(isControlled() .and. hasToReturn()) then
      call checkOut(csubid)
      return
    end if

    this%date_start%mjd = 37300.d0  ! Jan 1st, 1961
    this%date_end%mjd = 71183.d0    ! Dec 31st, 2051

    !** add three 'previous days' (for NRLMSISE-00) and an additional day at the end (for interpolation)
    mjd_start = this%date_start%mjd - 3.d0
    mjd_end   = this%date_end%mjd + 1.d0

    if(sgaDataType == FILE_FAPDAY .and. this%nmodel == MSIS2000) then           ! F10.7 81-day centered will be required, but ESA daily data provides last 81-day data only,
                                                                                ! so 41 additional days are required for a shift later on
      mjd_start = mjd_start - 81.d0                                             ! in order to compute F10.7 last 81-day averages
      mjd_end   = mjd_end   + 41.d0

    else if (sgaDataType == FILE_JB08 .and. this%nmodel == MSIS2000) then

      cmess = "Wrong input file (JB08) specified for the selected NRLMSISE-00 atmospheric model."
      call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
      return

    end if

    if(sgaDataType == FILE_JB08 .and. this%nmodel == JB08) then

      !mjd_start = 50449.5   ! Jan 1st, 1997
      call this%get_jb2008_date_interval(ich, this%date_start%mjd, this%date_end%mjd)
      mjd_start = this%date_start%mjd
      mjd_end   = this%date_end%mjd

    else if (sgaDataType == FILE_FAPDAY .and. this%nmodel == JB08) then

      cmess = "Wrong input file (FAPDAY) specified for the selected JB2008 atmospheric model."
      call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
      return

    else if (sgaDataType == FILE_CSSI .and. this%nmodel == JB08) then

      cmess = "Wrong input file (CSSI) specified for the selected JB2008 atmospheric model."
      call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
      return

    end if

    !** allocate SGA array according to number of days... (+60 to account mainly for monthly predictions, requiring more input to read for interpolation
    !---------------------------------------------------
    idays = int(mjd_end) - int(mjd_start) + 60

    if(allocated(sga_data)) then
      deallocate(sga_data)
    end if
    allocate(sga_data(idays))

    !---------------------------------------------------

    !** now start reading data into the allocated sga_data array depending on type of file
    !----------------------------------------------------------------------------------------
    select case(sgaDataType)
      case(FILE_CSSI)
        call this%readCssiDataFile(ich, mjd_start, mjd_end)
        if(hasFailed()) return
      case(FILE_FAPDAY)
        call this%readFapDayDataFile(ich, mjd_start, mjd_end)
        if(hasFailed()) return
      case(FILE_JB08)
        call this%read_jb2008_files(ich, mjd_start, mjd_end)
        if(hasFailed()) return
    end select

    !** IF wind is considered, open horizontal wind data files
    !----------------------------------------------------------------------------------
    if(this%considerWind) then

      !** Horizontal Wind Model data file
      call message(' - Reading horizontal wind model data...', LOG_AND_STDOUT)
      call this%hwm07_model%loadmodel(trim(adjustl(cpath))//cdelimit//trim(adjustl(this%hwmDataFile)))
      if(hasFailed()) return

      !** Pre-computed interpolation grid for quasi-dipole coordinates
      call message(' - Reading QD interpolation grid...', LOG_AND_STDOUT)
      call this%hwm07_model%apxrda(trim(adjustl(cpath))//cdelimit//trim(adjustl(this%qdGridDataFile)))
      if(hasFailed()) return

      !** Disturbance winds
      call message(' - Reading disturbance wind model data...', LOG_AND_STDOUT)
      call this%hwm07_model%loaddwm(trim(adjustl(cpath))//cdelimit//trim(adjustl(this%distWindDataFile)))
      if(hasFailed()) return

    end if

    this%atmosphereInitialized = .true.

    ich = closeFile(ich)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine initAtmosphere

!=============================================================================
!
!> @anchor      getAtmosphereAcceleration
!!
!> @brief       Computes acceleration due to drag
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 19.01.2013 (initial design)</li>
!!                <li> 04.02.2014 (re-structured in order to harmonize with covariance propagation routine)</li>
!!              </ul>
!!
!> @param[in]   gravity_model     the gravity model
!> @param[in]   satellite_model   the satellite model
!> @param[in]   solarsystem_model the solar system model
!> @param[in]   r_gcrf            radius vector in GCRF / km
!> @param[in]   v_gcrf            velocity vector in GCRF / km/s
!> @param[in]   r_itrf            radius vector in ITRF / km
!> @param[in]   time_mjd          current time MJD
!!
!> @param[out]  acc_atmosphere    acceleration vector in inertial frame / km/s**2
!!
!> @details     This routine computes the acceleration in the GCRF due to atmospheric drag.
!!              It uses the NRLMSISE-00 model to estimate the atmospheric density.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Finish. </li>
!!              </ol>
!!
!> @todo        Ap(1) in NRLMSISE-00 - rounded or not?!
!!
!!------------------------------------------------------------------------------------------------
  subroutine getAtmosphereAcceleration(         &
                             this,              &
                             gravity_model,     &  ! <-> TYPE gravity model
                             satellite_model,   &  ! <-> TYPE satellite model
                             solarsystem_model, &  ! <-> TYPE solarsystem model
                             reduction,         &  ! <-> TYPE reduction
                             r_gcrf,            &  ! <-- DBL(3) radius vector in GCRF
                             v_gcrf,            &  ! <-- DBL(3) velocity vector in GCRF
                             r_itrf,            &  ! <-- DBL(3) radius vector in ITRF
                             v_itrf,            & ! <-- DBL(3) velocity vector in ITRF
                             time_mjd,          &  ! <-- DBL    current time MJD
                             acc_atmosphere     &  ! --> DBL(3) acceleration vector in inertial frame
                                      )


    implicit none

    !** interface
    !---------------------------------------------------
    class(Atmosphere_class)                 :: this
    type(Gravity_class),intent(inout)       :: gravity_model
    type(Satellite_class),intent(inout)     :: satellite_model
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)     :: reduction
    real(dp), dimension(3), intent(in)      :: r_gcrf
    real(dp), dimension(3), intent(in)      :: v_gcrf
    real(dp), dimension(3), intent(in)      :: r_itrf
    real(dp), dimension(3), intent(in)      :: v_itrf
    real(dp),               intent(in)      :: time_mjd

    real(dp), dimension(3), intent(out) :: acc_atmosphere
    !---------------------------------------------------

    character(len=*), parameter :: csubid = "getAtmosphereAcceleration"

    real(dp)               :: crossSection    ! cross-section
    real(dp)               :: rho             ! atmospheric density from MSIS model
    real(dp), dimension(3) :: v_rel           ! relative velocity (due to atmosphere rotation and wind)


    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether initialization has been done already...
    if(.not. this%atmosphereInitialized) then
      call setNeptuneError(E_ATMOSPHERE_INIT, FATAL)
      return
    end if

    !** get atmospheric density
    rho = this%getAtmosphericDensity(gravity_model, r_gcrf, v_gcrf, r_itrf, v_itrf, time_mjd)

    if(hasFailed()) return

    !** get drag coefficient and mass
    if(satellite_model%hasChangedSatellite(ID_ATMOSPHERE) .or. this%first)  then                ! only required, if configuration has changed or called for first time

      this%cdrag = satellite_model%getObjectDragCoefficient()
      if(hasFailed()) return
      this%dmass = satellite_model%getObjectMass()
      if(hasFailed()) return
      this%cdom  = this%cdrag/this%dmass
      this%first = .false.

    end if

    !** get relative velocity (considering wind IF selected)
    v_rel = this%getRelativeVelocity(reduction, r_gcrf, v_gcrf, time_mjd)
    if(hasFailed()) return

    !** get cross-section
    crossSection = satellite_model%getObjectCrossSection(solarsystem_model, reduction, time_mjd, r_gcrf, v_rel)
    if(hasFailed()) return

    !write(*,*) "crossSection = ", crossSection
    !write(*,*) "cdom = ", cdom
    !write(*,*) "rho = ", rho
    !write(*,*) "vrel = ", v_rel
    !write(*,*) "vabs = ", mag(v_rel)

    !write(67,'(f12.6,x,3(e12.5e2,x),f8.4)') time_mjd, crossSection, cdom, rho, mag(v_rel)

    !** compute acceleration
    acc_atmosphere(1:3) = -0.5d0*crossSection*this%cdom*rho*v_rel(1:3)*mag(v_rel)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine getAtmosphereAcceleration

!!------------------------------------------------------------------------------------------------
!
!> @anchor      getAtmosphericDensity
!!
!> @brief       Gets the density from the atmosphere model using the radius vector (ITRF) and the MJD.
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial design)</li>
!!                <li> 23.08.2015 (implemented exponential atmosphere)</li>
!!              </ul>
!!
!> @param[in]   gravity_model     the gravity model
!> @param[in]   r_itrf            radius vector of satellite in GCRF / km
!> @param[in]   time_mjd          MJD
!!
!> @return      Rho - the atmospheric density / kg/km**3
!!
!> @details     This routine returns the density (in kg/km**3) from the selected atmosphere model
!!              for the passed radius vector and the MJD.
!!
!> @todo        JB2008 data only go up to 03.11.2020 after that NRLMSISE-00 is used (this is hardcoded)
!!
!!------------------------------------------------------------------------------------------------
  real(dp) function getAtmosphericDensity(this, gravity_model, r_gcrf, v_gcrf, r_itrf, v_itrf, time_mjd) result(rho)
    use slam_strings, only: toString

    implicit none

    class(Atmosphere_class)            :: this
    type(Gravity_class),intent(inout)  :: gravity_model
    real(dp), dimension(3), intent(in) :: r_gcrf
    real(dp), dimension(3), intent(in) :: v_gcrf
    real(dp), dimension(3), intent(in) :: r_itrf
    real(dp), dimension(3), intent(in) :: v_itrf
    real(dp),               intent(in) :: time_mjd

    character(len=*), parameter :: csubid =  "getAtmosphericDensity"

    integer  :: ierr                                                            ! error flag
    real(dp) :: altitude                                                        ! geodetic altitude (km)
    real(dp) :: height                                                          ! geocentric altitude (km)

    real(dp) :: lat_gd                                                          ! geodetic latitude (rad)
    real(dp) :: lon_gd                                                          ! geodetic longitude (rad)

    real(dp) :: lat_gc                                                          ! geocentric latitude (rad)
    real(dp) :: lon_gc                                                          ! geocentric longitude (rad)
    real(dp) :: right_ascension                                                 ! right ascension of satellite position (not RAAN!)
    real(dp) :: rabs                                                            ! magnitude of radius (km)

    character(len=255) :: cmess

    real(dp), dimension(3)       :: tmp_r
    this%last_r = -1.d0
    this%last_time_mjd = -1.
    this%last_rho = -1.d0
    tmp_r = -1.d0

    rho = -1.d0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    tmp_r = r_itrf - this%last_r
    if(mag(tmp_r) < eps9 .and. time_mjd == this%last_time_mjd) then   ! do not compute anew, as already available

      rho = this%last_rho
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return

    end if

    !** get current geocentric latitude and longitude
    call getGeodeticLatLon(r_itrf, altitude, lat_gd, lon_gd)
    call getRadiusLatLon(r_itrf, v_itrf, rabs, lat_gc, lon_gc)

    ! Calculate Right Ascension of Position (Vallado Algorithm 25):
    if(r_gcrf(2) > 0.0 ) then
      right_ascension = acos(r_gcrf(1)/sqrt(r_gcrf(1)**2+r_gcrf(2)**2))
    else
      right_ascension = 2*pi - acos(r_gcrf(1)/sqrt(r_gcrf(1)**2+r_gcrf(2)**2))
    end if

    ! height = rabs - 6371.0 !height not clearly defined in jb2008 this is for geocentric height



    !** check whether altitude is below re-entry altitude limit
    if(altitude < this%p_altitude_min) then
      call setNeptuneError(E_MIN_ALTITUDE, FATAL)
      return
    end if

    ierr = gravity_model%setCurrentAltitude(altitude)

    if(this%nmodel == EXPONENTIAL) then
      rho = this%getDensityExponential(altitude)
    else if(this%nmodel == MSIS2000) then
      rho = this%getDensityMSIS2000(altitude, lat_gd, lon_gd, time_mjd)
    else if(this%nmodel == JB08) then
      if(this%date_start%mjd <= time_mjd .and. time_mjd <= this%date_end%mjd) then                                                       ! JB2008 data currently only available until 03.11.2020
        rho = this%getDensityJB2008(altitude, lat_gc, right_ascension, time_mjd)
      else
        cmess = "No solar and geomagnetic data available for the given time: "//toString(time_mjd)
        call setNeptuneError(E_SPECIAL, FATAL, (/cmess/))
        return
      end if
    end if

    this%last_rho      = rho
    this%last_time_mjd = time_mjd
    this%last_r        = r_itrf

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getAtmosphericDensity

  !!------------------------------------------------------------------------------------------------
  !
  !> @anchor      getDensityMSIS2000
  !!
  !> @brief       Gets the density from the NRLMSISE-00 atmosphere model
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 23.08.2015 (initial design, after re-structuring to introduce a new atmosphere model)</li>
  !!                <li> 24.01.2016 (fix: latitude and longitude are essential input for this routine)</li>
  !!              </ul>
  !!
  !> @param[in]   altitude          geodetic altitude / km
  !> @param[in]   latitude          geodetic latitude / rad
  !> @param[in]   longitude         longitude / rad
  !> @param[in]   time_mjd          MJD
  !!
  !> @return      Rho - the atmospheric density / kg/km**3
  !!
  !> @details     This routine returns the density (in kg/km**3) from the NRLMSISE-00 model
  !!              for the passed altitude and the MJD.
  !!
  !!------------------------------------------------------------------------------------------------
  real(dp) function getDensityMSIS2000(this, altitude, latitude, longitude, time_mjd) result(rho)
    use nrlmsise00Class,        only: nrlmsise_flags,nrlmsise_input,nrlmsise_output

    implicit none

    class(Atmosphere_class) :: this
    real(dp), intent(in)    :: altitude
    real(dp), intent(in)    :: time_mjd
    real(dp), intent(in)    :: latitude
    real(dp), intent(in)    :: longitude

    real(dp)                :: lat_gd
    real(dp)                :: lon
    type(nrlmsise_input)    :: input
    type(nrlmsise_flags)    :: flags
    type(nrlmsise_output)   :: output

    integer                 :: ap_pointer                                       ! Ap pointer
    integer                 :: hour                                             ! hour
    integer                 :: i                                                ! loop counter
    integer                 :: idate                                            ! date in YYDDD format for MSIS atmosphere model
    integer                 :: year                                             ! date in YY format for MSIS atmosphere model
    integer                 :: doy                                              ! date in DDD format for MSIS atmosphere model
    integer                 :: idx
    integer                 :: mass                                             ! mass number
    integer                 :: offset                                           ! index offset in space weather data array

    real(dp), dimension(7)  :: ap                                               ! planetary amplitude Ap
    real(dp)                :: loc_solar_time                                   ! local solar time (hrs.)
    real(dp)                :: sec                                              ! seconds of day

    !** convert to time format YYDDD with DDD being the day of the year
    idate = mjd2yyddd(time_mjd)
    ! Get only day of year, as years are ignored anyhow
    year = int(idate / 1d3)
    doy = idate - (year * 1d3)

    lat_gd = latitude*rad2deg
    lon    = longitude*rad2deg

    !** get seconds of day
    sec  = mjd2daySeconds(time_mjd)
    hour = int(sec/3600.d0)

    !** get local solar time (in hrs.)
    loc_solar_time = getLocalSolarTime(time_mjd, lon)

    !** find index in space weather data for current date
    idx = int(time_mjd - int(sga_data(1)%mjd)) + 1

    !** build ap array
    !------------------------------------------
    ap_pointer = int(hour/3.d0) + 1
    offset     = 0

    ! (1) DAILY AP
    ap(1) = 0.d0

    do i = 1,8
      ap(1) = ap(1) + sga_data(idx)%ap(i)
    end do

    ap(1) = ap(1)/8.d0

    ! (2) 3 HR AP INDEX FOR CURRENT TIME
    ap(2) = sga_data(idx)%ap(ap_pointer)

    do i=3,5  ! (3,4,5) 3 HR AP INDEX FOR 3 HR, 6 HR and 9 HR BEFORE CURRENT TIME
      ap_pointer = ap_pointer - 1
      if(ap_pointer < 1) then
        ap_pointer = 8
        offset     = offset + 1
      end if
      ap(i) = sga_data(idx - offset)%ap(ap_pointer)
    end do

    ! (6) AVERAGE OF EIGHT 3 HR AP INDICES FROM 12 TO 33 HRS PRIOR
    !     TO CURRENT TIME
    ap(6) = 0.d0

    do i=1,8
      ap_pointer = ap_pointer - 1
      if(ap_pointer < 1) then
        ap_pointer = 8
        offset     = offset + 1
      end if
      ap(6) = ap(6) + sga_data(idx-offset)%ap(ap_pointer)
    end do

    ap(6) = ap(6)/8.d0

    ! (7) AVERAGE OF EIGHT 3 HR AP INDICES FROM 36 TO 57 HRS PRIOR
    !     TO CURRENT TIME
    ap(7) = 0.d0

    do i=1,8
      ap_pointer = ap_pointer - 1
      if(ap_pointer < 1) then
        ap_pointer = 8
        offset     = offset + 1
      end if
      ap(7) = ap(7) + sga_data(idx-offset)%ap(ap_pointer)
    end do

    ap(7) = ap(7)/8.d0

    !** set switch sw(9) = 1 in order to be compliant with the MSIS
    !   description
    !----------------------
    !sw(:)   = 1.d0
    !sw(9)  = -1.d0
    !sw(15) = 1.d0
    !call tselec(sw)
    flags%switches(0:23) = 1.d0
    flags%switches(0) =-1.d0
    flags%switches(9) =-1.d0
    call this%nrlmsise00_model%tselec(flags)
    !---------------------

    mass = 48 ! all components...

    if(altitude <= 500.d0) then
       !call gtd7 (idate,real(sec),real(altitude),real(lat_gd),real(lon),real(loc_solar_time), &
       !           real(sga_data(idx)%f81ctr),real(sga_data(idx-1)%f107),real(ap), &
       !           mass,sp_density,temp)

        input%year      = 0                             ! year, currently ignored
        input%doy       = doy                           ! day of year
        input%sec       = sec                           ! seconds in day (UT)
        input%alt       = altitude                      ! altitude in kilometers
        input%g_lat     = lat_gd                        ! geodetic latitude
        input%g_long    = lon                           ! geodetic longitude
        input%lst       = loc_solar_time                ! local apparent solar time (hours), see note below
        input%f107A     = sga_data(idx)%f81ctr     ! 81 day average of F10.7 flux (centered on doy)
        input%f107      = sga_data(idx-1)%f107     ! daily F10.7 flux for previous day
        input%ap        = ap(1)                         ! magnetic index(daily)
        input%ap_a%a(0:6) = ap(1:7)

        call this%nrlmsise00_model%gtd7(input,flags,output)

    else
      ! For atmospheric drag calculations at altitudes above 500 km,
      ! call SUBROUTINE GTD7D to compute the "effective total mass
      ! density" by including contributions from "anomalous oxygen."
      !call gtd7d (idate,real(sec),real(altitude),real(lat_gd),real(lon),real(loc_solar_time), &
      !             real(sga_data(idx)%f81ctr),real(sga_data(idx-1)%f107),real(ap), &
      !             mass,sp_density,temp)

        input%year      = 0                                     ! year, currently ignored
        input%doy       = doy                                   ! day of year
        input%sec       = real(sec)                             ! seconds in day (UT)
        input%alt       = real(altitude)                        ! altitude in kilometers
        input%g_lat     = real(lat_gd)                          ! geodetic latitude
        input%g_long    = real(lon)                             ! geodetic longitude
        input%lst       = real(loc_solar_time)                  ! local apparent solar time (hours), see note below
        input%f107A     = real(sga_data(idx)%f81ctr)       ! 81 day average of F10.7 flux (centered on doy)
        input%f107      = real(sga_data(idx-1)%f107)       ! daily F10.7 flux for previous day
        input%ap        = real(ap(1))                           ! magnetic index(daily)
        input%ap_a%a(0:6) = real(ap(1:7))

        call this%nrlmsise00_model%gtd7d(input,flags,output)

    end if
    !** compute density in kg/km**3
    !rho = dble(sp_density(6))*1.d12
    rho = dble(output%d(5))*1.d12
    !write(*,*) time_mjd, longitude, latitude, altitude


    return

  end function getDensityMSIS2000

  !!------------------------------------------------------------------------------------------------
  !
  !> @anchor      getDensityExponential
  !!
  !> @brief       Gets the density from the exponential atmosphere model
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 23.08.2015 (initial design)</li>
  !!              </ul>
  !!
  !> @param[in]   altitude            geodetic altitude / km
  !!
  !> @return      Rho - the atmospheric density / kg/km**3
  !!
  !> @details     This routine returns the density (in kg/km**3) from the exponential model, which
  !!              is given by Wertz (1978). It is a combination of the USSA (1976) for 0 km,
  !!              CIRA-72 between 25 km and 500 km and CIRA-72 with exospheric temperature, T_inf = 1000 K
  !!              for 500 - 1000 km. Above 1000 km, the density is considered zero.
  !!
  !!------------------------------------------------------------------------------------------------
  real(dp) function getDensityExponential(this,altitude) result(rho)

    class(Atmosphere_class) :: this
    real(dp), intent(in)    :: altitude

    real(dp), dimension(27), parameter :: expAltitudes = (/ 25.d0, 30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 90.d0,100.d0, &
                     110.d0,120.d0,130.d0,140.d0,150.d0,180.d0,200.d0,250.d0,300.d0, &
                     350.d0,400.d0,450.d0,500.d0,600.d0,700.d0,800.d0,900.d0,1000.d0/)
    real(dp), dimension(27), parameter :: exph0 = (/ 0.d0, 25.d0, 30.d0, 40.d0, 50.d0, 60.d0, 70.d0, 80.d0, 90.d0,100.d0, &
                     110.d0,120.d0,130.d0,140.d0,150.d0,180.d0,200.d0,250.d0,300.d0, &
                     350.d0,400.d0,450.d0,500.d0,600.d0,700.d0,800.d0,900.d0/)
    real(dp), dimension(27), parameter :: expDensity = (/1.225d0,3.899d-2,1.774d-2,3.972d-3,1.057d-3,3.206d-4,8.770d-5, &
                     1.905d-5, 3.396d-6, 5.297d-7, 9.661d-8, 2.438d-8, 8.484d-9, 3.845d-9, &
                     2.070d-9, 5.464d-10,2.789d-10,7.248d-11,2.418d-11,9.518d-12,3.725d-12,&
                     1.585d-12,6.967d-13,1.454d-13,3.614d-14,1.170d-14,5.245d-15/)
    real(dp), dimension(27), parameter :: expHeight = (/ 7.249d0, 6.349d0, 6.682d0, 7.554d0, 8.382d0, 7.714d0, 6.549d0, 5.799d0, 5.382d0, &
                                                         5.877d0, 7.263d0, 9.473d0,12.636d0,16.149d0,22.523d0,29.740d0,37.105d0,45.546d0, &
                                                        53.628d0,53.298d0,58.515d0,60.828d0,63.822d0,71.835d0,88.667d0,124.64d0,181.05d0/)
    real(dp) :: maxAltitude = 1000.d0
    integer :: i,idx

    if(altitude > maxAltitude) then
      rho = 0.d0
      return
    end if

    ! find index in altitudes array
    do i=1,size(expAltitudes)
      if(altitude < expAltitudes(i)) then
        idx = i
        exit
      end if
    end do

    rho = expDensity(idx)*exp(-(altitude - exph0(idx))/expHeight(idx))

    !** convert to kg/km**3
    rho = rho*1.d9

    return

  end function getDensityExponential

  !!------------------------------------------------------------------------------------------------
  !
  !> @anchor      getDensityJB2008
  !!
  !> @brief       Gets the density from the JB2008 atmosphere model
  !> @author      Daniel Lück
  !!
  !> @date        <ul>
  !!                <li> 18.12.2020 (initial design, to implement new atmosphere model)</li>
  !!              </ul>
  !!
  !> @param[in]   altitude          geodetic altitude / km    (althought this might also be geocentric)
  !> @param[in]   latitude          geocenric latitude / rad
  !> @param[in]   right ascension   right ascension / rad
  !> @param[in]   time_mjd          MJD
  !!
  !> @return      Rho - the atmospheric density / kg/km**3
  !!
  !> @details     This routine returns the density (in kg/km**3) from the JB2008 model
  !!              for the passed altitude and the MJD.
  !!
  !!------------------------------------------------------------------------------------------------
  real(dp) function getDensityJB2008(this, altitude, latitude, right_ascension, time_mjd) result(rho)
    use nrlmsise00Class,        only: nrlmsise_flags,nrlmsise_input,nrlmsise_output
    use JB2008module


    implicit none

    class(Atmosphere_class) :: this
    real(dp), intent(in)    :: altitude
    real(dp), intent(in)    :: time_mjd
    real(dp), intent(in)    :: latitude
    real(dp), intent(in)    :: right_ascension

    integer                 :: idx

    real(dp), dimension(2)  :: sun              ! position of the sun
    real(dp), dimension(3)  :: sat              ! position of the sun
    real(dp)                :: f10              ! f10 solar value
    real(dp)                :: f10b             ! f10 solar value averaged over 81 days
    real(dp)                :: s10              ! see f10/f10b
    real(dp)                :: s10b             ! see f10/f10b
    real(dp)                :: m10              ! see f10/f10b
    real(dp)                :: m10b             ! see f10/f10b
    real(dp)                :: y10              ! see f10/f10b
    real(dp)                :: y10b             ! see f10/f10b
    real(dp)                :: dstdtc           ! temperature change through geomagnetic activity
    real(dp), dimension(2)  :: temp             ! jb2008 output temperature

    ! Parameters to calculate sun position
    real(dp)                :: solras           ! right ascension of the sun
    real(dp)                :: soldec           ! declination of the sun
    real(dp)                :: solan
    real(dp)                :: sin1l
    real(dp)                :: sin2l
    real(dp)                :: sin4l
    real(dp)                :: eps
    real(dp)                :: eclon
    real(dp)                :: tanhalfeps1
    real(dp)                :: tanhalfeps2
    real(dp)                :: tanhalfeps4
    real(dp)                :: solon
    real(dp)                :: d2000
    real(dp)                :: hourofday
    real(dp)                :: fractionofhour

    ! calculate position of the sun (taken from JB08DRVY2k.for)

    d2000 = time_mjd - 51544.5D0
    solan = 357.528D0 + 0.9856003 * d2000
    solan = solan*deg2rad
    solon = 280.460D0 + 0.9856474 * d2000
    solon = dmod(solon,360.D0)

    if (solon .LT. 0.D0) then
      solon = solon + 360.D0
    endif

    eclon = solon + 1.915D0 * dsin(solan) + 0.02D0 * dsin(2.D0*solan)
    eclon = eclon * deg2rad


    eps   = 23.439D0 - 0.0000004 * D2000
    eps   = eps * deg2rad


    sin1l = dsin(1.D0 * eclon)
    sin2l = dsin(2.D0 * eclon)
    sin4l = dsin(4.D0 * eclon)


    tanhalfeps1 = dtan(0.5 * eps)
    tanhalfeps2 = tanhalfeps1 * tanhalfeps1
    tanhalfeps4 = tanhalfeps2 * tanhalfeps2



    solras = eclon - tanhalfeps2 * sin2l + 0.5 * tanhalfeps4 * sin4l

    if (solras .LT. 0.D0) then
      solras = solras + 2*pi
    elseif (solras .GT. 2*pi) then
      solras = solras - 2*pi
    end if

    soldec = dasin(dsin(eps) * sin1l)

    sun(1) = solras
    sun(2) = soldec


    !** find index in space weather data for current date
    idx = int(time_mjd - int(sga_data(1)%mjd)) + 1

    ! set solar parameters
    f10 = sga_data(idx-1)%f10_jb
    f10b = sga_data(idx-1)%f10b_jb
    s10 = sga_data(idx-1)%s10_jb
    s10b = sga_data(idx-1)%s10b_jb
    m10 = sga_data(idx-2)%m10_jb
    m10b = sga_data(idx-2)%m10b_jb
    y10 = sga_data(idx-5)%y10_jb
    y10b = sga_data(idx-5)%y10b_jb

    if (f10.LT.40..OR.f10b.LT.40.) then
      f10  = 0.
      f10b = 0.
    endif
    if (s10.LT.40..OR.s10b.LT.40.) then
      s10  = 0.
      s10b = 0.
    endif
    if (m10.LT.40..OR.s10b.LT.40.) then
      m10  = 0.
      m10b = 0.
    endif
    if (y10.LT.40..OR.y10b.LT.40.) then
      y10  = 0.
      y10b = 0.
    endif

    ! set satellite position
    sat(1) = right_ascension
    sat(2) = latitude
    sat(3) = altitude

    ! set delta T caused by geomagnetic activity
    hourofday = 1+int(24*(time_mjd-int(time_mjd)))
    fractionofhour = 24*(time_mjd-int(time_mjd))-int(24*(time_mjd-int(time_mjd)))
    if(hourofday<24) then
      dstdtc = (sga_data(idx)%dtc(hourofday+1)-sga_data(idx)%dtc(hourofday))*fractionofhour + sga_data(idx)%dtc(hourofday)
    else
      dstdtc = (sga_data(idx+1)%dtc(1)-sga_data(idx)%dtc(hourofday))*fractionofhour + sga_data(idx)%dtc(hourofday)
    end if



    ! call atmosphere model
    call JB2008(time_mjd, &
                sun,      &
                sat,      &
                f10,      &
                f10b,     &
                s10,      &
                s10b,     &
                m10,     &
                m10b,    &
                y10,      &
                y10b,     &
                dstdtc,   &
                temp,     &
                rho)
    rho = rho * 1.d9

    return

  end function getDensityJB2008

  !==============================================================
  !
  !> @brief     Returns a the daily F10.7 value for a given date
  !!
  !> @author    Vitali Braun
  !!
  !> @date      <ul>
  !!              <li> 03.06.2014 (initial design) </li>
  !!            </ul>
  !!
  !> @param[in] mjd   Modified julian day
  !!
  !> @return    f107 value / sfu
  !!
  !> @anchor    getSolarActivityForDate
  !!
  !------------------------------------
  real(dp) function getSolarActivityForDate(this,mjd) result(f107)

    class(Atmosphere_class) :: this
    real(dp), intent(in)    :: mjd

    character(len=*), parameter :: csubid = 'getSolarActivityForDate'

    integer :: i

    real(dp) :: diff      ! difference between array date and variable 'mjd'
    real(dp) :: diff_min  ! min. difference between array date and variable 'mjd'

    f107 = undefined

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether initialization has been done already...
    if(.not. this%atmosphereInitialized) then
      call setNeptuneError(E_ATMOSPHERE_INIT, FATAL)
      return
    end if

    diff_min = infinite ! arbitrary high value..

    !** loop solmag array
    do i = 1, size(sga_data)

      diff = abs(sga_data(i)%mjd - mjd)

      if(diff < eps15) then ! return immediately - near perfect match!
        f107 = sga_data(i)%f107
        exit
      else if(diff < diff_min) then
        diff_min = diff
        f107     = sga_data(i)%f107
      end if

    end do

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSolarActivityForDate

  !==============================================================
  !
  !> @brief     Returns a date for a given solar activity
  !!
  !> @author    Vitali Braun
  !!
  !> @date      <ul>
  !!              <li> 09.05.2014 (initial design) </li>
  !!            </ul>
  !!
  !> @param[in] sfu       Solar activity (sfu) the date is searched for
  !< @param[in] minDate   Minimum MJD (no epoch earlier than this date is provided)
  !!
  !> @return    date as time_t
  !!
  !> @anchor    getDateForSolarActivity
  !!
  !------------------------------------
  type(time_t) function getDateForSolarActivity(this, sfu, minDate) result(date)

    class(Atmosphere_class)             :: this
    real(dp), intent(in)                :: sfu
    real(dp), intent(in)                :: minDate

    character(len=*), parameter :: csubid = 'getDateForSolarActivity'

    integer :: i

    real(dp) :: diff      ! difference between array data and variable 'sfu'
    real(dp) :: diff_min  ! min. difference between array data and variable 'sfu'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether initialization has been done already...
    if(.not. this%atmosphereInitialized) then
      call setNeptuneError(E_ATMOSPHERE_INIT, FATAL)
      return
    end if

    diff_min = infinite ! give arbitrary but high number...

    !** loop solmag array
    do i = 1, size(sga_data)

      if(sga_data(i)%mjd < minDate) cycle

      diff = abs(sga_data(i)%f107 - sfu)

      if(diff < eps15) then ! return immediately - near perfect match!
        date%mjd = sga_data(i)%mjd
        exit
      else if(diff < diff_min) then
        diff_min = diff
        date%mjd = sga_data(i)%mjd
      end if

    end do

    !** compute missing data...
    call mjd2gd(date)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if


    return

  end function getDateForSolarActivity

  !==========================================================
  !
  !> @brief     Returns minimum altitude
  !!
  !> @author    Vitali Braun
  !!
  !> @date      <ul>
  !!              <li> 09.10.2013 (initial design) </li>
  !!              <li> 17.11.2013 (transferred to 'derivatives' module)</li>
  !!              <li> 04.02.2014 (transferred to 'atmosphere' module)</li>
  !!            </ul>
  !!
  !> @return    minimum altitude / km
  !!
  !> @anchor    getMinAltitude
  !!
  !------------------------------------
  real(dp) function getMinAltitude(this)

    class(Atmosphere_class)  :: this
    getMinAltitude = this%p_altitude_min
    return

  end function getMinAltitude

  !==========================================================
  !
  !> @brief     Sets minimum altitude
  !!
  !> @author    Vitali Braun
  !!
  !> @date      <ul>
  !!              <li> 09.10.2013 (initial design) </li>
  !!              <li> 17.11.2013 (transferred to 'derivatives' module and added minimum value of 0.d0 km)</li>
  !!              <li> 04.02.2014 (transferred to 'atmosphere' module)</li>
  !!            </ul>
  !!
  !> @anchor  setMinAltitude
  !!
  !------------------------------------
  subroutine setMinAltitude(this,minalt)
    class(Atmosphere_class) :: this
    real(dp), intent(in)    :: minalt

    this%p_altitude_min = max(minalt, 0.d0)
    return

  end subroutine setMinAltitude


!!------------------------------------------------------------------------------------------------
!
!> @anchor      getDragCovariance
!!
!> @brief       Gets the contributions to the variational equations due to atmospheric drag
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial design)</li>
!!              </ul>
!!
!> @param[in]   gravity_model     the gravity model
!> @param[in]   satellite_model   the satellite model
!> @param[in]   solarsystem_model the solar system model
!> @param[in]   reduction         reduction handler
!> @param[in]   r_gcrf            radius vector of satellite in GCRF / km
!> @param[in]   v_gcrf            velocity vector of satellite in GCRF / km/s
!> @param[in]   time_mjd          MJD
!!
!> @return      covariance as 3x6 matrix
!!
!!------------------------------------------------------------------------------------------------
  function getDragCovariance(this,              &
                            gravity_model,      &
                            satellite_model,    &
                            solarsystem_model,  &
                            reduction,          &
                            r_gcrf,             &
                            v_gcrf,             &
                            time_mjd)           &
                            result(cov)

    implicit none

    !** interface
    !---------------------------------------------------
    class(Atmosphere_class)                 :: this
    type(Gravity_class),intent(inout)       :: gravity_model
    type(Satellite_class),intent(inout)     :: satellite_model
    type(Solarsystem_class),intent(inout)   :: solarsystem_model
    type(Reduction_type),intent(inout)      :: reduction
    real(dp), dimension(3), intent(in)      :: r_gcrf
    real(dp), dimension(3), intent(in)      :: v_gcrf
    real(dp),               intent(in)      :: time_mjd

    real(dp), dimension(3,6) :: cov     ! result in 3x6 matrix
    !---------------------------------------------------

    character(len=*), parameter :: csubid =  "getDragCovariance"

    real(dp)                 :: crossSection  ! satellite cross section
    real(dp), dimension(3,3) :: mati          ! identity matrix
    real(dp), dimension(3)   :: r_itrf        ! radius vector in ITRF
    real(dp), dimension(3)   :: v_itrf        ! radius vector in ITRF
    real(dp), dimension(3,3) :: rotTIRS       ! rotation matrix to convert from GCRF to TIRS
    real(dp)                 :: vabs          ! magnitude of v_rel
    real(dp), dimension(3)   :: v_rel         ! relative velocity in GCRF

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether initialization has been done already...
    if(.not. this%atmosphereInitialized) then
      call setNeptuneError(E_ATMOSPHERE_INIT, FATAL)
      return
    end if

    !** get drag coefficient and mass
    if(satellite_model%hasChangedSatellite(ID_ATMOSPHERE) .or. this%first)  then                ! only required, if configuration has changed or called for first time

      this%cdrag = satellite_model%getObjectDragCoefficient()
      if(hasFailed()) return
      this%dmass = satellite_model%getObjectMass()
      if(hasFailed()) return
      this%cdom  = this%cdrag/this%dmass
      this%first = .false.

    end if

    !** get relative velocity (without wind!)
    v_rel = this%getRelativeVelocity(reduction, r_gcrf, v_gcrf, time_mjd)

    !** get cross-section
    crossSection = satellite_model%getObjectCrossSection(solarsystem_model, reduction, time_mjd, r_gcrf, v_rel)

    !** get r_itrf
    call reduction%inertial2earthFixed_rv(r_gcrf, v_gcrf, time_mjd, r_itrf, v_itrf)

    !** get density
    this%rho = this%getAtmosphericDensity(gravity_model, r_gcrf, v_gcrf, r_itrf, v_itrf, time_mjd)
    if(hasFailed()) return

    call identity_matrix(mati)
    vabs = mag(v_rel)

   !write(*,*) "cross = ", crossSection
   !write(*,*) "cdom  = ", cdom
   !write(*,*) "rho   = ", this%rho
   !write(*,*) "mati = ", mati
   !write(*,*) "v_rel = ", v_rel
   !write(*,*) "vabs = ", vabs
   !
   !write(66,'(f12.6,x,3(e12.5e2,x),f8.4)') time_mjd, crossSection, cdom, this%rho, vabs

    !** build part of matrix, which contains derivatives wrt. the velocity
    cov(1:3,4:6) = -0.5d0*crossSection*this%cdom*this%rho*(vabs*mati + outerproduct(v_rel, v_rel)/vabs)

    !** correct for dv_r/dv (equivalent to transformation between TIRS and GCRF)
    rotTIRS      = reduction%getRotationMatrixRC2TI(time_mjd)
    !write(*,*) "rotTIRS = "
    !write(*,*) rotTIRS(1,1:3)
    !write(*,*) rotTIRS(2,1:3)
    !write(*,*) rotTIRS(3,1:3)
!
!    cov(1:3,4:6) = matmul(cov(1:3,4:6), transpose(rotTIRS))

    !write(*,*) "out/vabs = ", outerproduct(v_rel, v_rel)/vabs
    !read(*,*)

    !** build part of matrix, which contains derivatives wrt. the position

    cov(1:3,1) = -cov(1:3,5)*getEarthRotation(0.d0, EGM96)
    cov(1:3,2) = cov(1:3,4)*getEarthRotation(0.d0, EGM96)
    cov(1:3,3) = 0.d0
    !cov(1:3,1:3) = -getEarthRotation(0.d0,EGM96)*matmul(cov(1:3,4:6), transpose(rotTIRS))

!   write(*,*) "cov(1:3,1:3) = "!, cov(1:3,1), cov(1:3,4), getEarthRotation(0.d0,EGM96)
!   write(*,*) cov(1,1:3)
!   write(*,*) cov(2,1:3)
!   write(*,*) cov(3,1:3)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getDragCovariance

!-----------------------------------------------------------------------------------------------
!
!> @anchor      setAtmosphereModel
!!
!> @brief       Set the atmosphere model to be used
!> @author      Vitali Braun
!!
!> @param[in]   imodel - either MSIS2000 or EXPONENTIAL
!!
!> @date        <ul>
!!                <li>23.08.2015: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  subroutine setAtmosphereModel(this,imodel)

    class(Atmosphere_class) :: this
    integer, intent(in)     :: imodel

    character(len=*), parameter :: csubid = 'setAtmosphereModel'
    character(len=3) :: cerr
    integer          :: currentModel

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    currentModel = this%nmodel ! save current model

    select case(imodel)
      case(JB08)
        this%nmodel = JB08
      case(MSIS2000)
        this%nmodel = MSIS2000
      case(EXPONENTIAL)
        this%nmodel = EXPONENTIAL
      case default
        write(cerr,'(i3)') imodel
        call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/cerr/))
        return
    end select

    !** after model change, a re-initialization is required
    if(this%nmodel /= currentModel) then
      this%atmosphereInitialized = .false.
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine setAtmosphereModel

!-----------------------------------------------------------------------------------------------
!
!> @anchor      getAtmosphereModelName
!!
!> @brief       Get the name of the atmosphere model used
!> @author      Vitali Braun
!!
!> @return      model name as string
!!
!> @date        <ul>
!!                <li> 21.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=len(modelName)) function getAtmosphereModelName(this)

    class(Atmosphere_class)  :: this
    getAtmosphereModelName = modelName(this%nmodel)
    return

  end function getAtmosphereModelName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      getWindModelName
!!
!> @brief       Get the name of the horizontal wind model used
!> @author      Vitali Braun
!!
!> @return      model name as string
!!
!> @date        <ul>
!!                <li> 25.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=len(windModelName)) function getWindModelName(this)

    class(Atmosphere_class)  :: this
    getWindModelName = windModelName(this%nwmodel)
    return

  end function getWindModelName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      getHorizontalWindFileName
!!
!> @brief       Get the name of the data file for the horizontal wind model used
!> @author      Vitali Braun
!!
!> @return      file name as string
!!
!> @date        <ul>
!!                <li> 29.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=255) function getHorizontalWindFileName(this)

    class(Atmosphere_class)  :: this
    getHorizontalWindFileName = this%hwmDataFile
    return

  end function getHorizontalWindFileName
!-----------------------------------------------------------------------------------------------
!
!> @anchor      setHorizontalWindFileName
!!
!! @brief       Set the name of the data file for the horizontal wind model used
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 29.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  subroutine setHorizontalWindFileName(this,val)

    class(Atmosphere_class)         :: this
    character(len=*), intent(in)    :: val

    this%hwmDataFile = val(1:min(len(val),len(this%hwmDataFile)))
    return

  end subroutine setHorizontalWindFileName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      getSolmagFileVersion
!!
!> @brief       Get the version and update date of the SOLMAG data files
!> @author      Vitali Braun
!!
!> @return      file version as string
!!
!> @param[in]   typ   type of data file (DAILY or MONTHLY)
!!
!> @date        <ul>
!!                <li>VB: 23.10.2016: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=24) function getSolmagFileVersion(this,typ) result(res)

    implicit none

    class(Atmosphere_class) :: this
    integer, intent(in)     :: typ

    if(typ == DAILY) then
        res = this%solmagFileVersion//' '//date2string(this%solmagFileDateDaily)
    else if(typ == MONTHLY) then
        if(this%monthlyDataRequired) then
            res = this%solmagFileVersion//' '//date2string(this%solmagFileDateMonthly)
        else
            res = ' - Not required - '
        end if
    else
        res = '---'
    end if
    return
  end function


!-----------------------------------------------------------------------------------------------
!
!> @anchor      getSolmagFileName
!!
!> @brief       Get the name of the data file for the solar and geomagnetic activity used
!> @author      Vitali Braun
!!
!> @param[in]   typ   type of data file (DAILY or MONTHLY)
!!
!> @return      file name as string
!!
!> @date        <ul>
!!                <li> 30.10.2013: initial design</li>
!!                <li> 02.11.2013: added DAILY/MONTHLY switch</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=255) function getSolmagFileName(this,typ)

    class(Atmosphere_class) :: this
    integer, intent(in)     :: typ                                              ! type of data file (DAILY or MONTHLY)

    select case(typ)
      case(DAILY)
        getSolMagFileName = this%sgaDataFile
      case(MONTHLY)
        getSolMagFileName = this%sgaMonthlyDataFile
      case default
        getSolMagFileName = this%sgaDataFile
    end select

    return

  end function getSolMagFileName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      setSolMagFileName
!!
!> @brief       Set the name of the data file for the solar and geomagnetic activity used
!> @author      Vitali Braun
!!
!> @param[in]   val   name of solmag data file
!> @param[in]   typ   type of data file (DAILY or MONTHLY)
!!
!> @date        <ul>
!!                <li> 30.10.2013: initial design</li>
!!                <li> 02.11.2013: added DAILY/MONTHLY switch</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  subroutine setSolMagFileName(this, val, typ)

    class(Atmosphere_class)         :: this
    integer,          intent(in)    :: typ
    character(len=*), intent(in)    :: val

    select case(typ)
      case(DAILY)
        this%sgaDataFile = val(1:min(len(val),len(this%sgaDataFile)))
      case(MONTHLY)
        this%sgaMonthlyDataFile = val(1:min(len(val),len(this%sgaMonthlyDataFile)))
      case default
        this%sgaDataFile = val(1:min(len(val),len(this%sgaDataFile)))
    end select

    return

  end subroutine setSolMagFileName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      getDisturbanceWindFileName
!!
!> @brief       Get the name of the data file for the disturbance wind model used
!> @author      Vitali Braun
!!
!> @return      file name as string
!!
!> @date        <ul>
!!                <li> 29.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=255) function getDisturbanceWindFileName(this)

    class(Atmosphere_class)  :: this
    getDisturbanceWindFileName = this%distWindDataFile
    return

  end function getDisturbanceWindFileName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      setDisturbanceWindFileName
!!
!> @brief       Set the name of the data file for the disturbance wind model used
!> @author      Vitali Braun
!!
!> @param[in]   val   name of disturbance wind file
!!
!> @date        <ul>
!!                <li> 29.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  subroutine setDisturbanceWindFileName(this,val)

    class(Atmosphere_class)     :: this
    character(len=*),intent(in) :: val

    this%distWindDataFile = val(1:min(len(val),len(this%distWindDataFile)))
    return

  end subroutine setDisturbanceWindFileName


!-----------------------------------------------------------------------------------------------
!
!> @anchor      getQDGridFileName
!!
!> @brief       Get the name of the data file for the quasi-dipole interpolation grid (HWM)
!> @author      Vitali Braun
!!
!> @return      file name as string
!!
!> @date        <ul>
!!                <li> 29.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  character(len=255) function getQDGridFileName(this)

    class(Atmosphere_class)  :: this
    getQDGridFileName = this%qdGridDataFile
    return

  end function getQDGridFileName

!-----------------------------------------------------------------------------------------------
!
!> @anchor      setQDGridFileName
!!
!> @brief       Set the name of the data file for the quasi-dipole interpolation grid (HWM)
!> @author      Vitali Braun
!!
!> @param[in]   val   name of QD grid file
!!
!> @date        <ul>
!!                <li> 29.10.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------------------------
  subroutine setQDGridFileName(this,val)

    class(Atmosphere_class)         :: this
    character(len=*), intent(in)    :: val

    this%qdGridDataFile = val(1:min(len(val),len(this%qdGridDataFile)))
    return

  end subroutine setQDGridFileName

!==============================================================================================
!
!> @anchor      getRelativeVelocity
!!
!> @brief       Get satellite's relative velocity wrt. the Earth's atmosphere
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 24.01.2013 (initial design)</li>
!!                <li> 25.10.2013 (added HWM07 wind model)</li>
!!                <li> 04.02.2014 (changed behaviour, so that for wind model, solmag quantities are computed instead of being passed) </li>
!!              </ul>
!!
!> @param[in]   reduction   reduction handler
!> @param[in]   r_inr       radius vector wrt. inertial frame / km
!> @param[in]   v_inr       velocity with respect to an inertial frame / km/s
!> @param[in]   mjd         current MJD
!!
!> @return      file name as string
!!
!> @details     This routine returns the relative velocity vector
!!              with respect to the Earth's atmosphere. The horizontal wind model HWM07 is
!!              used besides the rotation of the atmosphere.
!!
!> @todo        Check influence of LOD on relative velocity
!!
!!------------------------------------------------------------------------------------------------
  function getRelativeVelocity(this, reduction, r_inr, v_inr, mjd) result(v_rel)

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    type(Reduction_type)               :: reduction
    real(dp), dimension(3), intent(in)  :: r_inr
    real(dp), dimension(3), intent(in)  :: v_inr
    real(dp),               intent(in)  :: mjd

    real(dp), dimension(3) :: v_rel

    !-------------------------------------------

    integer  :: hour
    integer  :: idx
    integer  :: iyd
    real(dp) :: sec
    real(dp) :: alt
    real(dp) :: glat
    real(dp) :: glon
    real(dp) :: stl

    character(len=*), parameter :: csubid = "getRelativeVelocity"

    real(dp), dimension(2) :: ap     ! Ap array, however, currently only ap2(2) is used by HWM07
    real(dp), dimension(2) :: horizontalWindSpeed   ! horizontal wind speed as resulting from HWM07 (m/s)
    real(dp), dimension(3) :: rdummy  ! dummy radius vector
    real(dp), dimension(3) :: r_itrf  ! radius in itrf frame
    real(dp), dimension(3) :: vw_sez  ! wind speed in SEZ frame
    real(dp), dimension(3) :: vw_gcrf ! wind speed in GCRF
    real(dp), dimension(3) :: vw_itrf ! wind speed in ITRF

    real(4) :: sec_4   ! real value of sec
    real(4) :: alt_4   ! real value of alt
    real(4) :: glat_4  ! real value of glat
    real(4) :: glon_4  ! real value of glon
    real(4) :: stl_4   ! real value of stl
    real(4) :: f107, f107a  ! auxiliaries, not used
    real(4), dimension(2) :: ap_4    ! real value of ap
    real(4), dimension(2) :: w_4     ! horizontal wind components - real value

    !real(4), dimension(2) :: apqt_4    ! real value of ap
    !real(4), dimension(2) :: qw,dw,w
    !integer :: iyd_tmp, i
    !real(4), external :: pershift
    !real(4) :: kp, mlt, mlat, mmpwind, mzpwind

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** Get radius in ITRF (NOTE: the reductions are only performed as long
    !   as the rotation matrices are not available for the given epoch.
    !   When the function is called here, the matrices should already be
    !   there.
    !-----------------------------------------------------------------------
    call reduction%inertial2earthFixed(r_inr, mjd, r_itrf)
    !** handle relative velocity contribution due to horizontal wind
    !-----------------------------------------------------------------
    if(this%considerWind) then

      ap(1) = 0.d0

      !** get ap(2) which is the current 3hr ap index
      idx   = int(mjd - int(sga_data(1)%mjd)) + 1
      hour  = int(mod(mjd, 1.d0)*24.d0)
      ap(2) = sga_data(idx)%ap(int(hour/3.d0)+1)

      !** get time in yyddd format
      iyd   = mjd2yyddd(mjd)

      !** seconds
      sec   = mjd2daySeconds(mjd)

      !** altitude, longitude and latitude
      call getGeodeticLatLon(r_itrf, alt, glat, glon)

      glat = glat*rad2deg
      glon = glon*rad2deg

      !** local solar time
      stl = getLocalSolarTime(mjd, glon)

      !** convert to real(kind=4) numbers
      sec_4  = real(sec, 4)
      alt_4  = real(alt, 4)
      glat_4 = real(glat,4)
      glon_4 = real(glon,4)
      stl_4  = real(stl, 4)
      ap_4   = real(ap,  4)

      call this%hwm07_model%HWM07(iyd, sec_4, alt_4, glat_4 ,glon_4, stl_4, f107a, f107, ap_4, w_4)

      horizontalWindSpeed = real(w_4, dp)

      !** convert horizontal wind speed to inertial frame (GCRF), in order to add it to v_rel
      !----------------------------------------------------------------------------------------
      vw_sez(1) = -horizontalWindSpeed(1)                                       ! positive in south direction
      vw_sez(2) =  horizontalWindSpeed(2)                                       ! positive in east direction
      vw_sez(3) =  0.d0                                                         ! only horizontal wind considered.

      vw_sez = vw_sez*1.d-3   ! in km/s

      !** convert from SEZ to ITRF
      call reduction%sez2itrf(vw_sez, glat*deg2rad, glon*deg2rad, vw_itrf)

    else

      vw_itrf = 0.d0

    end if

    !** convert from ITRF to GCRF
    call reduction%earthFixed2inertial(r_itrf, vw_itrf, mjd, rdummy, vw_gcrf)

    !** in order to obtain relative velocity, subtract wind in GCRF from inertial velocity vector
    v_rel = v_inr - vw_gcrf

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getRelativeVelocity

!==============================================================================================
!
!> @anchor      getSolMagTime
!!
!> @brief       Get time of solar and geomagnetic activity observations
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 18.01.2013 (initial design)</li>
!!              </ul>
!!
!> @param[in]   year        Year of given date
!> @param[in]   month       Month of given date
!!
!> @return      20 if beyond 1991 or 17 for any other case
!!
!> @details     This function returns the time in UT for which solar and geomagnetic activity observations were
!!              recorded. Observations were made at 1700 UT until May 31, 1991 and at
!!              2000 UT from that time on.
!!              (Vallado, D., Using EOP and Space Weather Data for Satellite Operations, AAS/AIAA, 2005)
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Check given date: </li>
!!                <li> Return '20' if it is beyond May 1991 </li>
!!                <li> Return '17' else </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  integer function getSolMagTime(year, month)

    integer, intent(in)                 :: year
    integer, intent(in)                 :: month

    character(len=*), parameter :: csubid = "getSolMagTime"

    getSolMagTime = - 1

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if


    if(year < 1991) then
      getSolMagTime = 17
    else
      if(year == 1991 .and. month < 6) then
        getSolMagTime = 17
      else
        getSolMagTime = 20
      end if
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSolMagTime

!========================================================================================
!
!> @anchor      setApForecast
!!
!> @brief       Set value of Ap for long-term forecast
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 19.01.2013 (initial design)</li>
!!                <li> 04.11.2013 (changed functionality: value of Ap is now set directly)</li>
!!              </ul>
!!
!> @param[in]   val   desired Ap level for long-term forecast
!!
!> @details     This routine sets the parameter 'ap_predict', which is required for long-term
!!              geomagnetic activity forecasts.
!!
!!------------------------------------------------------------------------------------------------
  subroutine setApForecast(this,val)

    class(Atmosphere_class) :: this
    integer,intent(in)      :: val

    this%ap_predict = max(val,0)
    return

  end subroutine setApForecast

!--------------------------------------------------------------------------------------------
!
!> @anchor      getApForecast
!!
!! @brief       Get value of Ap for long-term forecast
!! @author      Vitali Braun
!!
!> @return      AP prediction (integer)
!!
!! @date        <ul>
!!                <li> 22.10.2013 (initial design)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  integer function getApForecast(this)

    class(Atmosphere_class)  :: this
    getApForecast = this%ap_predict
    return

  end function getApForecast

!========================================================================================
!
!> @anchor      setSolarFluxForecast
!!
!> @brief       Set value of solar flux for long-term forecast
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.11.2013 (initial design)</li>
!!              </ul>
!!
!> @param[in]   val   Solar flux in sfu
!!
!> @details     This function sets the parameter 'sol_predict', which is required for long-term
!!              solar activity forecasts beyond available data.
!!
!!------------------------------------------------------------------------------------------------
  subroutine setSolarFluxForecast(this,val)

    class(Atmosphere_class) :: this
    real(dp), intent(in)    :: val

    this%sol_predict = max(val,0.d0)

  end subroutine setSolarFluxForecast

!========================================================================================
!
!> @anchor      getSolarFluxForecast
!!
!> @brief       Get value of solar flux for long-term forecast
!> @author      Vitali Braun
!!
!> @return      f10.7 forecast / sfu
!!
!> @date        <ul>
!!                <li> 04.11.2013 (initial design)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  real(dp) function getSolarFluxForecast(this)

    class(Atmosphere_class)  :: this
    getSolarFluxForecast = this%sol_predict
    return

  end function getSolarFluxForecast

!!------------------------------------------------------------------------------------------------
!> @anchor      getDragCovFlag
!!
!> @brief       Gets the appropriate flag covDrag
!> @author      Vitali Braun
!!
!> @return      .true. or .false.
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial design)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  logical function getDragCovFlag(this)

    class(Atmosphere_class)  :: this
    getDragCovFlag = this%covDrag
    return

  end function getDragCovFlag

!!------------------------------------------------------------------------------------------------
!> @anchor      setDragCovFlag
!!
!> @brief       Sets the appropriate flag covDrag
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2014 (initial design)</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  subroutine setDragCovFlag(this,val)

    class(Atmosphere_class)  :: this
    logical, intent(in)      :: val

    this%covDrag = val
    return

  end subroutine setDragCovFlag


!===============================================================================
!
!> @anchor      kp2ap
!!
!> @brief       Convert from planetary index Kp to planetary amplitude Ap
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 31.10.2013 (initial design)</li>
!!              </ul>
!!
!> @param[in]   kp   value for kp
!!
!> @return      ap value
!!
!> @details     This function computes for a given geomagnetic activity value of kp
!!              (planetary index) the equivalent planetary amplitude (Ap) value,
!!              based on the algorithm as given in Vallado,D., "Fundamentals of
!!              Astrodynamics and Applications", 4th Edition, p.558, Microcosm Press, 2013.
!!
!!------------------------------------------------------------------------------------------------
  real(dp) elemental function kp2ap(this,kp)

    class(Atmosphere_class),intent(in)  :: this
    real(dp), intent(in)                :: kp

    integer :: idx
    real(dp), parameter, dimension(28) :: convArray = (/0.d0, 2.d0,   3.d0,   4.d0,   5.d0,   6.d0,   7.d0, &
                                             9.d0,  12.d0,  15.d0,  18.d0,  22.d0,  27.d0,  32.d0, &
                                            39.d0,  48.d0,  56.d0,  67.d0,  80.d0,  94.d0, 111.d0, &
                                           132.d0, 154.d0, 179.d0, 207.d0, 236.d0, 300.d0, 400.d0/)
    !** compute array index
    idx = int(kp*3.d0) + 1

    kp2ap = convArray(idx)

  end function kp2ap

  !=============================================================================
  !
  !> @anchor      get_jb2008_date_interval
  !!
  !! @brief       Getting first and last MJD entry in JB2008 solmag file
  !! @author      Christopher Kebschull (ChK)
  !!
  !! @date        <ul>
  !!                <li> 14.01.2021 (ChK: initial design)</li>
  !!              </ul>
  !!
  !!
  !! @details     This routine reads solar and geomagnetic and returns the last MJD entry
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine get_jb2008_date_interval( this, &
                                  ich,       &  ! <-- INT     input channel to read data from
                                  mjd_start, &  ! <--> DBL    MJD of first SGA entry
                                  mjd_end    &  ! <--> DBL    MJD of last SGA entry
                                )

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    integer,  intent(in)                :: ich
    real(dp), intent(inout)             :: mjd_start
    real(dp), intent(inout)             :: mjd_end
    character(len=255)                  :: cbuf,cbuf_dtc
    !-------------------------------------------

    real(dp)        :: dummy
    integer         :: year_sol
    integer         :: doy_sol
    integer         :: ios                      ! I/O status
    logical         :: first
    integer         :: dtc_ich

    character(len=*), parameter :: csubid = "get_jb2008_date_interval"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    rewind(ich)

    dtc_ich = openFile(trim(adjustl(this%dataPath))//cdelimit//"DTCFILE.TXT", SEQUENTIAL, IN_FORMATTED)

    first = .true.
    do
      read(ich,'(A)',iostat=ios) cbuf
      if (ios < 0) exit
      if(index(cbuf,'#') /= 0) cycle

      ! Now also start reading the DTC file
      read(dtc_ich,'(A)',iostat=ios) cbuf_dtc
      if (ios < 0) exit

      read (cbuf,*) year_sol, doy_sol, mjd_end
      mjd_end = mjd_end - 2400000.5
      if (first) then
        first= .false.
        mjd_start = mjd_end
      end if
    end do

    dtc_ich = closeFile(dtc_ich)

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine get_jb2008_date_interval

!=============================================================================
!
!> @anchor      read_jb2008_files
!!
!! @brief       Reading solar and geomagnetic activity data from CSSI data file
!! @author      Daniel Lück (DLU)
!!
!! @date        <ul>
!!                <li> 07.01.2021 (DLU: initial design)</li>
!!              </ul>
!!
!!
!! @details     This routine reads solar and geomagnetic activity files for JB2008
!!
!> @todo        This will crash if a wrong date is used
!!------------------------------------------------------------------------------------------------

  subroutine read_jb2008_files(  this,      &
                              ich,       &  ! <-- INT     input channel to read data from
                              mjd_start, &  ! --> DBL     MJD of first SGA data entry
                              mjd_end    &  ! --> DBL     MJD of last SGA entry
                            )

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    integer,  intent(in)                :: ich
    real(dp), intent(in)                :: mjd_start
    real(dp), intent(in)                :: mjd_end
    !-------------------------------------------

    real(dp)        :: JD_start
    real(dp)        :: JD_end
    real(dp)        :: JD_temp
    integer         :: year_sol, year_mag
    integer         :: doy_sol, doy_mag
    integer         :: i
    character(len=3):: buffer
    integer         :: dtc_ich

    character(len=*), parameter :: csubid = "read_jb2008_files"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    write(*,*) '- Reading JB2008 files...'

    rewind(ich)

    dtc_ich = openFile(trim(adjustl(this%dataPath))//cdelimit//"DTCFILE.TXT", SEQUENTIAL, IN_FORMATTED)

    do i = 1, 4
      read (ich,*)
    end do

    JD_temp = 0.0
    do while(JD_temp < mjd_start + 2399999.5D0)
      read (ich,*) year_sol,doy_sol, JD_temp
    end do
    year_mag = 0
    doy_mag = 0
    do while(.not.(year_mag==year_sol .and. doy_mag==doy_sol))
      read (dtc_ich,*) buffer, year_mag, doy_mag
    end do

    i = 1
    do while(JD_temp < mjd_end + 2400000.5)

      read (ich,*) year_sol, doy_sol, JD_temp, sga_data(i)%f10_jb, sga_data(i)%f10b_jb, sga_data(i)%s10_jb, &
      sga_data(i)%s10b_jb, sga_data(i)%m10_jb,sga_data(i)%m10b_jb,sga_data(i)%y10_jb,sga_data(i)%y10b_jb

      read (dtc_ich,*) buffer, year_mag, doy_mag, sga_data(i)%dtc

      sga_data(i)%mjd = JD_temp - 2400000.5
      i = i+1
    end do

    dtc_ich = closeFile(dtc_ich)

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine read_jb2008_files

!===============================================================================
!
!> @anchor      getSgaDataFileType
!!
!> @brief       Determine the type of solar & geomagnetic data file which is used for the propagation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 30.10.2013 (initial design)</li>
!!              </ul>
!!
!> @param[in]   ich   input channel to read from
!!
!> @details     This function determines the type of solar and geomagnetic data file
!!              by checking the first few lines for typical keywords.
!!
!!------------------------------------------------------------------------------------------------
  integer function getSgaDataFileType(this,ich)

    class(Atmosphere_class) :: this
    integer, intent(in)     :: ich                                              ! input channel

    character(len=*), parameter :: csubid = "getSgaDataFileType"
    character(len=255) :: cbuf

    getSgaDataFileType = -1

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** read first line of file
    read(ich,'(a)') cbuf

    if(index(cbuf,fileIdent(FILE_CSSI)) /= 0) then
      getSgaDataFileType = FILE_CSSI
    else if(index(cbuf,fileIdent(FILE_FAPDAY)) /= 0) then
      getSgaDataFileType = FILE_FAPDAY
    else if(index(cbuf,fileIdent(FILE_JB08)) /= 0) then
      getSgaDataFileType = FILE_JB08
    else
      call setNeptuneError(E_SGA_UNSUPPORTED_TYPE, FATAL, suppSgaFiles)
      getSgaDataFileType = -1
      return
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end function getSgaDataFileType

!=============================================================================
!
!> @anchor      readCssiDataFile
!!
!! @brief       Reading solar and geomagnetic activity data from CSSI data file
!! @author      Vitali Braun (VB)
!!
!! @date        <ul>
!!                <li> 31.10.2013 (VB: initial design)</li>
!!                <li> 18.07.2016 (VB: fixed several issues, where especially predicted information was incorrectly parsed)</li>
!!              </ul>
!!
!! @param[in]   ich         channel to read solmag data from
!!
!! @details     This routine parses the CSSI data file for the supported versions. It
!!              finds the required data which is required for the propagation and stores
!!              it within the sga data array.
!!
!!------------------------------------------------------------------------------------------------

  subroutine readCssiDataFile(  this,      &
                                ich,       &  ! <-- CHR()   input channel to read data from
                                mjd_start, &  ! --> DBL     MJD of first SGA data entry
                                mjd_end    &  ! --> DBL     MJD of last SGA entry
                             )

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    integer,  intent(in)                :: ich
    real(dp), intent(in)                :: mjd_start
    real(dp), intent(in)                :: mjd_end
    !-------------------------------------------

    character(len=*), dimension(1), parameter :: c_sga_version = (/"1.2"/)      ! supported CSSI data file versions
    character(len=*), parameter :: csubid = "readCssiDataFile"
    character(len=255) :: cbuf          ! character buffer
    character(len=255) :: cmess         ! message string
    character(len=10)  :: ctemp10   ! temporary string
    character(len=3)   :: ctemp3    ! temporary string
    character(len=3)   :: cversion      ! version as read from data file

    integer :: day                  ! day
    integer :: ios                  ! I/O status
    integer :: i,j,k                ! loop counter
    integer :: idx_start            ! interpolation array index (lower boundary)
    integer :: idx_end              ! interpolation array index (upper boundary)
    integer :: month                ! month
    integer :: year                 ! year
    integer :: m_pred_points        ! number of monthly predicted data points
    integer :: firstMonthlyIndex    ! index in monthly array which provides a starting point for the start date being within the monthly predicted values

    real(dp)  :: dt                       ! time offset for solar activity fit
    real(dp)  :: dtemp                    ! temporary
    real(dp)  :: c1,c2,c3,c4,c5,c6        ! monthly fit coefficients
    real(dp), dimension(:,:), allocatable :: monthlies    ! array containing monthly values for F10.7 and F10.7 81ctr (obs)

    logical :: flag_begin_obs           ! BEGIN OBSERVED tag has been found
    logical :: flag_begin_daily_pred    ! BEGIN DAILY PREDICTED tag has been found
    logical :: flag_begin_monthly_pred  ! BEGIN MONTHLY PREDICTED tag has been found
    logical :: flag_end_obs             ! END OBSERVED tag has been found
    logical :: flag_end_daily_pred      ! END DAILY PREDICTED tag has been found
    logical :: flag_end_monthly_pred    ! END MONTHLY PREDICTED tag has been found
    logical :: flag_obs_only            ! Solar and geomagnetic activity completely through available observation data
    logical :: flag_obs_pred            ! Solar and geomagnetic activity covered by observed and daily predicted data
    logical :: flag_obs_m_pred          ! Solar and geomagnetic activity covered by observed and monthly predicted data
    logical :: flag_obs_fit             ! Solar and geomagnetic activity covered by observed, monthly predicted and monthly fit data
    logical :: flag_pred_only           ! Solar and geomagnetic activity covered by daily predicted data
    logical :: flag_pred_m_pred         ! Solar and geomagnetic activity covered by daily and monthly predicted data
    logical :: flag_pred_fit            ! Solar and geomagnetic activity covered by daily and monthly predicted data, as well as monthly fit
    logical :: flag_m_pred_only         ! Solar and geomagnetic activity completely covered by monthly predictions
    logical :: flag_m_pred_fit          ! Solar and geomagnetic activity covered by monthly predictions and monthly fit
    logical :: flag_fit_only            ! monthly fit forecast for complete interval required

    type(time_t) :: sga_time_tag     ! time tag of SGA data file
    type(sga_t)  :: sga_obs_start    ! time of first observation
    type(sga_t)  :: sga_obs_end      ! time of last observation
    type(sga_t)  :: sga_pred_start   ! time of first daily predicted observation
    type(sga_t)  :: sga_pred_end     ! time of last daily predicted observation
    type(sga_t)  :: sga_m_pred_start ! time of first daily predicted observation
    type(sga_t)  :: sga_m_pred_end   ! time of last daily predicted observation
    type(sga_t)  :: spwTemp          ! auxiliary
    type(time_t) :: system_time      ! system time
    type(time_t) :: tempDate         ! auxiliary

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** make sure that file is positioned at the beginning
    !-----------------------------------------------------
    rewind(ich)

    !** initialise flags
    !------------------------------------
    flag_begin_obs           = .false.
    flag_end_obs             = .false.
    flag_begin_daily_pred    = .false.
    flag_end_daily_pred      = .false.
    flag_begin_monthly_pred  = .false.
    flag_end_monthly_pred    = .false.

    flag_obs_only            = .false.
    flag_obs_pred            = .false.
    flag_obs_m_pred          = .false.
    flag_obs_fit             = .false.
    flag_pred_only           = .false.
    flag_pred_m_pred         = .false.
    flag_pred_fit            = .false.
    flag_m_pred_only         = .false.
    flag_m_pred_fit          = .false.
    flag_fit_only            = .false.
    !-------------------------------------

    !** find 'VERSION' key word
    !------------------------------------------
    do
        read(ich,'(A)',iostat=ios) cbuf
        if(ios /= 0) then
            call setNeptuneError(E_SGA_NO_VERSION, FATAL, (/c_sga_version(1)/))
            return
        else if(index(cbuf,'VERSION') /= 0) then
            read(cbuf(index(cbuf,'ION')+4:index(cbuf,'ION')+6),'(a)') cversion
            if(.not. any(c_sga_version == cversion)) then
                call setNeptuneError(E_SGA_UNSUPPORTED_VERSION, FATAL)
                return
            else
                ! store global variable
                this%solmagFileVersion = cversion
                exit
            end if
        end if
    end do

    rewind(ich)

    !** find 'UPDATED' key word
    !------------------------------------------
    do
        read(ich,'(A)',iostat=ios) cbuf
        if(ios /= 0) then
            call setNeptuneError(E_SGA_TIME_TAG, WARNING, (/"N/A"/))
        else if(index(cbuf,'UPDATED') /= 0) then
            exit
        end if
    end do

    !** read epoch
    read(cbuf(index(cbuf,'UPDATED')+ 8:index(cbuf,'UPDATED')+11),*) sga_time_tag%year
    read(cbuf(index(cbuf,'UPDATED')+13:index(cbuf,'UPDATED')+15),*) ctemp3

    call getMonthNumber(ctemp3,sga_time_tag%month)

    if(hasFailed() .and. isControlled()) then
        call checkOut(csubid)
        return
    end if

    read(cbuf(index(cbuf,'UPDATED')+17:index(cbuf,'UPDATED')+18),*) sga_time_tag%day
    read(cbuf(index(cbuf,'UPDATED')+20:index(cbuf,'UPDATED')+21),*) sga_time_tag%hour
    read(cbuf(index(cbuf,'UPDATED')+23:index(cbuf,'UPDATED')+24),*) sga_time_tag%minute
    read(cbuf(index(cbuf,'UPDATED')+26:index(cbuf,'UPDATED')+27),*) sga_time_tag%second

    call gd2mjd(sga_time_tag)
    if(hasFailed()) return

    ! save global variable
    this%solmagFileDateDaily = sga_time_tag

    !** compare to system time
    system_time = getDateTimeNow()
    if(hasFailed()) return

    if((system_time%mjd - sga_time_tag%mjd) > 60.d0) then
      cmess = "Solar and geomagnetic activity data file has to be updated."
      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
    end if

    cmess = "Using solar and geomagnetic activity data file from '"//date2string(sga_time_tag)//"'."
    call setNeptuneError(E_SPECIAL, REMARK, (/cmess/))

    rewind(ich)

    !** read first and last obs date as well as first and last
    !   daily and monthly predicted date in data file
    !------------------------------------------
    do

      read(ich,'(A)',iostat=ios) cbuf
      if((index(cbuf,'BEGIN OBSERVED') /= 0) .and. (.not. flag_begin_obs)) then

        read(ich,'(A)') cbuf
        sga_obs_start  = this%parseCSSIDataLine(cbuf)
        flag_begin_obs = .true.
        if(hasFailed()) return

      else if((index(cbuf,'END OBSERVED')/= 0) .and. (.not. flag_end_obs)) then

        backspace(ich)
        backspace(ich)

        read(ich,'(A)') cbuf
        sga_obs_end  = this%parseCSSIDataLine(cbuf)
        flag_end_obs = .true.
        if(hasFailed()) return

      else if((index(cbuf,'BEGIN DAILY_PREDICTED') /= 0) .and. (.not. flag_begin_daily_pred)) then

        read(ich,'(A)') cbuf
        sga_pred_start        = this%parseCSSIDataLine(cbuf)
        flag_begin_daily_pred = .true.
        if(hasFailed()) return

      else if((index(cbuf,'END DAILY_PREDICTED') /= 0) .and. (.not. flag_end_daily_pred)) then

        backspace(ich)
        backspace(ich)

        read(ich,'(A)') cbuf    ! last predicted date now in cbuf
        sga_pred_end        = this%parseCSSIDataLine(cbuf)
        flag_end_daily_pred = .true.
        if(hasFailed()) return

      else if((index(cbuf,'BEGIN MONTHLY_PREDICTED') /= 0) .and. (.not. flag_begin_monthly_pred)) then

        read(ich,'(A)') cbuf    ! first daily predicted date now in cbuf
        sga_m_pred_start        = this%parseCSSIDataLine(cbuf)
        flag_begin_monthly_pred = .true.
        if(hasFailed()) return

      else if((index(cbuf,'END MONTHLY_PREDICTED') /= 0) .and. (.not. flag_end_monthly_pred)) then

        backspace(ich)
        backspace(ich)

        read(ich,'(A)') cbuf    ! last predicted date now in cbuf
        sga_m_pred_end        = this%parseCSSIDataLine(cbuf)
        flag_end_monthly_pred = .true.
        if(hasFailed()) return

      else if((index(cbuf,'MONTHLY_FIT_TYPE') /= 0)) then

        !** read equation for monthly fit
        !------------------------------------------

        !** get constant term
        read(cbuf(index(cbuf,'=')+1:index(cbuf,'+')-1),*) c1
        !** get amplitude
        read(cbuf(index(cbuf,'+')+1:index(cbuf,'*')-1),*) c2
        !** get linear term in cos argument
        read(cbuf(index(cbuf,'COS(')+4:index(cbuf,'*T')-1),*) c3
        !** get amplitude of sin term
        read(cbuf(index(cbuf,'*T')+5:index(cbuf,'*SIN')-1),*) c4
        !** get sin argument
        read(cbuf(index(cbuf,'SIN(')+4:index(cbuf,'*T))')-1),*) c5
        !** get offset
        read(cbuf(index(cbuf,'MJD -')+5:len_trim(cbuf)),*) c6

      else if(index(cbuf,'MONTHLY_PREDICTED_POINTS') /= 0) then ! number of predicted points, required for later interpolation

        read(cbuf(index(cbuf,'NTS')+3:len_trim(cbuf)),*) m_pred_points
        allocate(monthlies(m_pred_points+1,3))  ! 1=date, 2=f107, 3=f81ctr

        !** put last predicted daily data into monthly data temporary array
        monthlies(1,1) = sga_pred_end%mjd
        monthlies(1,2) = sga_pred_end%f107
        monthlies(1,3) = sga_pred_end%f81ctr

      else if(ios /= 0) then
        exit
      end if

    end do

    rewind(ich)

    !** now check for begin and end date within the observed/predicted
    !   interval
    !-------------------------------------------
    if(sga_obs_start%mjd .gt. mjd_start) then

      call setNeptuneError(E_SGA_MISSING, FATAL)
      return

    else if((mjd_start >= sga_obs_start%mjd) .and. &
            (mjd_end   <= sga_obs_end%mjd)) then ! interval completely covered by observed data
      flag_obs_only = .true.
    else if((mjd_start >= sga_obs_start%mjd) .and. &
             (mjd_end  <= sga_pred_end%mjd)) then  ! interval covered by observed and daily predicted data
      flag_obs_pred = .true.
    else if((mjd_start >= sga_obs_start%mjd) .and. &
            (mjd_start <= sga_obs_end%mjd)   .and. &
            (mjd_end   <= sga_m_pred_end%mjd)) then  ! interval covered by observed and monthly predicted data
      flag_obs_m_pred = .true.
    else if((mjd_start >= sga_obs_start%mjd) .and. &
            (mjd_start <= sga_obs_end%mjd)   .and. &
            (mjd_end   >  sga_m_pred_end%mjd)) then  ! interval covered by observed, monthly predicted and monthly fit data
      flag_obs_fit = .true.
    else if((mjd_start >= sga_pred_start%mjd) .and. &
            (mjd_end   <= sga_pred_end%mjd)) then  ! interval covered by daily predicted data
      flag_pred_only = .true.
    else if((mjd_start >= sga_pred_start%mjd) .and. &
            (mjd_start <= sga_pred_end%mjd)   .and. &
            (mjd_end   <= sga_m_pred_end%mjd)) then  ! interval covered by daily and monthly predicted data
      flag_pred_m_pred = .true.
    else if((mjd_start >= sga_pred_start%mjd) .and. &
            (mjd_start <= sga_pred_end%mjd)   .and. &
            (mjd_end   >  sga_m_pred_end%mjd)) then  ! interval covered by daily and monthly predicted data, as well as monthly fit
      flag_pred_fit = .true.
    else if((mjd_start >  sga_pred_end%mjd) .and. &
            (mjd_end   <= sga_m_pred_end%mjd)) then  ! interval completely covered by monthly predictions
      flag_m_pred_only = .true.
    else if((mjd_start >  sga_pred_end%mjd) .and. &
            (mjd_start <= sga_m_pred_end%mjd)   .and. &
            (mjd_end   >  sga_m_pred_end%mjd)) then  ! interval covered by monthly predictions and monthly fit
      flag_m_pred_fit = .true.
    else if(mjd_start  > sga_m_pred_end%mjd) then  ! monthly fit forecast for complete interval required
      flag_fit_only = .true.
    end if

    if(flag_obs_only .or. flag_obs_pred  .or. flag_obs_m_pred  .or.  &
       flag_obs_fit  .or. flag_pred_only .or. flag_pred_m_pred .or.  &
       flag_pred_fit) then

      call mjd2gd(mjd_start,year,month,day,dtemp)
      if(hasFailed()) return
      write(ctemp10,'(i4,x,i2.2,x,i2.2)') year, month, day

      !** read until begin date has been reached
      do
        read(ich,'(A)',iostat=ios) cbuf
        if(ios /= 0) then
            call setNeptuneError(E_SGA_MISSING, FATAL)
            return
        else if(index(cbuf(1:10),ctemp10) /= 0) then
            exit
        end if
      end do

      !** set counter
      i = 1

      call mjd2gd(mjd_end,year,month,day,dtemp)
      if(hasFailed()) return
      !** set end date
      write(ctemp10,'(i4,x,i2.2,x,i2.2)') year, month, day

      !** start reading data
      do

        sga_data(i) = this%parseCSSIDataLine(cbuf)
        if(hasFailed()) return

        if(index(cbuf(1:10),ctemp10) /= 0) exit ! observed data only

        read(ich,'(A)',iostat=ios) cbuf

        if(index(cbuf,'END OBSERVED') /= 0) then  ! skip comment lines and then read on in daily predicted data
          do
            read(ich,'(A)',iostat=ios) cbuf
            if(index(cbuf,'BEGIN DAILY_PREDICTED') /= 0) then
              read(ich,'(A)',iostat=ios) cbuf ! next line right after tag
              exit
            end if
          end do
        else if(index(cbuf,'END DAILY_PREDICTED') /= 0) then  ! leave for monthly predicted and/or monthly fit
          i = i + 1
          exit
        end if
        i = i + 1
      end do
    end if

    !** now read monthly predicted and do interpolation to get daily values
    if(flag_obs_m_pred  .or. flag_obs_fit  .or.     &
       flag_pred_m_pred .or. flag_pred_fit .or.     &
       flag_m_pred_only .or. flag_m_pred_fit) then

        !** generate warning that geomagnetic activity forecast is not
        !   available anymore, setting to a level specified via input....
        write(ctemp3,'(i3)') this%ap_predict

        tempDate%mjd = sga_pred_end%mjd
        call mjd2gd(tempDate)

        cmess = "Geomagnetic activity (Ap) data not available for part of the specified time interval. Note that pred"//  &
                  "ictions are available only until "//date2string(tempDate)//". For dates beyond this epoch, "//  &
                  "geomagnetic activity is set according to specified input: "//ctemp3
        call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))

        do
            read(ich,'(A)',iostat=ios) cbuf
            if(index(cbuf,'BEGIN MONTHLY_PREDICTED') /= 0) exit
        end do

        do j = 2, m_pred_points+1

            read(ich,'(A)',iostat=ios) cbuf
            spwTemp = this%parseCSSIDataLine(cbuf)
            if(hasFailed()) return

            monthlies(j,1) = spwTemp%mjd
            monthlies(j,2) = spwTemp%f107
            monthlies(j,3) = spwTemp%f81ctr

        end do

        ! search for start date in monthlies array - it has to be in between two array elements
        if(i > 0) then
            ! daily data already available, so continue from the first monthly data
            firstMonthlyIndex = 1
        else
            ! start epoch is somewhere in the monthly predicted data
            do j = 2, m_pred_points + 1  ! as first index is last line in predicted daily data
                if(monthlies(j,1) > mjd_start) then
                    firstMonthlyIndex = j - 1
                    exit
                end if
            end do
            ! finally set the start date as the first entry (will be required in subsequent loop)
            sga_data(1)%mjd = mjd_start
            i = 1
        end if

        idx_start = firstMonthlyIndex
        idx_end   = firstMonthlyIndex + 3
        if(idx_end > size(monthlies,1)) then
            idx_end   = size(monthlies,1)
            idx_start = idx_end - 3
        end if

        !** now do interpolation (lagrange, cubic)
        do

            if(i > 1) then
                sga_data(i)%mjd = sga_data(i-1)%mjd + 1.d0
            end if

            !** find index in monthly predicted array
            if(sga_data(i)%mjd > monthlies(min(idx_start+2,size(monthlies,1)),1)) then

                idx_start = idx_start + 1
                idx_end   = idx_start + 3

                if(idx_end > size(monthlies,1)) then

                    idx_end   = size(monthlies,1)
                    idx_start = idx_end - 3

                end if

            end if

            !** F10.7 obs
            call lagrange_interpolation(monthlies(idx_start:idx_end,1), &
                                        monthlies(idx_start:idx_end,2), &
                                        4, &
                                        sga_data(i)%mjd, &
                                        sga_data(i)%f107)

            !** F10.7 81ctr
            call lagrange_interpolation(monthlies(idx_start:idx_end,1), &
                                        monthlies(idx_start:idx_end,3), &
                                        4, &
                                        sga_data(i)%mjd, &
                                        sga_data(i)%f81ctr)

            !** Ap
            do j=1,8
                sga_data(i)%ap(j) = this%ap_predict
            end do

            if((sga_data(i)%mjd > sga_m_pred_end%mjd) .or. sga_data(i)%mjd >= mjd_end) exit
            i = i + 1

        end do

    end if

    !** last step: do monthly fit, if required
    if(flag_obs_fit .or. flag_pred_fit .or. flag_m_pred_fit .or. flag_fit_only) then

        !** set counter
        if(flag_fit_only) then
            i = 1
        else
            i = i + 1
        end if

        do
            if(i == 1) then
                sga_data(i)%mjd = mjd_start
            else
                sga_data(i)%mjd = sga_data(i-1)%mjd + 1.d0
            end if

            !** use fit equation
            dt = sga_data(i)%mjd - c6

            sga_data(i)%f107   = c1 + c2*cos(c3*dt + c4*sin(c5*dt))
            sga_data(i)%f81ctr = sga_data(i)%f107

            do k=1,8
                sga_data(i)%ap(k) = this%ap_predict
            end do

            if(sga_data(i)%mjd >= mjd_end) exit
            i = i + 1

        end do

    end if

    if (allocated(monthlies)) deallocate(monthlies)

    !open(unit=99, file="test_sga.dat")
    !open(unit=100, file="test_mon.dat")

    !do j=1,i

      !write(99,200) sga_data(j)%mjd, sga_data(j)%f107, sga_data(j)%f81ctr, (sga_data(j)%ap(k), k=1,8)

    !end do

    !do j=1,m_pred_points

      !write(100,*) (monthlies(j,k), k=1,3)

    !end do

!200 format(f10.4,1x,f6.2,1x,f6.2,1x,8(i3,1x))
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine readCssiDataFile

!=============================================================================
!
!> @anchor      parseCssiDataLine
!!
!! @brief       Parsing a line of data from CSSI space weather file
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 18.07.2016 (initial design)</li>
!!              </ul>
!!
!! @param[in]   line         line to be parsed
!!
!! @details     It returns an sga_data type containing all information from one
!!              line of the CSSI space weather data file
!!
!!------------------------------------------------------------------------------------------------

  function parseCSSIDataLine(this,line) result(spwData)

    class(Atmosphere_class)                 :: this
    character(len=*), intent(in)            :: line
    type(sga_t) :: spwData

    character(len=*), parameter :: csubid = 'parseCSSIDataLine'
    integer  :: year, month, day, hour, j
    real(dp) :: dtemp

    if(isControlled()) then
        if(hasToReturn()) return
        call checkIn(csubid)
    end if

    read(line(1:4), *) year
    read(line(6:7), *) month
    read(line(9:10),*) day
    hour = getSolMagTime(year, month)
    if(hasFailed()) return

    call gd2mjd(year, month, day, hour, 0, 0.d0, dtemp, spwData%mjd)
    if(hasFailed()) return

    !** read 3hr Ap indices
    do j=0,28,4
        if(line(50+j:50+j) == ' ') then
            spwData%ap(j/4+1) = this%ap_predict
        else
            read(line(47+j:50+j),*) spwData%ap(j/4+1)
        end if
    end do

    read(line(119:124),*) spwData%f81ctr
    read(line(113:118),*) spwData%f107

    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end function

!=============================================================================
!
!> @anchor      readFapDayDataFile
!!
!> @brief       Reading daily solar and geomagnetic activity data from ESA data file
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 31.10.2013 (initial design)</li>
!!              </ul>
!!
!> @param[in]   ich         channel to read solmag data from
!> @param[in]   mjd_start   MJD of first entry SGA data array
!> @param[in]   mjd_end     MJD of last entry in SGA data array
!!
!> @details     This routine parses the ESA daily data file. It finds the required
!!              data which is required for the propagation and stores
!!              it within the sga data array.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Check update status </li>
!!                <li> Read data           </li>
!!                <li> If propagation span is not covered, data from fap_mon is taken (by calling readFapMonDataFile)</li>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  subroutine readFapDayDataFile(  this,      &
                                  ich,       &  ! <-- INT     input channel to read data from
                                  mjd_start, &  ! <-- DBL     MJD of first SGA data entry
                                  mjd_end    &  ! <-- DBL     MJD of last SGA entry
                               )

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    integer,  intent(in)                :: ich
    real(dp), intent(in)                :: mjd_start
    real(dp), intent(in)                :: mjd_end
    !-------------------------------------------

    character(len=*), parameter :: csubid = "readFapDayDataFile"
    character(len=255) :: cbuf          ! character buffer
    character(len=255) :: cbuf_save     ! saving previous character buffer
    character(len=3)   :: cm            ! month title
    character(len=255) :: cmess         ! message string
    character(len=10)  :: ctemp10       ! auxiliary
    character(len=2)   :: kp_temp       ! planetary index as read from file

    integer :: apDaily                  ! daily Ap coefficient
    integer :: day                      ! last update day from fap_day data file
    integer :: ios                      ! I/O status
    integer :: i,j,k                    ! loop counter
    integer :: hour
    integer :: month
    integer :: year                     ! last update year from fap_day data file

    logical :: first = .true.           ! required to find first data set in data file
    logical :: flag_daily_only = .false.         ! only data from fap_day is required
    logical :: flag_daily_monthly  = .false.       ! data from fap_day and fap_mon required (monthlies)
    logical :: flag_monthly_only  = .false.        ! data from fap_mon only required

    real(dp) :: dtemp                   ! auxiliary
    real(dp), dimension(8) :: kp        ! 3hr planetary index for given day

    type(time_t) :: sga_date_start        ! first data set in fap_day
    type(time_t) :: sga_date_end          ! last data set in fap_day
    type(time_t) :: sga_time_tag          ! time tag of SGA data file
    type(time_t) :: system_time           ! system time

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if


    rewind(ich)

    !===================================================================================================
    !
    !   Find 'last update' key word in daily file and check whether update is
    !   required or not
    !
    !--------------------------------------------------------------------------
    do

      read(ich,'(a255)',iostat=ios) cbuf

      if(ios /= 0) then
        call setNeptuneError(E_SGA_TIME_TAG, WARNING, (/"N/A"/))
      else if(index(cbuf,'last update') /= 0) then
        exit
      end if

    end do

    !** read date string from buffer
    read(cbuf(index(cbuf,'date:')+ 6:index(cbuf,'date:')+ 7), *, iostat=ios) sga_time_tag%day
    read(cbuf(index(cbuf,'date:')+ 9:index(cbuf,'date:')+11), *, iostat=ios) cm
    read(cbuf(index(cbuf,'date:')+13:index(cbuf,'date:')+16), *, iostat=ios) sga_time_tag%year

    if(ios /= 0) then ! must be unsupported version...
      call setNeptuneError(E_SGA_UNSUPPORTED_VERSION, FATAL, (/"N/A"/))
    else

      call getMonthNumber(cm,sga_time_tag%month)

      sga_time_tag%hour   = 0
      sga_time_tag%minute = 0
      sga_time_tag%second = 0.d0

      call gd2mjd(sga_time_tag)

      ! save global variable
      this%solmagFileDateDaily = sga_time_tag

      !** compare to system time
      system_time = getDateTimeNow()
      if(hasFailed()) return

      if((system_time%mjd - sga_time_tag%mjd) > 60.d0) then
        cmess = "Solar and geomagnetic activity data file has to be updated."
        call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
      end if

      cmess = "Using solar and geomagnetic activity data file from '"//date2string(sga_time_tag)//"'."
      call setNeptuneError(E_SPECIAL, REMARK, (/cmess/))

    end if
    !-------------------------------------------------------------------------------------------------------

    rewind(ich)

    !=====================================================================
    !
    !   Find epoch of first and last data set
    !
    !------------------------------------------
    do

      read(ich,'(A)',iostat=ios) cbuf

      if(index(cbuf,'#') == 0 .and. first) then

        read(cbuf(7:10), *) sga_date_start%year
        read(cbuf(4:5),  *) sga_date_start%month
        read(cbuf(1:2),  *) sga_date_start%day

        sga_date_start%hour = getSolMagTime(sga_date_start%year, sga_date_start%month)
        if(hasFailed()) return

        sga_date_start%minute = 0
        sga_date_start%second = 0.d0

        call gd2mjd(sga_date_start)

        first = .false.

      else if(ios < 0) then ! end of file found, get date of last data set from previous cbuf (cbuf_save)

        read(cbuf_save(7:10), *) sga_date_end%year
        read(cbuf_save(4:5),  *) sga_date_end%month
        read(cbuf_save(1:2),  *) sga_date_end%day

        sga_date_end%hour = getSolMagTime(sga_date_end%year, sga_date_end%month)
        if(hasFailed()) return

        sga_date_end%minute = 0
        sga_date_end%second = 0.d0

        call gd2mjd(sga_date_end)
        exit

      end if

      cbuf_save  = cbuf ! required to find last data set

    end do

    rewind(ich)

    !=======================================================================
    !
    !   Now check for begin and end date within the data
    !
    !-------------------------------------------
    if(sga_date_start%mjd > mjd_start) then

      call setNeptuneError(E_SGA_MISSING, FATAL)
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return

    else if((mjd_start >= sga_date_start%mjd) .and. (mjd_end <= sga_date_end%mjd)) then ! interval completely covered by fap_day data
      flag_daily_only    = .true.
    else if((mjd_start >= sga_date_start%mjd) .and. (mjd_start <= sga_date_end%mjd) .and. (mjd_end > sga_date_end%mjd)) then  ! interval covered by daily and monthly data
      flag_daily_monthly = .true.
    else
      flag_monthly_only  = .true.
    end if

    !** set counter
    i = 0

    if(flag_daily_only .or. flag_daily_monthly) then

      call mjd2gd(mjd_start,year,month,day,dtemp)
      write(ctemp10,"(i2.2,'/',i2.2,'/',i4.4)") day, month, year

      !** read until begin date has been reached
      do

        read(ich,'(A)',iostat=ios) cbuf
        if(index(cbuf(1:10),ctemp10).ne.0) exit

      end do

      call mjd2gd(mjd_end,year,month,day,dtemp)
      !** set end date
      write(ctemp10,"(i2.2,'/',i2.2,'/',i4.4)") day, month, year

      !** start reading data
      do

        !** get date...
        read(cbuf(7:10),*) year
        read(cbuf(4:5), *) month
        read(cbuf(1:2), *) day

        i = i + 1

        hour = getSolMagTime(year, month)

        call gd2mjd(year, month, day, hour, 0, 0.d0, dtemp, sga_data(i)%mjd)

        !** read daily Ap
        read(cbuf(24:26),*) apDaily

        !** read 3hr Kp indices
        do j=1,15,2

          read(cbuf(27+j:29+j),*) kp_temp

          read(kp_temp(1:1),*) dtemp

          select case(kp_temp(2:2))
            case('+')
              dtemp = dtemp + 1/3.d0
            case('-')
              dtemp = dtemp - 1/3.d0
          end select

          kp((j+1)/2) = dtemp

        end do

        !** convert kp to ap, if kp > 0 (else, there is no actual value '0-0-...' - use same value for all Ap's from Ap column)
        do k=1,8
          if(kp(k) > 0.d0) then
            sga_data(i)%ap(k) = nint(this%kp2ap(kp(k)))
          else
            sga_data(i)%ap(k) = apDaily
          end if
        end do

        !** F10.7 obs, 81ctr - note that in fap_day the average value is last-81 - if NRLMSISE-00 is used
        !   as atmosphere model, than a correction has to be performed (shift) later on to get ctr-81!
        read(cbuf(16:18),*) sga_data(i)%f81ctr

        !** F10.7 obs
        read(cbuf(12:14),*) sga_data(i)%f107

        !** check whether end date has been reached (if only daily data used)
        if(flag_daily_only) then
          if(index(cbuf(1:10),ctemp10) /= 0) exit ! observed data only
        end if

        !** read next line
        read(ich,'(A)',iostat=ios) cbuf

        !** check for end of file, if also monthly data has to be read
        if(ios < 0) then

          !** The next call also accounts for filling the remaining array elements, even if no
          !   monthly data is available - the routine then falls back to a constant value for
          !   solar and geomagnetic activity
          this%monthlyDataRequired = .true.
          call this%readFapMonDataFile(i, mjd_start, mjd_end)
          exit

        end if

      end do

    else if(flag_monthly_only) then ! read monthly data file with i being zeroth index (= 0)
      !write(*,*) "reading FAP MON", i, mjd_start, mjd_end
      this%monthlyDataRequired = .true.
      call this%readFapMonDataFile(i, mjd_start, mjd_end)                       ! called with i = 0

    end if

    !=================================================================
    !
    ! Shift F10.7 last 81-day average by 40 days, in order to have
    ! centered 81-day for NRLMSISE-00
    !
    !-------------------------------------------------
    if(this%nmodel == MSIS2000) then

      sga_data(:)%f81ctr = eoshift(sga_data(:)%f81ctr, 40)

    end if

    !open(unit=99, file="test_sga.dat")
    !open(unit=100, file="test_mon.dat")

    !do j=1,size(sga_data)
!
!      write(99,200) sga_data(j)%mjd, sga_data(j)%f107, sga_data(j)%f81ctr, (sga_data(j)%ap(k), k=1,8)
!
!    end do
!
!    close(99)
!    stop

    !do j=1,m_pred_points

      !write(100,*) (monthlies(j,k), k=1,3)

    !end do

!200 format(f10.4,1x,f6.2,1x,f6.2,1x,8(i3,1x))

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine readFapDayDataFile

!=============================================================================
!
!> @anchor      readFapMonDataFile
!!
!> @brief       Reading monthly solar and geomagnetic activity data from ESA data file
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 02.11.2013 (initial design)</li>
!!              </ul>
!!
!> @param[in]   idx         solar and geomagnetic activity data array index where monthly data sets in
!> @param[in]   mjd_start   MJD of first SGA data entry
!> @param[in]   mjd_end     MJD of last SGA entry
!!
!> @details     This routine parses the ESA monthly data file. It finds the required
!!              data which is required for the propagation and stores
!!              it within the sga data array. The daily data is obtained by interpolation
!!              using four data points. If the monthly data is not enough to cover the
!!              propagation span, the remaining array elements are filled with a constant
!!              solar and geomagnetic activity, as specified via the input.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Read data           </li>
!!                <li> If propagation span is not covered, data a constant activity according to input file is assumed.</li>
!!                <li> Perform Lagrange interpolation to get daily data when monthly data is used </li>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  subroutine readFapMonDataFile(  this,      &
                                  idx,       &  ! <-- INT     current array index in sga_data array
                                  mjd_start, &  ! <-- DBL     MJD of first SGA data entry
                                  mjd_end    &  ! <-- DBL     MJD of last SGA entry
                               )

    !** interface
    !-------------------------------------------
    class(Atmosphere_class)             :: this
    integer,  intent(in)                :: idx
    real(dp), intent(in)                :: mjd_start
    real(dp), intent(in)                :: mjd_end
    !-------------------------------------------

    character(len=*), parameter :: csubid = "readFapMonDataFile"
    character(len=255) :: cbuf      ! character buffer
    character(len=255) :: cbuf_save ! character buffer
    character(len=255) :: cmess     ! message string
    character(len=3)   :: cm        ! month name of time tag
    character(len=255) :: fileName  ! name (+path) of monthly data file

    integer :: ich                ! input channel
    integer :: ios                ! I/O status
    integer :: i,j                ! loop counter
    integer :: idx_start          ! start index
    integer :: idx_end            ! end index
    integer :: m, y               ! month and year

    logical :: flag_monthly = .true.    ! indicating that monthly data is available
    real(dp) :: dtemp
    real(dp) :: mjd_temp
    real(dp) :: mjd_start_monthly ! MJD of first data set
    real(dp) :: mjd_end_monthly   ! MJD of last data set
    real(dp), dimension(:,:), allocatable :: monthlyTemp    ! holding information from monthly data file

    type(time_t) :: sga_time_tag  ! time tag of data file
    type(time_t) :: system_time

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** open monthly data output file
    !-------------------------------------------------------------------------------------
    fileName = trim(adjustl(this%dataPath))//cdelimit//trim(adjustl(this%sgaMonthlyDataFile))
    ich = openFile(fileName, SEQUENTIAL, IN_FORMATTED)

    !===================================================================================================
    !
    !   Find 'Last Update' keyword in monthly file and check whether update is
    !   required or not
    !
    !--------------------------------------------------------------------------
    do

      read(ich,'(a255)',iostat=ios) cbuf

      if(ios /= 0) then
        call setNeptuneError(E_SGA_TIME_TAG, WARNING, (/fileName/))
      else if(index(cbuf,'Last Update') /= 0) then
        exit
      end if

    end do

    !** read date string from buffer
    read(cbuf(index(cbuf,'date:')+ 6:index(cbuf,'date:')+ 9), *, iostat=ios) sga_time_tag%year
    read(cbuf(index(cbuf,'date:')+11:index(cbuf,'date:')+13), *, iostat=ios) cm
    read(cbuf(index(cbuf,'date:')+15:index(cbuf,'date:')+16), *, iostat=ios) sga_time_tag%day

    if(ios /= 0) then ! must be unsupported version...
      call setNeptuneError(E_SGA_UNSUPPORTED_VERSION, FATAL, (/"N/A"/))
    else
      call getMonthNumber(cm,sga_time_tag%month)

      sga_time_tag%hour   = 0
      sga_time_tag%minute = 0
      sga_time_tag%second = 0.d0

      call gd2mjd(sga_time_tag)

      ! save global variable
      this%solmagFileDateMonthly = sga_time_tag

      !** compare to system time
      system_time = getDateTimeNow()
      if(hasFailed()) return

      if((system_time%mjd - sga_time_tag%mjd) > 360.d0) then
        cmess = "Monthly solar and geomagnetic activity data file has to be updated."
        call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
      end if

      cmess = "Using monthly solar and geomagnetic activity data file from '"//date2string(sga_time_tag)//"'."
      call setNeptuneError(E_SPECIAL, REMARK, (/cmess/))

    end if
    !-------------------------------------------------------------------------------------------------------

    rewind(ich)

    !** find start date of MONTHLY data (first uncommented line)
    !-----------------------------------------------------------------
    do

      read(ich,'(A)',iostat=ios) cbuf

      if(index(cbuf,'mm_yyyy') /= 0) then
        read(ich,'(A)',iostat=ios) cbuf
        exit
      end if

    end do

    !** convert date to MJD (1950.0)
    read(cbuf,*) m, y

    call gd2mjd (y, m, 15, 12, 0, 0.d0, dtemp, mjd_start_monthly)

    !** find end date of MONTHLY data (last uncommented line)
    !---------------------------------------------------------------
    do

      read(ich,'(A)',iostat=ios) cbuf

      if(index(cbuf,'#') /= 0 .or. ios /= 0) exit
      cbuf_save = cbuf

    end do

    !** convert date to MJD (1950.0)
    read(cbuf_save,*) m, y

    call gd2mjd(y, m, 15, 12, 0, 0.d0, dtemp, mjd_end_monthly)
    !--------------------------------------------------------------

    !======================================================================
    !
    ! Now start reading data file
    !
    !---------------------------------------------------
    rewind(ich)

    do    ! find first entry
      read(ich,'(A)',iostat=ios) cbuf
      if(index(cbuf,'mm_yyyy') /= 0) exit
    end do

    if(idx == 0) then  ! no daily data available

      do    ! find first data point which is prior to start date...

        read(ich,'(a)',iostat=ios) cbuf
        read(cbuf,*) m, y
        call gd2mjd(y,m,15,12,0,0.d0,dtemp,mjd_temp)
        if(mjd_temp > mjd_start) exit
        cbuf_save = cbuf  ! storing character buffer which corresponds to earlier data set, which
                          ! is actually what is required after the loop has been left
      end do

      i = 1
      read(cbuf_save,*) m, y
      call gd2mjd(y,m,15,12,0,0.d0,dtemp,sga_data(i)%mjd)

      read(cbuf(31:36),*) sga_data(i)%f81ctr
      read(cbuf(67:71),*) dtemp

      sga_data(i)%ap   = nint(dtemp)
      sga_data(i)%f107 = sga_data(i)%f81ctr

    else

      do    ! find first data point which is after available daily data...

        read(ich,'(a)',iostat=ios) cbuf

        if(ios /= 0) then ! no monthly data available!
          flag_monthly = .false.
          exit
        end if

        read(cbuf,*) m, y
        call gd2mjd(y,m,15,12,0,0.d0,dtemp,mjd_temp)
        if(mjd_temp > sga_data(idx)%mjd) exit

      end do

      i = idx

    end if

    if(flag_monthly) then

      !** allocate temporary array to hold monthly data
      allocate(monthlyTemp(int((min(mjd_end,mjd_end_monthly) - sga_data(i)%mjd)/30.d0)+4,3)) ! '+4' for interpolation purposes

      j = 1

      !** now read data into temporary array, passing values of last daily data to first array element
      monthlyTemp(j,1) = sga_data(i)%mjd
      monthlyTemp(j,2) = sga_data(i)%f81ctr
      monthlyTemp(j,3) = mean(sga_data(i)%ap)

      do

        j = j + 1
        read(cbuf(1:3),*) m
        read(cbuf(5:8),*) y

        call gd2mjd(y,m,15,12,0,0.d0, dtemp, monthlyTemp(j,1))

        read(cbuf(31:36),*) monthlyTemp(j,2)
        read(cbuf(67:71),*) monthlyTemp(j,3)

        read(ich,'(a)',iostat=ios) cbuf
        if(ios /= 0 .or. index(cbuf,'#') /= 0 .or. j == size(monthlyTemp,1)) exit

      end do

      !** now do interpolation (lagrange, cubic)

      idx_start = 1
      idx_end   = 4

      do

        i = i + 1
        sga_data(i)%mjd = sga_data(i-1)%mjd + 1.d0

        !** find index in monthly predicted array
        if(sga_data(i)%mjd > monthlyTemp(min(idx_start+2, size(monthlyTemp,1)),1)) then

          idx_start = idx_start + 1
          idx_end   = idx_start + 3

          if(idx_end > size(monthlyTemp,1)) then

            idx_end   = size(monthlyTemp,1)
            idx_start = idx_end - 3

          end if

        end if

        !** F10.7 obs
        call lagrange_interpolation(monthlyTemp(idx_start:idx_end,1),  &
                                    monthlyTemp(idx_start:idx_end,2),  &
                                    4,                                 &
                                    sga_data(i)%mjd,              &
                                    sga_data(i)%f107)

        !** Ap
        call lagrange_interpolation(monthlyTemp(idx_start:idx_end,1),  &
                                    monthlyTemp(idx_start:idx_end,3),  &
                                    4,                                 &
                                    sga_data(i)%mjd,              &
                                    dtemp)

        sga_data(i)%ap = nint(dtemp)

        !sga_data(i)%f107 = this%sol_predict
        !sga_data(i)%ap = this%ap_predict

        if(sga_data(i)%mjd > mjd_end_monthly .or. i == size(sga_data)) exit

      end do

    end if

    !==================================================================
    !
    ! Fill remaining array elements with constant activity data
    !
    !-----------------------------------------------------------
    if(i < size(sga_data)) then

      do

        i = i + 1
        sga_data(i)%mjd  = sga_data(i-1)%mjd + 1.d0
        sga_data(i)%f107 = this%sol_predict
        sga_data(i)%ap   = this%ap_predict

        if(i == size(sga_data)) exit

      end do

    end if

    !==================================================================
    !
    ! Process F10.7 last 81-day averages
    !
    !----------------------------------------------------
    do i = max(82,idx+1), size(sga_data)

      sga_data(i)%f81ctr = mean(sga_data(i-81:i-1)%f107)

    end do

    ich = closeFile(ich)

    if (allocated(monthlyTemp)) deallocate(monthlyTemp)

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

 end subroutine readFapMonDataFile

end module atmosphere
