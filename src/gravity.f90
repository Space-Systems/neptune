!===============================================================================================
!
!> @anchor      gravity
!!
!> @brief       Geopotential modeling
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!> @author      Daniel Lück (DLU)
!!
!> @date        <ul>
!!                <li>VB:  02.10.2012 (initial design)                        </li>
!!                <li>VB:  31.05.2013 (code optimization)                     </li>
!!                <li>VB:  16.07.2013 (modifications to support tides)        </li>
!!                <li>VB:  29.07.2013 (added function 'getLegendrePolynomial')</li>
!!                <li>VB:  13.01.2014 (added interface to set harmonic coefficients)</li>
!!                <li>VB:  28.04.2014 (added module parameters: 'MODEL_EIGEN_GL04C', ... </li>
!!                <li>VB:  10.06.2014 (changed handling of Earth's radius, implemented new function 'getEarthGeopotentialRadius')</li>
!!                <li>CHK: 13.11.2015 (updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>DLU: 13.11.2015 (added calculation of legendre polynomials to gravity covariance)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for Earth's
!!              geopotential modeling, specifically in the context of numerical
!!              integration of a satellite trajectory.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      derivatives
!!
!!------------------------------------------------------------------------------------------------
module gravity

  use slam_astro,             only: EGM96, EGM08, EIGEN_GL04C, setEarthGravity, setEarthGeopotentialRadius, &
                                    getEarthGeopotentialRadius, getEarthGravity
  use slam_astro_conversions, only: getRadiusLatLon
  use neptune_error_handling, only: E_GRAVITY_COEFFICIENTS, E_GEO_MAX_DEGREE, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, &
                                    E_SPECIAL, E_UNKNOWN_PARAMETER, E_OUT_OF_BOUNDS, checkIn, checkOut
  use slam_io,                only: openFile, closeFile, SEQUENTIAL, IN_FORMATTED, LOG_AND_STDOUT, cdelimit, message
  use slam_math,              only: UNDEFINED, identity_matrix, rad2deg, redang, eps9, mag, outerproduct, uvec1, uvec2, uvec3
  use slam_reduction_class,   only: Reduction_type
  use slam_time,              only: time_t, gd2jd
  use slam_types,             only: dp

  implicit none

  private

  character(len=*), parameter :: DEFAULT_DATA_PATH = 'data' !< default path (relative to executable) where geopotential coefficients are found

  integer, parameter, public :: COEFF_C  = 1
  integer, parameter, public :: COEFF_DC = 10
  integer, parameter, public :: COEFF_S  = 2
  integer, parameter, public :: COEFF_DS = 20

  integer, parameter, public :: MODEL_EGM96       = 1
  integer, parameter, public :: MODEL_EGM2008     = 2
  integer, parameter, public :: MODEL_EIGEN_GL04C = 3
  integer, parameter         :: DEFAULT_MODEL     = 3 !< the EIGEN-GL04C model is used as default (according to ECSS)

  integer, parameter :: n_supported_models = 3    !< number of supported geopotential models

  ! model names according to the integer numbers associated through astro module
  character(len=*), dimension(n_supported_models), parameter :: modelName = (/"EGM96      ", &
                                                                              "EGM2008    ", &
                                                                              "EIGEN-GL04C"/)
  ! model names according to the integer numbers associated through astro module
  character(len=*), dimension(n_supported_models), parameter :: dataFileName = (/"egm96.dat      ", &
                                                                                 "egm2008.dat    ", &
                                                                                 "eigen_gl04c.dat"/)
  ! data file formats for each supported model
  character(len=37), dimension(n_supported_models,2), parameter :: cDataFormat = reshape((/"(a1,2(i4),2(e20.12e2),2(e16.8e2),a9)",  &
                                                                                          "(a1,2(i5),2(e25.15e2),2(e20.10e2))", &
                                                                                          "(a4,2(i5),2(e19.12e2),2(e11.4e2),a9),", &
                                                                                          "(a1,2(i4),2(e20.12e2),2(e16.8e2))",  &
                                                                                          "(a1,2(i5),2(e25.15e2),2(e20.10e2))", &
                                                                                          "(a3,2(i5),2(e19.12e2),2(e11.4e2))"/), (/3,2/))

  integer, parameter :: maxdegree = 85                                          ! maximum degree of geopotential
  integer, parameter :: maxDistinctHarmonics = 80                               ! maximum number of individual harmonics to be analysed

  ! The gravity model class (type defintion)
  type, public :: Gravity_class

      integer :: nmodel                                                       !< geopotential model to be used

      logical :: analyseDistinctHarmonics                                     ! tells if only distinct spherical harmonics are to be used in the computation
      logical :: areAvailableDistinctHarmonics                                ! Is set to .true. if distinct harmonics are available for analysis
      logical :: coeffInitialized                                             ! tells if coefficients have been already read from input file
      logical :: supportingTides                                              ! providing support for tides module

      integer :: covDegree                                                    ! geopotential degree for covariance matrix
      integer :: degree                                                       ! geopotential degree used in propagation
      integer :: degree_data                                                  ! geopotential degree for which data has been read

      integer :: numDistinctHarmonics                                         ! number of distinct harmonics to be analysed
      integer, dimension(maxDistinctHarmonics) :: distinctHarmonics           ! containing for each array element a combination of degree and order to be analysed, e.g.:
                                                                              ! arr(1) = 20; means that C_{2,0} and S_{2,0} will be considered only
                                                                              ! arr(2) = 31; means that besides the above harmonics, also C_{3,1} and S_{3,1} are considered, etc.

      !** variables which are potentially shared between state vector geop. and covariance geop.
      !--------------------------------------------------------------------------------------------
      real(dp), dimension(:,:), allocatable :: coeffC                         ! gravity potential coefficients 'C'
      real(dp), dimension(:,:), allocatable :: coeffS                         ! gravity potential coefficients 'S'
      real(dp), dimension(:),   allocatable :: costerm                        ! cos(m*lambda)
      real(dp)                              :: current_alt                    ! current altitude in km
      real(dp)                              :: dudlambda                      ! dU/d(lambda)
      real(dp)                              :: dudphi                         ! dU/d(phi_gc)
      real(dp)                              :: dudr                           ! dU/d(rabs)
      real(dp)                              :: lambda                         ! current geocentric latitude
      real(dp), dimension(:,:), allocatable :: lp                             ! legendre polynomials
      real(dp)                              :: mu                             ! Earth's gravity constant
      real(dp)                              :: oorabs                         ! 1/rabs
      real(dp)                              :: oorabs2                        ! 1/rabs2
      real(dp)                              :: oorabs3                        ! 1/rabs3
      real(dp)                              :: rabs                           ! magnitude of radius vector
      real(dp)                              :: rabs2                          ! squared radius magnitude
      real(dp)                              :: rabs3                          ! cubed radius magnitude
      real(dp)                              :: phi_gc                         ! current geocentric longitude
      real(dp)                              :: rekm                           ! Earth's radius in km
      real(dp)                              :: r1r2                           ! radius(1)**2 + radius(2)**2
      real(dp), dimension(:),   allocatable :: rrfac                          ! temporary
      real(dp), dimension(:,:), allocatable :: sigmaC                         ! error in coefficient C
      real(dp), dimension(:,:), allocatable :: sigmaS                         ! error in coefficient S
      real(dp), dimension(:,:), allocatable :: sigmaCdenorm                   ! error in coefficient C (de-normalized)
      real(dp), dimension(:,:), allocatable :: sigmaSdenorm                   ! error in coefficient S (de-normalized)
      real(dp), dimension(:),   allocatable :: sinterm                        ! sin(m*lambda)
      real(dp)                              :: sqrt_r1r2                      ! sqrt(r1r2)
      real(dp)                              :: tanphi                         ! tan(phi_gc)

      real(dp), dimension(3)                :: current_radius                 ! saving the state for which accelerations are computed in ITRF

      logical, dimension(:,:), allocatable :: useCoeffC                       !< specific coefficients only used, if value set to .true.
      logical, dimension(:,:), allocatable :: useCoeffS                       !< specific coefficients only used, if value set to .true.
      !----------------------------------------------------------------------------------------------

      real(dp), dimension(:),   allocatable :: dCdt                           ! time derivatives of C coefficients
      real(dp), dimension(:),   allocatable :: dSdt                           ! time derivatives of S coefficients
      !real(dp), dimension(2:6,0:6,1:5)      :: tidepars_lm                   ! parameters which can be used by tide module, f(l,m)
      type(time_t), dimension(:), allocatable :: cepoch                       ! epoch of coefficients C
                                                                              ! this state is, e.g. required by the tides module
  contains

      procedure :: initGravityPotential

      !** getter
      procedure :: getGravityAcceleration
      procedure :: getGravityCovariance
      procedure :: getCoefficientSigma
      procedure :: getGeoCovDegree
      procedure :: getGeopotentialDegree
      procedure :: getGeopotentialModelName
      procedure :: getHarmonicCoefficient
      procedure :: getCurrentLatLon
      procedure :: getCurrentRadius
      procedure :: getLegendrePolynomial
      procedure :: getDistinctHarmonicsString
      procedure :: getGeoModelFileName
      procedure,private :: getPotentialDerivatives

      procedure :: usedDistinctHarmonics
      procedure :: isInitializedGeopotential

      !** setter
      procedure :: setCurrentAltitude
      procedure :: setDistinctHarmonics
      procedure :: setDistinctHarmonicsArray
      procedure :: setTideSupport
      procedure :: setGeoDegree
      procedure :: setGeoCovDegree
      procedure :: setHarmonicCoefficient
      procedure :: toggleHarmonicCoefficient
      procedure :: destroy

  end type Gravity_class

      ! Constructor
  interface Gravity_class
      module procedure constructor
  end interface Gravity_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !>  @author     Christopher Kebschull
    !>  @date       <ul>
    !!                  <li>ChK: 03.01.2018 (initial implementation)</li>
    !!              </ul>
    !>  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Gravity_class) function constructor()
        constructor%analyseDistinctHarmonics      = .false.                     ! tells if only distinct spherical harmonics are to be used in the computation
        constructor%areAvailableDistinctHarmonics = .false.                     ! Is set to .true. if distinct harmonics are available for analysis
        constructor%coeffInitialized              = .false.                     ! tells if coefficients have been already read from input file
        constructor%supportingTides               = .false.                     ! providing support for tides module

        constructor%covDegree   = -1                                            ! geopotential degree for covariance matrix
        constructor%degree      = -1                                            ! geopotential degree used in propagation
        constructor%degree_data = -1                                            ! geopotential degree for which data has been read
    end function constructor

    !===========================================================================
    !!
    !>  @anchor     destroy
    !!
    !>  @brief      Destroys all that needs destruction
    !>  @author     Christopher Kebschull
    !!
    !!
    !>  @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Gravity_class)    :: this

        if (allocated(this%coeffC)) deallocate(this%coeffC)
        if (allocated(this%coeffS)) deallocate(this%coeffS)
        if (allocated(this%costerm)) deallocate(this%costerm)
        if (allocated(this%lp)) deallocate(this%lp)
        if (allocated(this%rrfac)) deallocate(this%rrfac)
        if (allocated(this%sigmaC)) deallocate(this%sigmaC)
        if (allocated(this%sigmaS)) deallocate(this%sigmaS)
        if (allocated(this%sigmaCdenorm)) deallocate(this%sigmaCdenorm)
        if (allocated(this%sigmaSdenorm)) deallocate(this%sigmaSdenorm)
        if (allocated(this%sinterm)) deallocate(this%sinterm)

        if (allocated(this%useCoeffC)) deallocate(this%useCoeffC)
        if (allocated(this%useCoeffS)) deallocate(this%useCoeffS)
        !----------------------------------------------------------------------------------------------

        if (allocated(this%dCdt)) deallocate(this%dCdt)
        if (allocated(this%dSdt)) deallocate(this%dSdt)
        if (allocated(this%cepoch)) deallocate(this%cepoch)
    end subroutine destroy

!--------------------------------------------------------------------------------------------------
!
!> @anchor      getGeoModelFileName
!!
!> @brief       Returns the name of the file the Stokes coefficients are read from
!!
!> @return      file name of the geo model
!!
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> VB: 21.10.2016 (initial design)    </li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  character(len=len(dataFileName(1))) function getGeoModelFileName(this) result(fname)

    class(Gravity_class)    :: this

    if(.not. this%coeffInitialized) then
        fname = ''
    else
        fname = dataFileName(this%nmodel)
    end if
    return
  end function
!--------------------------------------------------------------------------------------------------
!
!> @anchor      initGravityPotential
!!
!> @brief       Initialization of geopotential module
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 02.10.2012 (initial design)    </li>
!!                <li> 31.05.2013 (code optimization) </li>
!!                <li> 19.08.2013 (added EGM96)       </li>
!!                <li> 04.02.2015 (now accounts for distinct harmonics) </li>
!!                <li> 17.01.2016 (path of data file and geo model are now optional: if not provided, default values are used)</li>
!!              </ul>
!!
!> @param[in]   cpath       (optional) Path to read geopotential data from. If not provided, default path is used
!> @param[in]   imodel      (optional) model to be used, e.g. EIGEN_GL04C. If not provided, default model is used
!!
!> @details     This routine initializes the gravity module. Global parameters are set and
!!              geopotential data is read from the data file according to the selected model
!!              into the specific arrays.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Open data file according to selected geopotential (default: EIGEN-GL04C) </li>
!!                <li> Find Earth's geocentric radius and gravity constant </li>
!!                <li> Allocate coefficients arrays </li>
!!                <li> Read data from file into arrays </li>
!!                <li> De-normalize coefficients </li>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  subroutine initGravityPotential( this,   &
                                   cpath,  & ! <-- CHR directory path containing gravity potential data
                                   imodel  & ! <-- INT geopotential model to be used
                                 )

    !** interface
    !------------------------------------------
    class(Gravity_class)                    :: this
    character(len=*), optional, intent(in)  :: cpath
    integer,          optional, intent(in)  :: imodel
    !------------------------------------------

    character(len=255) :: cbuf      ! character buffer
    character(len=255) :: cmess     ! message string
    character(len=*), parameter  :: csubid = "initGravityPotential"
    character :: ctemp              ! temporary character
    character(len=4)   :: ctemp4    ! temporary string
    character(len=255) :: fileName
    character(len=9)   :: tempEpoch ! temporary epoch
    integer :: i                    ! loop counter
    integer :: ich                  ! input channel
    integer, parameter :: imax = 1000 ! maximum number of lines to search for end_of_head in input file
    integer :: ios                ! I/O status
    integer :: k
    integer :: l                  ! degree of gravity potential
    integer :: lm                 ! combined number of l and m: lm = 10*l + m
    integer :: m                  ! order of gravity potential
    integer :: maxd               ! maximum degree found for distinct harmonics

    logical :: flag_mu            ! flag for mu being found
    logical :: flag_rekm          ! flag for rekm being found

    real(dp)  :: fac                ! multiplication factor
    real(dp), dimension(:), allocatable :: factorial    ! factorials of n (n!)
    real(dp)  :: tempC              ! temporary coefficient C
    real(dp)  :: tempdCdt           ! temporary derivative dCdt
    real(dp)  :: tempdSdt           ! temporary derivative dSdt
    real(dp)  :: tempS              ! temporary coefficient S
    real(dp)  :: tempSigmaC         ! temporary Sigma C
    real(dp)  :: tempSigmaS         ! temporary Sigma S

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%degree < 0) then   ! has not been set yet probably?!
      call this%setGeoDegree(this%degree)
    else if(this%analyseDistinctHarmonics .and. this%areAvailableDistinctHarmonics) then ! set degree according to passed harmonics

      !** find maximum degree value in array
      maxd = 0

      do i = 1, this%numDistinctHarmonics
        if(this%distinctHarmonics(i)/10 > maxd) then
          maxd = this%distinctHarmonics(i)/10
        end if
      end do

      call this%setGeoDegree(maxd)

    end if

    if(this%coeffInitialized .or. (this%degree < 2)) then   ! return if already initialized, or no geopotential is required.

      if(this%degree < 2) this%nmodel = 0

      if(isControlled()) then
        call checkOut(csubid)
      end if
      return

    end if

    this%coeffInitialized = .false.    ! as a new initialization is started...

    !=====================================================================================
    !
    ! Decide on which model to use (default: EIGEN-GL04C)
    !
    !---------------------------------------------------------

    !** check imodel validity
    if(present(imodel)) then
      if(imodel == EGM96 .or. imodel == EGM08 .or. imodel == EIGEN_GL04C) then
        this%nmodel = imodel
      else
        this%nmodel = DEFAULT_MODEL
      end if
    else
      this%nmodel = DEFAULT_MODEL
    end if

    !** open data file
    if(present(cpath)) then
      fileName = trim(adjustl(cpath))//cdelimit//trim(dataFileName(this%nmodel))
    else
      fileName = trim(adjustl(DEFAULT_DATA_PATH))//cdelimit//trim(dataFileName(this%nmodel))
    end if

    ich = openFile(fileName, SEQUENTIAL, IN_FORMATTED)

    call message(' - Reading '//trim(modelName(this%nmodel))//' geopotential coefficients...', LOG_AND_STDOUT)

    !====================================================================================================
    !
    !   Read earth radius and gravity constant and replace in
    !   astronomical constants module (it may be possible that it was not initialized appropriately)...
    !
    !--------------------------------------------------------------------------

    flag_mu   = .false.
    flag_rekm = .false.

    do i = 1,imax

      read(ich,'(a)',iostat=ios) cbuf

      !** earth gravity constant
      if(index(cbuf, "earth_gravity_constant") /= 0) then
        read(cbuf,*) ctemp, this%mu
        this%mu      = this%mu*1.d-9
        call setEarthGravity(this%mu)
        flag_mu = .true.
      else if(index(cbuf, "radius") /= 0) then
        read(cbuf,*) ctemp, this%rekm
        this%rekm      = this%rekm*1.d-3
        call setEarthGeopotentialRadius(this%rekm)
        flag_rekm = .true.
      else if(index(cbuf, "end_of_head") /= 0) then
        exit
      end if

      if(i == imax) then
        cmess = "Marker 'end_of_head' not found in file '"//trim(adjustl(fileName))//"'."
        call setNeptuneError(E_GRAVITY_COEFFICIENTS, FATAL, (/cmess/))
        ich  = closeFile(ich)
        return
      end if

    end do

    if(.not. flag_mu) then
      cmess ="Earth gravity constant not found for "//trim(modelName(this%nmodel))//" model. This may lead to inconsistencies when using the"// &
             " geopotential coefficients."
      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
    end if

    if(.not. flag_rekm) then
      cmess ="Earth radius not found for "//trim(modelName(this%nmodel))//" model. This may lead to inconsistencies when using the"// &
            " geopotential coefficients."
      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
    end if
    !-----------------------------------------------------------------

    !** check for maximum degree
    if(this%degree > maxdegree) then
      write(ctemp4,'(i2)') maxdegree
      call setNeptuneError(E_GEO_MAX_DEGREE, WARNING, (/ctemp4/))
    end if

    !** allocate coefficients
    !----------------------------------------
    if (allocated(this%coeffC)) deallocate(this%coeffC)
    allocate(this%coeffC(0:this%degree, 0:this%degree))
    if (allocated(this%coeffS)) deallocate(this%coeffS)
    allocate(this%coeffS(0:this%degree, 0:this%degree))
    if (allocated(this%useCoeffC)) deallocate(this%useCoeffC)
    allocate(this%useCoeffC(0:this%degree, 0:this%degree))
    if (allocated(this%useCoeffS)) deallocate(this%useCoeffS)
    allocate(this%useCoeffS(0:this%degree, 0:this%degree))
    if (allocated(this%sigmaC)) deallocate(this%sigmaC)
    allocate(this%sigmaC(0:this%degree, 0:this%degree))
    if (allocated(this%sigmaS)) deallocate(this%sigmaS)
    allocate(this%sigmaS(0:this%degree, 0:this%degree))
    if (allocated(this%sigmaCdenorm)) deallocate(this%sigmaCdenorm)
    allocate(this%sigmaCdenorm(0:this%degree, 0:this%degree))
    if (allocated(this%sigmaSdenorm)) deallocate(this%sigmaSdenorm)
    allocate(this%sigmaSdenorm(0:this%degree, 0:this%degree))


    if(this%degree >= 2) then  !** allocate time derivative of coefficients C and S
      if (allocated(this%dCdt)) deallocate(this%dCdt)
      allocate(this%dCdt(2:min(this%degree,4)))
      if (allocated(this%dSdt)) deallocate(this%dSdt)
      allocate(this%dSdt(2:min(this%degree,4)))
      if (allocated(this%cepoch)) deallocate(this%cepoch)
      allocate(this%cepoch(2:min(this%degree,4)))
    end if
    !-----------------------------------------

    l = 0

    !** read data
    !---------------------------------------------------
    do
      read(ich,'(a)', iostat=ios) cbuf
      if(ios /= 0) then
        if(this%degree > l) then
          write(ctemp4,'(i4)') l
          cmess = "Maximum degree of geopotential is "//trim(adjustl(ctemp4))// &
                  ". Maximum degree and order have been adapted accordingly."

          call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
          this%degree = l
        end if
        exit    !** end of file
      else if( (this%nmodel == EIGEN_GL04C .and. index(cbuf,'gfct') /= 0) .or. &
               (this%nmodel == EGM96       .and. cbuf(1:1) == 't') ) then
        ! write (*,*) cDataFormat(this%nmodel,1)
        read(cbuf,cDataFormat(this%nmodel,1)) ctemp, l, m, tempC, tempS, tempSigmaC, tempSigmaS, tempEpoch

        if(l > this%degree) cycle
        if(m > this%degree) exit

        this%coeffC(l,m) = tempC
        this%coeffS(l,m) = tempS
        this%sigmaC(l,m) = tempSigmaC
        this%sigmaS(l,m) = tempSigmaS

        read(tempEpoch(2:5),*) this%cepoch(l)%year
        read(tempEpoch(6:7),*) this%cepoch(l)%month
        read(tempEpoch(8:9),*) this%cepoch(l)%day

        this%cepoch(l)%hour   = 0
        this%cepoch(l)%minute = 0
        this%cepoch(l)%second = 0.d0

        call gd2jd(this%cepoch(l))

      else if( (this%nmodel == EIGEN_GL04C .and. index(cbuf,'gfc') /= 0) .or. &
               (this%nmodel == EGM08                                   ) .or. &
               (this%nmodel == EGM96       .and. cbuf(1:1) == ' ') ) then
        ! write (*,*) cDataFormat(this%nmodel,2)
        read(cbuf,cDataFormat(this%nmodel,2)) ctemp, l, m, tempC, tempS, tempSigmaC, tempSigmaS

        if(l > this%degree) cycle
        if(m > this%degree) exit

        this%coeffC(l,m) = tempC
        this%coeffS(l,m) = tempS
        this%sigmaC(l,m) = tempSigmaC
        this%sigmaS(l,m) = tempSigmaS

      else if( (this%nmodel == EIGEN_GL04C .and. index(cbuf,'dot') /= 0) .or. &
               (this%nmodel == EGM96       .and. cbuf(1:1) == 'd') ) then
        ! write (*,*) cDataFormat(this%nmodel,2)
        read(cbuf,cDataFormat(this%nmodel,2)) ctemp, l, m, tempdCdt, tempdSdt, tempSigmaC, tempSigmaS

        if(l > this%degree) cycle
        if(m > this%degree) exit

        this%dCdt(l) = tempdCdt
        this%dSdt(l) = tempdSdt    ! works only if there is ONE dSdt for each l --> S(2,1) in the case of EGM96

      else
        cmess = "Unknown identifier in file '"//trim(adjustl(fileName))
        call setNeptuneError(E_GRAVITY_COEFFICIENTS, FATAL, (/cmess/))
        ich  = closeFile(ich)
        return
      end if
    end do
    !-------------------------------------------------------

    !========================================================================
    !
    ! FEATURE/REMOVE: Output of coefficient sigma for certain degree
    !
    !------------------------------------------------------------------
!  open(unit=111, file="coeff_"//trim(dataFileName(nmodel)))
!  open(unit=112, file="csig_"//trim(dataFileName(nmodel)))
!
!  do l=0,degree
!
!    tempC = 0.d0
!    tempS = 0.d0

!    do m=0,l

!      tempS = tempS + sigmaC(l,m) + sigmaS(l,m)
!      tempC = tempC + coeffC(l,m)**2.d0 + coeffS(l,m)**2.d0

!    end do

!    tempS = tempS/(2.d0*l+1.d0)
!    tempC = sqrt(tempC/(l+1.d0))

!    write(111,'(i3,x,2(e20.12e2))') l, tempC
!    write(112,'(i3,x,e20.12e2)') l, tempS

!  end do
!
!  close(111)
!  close(112)

!  stop

    !========================================================================
    !
    ! De-normalize coefficients
    !
    !---------------------------------------------

    !** get factorials
    allocate(factorial(0:2*this%degree))

    factorial(0) = 1.d0
    do i = 1, 2*this%degree
      factorial(i) = factorial(i-1)*i
    end do

    do l = 0,this%degree
      do m = 0,l
        if(m == 0) then
          k = 1
        else
          k = 2
        end if
        fac = sqrt(factorial(l-m)*k*(2*l+1)/factorial(l+m))
        this%coeffC(l,m) = this%coeffC(l,m)*fac
        this%coeffS(l,m) = this%coeffS(l,m)*fac
        !sigmaC(l,m) = sigmaC(l,m)*fac  ! no de-normalization of sigma values, as these are used in normalized state by
        !sigmaS(l,m) = sigmaS(l,m)*fac  ! correlation routine.
        this%sigmaCdenorm(l,m) = this%sigmaC(l,m)*fac
        this%sigmaSdenorm(l,m) = this%sigmaS(l,m)*fac
      end do

      if(l < size(this%dCdt) .and. l >= 2) then
        this%dCdt(l) = this%dCdt(l)*fac
        this%dSdt(l) = this%dSdt(l)*fac
      end if
    end do

    deallocate(factorial)
    !---------------------------------------------------

    !============================================================
    !
    ! Set those harmonics to .false., which are not required, if
    ! a distinct harmonics analysis is requested.
    !
    !-------------------------------------------------------------
    if(this%analyseDistinctHarmonics .and. this%areAvailableDistinctHarmonics) then
      this%useCoeffC = .false.
      this%useCoeffS = .false.
      do l = 2, this%degree
        do m = 0, l
          !** search in distinct harmonics array for the combined number 'lm'
          lm = 10*l + m
          if(any(this%distinctHarmonics(:this%numDistinctHarmonics) == lm)) then
            this%useCoeffC(l,m) = .true.
            this%useCoeffS(l,m) = .true.
          end if
        end do
      end do
    else  ! otherwise, set all harmonics to true
      this%useCoeffC = .true.
      this%useCoeffS = .true.
    end if

    !============================================================
    !
    ! Allocate legendre polynomials
    !
    !-------------------------------------------------------------
    if (allocated(this%lp)) deallocate(this%lp)
    allocate(this%lp(0:this%degree+1, 0:this%degree+1))
    if (allocated(this%costerm)) deallocate(this%costerm)
    allocate(this%costerm(0:this%degree+1))
    if (allocated(this%sinterm)) deallocate(this%sinterm)
    allocate(this%sinterm(0:this%degree+1))
    if (allocated(this%rrfac)) deallocate(this%rrfac)
    allocate(this%rrfac(2:this%degree))
    !-------------------------------------------------------------

    ich = closeFile(ich)

    this%coeffInitialized = .true.
    this%degree_data      = this%degree

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine initGravityPotential
  !--------------------------------------------------------------------------------------------------
  !
  !> @anchor      getGravityAcceleration
  !!
  !> @brief       Providing accelerations due to geopotential
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 02.10.2012 (initial design)          </li>
  !!                <li> 31.05.2013 (code optimization)       </li>
  !!                <li> 17.07.2013 (providing tides support) </li>
  !!                <li> 17.12.2013 (added time_mjd in parameter list) </li>
  !!                <li> 14.03.2017 (removed unused r_gcrf in parameter list) </li>
  !!              </ul>
  !!
  !> @param[in]   reduction   reduction handler
  !> @param[in]   r_itrf      radius vector in ITRF / km
  !> @param[in]   v_itrf      velocity vector in ITRF / km/s
  !> @param[in]   time_mjd    date in MJD
  !!
  !> @param[out]  accel       acceleration vector in GCRF / km/s**2
  !!
  !> @details     This routine computes the acceleration due to Earth's geopotential
  !!              in an inertial frame (ECI).
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
  subroutine getGravityAcceleration(  this,     &
                                      reduction,&  ! <-- type responsible for coordinate transformations
                                      r_itrf,   &  ! <-- DBL(3) radius vector in ITRF / km
                                      v_itrf,   &  ! <-- DBL(3) velocity vector in ITRF / km/s
                                      time_mjd, &  ! <-- DBL    current time in MJD
                                      accel     &  ! --> DBL(3) acceleration vector in GCRF / km/s**2
                                   )

    implicit none

    !** interface
    !----------------------------------------------
    class(Gravity_class)                :: this
    type(Reduction_type)                :: reduction
    real(dp), dimension(3), intent(in)  :: r_itrf
    real(dp), dimension(3), intent(in)  :: v_itrf
    real(dp),               intent(in)  :: time_mjd

    real(dp), dimension(3), intent(out) :: accel
    !----------------------------------------------

    character(len=*), parameter :: csubid = "getGravityAcceleration" ! subroutine id

    integer :: l    ! loop counter
    integer :: m    ! loop counter

    real(dp), dimension(3) :: accel_temp    ! temporary
    real(dp) :: oosqrt_r1r2         ! 1/sqrt(r1r2)
    real(dp) :: temp, temp2         ! temporary
    real(dp) :: temp_t1             ! temporary
    real(dp) :: temp_t3             ! temporary
    real(dp) :: temp_t4             ! temporary
    real(dp) :: temp_t5             ! temporary

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether coefficients are already available
    if((this%degree > 1) .and. (.not. this%coeffInitialized)) then
      ! try to initialise
      call this%initGravityPotential()
      if(hasFailed()) return
    end if

    !** get radius, geocentric longitude and latitude
    call getRadiusLatLon(r_itrf, v_itrf, this%rabs, this%phi_gc, this%lambda)

    !** orbital radius
    this%rabs2        = this%rabs*this%rabs
    this%oorabs2      = 1.d0/this%rabs2
    this%rabs3        = this%rabs2*this%rabs
    this%oorabs       = this%rabs**(-1.d0) !1.d0/rabs
    this%oorabs3      = this%oorabs*this%oorabs2
    this%rekm = getEarthGeopotentialRadius()
    this%mu   = getEarthGravity()

    if(this%degree >= 2) then

      this%lp(0,0) = 1.d0
      this%lp(0,1) = 0.d0
      this%lp(1,0) = sin(this%phi_gc)
      this%lp(1,1) = cos(this%phi_gc)

      !** determine legendre polynomials recursively
      !if(.true.) then

      ! STABLE
      do m = 0, this%degree

        do l = max(2,m), this%degree

          if(l == m) then
            this%lp(m,m) = (2*m-1)*this%lp(1,1)*this%lp(m-1,m-1)
          else if(l == m + 1) then
            this%lp(l,m) = (2*m+1)*this%lp(1,0)*this%lp(l-1,m)
          else
            this%lp(l,m) = ((2*l-1)*this%lp(1,0)*this%lp(l-1,m) - (l+m-1)*this%lp(l-2,m))/(l-m)
          end if

        end do

      end do
      ! DEBUG: dump legendre polynomials
!      do l = 2, degree
!        do m = 0, l
!          write(34,*) l, m, lp(l,m)
!        end do
!      end do

!      else
      ! UNSTABLE??
      !--------------------------------------------------------------------
!      do l = 2,degree
!
!         lp(0,l-1) = 0.d0
!
!         do m = 0,l
!
!           if (m == 0) then
!
!             lp(l,0) = ((2*l-1)*lp(1,0) * lp(l-1,0) - (l-1)*lp(l-2,0) )/l
!
!           else
!
!             if(m == l .or. m == (l - 1)) then
!
!               lp(l,m) = (2*l-1) * lp(1,1) * lp(l-1,m-1)
!
!             else
!
!               lp(l,m) = lp(l-2,m) + (2*l-1) * lp(1,1) * lp(l-1,m-1)
!
!             end if
!
!           end if
!
!         end do ! DO m
!
!       end do ! DO L
!
!       ! DEBUG: dump legendre polynomials
!       do l = 2, degree
!         do m = 0, l
!           write(33,*) l, m, lp(l,m)
!         end do
!       end do
!
!       end if

!       stop
      !----------------------------------------------------------------------

      ! determine partial derivatives of the disturbing potential

      this%costerm(0) = 1.d0
      this%costerm(1) = cos(this%lambda)

      this%sinterm(0) = 0.d0
      this%sinterm(1) = sin(this%lambda)

      this%tanphi = tan(this%phi_gc)

      !** compute derivatives of potential
      call this%getPotentialDerivatives(this%oorabs, this%oorabs2, this%dudr, this%dudphi, this%dudlambda)

      ! pre-compute terms which are used for the acceleration components
      this%r1r2        = r_itrf(1)*r_itrf(1) + r_itrf(2)*r_itrf(2)
      temp_t3     = 1.d0/this%r1r2
      this%sqrt_r1r2   = sqrt(this%r1r2)
      oosqrt_r1r2 = 1.d0/this%sqrt_r1r2
      temp_t4     = r_itrf(3)*this%oorabs2*oosqrt_r1r2
      temp_t5     = this%oorabs2*this%sqrt_r1r2

    else    ! two body problem only

      this%dudlambda = 0.d0
      this%dudphi    = 0.d0
      this%dudr      = 0.d0

      temp_t3   = 0.d0
      temp_t4   = 0.d0

    end if

    !** add central body term to dU/dr
    temp_t1     = -this%mu*this%oorabs2 + this%dudr
    temp2       = this%dudlambda*temp_t3
    temp        = this%oorabs*temp_t1 - temp_t4*this%dudphi

    !=========================================================================
    !
    ! Compute the non-spherical perturbative accelerations in ITRF and convert
    !
    !-------------------------------------------------------------------------

    ! i-direction [km/s^2]
    accel_temp(1) = temp * r_itrf(1) - temp2 * r_itrf(2)
    ! j-direction [km/s^2]
    accel_temp(2) = temp * r_itrf(2) + temp2 * r_itrf(1)
    ! k-direction [km/s^2]
    accel_temp(3) = this%oorabs*temp_t1*r_itrf(3) + temp_t5 * this%dudphi

    call reduction%earthFixed2inertial(accel_temp, time_mjd, accel)

    !-------------------------------------------------------------------------

!if((time_mjd - 55276.75001d0) < 1.d-15) then
!    write(*,*) "acc_temp = ", accel_temp
!    write(*,*) "temp = ", temp, temp2
!    write(*,*) "r_itrf = ", r_itrf
!endif
    !write(*,*) "geop, accel = ", accel

    !** save state vector
    this%current_radius = r_itrf

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

  end subroutine getGravityAcceleration

  !--------------------------------------------------------------------------------------------------
  !
  !> @anchor      getPotentialDerivatives
  !!
  !> @brief       Computes dU/dr, dU/dphi and dU/dlambda for the geopotential
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 14.11.2013 (initial design)            </li>
  !!                <li> 17.01.2016 (implemented the option of switching off specific harmonics) </li>
  !!              </ul>
  !!
  !> @param[in]   oorabs   1/r
  !> @param[in]   oorabs2  1/r**2
  !> @param[out]  dur      dU/dr
  !> @param[out]  dup      dU/dphi
  !> @param[out]  dul      dU/dlambda
  !!
  !> @details     This routine computes the derivatives of the geopotential wrt. the radius, the
  !!              geocentric latitude and longitude in the ITRF. It is only to be called within the
  !!              gravity module, as it requires many parameters from within this module.
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine getPotentialDerivatives(this,oorabs, oorabs2, dur, dup, dul)

    class(Gravity_class)  :: this
    real(dp), intent(in)  :: oorabs
    real(dp), intent(in)  :: oorabs2

    real(dp), intent(out) :: dur
    real(dp), intent(out) :: dup
    real(dp), intent(out) :: dul

    integer  :: l, m   ! loop counter, being degree and order of geopotential
    real(dp) :: temp_t1
    real(dp) :: temp_t2
    real(dp) :: temp_t3
    real(dp) :: cfac, sfac  ! auxiliaries
    real(dp) :: insig1, insig2, insig3

    insig1 = 0.d0
    insig2 = 0.d0
    insig3 = 0.d0

    do l = 2,this%degree

      !** recursive computation of sin(ml) and cos(ml) (cos(0) and cos(1) as well as sin(0) and sin(1) have
      !   to be there already prior to call of this routine!
      !----------------------------------------------------------------------------------------------
      this%costerm(l) = 2.d0*this%costerm(1)*this%costerm(l-1)-this%costerm(l-2)
      this%sinterm(l) = 2.d0*this%costerm(1)*this%sinterm(l-1)-this%sinterm(l-2)

      ! determine pre-factor for expression inside the sigmas
      this%rrfac(l)   = (this%rekm*this%oorabs)**l
      this%lp(l,l+1) = 0.d0

      do m = 0,l !MIN(l,maxord)

        ! skip, if not requested for given harmonic
        if(.not. this%useCoeffC(l,m)) then
          cfac = 0.d0
        else
          cfac = 1.d0
        end if
        if(.not. this%useCoeffS(l,m)) then
          sfac = 0.d0
        else
          sfac = 1.d0
        end if

        temp_t1 = this%rrfac(l)*(l+1)*this%lp(l,m)
        temp_t2 = this%rrfac(l)*(this%lp(l,m+1) - (m * this%tanphi * this%lp(l,m)))
        temp_t3 = this%rrfac(l)* m * this%lp(l,m)

        ! radial partial derivative of the potential
        insig1 = insig1 + temp_t1 * (cfac*this%coeffC(l,m)*this%costerm(m) + sfac*this%coeffS(l,m)*this%sinterm(m))
        ! phi (latitudal) derivative of the potential
        insig2 = insig2 + temp_t2 * (cfac*this%coeffC(l,m)*this%costerm(m) + sfac*this%coeffS(l,m)*this%sinterm(m))
        ! lambda (longitudal) derivative of the potential
        insig3 = insig3 + temp_t3 * (sfac*this%coeffS(l,m)*this%costerm(m) - cfac*this%coeffC(l,m)*this%sinterm(m))

      end do
    end do

    !==========================================================================
    !
    ! Partial derivatives in the ITRF
    !
    !--------------------------------------------------------------------------

    temp_t1 = -this%mu*this%oorabs2
    temp_t2 = this%mu*this%oorabs

    ! compute the radial partial derivative of the potential
    dur = temp_t1 * insig1

    ! compute the latitutal partial derivative of the potential
    dup = temp_t2 * insig2

    ! compute the longitudal partial derivative of the potential
    dul = temp_t2 * insig3

    return

  end subroutine getPotentialDerivatives

  !================================================================================================
  !
  !> @anchor      getGravityCovariance
  !!
  !> @brief       Providing contributions to derivative of state transition matrix due to geopotential
  !> @author      Vitali Braun (VB)
  !> @author      Daniel Lück (DLU)
  !!
  !> @date        <ul>
  !!                <li>VB:  11.11.2013 (initial design)          </li>
  !!                <li>VB:  08.08.2014 (fixed a bug for d2U/dphi2) </li>
  !!                <li>VB:  09.08.2014 (fixed a bug for d2phi/dr2) </li>
  !!                <li>VB:  01.05.2015 (changed input to body-fixed) </li>
  !!                <li>VB:  14.03.2017 (removed unused dummy mjd) </li>
  !!                <li>DLU: 02.03.2021 (added calculation of legendre polynomials) </li>
  !!              </ul>
  !!
  !> @param[in]   r_itrf      radius vector in body-fixed frame
  !> @param[in]   v_itrf      velocity vector in body-fixed frame
  !!
  !> @details     This routine computes the contributions to the derivative of the state transition matrix
  !!              due to the geopotential. The equations from Vallado, D., "Fundamentals of Astrodynamics
  !!              and Applications", 4th Edition, Microcosm Press, 2013, pp.805-808
  !!
  !!------------------------------------------------------------------------------------------------
  function getGravityCovariance(this,r_itrf,v_itrf) result(cov)

    implicit none

    class(Gravity_class)               :: this
    real(dp), dimension(3), intent(in) :: r_itrf         ! radius vector (body-fixed)
    real(dp), dimension(3), intent(in) :: v_itrf         ! velocity vector (body-fixed)

    real(dp), dimension(3,3) :: cov           ! result in 3x3 matrix


    integer :: l,m                        ! loop counter
!    logical :: isAvailableFactors = .false.   ! flag to allow for faster execution, if some terms have already been computed by gravity acceleration routine

    real(dp), dimension(3,3) :: auxMat = reshape((/0,1,0,-1,0,0,0,0,0/), (/3,3/)) ! auxiliary matrix
    real(dp), dimension(3,3) :: ddr           ! first derivatives of r, phi and lambda wrt. r
    real(dp), dimension(3,3) :: ddr2          ! second derivative of radius wrt. radius vector
    real(dp), dimension(3,3) :: ddphi2        ! second derivative of latitude wrt. radius vector
    real(dp), dimension(3,3) :: ddlmb2        ! second derivative of longitude wrt. radius vector
    real(dp), dimension(3,3) :: imat          ! identity matrix
    real(dp) :: insig                         ! sum temporary
    real(dp) :: oor1r2                        ! 1/r1r2
    real(dp), dimension(3,3) :: sd            ! second derivatives matrix for geopotential
    real(dp) :: sec2                          ! secans of phi_gc squared
    real(dp) :: tanphi2                       ! tan(phi_gc) squared
    real(dp) :: temp

    if(.not. mag(r_itrf-this%current_radius) < eps9) then
      this%rabs    = mag(r_itrf)
      this%rabs2   = this%rabs*this%rabs
      this%oorabs  = 1.d0/this%rabs
      this%oorabs2 = this%oorabs/this%rabs
      this%oorabs3 = this%oorabs2/this%rabs

      !** get radius, geocentric longitude and latitude
      call getRadiusLatLon(r_itrf, v_itrf, this%rabs, this%phi_gc, this%lambda)

      this%costerm(0) = 1.d0
      this%costerm(1) = cos(this%lambda)

      this%sinterm(0) = 0.d0
      this%sinterm(1) = sin(this%lambda)

      this%tanphi    = tan(this%phi_gc)
      this%r1r2      = r_itrf(1)*r_itrf(1) + r_itrf(2)*r_itrf(2)
      this%sqrt_r1r2 = sqrt(this%r1r2)

      ! calculate legendre polynomials
      this%lp(0,0) = 1.d0
      this%lp(0,1) = 0.d0
      this%lp(1,0) = sin(this%phi_gc)
      this%lp(1,1) = cos(this%phi_gc)

      !** determine legendre polynomials recursively
      !if(.true.) then

      ! STABLE
      do m = 0, this%degree

        do l = max(2,m), this%degree

          if(l == m) then
            this%lp(m,m) = (2*m-1)*this%lp(1,1)*this%lp(m-1,m-1)
          else if(l == m + 1) then
            this%lp(l,m) = (2*m+1)*this%lp(1,0)*this%lp(l-1,m)
          else
            this%lp(l,m) = ((2*l-1)*this%lp(1,0)*this%lp(l-1,m) - (l+m-1)*this%lp(l-2,m))/(l-m)
          end if

        end do

      end do

      !** compute derivatives du/dr, du/dphi and du/dlmb
      !----------------------------------------------------------------------
      call this%getPotentialDerivatives(this%oorabs, this%oorabs2, this%dudr, this%dudphi, this%dudlambda)
      !----------------------------------------------------------------------

    end if

    !** secans
    sec2 = (1.d0/cos(this%phi_gc))**2.d0

    !** tanphi squared
    tanphi2 = this%tanphi*this%tanphi

    oor1r2 = 1.d0/this%r1r2

    !=============================================================
    !
    ! Compute 'matrix' components, representing the first
    ! part of the right hand side of the acceleration derivative
    ! wrt. the radius vector
    !
    !-------------------------------------------------------------

    !---------------------------------------------------
    !
    ! d2U/dr2
    !
    !---------------------------------------------------
    insig = 0.d0

    do l = 2, this%covDegree
      temp = this%rrfac(l)*(l+1)*(l+2)
      do m = 0,l !MIN(l,maxord)
        ! radial second partial derivative of the potential
        insig = insig + temp * this%lp(l,m)*(this%coeffC(l,m)*this%costerm(m) + this%coeffS(l,m)*this%sinterm(m))
      end do
    end do

    sd(1,1) = this%mu*this%oorabs3*insig

    !---------------------------------------------------
    !
    ! d2U/dphi2
    !
    !---------------------------------------------------
    insig = 0.d0

    do l = 2, this%covDegree
      do m = 0,l !MIN(l,maxord)
        temp = this%rrfac(l)*(this%tanphi*this%lp(l,m+1) + (m*m*sec2 - m*tanphi2- l*(l+1))*this%lp(l,m))
        ! radial second partial derivative of the potential
        insig = insig + temp*(this%coeffC(l,m)*this%costerm(m) + this%coeffS(l,m)*this%sinterm(m))
      end do
    end do

    sd(2,2) = this%mu*this%oorabs*insig

    !---------------------------------------------------
    !
    ! d2U/dlmbd2
    !
    !---------------------------------------------------
    insig = 0.d0

    do l = 2, this%covDegree
      do m = 0,l !MIN(l,maxord)
        temp  = m*m*this%lp(l,m)
        insig = insig + temp*this%rrfac(l)*(this%coeffC(l,m)*this%costerm(m) + this%coeffS(l,m)*this%sinterm(m))
      end do
    end do

    sd(3,3) = -this%mu*this%oorabs*insig

    !---------------------------------------------------
    !
    ! d2U/drdphi and d2U/dphidr
    !
    !---------------------------------------------------
    insig = 0.d0

    do l = 2, this%covDegree
     temp = this%rrfac(l)*(l+1)
     do m = 0,l !MIN(l,maxord)
        insig = insig + temp*(this%lp(l,m+1) - m*this%tanphi*this%lp(l,m))*(this%coeffC(l,m)*this%costerm(m) + this%coeffS(l,m)*this%sinterm(m))
      end do
    end do

    sd(1,2) = -this%mu*this%oorabs2*insig
    sd(2,1) = sd(1,2)

    !---------------------------------------------------
    !
    ! d2U/drdlmb and d2U/dlmbdr
    !
    !---------------------------------------------------
    insig = 0.d0

    do l = 2, this%covDegree
     temp = this%rrfac(l)*(l+1)
     do m = 0,l !MIN(l,maxord)
        insig = insig + temp*m*this%lp(l,m)*(this%coeffS(l,m)*this%costerm(m) - this%coeffC(l,m)*this%sinterm(m))
      end do
    end do

    sd(1,3) = -this%mu*this%oorabs2*insig
    sd(3,1) = sd(1,3)

    !---------------------------------------------------
    !
    ! d2U/dphidlmb and d2U/dlmbdphi
    !
    !---------------------------------------------------
    insig = 0.d0

    do l = 2, this%covDegree
     do m = 0,l !MIN(l,maxord)
        insig = insig + this%rrfac(l)*m*(this%lp(l,m+1) - m*this%tanphi*this%lp(l,m))*(this%coeffS(l,m)*this%costerm(m) - this%coeffC(l,m)*this%sinterm(m))
      end do
    end do

    sd(2,3) = this%mu*this%oorabs*insig
    sd(3,2) = sd(2,3)

    !=======================================================
    !
    ! Compute 'matrix' ddr components which multiply with
    ! the computed 'matrix' sd
    !
    !----------------------------------------------------

    !** dr/dr
    ddr(1,1:3) = this%oorabs*r_itrf(1:3)

    !** dphi/dr
    temp = 1.d0/this%sqrt_r1r2
    ddr(2,1:3) = -this%oorabs2*r_itrf(3)*r_itrf(1:3)*temp
    ddr(2,3)   = ddr(2,3) + temp

    !** dlamb/dr
    ddr(3,1) = -r_itrf(2)*oor1r2
    ddr(3,2) = r_itrf(1)*oor1r2
    ddr(3,3) = 0.d0
    !----------------------------------------------------

    !======================================================
    !
    ! Compute 'matrix' ddr2 components, which multiply with
    ! the derivatives of the geopotential
    !
    !------------------------------------------------------
    call identity_matrix(imat)

    ddr2   = this%oorabs*(imat - this%oorabs2*outerproduct(r_itrf,r_itrf))
    !write(*,*) " + + + DDR2 + + + "
    !do i=1,3
    !  write(*,*) (ddr2(i,j), j=1,3)
    !end do

    !write(*,*) " + + + Outerproduct + + + "
    !tempMat = outerproduct(r,r)
    !do i=1,3
    !  write(*,*) (1.d0 - tempMat(i,j)*oorabs2, j=1,3)
    !end do


    ddphi2 = -outerproduct((uvec3 - r_itrf(3)*this%oorabs2*r_itrf),(r_itrf(1)*uvec1 + r_itrf(2)*uvec2))/(this%sqrt_r1r2**3.d0) &
             -(outerproduct(r_itrf,uvec3) + r_itrf(3)*imat - 2.d0*this%oorabs2*r_itrf(3)* outerproduct(r_itrf,r_itrf))*this%oorabs2/this%sqrt_r1r2
    ddlmb2 = -2.d0*oor1r2*oor1r2*outerproduct((r_itrf(1)*uvec2 - r_itrf(2)*uvec1),(r_itrf(1)*uvec1 + r_itrf(2)*uvec2)) + oor1r2*auxMat

    !=================================================================
    !
    ! Now put everything together...
    !
    !-----------------------------------------------------------------
    cov = ddr2*this%dudr + ddphi2*this%dudphi + ddlmb2*this%dudlambda + matmul(transpose(matmul(sd,ddr)),ddr)

     !write(*,*) "leaving gravityCov"
    !do i=1,3
    !  write(*,*) (cov(i,j), j=1,3)
    !end do

    !write(*,*) " + + + DDR2 + + + "
    !do i=1,3
     ! write(*,*) (dudr*ddr2(i,j), j=1,3)
    !end do

    !write(*,*) " + + + DDPHI2 + + + "
    !do i=1,3
   !   write(*,*) (dudphi*ddphi2(i,j), j=1,3)
   ! end do

   ! write(*,*) " + + + DDLMB2 + + + "
   ! do i=1,3
   !   write(*,*) (dudlambda*ddlmb2(i,j), j=1,3)
   ! end do

    return

  end function getGravityCovariance


  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      getCoefficientSigma
  !!
  !> @brief       Returns the standard deviation for a harmonic coefficient C(l,m) or S(l,m)
  !> @author      Vitali Braun
  !!
  !> @param[in]   coeff     ! coefficient 'C' or 'S'
  !> @param[in]   l         ! degree
  !> @param[in]   m         ! order
  !!
  !> @return      sigma
  !!
  !> @date        <ul>
  !!                <li> 29.07.2013: initial design</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------------------------
  real(dp) function getCoefficientSigma(this, coeff, l, m)

    class(Gravity_class)  :: this
    character, intent(in) :: coeff
    integer,   intent(in) :: l
    integer,   intent(in) :: m

    if(l < 0 .or. m > l) then
      getCoefficientSigma = 0.d0
    else
      select case(coeff)
        case('C', 'c')
          getCoefficientSigma = this%sigmaC(l,m)
        case('S', 's')
          getCoefficientSigma = this%sigmaS(l,m)
        case default
          getCoefficientSigma = 0.d0
      end select
    end if

    return

  end function getCoefficientSigma

  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      getCurrentRadius
  !!
  !> @brief       Get the current radius vector in GCRF
  !> @author      Vitali Braun
  !!
  !> @return      position vector / km
  !!
  !> @date        <ul>
  !!                <li> 18.07.2013: initial design</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------------------------
  function getCurrentRadius(this)

    class(Gravity_class)    :: this
    real(dp), dimension(3)  :: getCurrentRadius

    getCurrentRadius = this%current_radius
    return

  end function getCurrentRadius

  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      getGeopotentialDegree
  !!
  !> @brief       Get the degree of the geopotential
  !< @author      Vitali Braun
  !!
  !> @return     geopotential degree
  !!
  !> @date        <ul>
  !!                <li> 17.07.2013: initial design</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------------------------
  integer function getGeopotentialDegree(this)

    class(Gravity_class)  :: this

    getGeopotentialDegree = this%degree
    return

  end function getGeopotentialDegree
  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      getGeoCovDegree
  !!
  !! @brief       Get the degree of the geopotential for the covariance matrix integration
  !! @author      Vitali Braun
  !!
  !> @return     geopotential degree covariance
  !!
  !! @date        <ul>
  !!                <li> 12.11.2013: initial design</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------------------------
  integer function getGeoCovDegree(this)

    class(Gravity_class)  :: this

    getGeoCovDegree = this%covDegree
    return

  end function getGeoCovDegree

  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      getGeopotentialModelName
  !!
  !> @brief       Get the name of the geopotential used
  !> @author      Vitali Braun
  !!
  !> @return     name of the geopotential as string
  !!
  !> @date        <ul>
  !!                <li> 21.10.2013: initial design</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------------------------
  character(len=len(modelName)) function getGeopotentialModelName(this)

    class(Gravity_class)  :: this

    if(this%nmodel == 0) then
      getGeopotentialModelName = "Kepler-only"
    else
      getGeopotentialModelName = modelName(this%nmodel)
    end if

    return

  end function getGeopotentialModelName

  !====================================================================================================
  !
  !> @anchor      getCurrentLatLon
  !!
  !> @brief       Get geocentric longitude and latitude, as well as altitude for current state vector
  !< @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 19.01.2013: initial design</li>
  !!                <li> 26.02.2013: added altitude output</li>
  !!                <li> 31.05.2013: code optimization </li>
  !!              </ul>
  !!
  !> @param[in]   imode   output format mode:
  !!                      0 = unformatted (rad)
  !!                      1 = formatted (deg, longitude [-180:180], latitude [-90:90])
  !> @param[out]  alt     altitude in km
  !> @param[out]  lat     geocentric latitude
  !> @param[out]  lon     geocentric longitude
  !!
  !> @details     This routine returns the geocentric longitude and latitude as well
  !!              as the altitude wrt. the reference ellipsoid for the current state vector
  !!
  !!------------------------------------------------------------------------------------------------
  subroutine getCurrentLatLon(this, imode, alt, lat, lon)

    class(Gravity_class)  :: this
    integer, intent(in)   :: imode  ! format of latitude and longitude
                                    !  0 = unformatted (rad)
                                    !  1 = formatted (deg,
                                    !      longitude: [-180:180]
                                    !      latitude:  [-90:90]
    real(dp), intent(out) :: alt    ! altitude in km
    real(dp), intent(out) :: lat    ! geocentric latitude
    real(dp), intent(out) :: lon    ! geocentric longitude

    character(len=*), parameter :: csubid = "getCurrentLatLon"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if


    if(imode == 0) then

      lat = this%phi_gc
      lon = this%lambda

    else if(imode == 1) then

      lat = this%phi_gc*rad2deg
      lon = this%lambda*rad2deg
      lon = redang(lon,1,1,.true.)

    end if

    alt = this%current_alt

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine getCurrentLatLon

  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      getLegendrePolynomial
  !!
  !> @brief       Returns zeroth order (m=0) Legendre polynomial using the recursion relation
  !> @author      Vitali Braun
  !!
  !> @date        <ul>
  !!                <li> 29.07.2013: initial design</li>
  !!              </ul>
  !!
  !> @param[in]   p2     legendre polynomial for step (t-2)
  !> @param[in]   p1     legendre polynomial for step (t-1)
  !> @param[in]   arg    argument of legendre polynomial
  !> @param[in]   n      degree of polynomial
  !!
  !> @details     The Legendre polynomial is computed using the recursion relation as defined in
  !!              Long et al. (1989, 4-13). Only the zero-th order term (m=0) is returned, using
  !!              the values for the polynomials for steps -1 and -2.
  !!
  !-----------------------------------------------------------------------------------------------
  real(dp) function getLegendrePolynomial(this, p2, p1, arg, n)

    class(Gravity_class) :: this
    real(dp), intent(in) :: p2    ! legendre polynomial for step -2
    real(dp), intent(in) :: p1    ! legendre polynomial for step -1
    real(dp), intent(in) :: arg   ! argument of legendre polynomial
    integer,  intent(in) :: n     ! degree of polynomial

    if(n == 0) then
      getLegendrePolynomial = 1.d0
    else if(n == 1) then
      getLegendrePolynomial = arg
    else
      getLegendrePolynomial = ((2.d0*n - 1.d0)*arg*p1 - (n - 1.d0)*p2)/n
    end if

    return

  end function getLegendrePolynomial

  !-----------------------------------------------------------------------------------------------
  !
  !> @anchor      isInitializedGeopotential
  !!
  !> @brief       Tells about the initialization state of geopotential module
  !> @author      Vitali Braun
  !!
  !> @return      .true. when the model is initialized
  !!
  !> @date        <ul>
  !!                <li> 29.07.2013: initial design</li>
  !!              </ul>
  !!
  !-----------------------------------------------------------------------------------------------
  logical function isInitializedGeopotential(this)

    class(Gravity_class) :: this

    isInitializedGeopotential = this%coeffInitialized
    return

  end function isInitializedGeopotential

!----------------------------------------------------------------------------------------------
!
!> @anchor      getDistinctHarmonicsString
!!
!> @brief       Get string containing all harmonics that are/were used during computation
!> @author      Vitali Braun
!!
!> @return      string
!!
!> @date        <ul>
!!                <li> 04.02.2015: initial design</li>
!!              </ul>
!!
!!------------------------------------------------------------------------------------------------
  character(len=3*maxDistinctHarmonics) function getDistinctHarmonicsString(this) result(str)

    class(Gravity_class) :: this
    integer              :: i,j                                                 ! loop counter

    str = ''
    j   = 1

    do i = 1, this%numDistinctHarmonics

      write(str(j:j+2),'(i2,a)') this%distinctHarmonics(i), ' '
      j = j + 3

    end do

    return

  end function getDistinctHarmonicsString

!----------------------------------------------------------------------------------------------
!
!> @anchor      setDistinctHarmonicsArray
!!
!> @brief       Set the array containing distinct harmonics for analysis
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 04.02.2015: initial design</li>
!!              </ul>
!!
!> @param[in]   arr   Array containing two-digit integers, which represent the individual spherical harmonics
!!
!> @detail      The spherical harmonics are provided as integers, e.g. '20' means that only
!!              C_{2,0} = -J_2 and S_{2,0} (= 0.0) will be considered. The first number thus stands for
!!              the degree, the second for the order. For instance, '31' means third degree harmonic
!!              of second order: C_{3,1}, S{3,1}.
!!
!!------------------------------------------------------------------------------------------------
  subroutine setDistinctHarmonicsArray(this,arr)

    class(Gravity_class)            :: this
    integer,dimension(:),intent(in) :: arr
    integer                         :: i, imax  ! loop counter

    !** init
    this%distinctHarmonics = 0

    imax = min(size(this%distinctHarmonics), size(arr))

    do i = 1, imax
      this%distinctHarmonics(i) = mod(arr(i), 100) ! modulo to make sure that only two-digit numbers are stored...
    end do

    !** now filter out those, where order > degree or where degree < 2
    do i = 1, imax

      if((this%distinctHarmonics(i)/10 < 2) .or. &
         (mod(this%distinctHarmonics(i), 10) > this%distinctHarmonics(i)/10)) then
        this%distinctHarmonics(i) = 0

      end if

    end do

    !** make sure that the array does not contain any zeroes "in-between"
    do i = 1, imax

      if(all(this%distinctHarmonics(i:imax) == 0)) then
        exit
      else if(this%distinctHarmonics(i) == 0) then
        do
          this%distinctHarmonics(i:imax) = eoshift(this%distinctHarmonics(i:imax),1)
          if(this%distinctHarmonics(i) /= 0) exit
        end do
      end if

    end do

    !** finally count number of harmonics
    this%numDistinctHarmonics = 0

    do i = 1, imax
      if(this%distinctHarmonics(i) > 0) then
        this%numDistinctHarmonics = this%numDistinctHarmonics + 1
      end if
    end do

    !** set flag
    if(this%numDistinctHarmonics > 0) then
      this%areAvailableDistinctHarmonics = .true.
    else
      this%areAvailableDistinctHarmonics = .false.
    end if

    return

  end subroutine setDistinctHarmonicsArray


!----------------------------------------------------------------------------------------------
!
!> @anchor      setDistinctHarmonics
!!
!> @brief       Set switch for analysing only distinct spherical harmonics
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 03.02.2015: initial design</li>
!!              </ul>
!!
!> @param[in]   flag
!!
!!------------------------------------------------------------------------------------------------
  subroutine setDistinctHarmonics(this,flag)

    class(Gravity_class)    :: this
    logical, intent(in)     :: flag

    this%analyseDistinctHarmonics = flag
    return

  end subroutine setDistinctHarmonics

!----------------------------------------------------------------------------------------------
!
!> @anchor      usedDistinctHarmonics
!!
!< @brief       Indicates if distinct harmonics are used
!> @author      Vitali Braun
!!
!> @return      .true. when distinct harmonics are used
!!
!> @date        <ul>
!!                <li> 04.02.2015: initial design</li>
!!              </ul>
!!
!> @param[in]   flag
!!
!!------------------------------------------------------------------------------------------------
  logical function usedDistinctHarmonics(this) result(flag)

    class(Gravity_class)    :: this

    if(this%analyseDistinctHarmonics .and. this%areAvailableDistinctHarmonics) then
      flag = .true.
    else
      flag = .false.
    end if

    return

  end function usedDistinctHarmonics

!----------------------------------------------------------------------------------------------
!
!> @anchor      setGeoDegree
!!
!> @brief       Set geopotential degree for state vector propagation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.11.2013: initial design</li>
!!                <li> 14.05.2015: handling input of 'ideg = 1', which equals to two-body solution </li>
!!              </ul>
!!
!> @param[in]   ideg     degree
!!
!!------------------------------------------------------------------------------------------------
  subroutine setGeoDegree(this,ideg)

    class(Gravity_class)    :: this
    integer, intent(in)     :: ideg

    character(len=255)      :: cmess
    character(len=*), parameter :: csubid = 'setGeoDegree'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(ideg < 0) then
      cmess = "Geopotential switch <0 (or undefined!). Set to '0', which is equivalent to the central body force."
      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
      this%degree = 0
    else if(ideg > this%degree_data) then
      this%coeffInitialized = .false.
      this%degree = ideg
    else if(ideg == 1) then
      this%degree = 0
    else
      this%degree = ideg
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine setGeoDegree

!----------------------------------------------------------------------------------------------
!
!> @anchor      setGeoCovDegree
!!
!> @brief       Set geopotential degree for covariance matrix integration
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 12.11.2013: initial design</li>
!!              </ul>
!!
!> @param[in]   ideg     degree
!!
!!------------------------------------------------------------------------------------------------
  subroutine setGeoCovDegree(this,ideg)

    class(Gravity_class)    :: this
    integer, intent(in)     :: ideg

    character(len=300)  :: cmess
    character(len=*), parameter :: csubid = 'setGeoCovDegree'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%degree < 0) then
      call this%setGeoDegree(this%degree)
    end if

    if(ideg > this%degree) then  ! covariance geopotential not allowed to be greater than geopotential for state vector

      cmess = "Geopotential switch for covariance matrix indicates a higher degree than used for state vector"// &
              " propagation. This is not supported. The degree of the covariance matrix geopotential has been"// &
              " set to the value of the geopotential as used for the state vector propagation."
      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
      this%covDegree = this%degree

    else if(ideg < 0) then

      cmess = "Geopotential switch <0 (or undefined!) for covariance matrix integration. Set to '0', which is"// &
              " equivalent to the central body force."
      call setNeptuneError(E_SPECIAL, WARNING, (/cmess/))
      this%covDegree = 0

    else

      this%covDegree = ideg

    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if


    return

  end subroutine setGeoCovDegree


!----------------------------------------------------------------------------------------------
!
!> @anchor      setCurrentAltitude
!!
!> @brief       Set altitude for current state vector
!> @author      Vitali Braun
!!
!< @date        <ul>
!!                <li> 26.02.2013: initial design</li>
!!                <li> 31.05.2013: code optimization</li>
!!              </ul>
!!
!> @param[in]   alt     current altitude
!!
!> @details     This function sets the current altitude to the specified value given via input
!!
!!------------------------------------------------------------------------------------------------
  integer function setCurrentAltitude(this,alt)

    class(Gravity_class)    :: this
    real(dp), intent(in)    :: alt    ! altitude in km

    this%current_alt = alt
    setCurrentAltitude = 0

    return

  end function setCurrentAltitude

!----------------------------------------------------------------------------------------------
!
!> @anchor      getHarmonicCoefficient
!!
!> @brief       Get value of a specific harmonic coefficient or its standard deviation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 13.01.2014 (initial design) </li>
!!                <li> 16.04.2017 (default return value HUGE) </li>
!!              </ul>
!!
!! @param[in]   coeff     Determines whether C (COEFF_C), S (COEFF_S) or the st.dev. dC (COEFF_DC) or dS (COEFF_DS) has to be returned
!! @param[in]   l         Degree of coefficient C_lm (dC_lm) or S_lm (dS_lm)
!! @param[in]   m         Order of coefficient C_lm (dC_lm) or S_lm (dS_lm)
!!
!> @return      specific harmonic coefficient
!!
!> @details     This routine returns the value of a specific harmonic coefficient or its standard
!!              deviation (being either C_lm, dC_lm, S_lm or dS_lm, depending on the three input
!!              variables) of the geopotential. Only the de-normalized are returned.
!!
!> @return      The value of the harmonic coefficient or its standard deviation
!!
!!------------------------------------------------------------------------------------------------
  real(dp) function getHarmonicCoefficient(this, coeff, l, m)

    class(Gravity_class)    :: this
    integer,  intent(in)    :: coeff
    integer,  intent(in)    :: l
    integer,  intent(in)    :: m

    character(len=*), parameter :: csubid = "getHarmonicCoefficient"
    character(len=200) :: cmess

    getHarmonicCoefficient = HUGE(1.d0)

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if


    !** check whether geopotential coefficients are already available
    if(.not. this%coeffInitialized) then
      ! try to initialise
      call this%initGravityPotential()
      if(hasFailed()) return
    end if

    !** check allowed ranges for l and m
    if((l < 2) .or. (l > this%degree) .or. (m < 0) .or. (m > l)) then
      call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
      return
    end if

    !** at this point, everything should be fine to finally return the coefficient
    select case(coeff)
      case(COEFF_C)
        getHarmonicCoefficient = this%coeffC(l,m)
      case(COEFF_DC)
        getHarmonicCoefficient = this%sigmaCdenorm(l,m)
      case(COEFF_S)
        getHarmonicCoefficient = this%coeffS(l,m)
      case(COEFF_DS)
        getHarmonicCoefficient = this%sigmaSdenorm(l,m)
      case default
        write(cmess,'(a,i3,4(a,i2),a)') "'", coeff, "', only possibilities are: '", COEFF_C, &
                                        "' (C_lm), '", COEFF_DC, "' (dC_lm), '", COEFF_S,    &
                                        "' (S_lm), '", COEFF_DS, "' (dS_lm)."
        call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/cmess/))
        return
    end select

    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end function getHarmonicCoefficient

!----------------------------------------------------------------------------------------------
!
!> @anchor      toggleHarmonicCoefficient
!!
!> @brief       Switch ON or OFF a specific harmonic (Stokes) coefficient
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 17.01.2016 (initial design) </li>
!!              </ul>
!!
!> @param[in]   coeff     Determines whether C (COEFF_C) or S (COEFF_S) are toggled
!> @param[in]   l         Degree of coefficient C_lm or S_lm
!< @param[in]   m         Order of coefficient C_lm or S_lm
!< @param[in]   switch    (optional) .true. switches on, .false. off
!!
!> @details     This routine allows to toggle individual Stokes coefficients. If the switch is
!!              provided, value of the useCoeffS/useCoeffC array will be set accordingly. Otherwise
!!              the state will be changed either from .true. to .false. or from .false. to .true..
!
!----------------------------------------------------------------------------------------------
  subroutine toggleHarmonicCoefficient(this, coeff, l, m, switch)

    class(Gravity_class)    :: this
    integer, intent(in)     :: coeff
    integer, intent(in)     :: l
    integer, intent(in)     :: m
    logical, optional, intent(in) :: switch

    character(len=*), parameter :: csubid = 'toggleHarmonicCoefficient'
    character(len=200) :: cmess

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether geopotential coefficients are already available
    if(.not. this%coeffInitialized) then
      ! try to initialise
      call this%initGravityPotential()
      if(hasFailed()) return
    end if

    !** check allowed ranges for l and m
    if((l < 2) .or. (l > this%degree) .or. (m < 0) .or. (m > l)) then
      call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
      return
    end if

    !** set according to switch, if provided
    if(present(switch)) then
      select case(coeff)
        case(COEFF_C)
          this%useCoeffC(l,m) = switch
        case(COEFF_S)
          this%useCoeffS(l,m) = switch
        case default
          write(cmess,'(a,i3,4(a,i2),a)') "'", coeff, "', only possibilities are: '", COEFF_C, &
                                          "' (C_lm), '", COEFF_S, "' (S_lm)."
          call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/cmess/))
          return
      end select
    else
      select case(coeff)
        case(COEFF_C)
          if(this%useCoeffC(l,m)) then
            this%useCoeffC(l,m) = .false.
          else
            this%useCoeffC(l,m) = .true.
          end if
        case(COEFF_S)
          if(this%useCoeffS(l,m)) then
            this%useCoeffS(l,m) = .false.
          else
            this%useCoeffS(l,m) = .true.
          end if
        case default
          write(cmess,'(a,i3,4(a,i2),a)') "'", coeff, "', only possibilities are: '", COEFF_C, &
                                          "' (C_lm), '", COEFF_S, "' (S_lm)."
          call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/cmess/))
          return
      end select
    end if

    if(isControlled()) then
      call checkOut(csubid)
    end if
    return

  end subroutine toggleHarmonicCoefficient
!----------------------------------------------------------------------------------------------
!
!> @anchor      setHarmonicCoefficient
!!
!> @brief       Set value of a specific harmonic coefficient or its standard deviation
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 13.01.2014 (initial design) </li>
!!              </ul>
!!
!> @param[in]   coeff     Determines whether C (COEFF_C), S (COEFF_S) or the st.dev. dC (COEFF_DC) or dS (COEFF_DS) has to be set
!> @param[in]   l         Degree of coefficient C_lm (dC_lm) or S_lm (dS_lm)
!> @param[in]   m         Order of coefficient C_lm (dC_lm) or S_lm (dS_lm)
!> @param[in]   val       Value the coefficient is set to
!!
!> @details     This routine sets the value of a specific harmonic coefficient or its standard
!!              deviation (being either C_lm, dC_lm, S_lm or dS_lm, depending on the first three input
!!              variables) of the geopotential according to the passed value. Only the de-normalized
!!              coefficients can be set.
!!
!!------------------------------------------------------------------------------------------------
  subroutine setHarmonicCoefficient(this,coeff, l, m, val)

    class(Gravity_class) :: this
    integer,  intent(in) :: coeff
    integer,  intent(in) :: l
    integer,  intent(in) :: m
    real(dp), intent(in) :: val

    character(len=*), parameter :: csubid = "setHarmonicCoefficient"
    character(len=200) :: cmess

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether geopotential coefficients are already available
    if(.not. this%coeffInitialized) then
      ! try to initialise
      call this%initGravityPotential()
      if(hasFailed()) return
    end if

    !** check allowed ranges for l and m
    if((l < 2) .or. (l > this%degree) .or. (m < 0) .or. (m > l)) then
      call setNeptuneError(E_OUT_OF_BOUNDS, FATAL)
      return
    end if

    !** at this point, everything should be fine to finally set the coefficient
    select case(coeff)
      case(COEFF_C)
        this%coeffC(l,m) = val
      case(COEFF_DC)
        this%sigmaCdenorm(l,m) = val
      case(COEFF_S)
        this%coeffS(l,m) = val
      case(COEFF_DS)
        this%sigmaSdenorm(l,m) = val
      case default
        write(cmess,'(a,i3,4(a,i2),a)') "'", coeff, "', only possibilities are: '", COEFF_C, &
                                        "' (C_lm), '", COEFF_DC, "' (dC_lm), '", COEFF_S,    &
                                        "' (S_lm), '", COEFF_DS, "' (dS_lm)."
        call setNeptuneError(E_UNKNOWN_PARAMETER, FATAL, (/cmess/))
        return
    end select

    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine setHarmonicCoefficient

!----------------------------------------------------------------------------------------------
!
!> @anchor      setTideSupport
!!
!> @brief       Provide support for the solid and oceand tides module
!> @author      Vitali Braun
!!
!> @date        <ul>
!!                <li> 16.07.2013 (initial design) </li>
!!              </ul>
!!
!< @param[in]   val     true/false
!!
!> @details     This routine activates the support for the tide module. The gravity module then
!!              will store computation parameters, which can be accessed by the tides module in
!!              order to save computation time, as similar parameters are required
!!
!!------------------------------------------------------------------------------------------------
  subroutine setTideSupport(this,val)

    class(Gravity_class)    :: this
    logical, intent(in)     :: val

    this%supportingTides = val
    return

  end subroutine setTideSupport

end module gravity
