!>-----------------------------------------------------------------------------------------------
!!
!> @brief       Reads NEPTUNE input file
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  16.12.2012 (Output file switches for groundtrack and geopotential accelerations)</li>
!!                <li>VB:  24.01.2013 (Switches for atmosphere handling)      </li>
!!                <li>VB:  28.02.2013 (Switches for Sun and Moon)             </li>
!!                <li>VB:  13.03.2013 (Satellite parameters for SRP)          </li>
!!                <li>VB:  14.03.2013 (Switches for Solar radiation pressure) </li>
!!                <li>VB:  20.04.2013 (Interface/API redesign)                </li>
!!                <li>VB:  04.06.2013 (Code optimization)                     </li>
!!                <li>VB:  04.07.2013 (Covariance matrix added)               </li>
!!                <li>VB:  06.11.2014 (Lookup table switch for PN theory)     </li>
!!                <li>VB:  01.02.2015 (added new output file: total atmospheric density) </li>
!!                <li>VB:  03.02.2015 (added new input for distinct spherical harmonics in geopotential) </li>
!!                <li>VB:  29.06.2015 (added TEME frame input) </li>
!!                <li>VB:  23.08.2015 (added atmosphere model switch) </li>
!!                <li>VB:  06.10.2015 (added all other planets of the solar system as third body gravitational attractors) </li>
!!                <li>CHK: 16.11.2015 (updated to use libslam) </li>
!!                <li>VB:  24.01.2016 (covariance reference frame is now being assigned, as NEPTUNE will check for it)</li>
!!                <li>VB:  20.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>VB:  14.06.2016 (removed autoconfiguration, as it is not used)</li>
!!                <li>CHK: 02.01.2018 (regards NEPTUNE class)
!!              </ul>
!!
!> @param[inout] neptune   NEPTUNE class
!> @param[in]    f_input   Input file name
!> @param[out]   state     Initial cartesian state vector (GCRF)
!> @param[out]   covar     Initial covariance matrix (GCRF)
!> @param[out]   epoch     Begin and end date for propagation
!!                        <ul>
!!                          <li>  index = 1: begin date </li>
!!                          <li>  index = 2: end date   </li>
!!                        </ul>
!
!> @details     This routine reads the NEPTUNE input file, which is \c input/neptune.inp
!!              in the default settings. After the input is read it performs some basic input checks
!!              and initialisation procedures required for NEPTUNE, e.g. the transformation of the input
!!              epoch from julian day to gregorian date and vice versa, the transformation
!!              from true to mean anomaly, etc.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      rdinp
!!
!!------------------------------------------------------------------------------------------------
subroutine rdinp(                &
                  neptune,       &   ! <->  CLASS NEPTUNE class
                  f_input,       &   ! <--  CHR() input file name
                  state,         &   ! -->  TYP   initial cartesian state vector
                  covar,         &   ! -->  TYP   initial covariance matrix
                  epoch          &   ! -->  TYP() begin and end date for propagation
                )                    !             1 = begin date
                                     !             2 = end date
  use slam_types,          only: dp
  use slam_astro,          only: initAstroConstants
  use slam_astro_conversions, only: coe2rv, true2mean
  use slam_error_handling, only: hasToReturn, isControlled, hasFailed, checkIn, checkOut
  use slam_io,             only: closeFile, openFile, SEQUENTIAL, IN_FORMATTED, nxtbuf
  !use neptuneInput,        only:  setNeptuneVar
  use neptuneClass,        only: Neptune_class
  use neptuneParameters,   only: INPUT_TEME, INPUT_ITRF, INPUT_ITRF_TEME, INPUT_OSCULATING, INPUT_COV_UVW, INPUT_COV_GCRF
  use slam_math,           only: deg2rad
  use slam_orbit_types,    only: state_t, covariance_t, kepler_t, idimcov
  use slam_rframes,        only: getFrameId, REF_FRAME_UVW, REF_FRAME_GCRF
  use solarsystem,         only: n_supported_bodies
  use slam_units,          only: UNIT_RAD, UNIT_KMPS, UNIT_KM
  use slam_time,           only: time_t, date2string, jd245, checkDate, jd2gd, gd2mjd
  use version,             only: version_t

  implicit none

  !** interface
  !-----------------------------------------------------
  type(Neptune_class),         intent(inout)    :: neptune
  character(len=*),            intent(in)       :: f_input
  type(state_t),               intent(out)      :: state
  type(covariance_t),          intent(out)      :: covar
  type(time_t), dimension(2),  intent(out)      :: epoch
  !-----------------------------------------------------

  character(len=20) :: cepoch
  character(len=512), dimension(11) :: paths

  type(kepler_t) :: oscKep

  !----------------------------------------------------

  character(len=255) :: cbuf    ! character buffer
  character(len=*), parameter :: csubid = "rdinp"   ! subroutine id
  character(len=*), dimension(n_supported_bodies), parameter :: bodyOutNames = (/ &
                                               "ACS",  &
                                               "ACM",  &
                                               "AME",  &
                                               "ACV",  &
                                               "AMA",  &
                                               "ACJ",  &
                                               "ASA",  &
                                               "ACU",  &
                                               "ACN"/)
  character(len=*), dimension(9), parameter :: bodyNames = (/"SUN    ", &
                                               "MOON   ", &
                                               "MERCURY", &
                                               "VENUS  ", &
                                               "MARS   ", &
                                               "JUPITER", &
                                               "SATURN ", &
                                               "URANUS ", &
                                               "NEPTUNE"/)
  character(len=50)  :: runid    ! temporary run id as read from input file
  character(len=50)  :: ctemp    ! temporary
  character(len=2)   :: cgeomod  ! geopotential model
  character(len=1)   :: catmomod ! atmosphere model

  integer :: i,j,k                ! loop counter
  integer :: ios                ! I/O status
  integer :: ich_inp            ! input file channel
  integer :: ierr               ! error flag
  integer :: igeomod            ! geopotential model option
  integer :: itemp              ! temporary
  integer :: input_type
  integer :: input_type_cov
  integer :: b_left, b_right    ! indices for version string parsing

  integer, dimension(100) :: distinctHarmonics

  real(dp), dimension(3) :: dtmp3
  real(dp), dimension(6,6) :: jac   ! jacobian to convert from ECI to UVW frame

  type(state_t)   :: stateTEME              ! required for TEME state input
  type(state_t)   :: stateITRF              ! required for ITRF state input

  if(isControlled()) then
    if(hasToReturn()) return
    call checkIn(csubid)
  end if

  !** open input file
  ich_inp = openFile(f_input, SEQUENTIAL, IN_FORMATTED)

  if(hasFailed()) return

  ! First step is to find out about the input file version - if it's not there, it is a pre 0.9.15-beta version
  ! => file version is required to have backward compatibility
  do
    read(ich_inp,'(a)',iostat=ios) cbuf
    if(ios /= 0) exit
    if(index(cbuf, '[v') /= 0) then
        b_left  = index(cbuf,'[v')
        b_right = index(cbuf,']')
        neptune%version_model%current_version = neptune%version_model%parse_version_string(cbuf(b_left+2:b_right-1))
        if(hasFailed()) return
    end if
  end do

  rewind(ich_inp)

  !** read input file
  !----------------------------------------------

  !==============================================
  !
  !   General parameters
  !
  !-------------------------------------

  !** run id
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) runid
  ierr =  neptune%setNeptuneVar("RUN_ID", runid)

  !** progress writing - only since version 0.9.15-beta
  if(neptune%version_model%current_version > version_t(0,9,14)) then
    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) itemp
    if(itemp == 0) then
        ierr =  neptune%setNeptuneVar('OPT_PROGRESS', 'OFF')
    else
        ierr =  neptune%setNeptuneVar('OPT_PROGRESS', 'ON')
    end if

    call nxtbuf('#', 0, ich_inp, cbuf)
    ierr =  neptune%setNeptuneVar('FILE_PROGRESS', trim(cbuf))
  end if

  !==============================================
  !
  !   Input
  !
  !-------------------------------------
  !** input type (state)
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) input_type
  call neptune%output%set_input_type(input_type)

  !** input type (covariance)
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) input_type_cov
  call neptune%output%set_input_type_cov(input_type_cov)

  !** begin date (epoch(1)) and end date (epoch(2))
  do i=1,2

    call nxtbuf('#', 0, ich_inp, cbuf)
    cbuf = adjustl(cbuf)

    if(cbuf(5:5) == ' ') then ! gregorian date

      read(cbuf(1:4),  '(i4)') epoch(i)%year
      read(cbuf(6:7),  '(i2)') epoch(i)%month
      read(cbuf(9:10), '(i2)') epoch(i)%day
      read(cbuf(12:13),'(i2)') epoch(i)%hour
      read(cbuf(15:16),'(i2)') epoch(i)%minute
      read(cbuf(18:),*)        epoch(i)%second

      !** check date
      call checkDate( epoch(i)%year, epoch(i)%month, epoch(i)%day, epoch(i)%hour, epoch(i)%minute, epoch(i)%second)

      !** compute modified julian day
      call gd2mjd(epoch(i))

      epoch(i)%jd = epoch(i)%mjd + jd245

    else  ! julian day

      read(cbuf,*) epoch(i)%jd

      !** compute gregorian date
      call jd2gd(epoch(i)%jd, epoch(i)%year, epoch(i)%month, epoch(i)%day, epoch(i)%hour, epoch(i)%minute, epoch(i)%second)

      epoch(i)%mjd = epoch(i)%jd - jd245

    end if

    cepoch = date2string(epoch(i))

    if(i == 1) then
      ierr =  neptune%setNeptuneVar("EPOCH_START_GD", cepoch)
    else if(i == 2) then
      ierr =  neptune%setNeptuneVar("EPOCH_END_GD", cepoch)
    end if

  end do

  !** check if start epoch < end epoch
  !if(epoch(2)%jd .lt. epoch(1)%jd) then
!
!     cmess = "Start epoch > end epoch."
!     call error(csubid, cmess, 1)
!     ierr = E_INPUT_EPOCH
!     return
!
!   end if

  !===========================================
  !
  ! State vector: cartesian (GCRF)
  !
  !-------------------------------------------

  !** radius
  do i = 1,3

    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) state%r(i)

  end do

  !** velocity
  do i = 1,3

    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) state%v(i)

  end do

  !-------------------------------------------
  !
  ! State vector: osculating kepler elements
  !
  !-------------------------------------------

  !** semi-major axis
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) oscKep%sma

  !** eccentricity
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) oscKep%ecc

  !** inclination
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) oscKep%inc
  oscKep%inc = oscKep%inc*deg2rad

  !** right ascension of ascending node
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*)  oscKep%raan
  oscKep%raan = oscKep%raan*deg2rad

  !** argument of perigee
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) oscKep%aop
  oscKep%aop = oscKep%aop*deg2rad

  !** true anomaly
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) oscKep%tran
  oscKep%tran = oscKep%tran*deg2rad

  if(input_type == INPUT_OSCULATING) then   ! do further computations

    !** compute mean and eccentric anomaly
    call true2mean(oscKep%ecc, oscKep%tran, oscKep%ecan, oscKep%man)

    oscKep%epoch%jd   = epoch(1)%jd

    call jd2gd(oscKep%epoch%jd, oscKep%epoch%year, oscKep%epoch%month, oscKep%epoch%day, oscKep%epoch%hour, &
               oscKep%epoch%minute, oscKep%epoch%second)

  end if

  !------------------------------------------------
  !
  ! covariance matrix
  !
  !-----------------------------------------
  do i = 1, idimcov
    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) covar%elem(i,i)
  end do

  do i = 1, idimcov-1
    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) (covar%elem(i+1,j), j=1,i)
  end do

  !** fill remaining elements
  do i=1,idimcov-1
    do j=i+1,idimcov
      covar%elem(i,j) = covar%elem(j,i)
    end do
  end do

  if(input_type_cov == INPUT_COV_UVW) then
    covar%frame = REF_FRAME_UVW
  else if(input_type_cov == INPUT_COV_GCRF) then
    covar%frame = REF_FRAME_GCRF
  end if

  !========================================================
  !
  ! perturbation switches
  !
  !---------------------------------------------------

  !** geopotential
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("GEOPOTENTIAL", ctemp)

  !** atmosphere
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"

    if(itemp == 2) then ! consider horizontal wind
      ierr =  neptune%setNeptuneVar("HORIZONTAL_WIND", "ON")
    else
      ierr =  neptune%setNeptuneVar("HORIZONTAL_WIND", "OFF")
    end if

  end if
  ierr =  neptune%setNeptuneVar("ATMOSPHERE", ctemp)

  ! ** third bodies (only Sun and Moon so far)
  ! -----------------------------------------
  do i = 1, 2 !n_supported_bodies
    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) itemp
    if(itemp == 0) then
      ctemp = "OFF"
    else
      ctemp = "ON"
    end if
    ierr =  neptune%setNeptuneVar(trim(bodyNames(i)), ctemp)
  end do
  ! -----------------------------------------

  !** solar radiation pressure
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("SRP", ctemp)

  !** albedo
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("ALBEDO", ctemp)

  !** solid earth tides
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("SOLID_TIDES", ctemp)

  !** ocean tides
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OCEAN_TIDES", ctemp)

  !** maneuvers
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("MANEUVERS", ctemp)

  !** geopotential model to be used
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) cgeomod
  ierr =  neptune%setNeptuneVar("OPT_GEO_MODEL", cgeomod)

  !** atmosphere model to be used
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) catmomod
  ierr =  neptune%setNeptuneVar("OPT_ATMOSPHERE_MODEL", catmomod)


  !** check whether only single harmonics are to be used
  call nxtbuf('#', 0, ich_inp, cbuf)

  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OPT_HARMONICS", ctemp)

  !** read single harmonics to be used and store in an array
  call nxtbuf('#', 0, ich_inp, cbuf)

  i = 0
  do
    i = i + 1
    if(cbuf(i:i) == ' ') cycle ! skip leading blanks
    exit
  end do

  k = 0
  j = i

  do

    if(cbuf(j:j) == ' ') then
      j = j + 1
      if(j > len(cbuf) - 2) exit
      cycle ! skip blanks
    end if

    !** read first two characters
    k = k + 1
    read(cbuf(j:j+1),'(i2)') distinctHarmonics(k)
    j = j + 2
    if(j > len(cbuf) - 2) exit

    !** skip all other stuff until the next blank comes...
    do
      if(j > len(cbuf) - 2) exit ! avoid too high indices
      if(cbuf(j:j) /= ' ') then
        j = j + 1
        cycle
      else
        exit
      end if
    end do

  end do

  ierr =  neptune%setNeptuneVar("HARMONICS", distinctHarmonics(1:k))

  !do j=1,k
  !  write(*,*) "distinct = ", j, distinctHarmonics(j)
  !end do

  !** shadow model to be used
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "NONE"
  else
    ctemp = "CONICAL"
  end if

  ierr =  neptune%setNeptuneVar("OPT_SHADOW", ctemp)

  !** SRP correction
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
 ierr =  neptune%setNeptuneVar("OPT_SRP_CORRECT", ctemp)

  !** covariance matrix propagation switch
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("COVARIANCE_PROPAGATION", ctemp)

  !** covariance matrix geopotential degree
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("COVARIANCE_GEOPOTENTIAL", ctemp)

  !** covariance matrix: variational equations due to atmospheric drag
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("COVARIANCE_DRAG", ctemp)

  !** covariance matrix: variational equations due to the Sun
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("COVARIANCE_SUN", ctemp)
  if(hasFailed()) return

  !** covariance matrix: variational equations due to the Moon
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("COVARIANCE_MOON", ctemp)

  !** covariance matrix: variational equations due to SRP
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("COVARIANCE_SRP", ctemp)

  !** correlation matrix computation switch
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("CORRELATION_MATRIX", ctemp)


  !** further options
  !----------------------------
  !** EOP switch
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OPT_EOP", ctemp)

  !** PN lookup switch
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OPT_PN_LOOKUP", ctemp)

  !** Satellite properties switch
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("OPT_SAT_PROPERTIES", ctemp)

  !** Mass
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_MASS", ctemp)

  !** Cross-section (spherical object)
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_CROSS_SECTION", ctemp)

  !** Drag coefficient
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_CDRAG", ctemp)

  !** Reflectivity coefficient
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_CREFL", ctemp)

  !** Geomagnetic activity long-term forecast
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("OPT_AP_FORECAST", ctemp)

  !** Solar activity long-term forecast
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("OPT_SOL_FORECAST", ctemp)

  !** Re-entry stop condition
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_REENTRY", ctemp)


  !==============================================
  !
  !   Output
  !
  !-------------------------------------
  ierr =  neptune%setNeptuneVar("OUTPUT_FILES", "ON")

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("OUTPUT_STEP", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("OPT_STORE_DATA", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_CSV", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_OSC", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_GLL", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ATM", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_VAR_ECI", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_VAR_UVW", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_COV_ECI", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_COV_UVW", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACC", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACG", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACD", ctemp)

  ! third bodies (only Sun and Moon so far...)
  do i = 1, 2 !n_supported_bodies
    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) itemp
    if(itemp == 0) then
      ctemp = "OFF"
    else
      ctemp = "ON"
    end if
    ierr =  neptune%setNeptuneVar("OUTPUT_"//trim(bodyOutNames(i)), ctemp)
  end do

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACR", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACA", ctemp)

  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACT", ctemp)

 call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_ACO", ctemp)

 call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OUTPUT_AMN", ctemp)


  !==============================================
  !
  !   File paths
  !
  !-------------------------------------
  do i=1,size(paths)

    call nxtbuf('#', 0, ich_inp, cbuf)
    read(cbuf,*) paths(i)

  end do

  ierr =  neptune%setNeptuneVar("PATH_DATA",           paths(1))
  ierr =  neptune%setNeptuneVar("PATH_INPUT",          paths(2))
  ierr =  neptune%setNeptuneVar("PATH_OUTPUT",         paths(3))
  ierr =  neptune%setNeptuneVar("FILE_MANEUVERS",      paths(4))
  ierr =  neptune%setNeptuneVar("FILE_SURFACES",       paths(5))
  ierr =  neptune%setNeptuneVar("FILE_SOLMAG",         paths(6))
  ierr =  neptune%setNeptuneVar("FILE_SOLMAG_MONTHLY", paths(7))
  ierr =  neptune%setNeptuneVar("FILE_QDGRID",         paths(8))
  ierr =  neptune%setNeptuneVar("FILE_DWIND",          paths(9))
  ierr =  neptune%setNeptuneVar("FILE_HWIND",          paths(10))
  ierr =  neptune%setNeptuneVar("FILE_EOP",            paths(11))

  !==============================================
  !
  !   Numerical integration parameters
  !
  !-------------------------------------

  !** integration method
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_INT_METHOD", ctemp)

  !** relative tolerance
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_INT_RELEPS", ctemp)

  !** absolute tolerance
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_INT_ABSEPS", ctemp)

  !** covariance matrix integration step size
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_INT_COV_STEP", ctemp)

  !** covariance matrix integration method
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) ctemp
  ierr =  neptune%setNeptuneVar("PAR_INT_COV_METHOD", ctemp)

  !** logfile
  call nxtbuf('#', 0, ich_inp, cbuf)
  read(cbuf,*) itemp
  if(itemp == 0) then
    ctemp = "OFF"
  else
    ctemp = "ON"
  end if
  ierr =  neptune%setNeptuneVar("OPT_INT_LOGFILE", ctemp)
  !-----------------------------------------------

  !** close input file
  ich_inp = closeFile(ich_inp)

  !==================================================
  !
  !   Conversions...
  !
  !-----------------------

  !** convert kepler parameters to state vector
  !   OR TEME state to GCRF state
  !-----------------------------------------------
  if(input_type == INPUT_OSCULATING) then

    oscKep%truelon = oscKep%aop  + oscKep%raan + oscKep%tran
    oscKep%arglat  = oscKep%aop  + oscKep%tran
    oscKep%lonper  = oscKep%raan + oscKep%aop

    !** initialize astronomical constants
    read(cgeomod,*) igeomod

    call initAstroConstants(igeomod)

    oscKep%angles_unit = UNIT_RAD
    oscKep%sma_unit    = UNIT_KM

    call coe2rv(          &
                 oscKep,  &
                 state    &
               )

  else if(input_type == INPUT_TEME) then

    ! At this point EOPs have not been initialized yet - let's do this now
    call neptune%reduction%initEOP(neptune%getDataPath())
    stateTEME = state
    call neptune%reduction%teme2eci(stateTEME%r, stateTEME%v, (/0.d0,0.d0,0.d0/), & ! no accelerations..
                  epoch(1)%mjd, state%r, state%v, dtmp3)

  else if(input_type == INPUT_ITRF) then

    ! At this point EOPs have not been initialized yet - let's do this now
    call neptune%reduction%initEOP(neptune%getDataPath())
    stateITRF = state
    call neptune%reduction%earthFixed2inertial(stateITRF%r, stateITRF%v, (/0.d0,0.d0,0.d0/), & ! no accelerations..
                  epoch(1)%mjd, state%r, state%v, dtmp3)

  else if(input_type == INPUT_ITRF_TEME) then

    ! At this point EOPs have not been initialized yet - let's do this now
    call neptune%reduction%initEOP(neptune%getDataPath())
    stateITRF = state
    call neptune%reduction%earthFixed2inertial(stateITRF%r, stateITRF%v, (/0.d0,0.d0,0.d0/), & ! no accelerations..
                  epoch(1)%mjd, state%r, state%v, dtmp3)
    stateTEME = state
    call neptune%reduction%teme2eci(stateTEME%r, stateTEME%v, (/0.d0,0.d0,0.d0/), & ! no accelerations..
                  epoch(1)%mjd, state%r, state%v, dtmp3)

  end if

  state%frame = getFrameId('GCRF')  ! set frame ID

  !** conversions for covariance matrix
  !----------------------------------------------

  if(input_type_cov == INPUT_COV_UVW) then
    !** convert to ECI frame
    call neptune%reduction%getJacobianEci2uvw(state%r, state%v, jac)
    covar%elem = matmul(matmul(transpose(jac),covar%elem),jac)
    covar%frame = getFrameId('GCRF')
  end if

  !** finally set initial state and covariance in NEPTUNE input module
  !--------------------------------------------------------------------
  state%epoch = epoch(1)
  if(input_type == INPUT_OSCULATING) then
    ierr =  neptune%setNeptuneVar("INITIAL_STATE", oscKep)
  else
    state%radius_unit   = UNIT_KM
    state%velocity_unit = UNIT_KMPS
    ierr =  neptune%setNeptuneVar("INITIAL_STATE", state)
  end if

  ierr =  neptune%setNeptuneVar("INITIAL_COVARIANCE", covar)

  !** done!
  if(isControlled()) then
    call checkOut(csubid)
  end if

  return

end subroutine rdinp
