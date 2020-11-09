!=================================================================================================
!
!> @anchor      m_maneuvers
!!
!> @brief       Orbital maneuver modeling
!> @author      Vitali Braun (VB)
!> @author      Christopher Kebschull (CHK)
!!
!> @date        <ul>
!!                <li>VB:  04.11.2013 (initial design)</li>
!!                <li>CHK: 13.11.2015 (updated to use libslam)</li>
!!                <li>VB:  15.05.2016 (updated to have only statements for all use statements)</li>
!!                <li>CHK: 03.01.2018 (created Manoeuvres_class)</li>
!!              </ul>
!!
!> @details     This module contains parameters, subroutines and functions required for the
!!              modeling of orbital maneuvers. The individual maneuvers are read from an input
!!              file which is specified in the general NEPTUNE input file.
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      maneuvers
!!------------------------------------------------------------------------------------------------
module maneuvers

  use slam_types,             only: dp
  use slam_io,                only: openFile, closeFile, SEQUENTIAL, IN_FORMATTED, cdelimit, LOG_AND_STDOUT, message
  use neptune_error_handling, only: E_NO_MANEUVERS, E_MAX_MANEUVERS, E_MANEUVER_DEFINITION, E_SGA_INPUT, &
                                    E_MANEUVER_INIT, setNeptuneError
  use slam_error_handling,    only: isControlled, hasToReturn, hasFailed, FATAL, WARNING, REMARK, &
                                    E_ALLOCATION, checkIn, checkOut, E_OUT_OF_BOUNDS
  use slam_math,              only: mag
  use slam_reduction_class,   only: Reduction_type
  use slam_time,              only: time_t, tokenizeDate, sec_per_day

  implicit none

  private

    !===================================
    !
    ! Maneuver type definition(s)
    !
    !-----------------------------------
    type mnvPhase_t

        real(dp) :: mjd_start                                                   ! MJD of phase start
        real(dp) :: mjd_end                                                     ! MJD of phase end
        real(dp), dimension(3) :: acc                                           ! acceleration vector in UVW frame / km/s**2
        real(dp) :: thrust_efficiency                                           ! Efficiency of the thrust

    end type mnvPhase_t

    type maneuver_t

        type(time_t) :: start_date                                              ! date of maneuver start
        type(time_t) :: end_date                                                ! date of maneuver end
        integer      :: nphases                                                 ! number of different phases within one maneuver
        type(mnvPhase_t), dimension(:), allocatable :: phase                    ! containing individual maneuver phases

    end type maneuver_t

    type mnv_sequence_t

        real(dp) :: mjd_start                                                   ! MJD start of phase this structure is pointing at
        real(dp) :: mjd_end                                                     ! MJD end of phase this structure is pointing at
        integer  :: midx                                                        ! maneuver index
        integer  :: pidx                                                        ! phase index
        integer  :: trans                                                       ! transition flag: 0 = no transition,
                                                                                !                  1 = transition interval at begin of maneuver
                                                                                !                  2 = transition between two phases
                                                                                !                  3 = transition after a maneuver
    end type mnv_sequence_t

    !** types
    public :: maneuver_t
    public :: mnvPhase_t
    !-----------------------------------

    integer, parameter :: MANEUVER_STACK_SIZE     = 4000                        ! number of maneuvers the initial array is capable of storing
    integer, parameter :: MANEUVER_STACK_SIZE_MAX = 20000                       ! maximum number of maneuvers which can be handled
    integer, parameter :: MAX_MANEUVER_PHASES     = 50                          ! maximum number of maneuver phases which can be handled per maneuver

    integer, parameter :: NO_TRANSITION = 0
    integer, parameter :: TRANS_START   = 1
    integer, parameter :: TRANS_PHASE   = 2
    integer, parameter :: TRANS_END     = 3

    real(dp), parameter :: transition_time = 10.d0                              ! transition time between two maneuver phases / s
    real(dp), parameter :: transDay       = transition_time/86400.d0            ! duration of thrust transition in days

    type, public :: Manoeuvres_class

        logical :: man_initialized                                              ! initialization flag
        logical :: flag_no_data                                                 ! .true. if no maneuver data available for processing

        character(len=255) :: man_file_name                                     ! maneuver specification file name

        type(maneuver_t),    dimension(:), allocatable :: manv                  ! maneuver data structure
        type(mnv_sequence_t), dimension(:), allocatable :: mnv_sequence         ! structure, which points at the first maneuver phase for index 1,
                                                                                ! the second one for index 2, etc.
    contains

        procedure :: init_maneuvers_file
        procedure :: init_maneuvers_api
        generic   :: init_maneuvers => init_maneuvers_file, init_maneuvers_api
        procedure :: destroy

        !** setter
        procedure :: set_maneuvers_init_flag
        procedure :: set_man_file_name

        !** getter
        procedure :: get_man_file_name
        procedure :: get_manoeuvre_data
        procedure :: get_number_of_manoeuvres
        procedure :: get_maneuver_acceleration
        procedure :: get_upcoming_manoeuvre_change_epoch

        procedure, private :: read_maneuvre_file
        procedure, private :: init_maneuver_sequence
        procedure, private :: get_current_acc
        procedure, private :: get_current_index

    end type Manoeuvres_class

    ! Constructor
    interface Manoeuvres_class
        module procedure constructor
    end interface Manoeuvres_class

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !!  @author     Christopher Kebschull
    !!  @date       <ul>
    !!                  <li>ChK: 03.01.2017 (initial implementation)</li>
    !!              </ul>
    !!
    !!  @anchor     constructor
    !!
    ! --------------------------------------------------------------------
    type(Manoeuvres_class) function constructor()
        constructor%man_initialized = .false.                                    ! initialization flag
        constructor%flag_no_data   = .false.
        constructor%man_file_name = "neptune.mnv"
    end function constructor

    !===========================================================================
    !!
    !>  @anchor     destroy
    !!
    !!  @brief      Destroys all that needs destruction
    !!  @author     Christopher Kebschull
    !!
    !!
    !!  @date       <ul>
    !!                <li>03.01.2018 (initial design)</li>
    !!              </ul>
    !!
    !---------------------------------------------------------------------------
    subroutine destroy(this)
        class(Manoeuvres_class)    :: this

        if (allocated(this%mnv_sequence)) deallocate(this%mnv_sequence)
        if (allocated(this%manv)) deallocate(this%manv)

    end subroutine destroy

!==========================================================================
!
!> @anchor      set_maneuvers_init_flag
!!
!! @brief       Set initialization flag to .false.
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 04.11.2013 (initial design) </li>
!!              </ul>
!!
!!-------------------------------------------------------------------------
  subroutine set_maneuvers_init_flag(this)

    class(Manoeuvres_class)     :: this

    this%man_initialized = .false.

  end subroutine set_maneuvers_init_flag


  !
  subroutine read_maneuvre_file(this,            &
                                cpath,           & ! <-- CHR()   path to read data from
                                number_phases    & ! --> INT(:)  number of phases from file
                                )
    !** interface
    !-------------------------------------------
    class(Manoeuvres_class)                         :: this
    character(len=*), intent(in)                    :: cpath
    integer, dimension(:), allocatable, intent(out) :: number_phases  ! counting the number of phases for individual maneuvers
    !------------------------------------------

    character(len=255) :: cbuf
    character(len=255) :: cfile
    character(len=3)   :: ctemp
    character(len=26)  :: dateString
    character(len=*), parameter :: csubid = "read_maneuvre_file"
    integer :: i,j    ! loop counter
    integer :: ich    ! input channel
    integer :: ierr   ! error flag
    integer :: iman   ! number of individual maneuvers
    integer :: ios    ! I/O status
    integer :: kphs   ! counter for maneuver phases
    integer :: new_size   ! auxiliary for the reallocation of 'number_phases'


    logical  :: flag_current_maneuver   ! indicator for current maneuver read

    real(dp) :: duration    ! maneuver duration / s

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%man_initialized) then   ! return if already initialized
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
    end if

    iman = 0

    !==========================================================
    !
    ! Start maneuver input file read
    !
    !------------------------------------------

    !** open maneuver data file
    !----------------------------------------------------------------------
    cfile = trim(adjustl(cpath))//cdelimit//trim(adjustl(this%man_file_name))
    ich = openFile(cfile, SEQUENTIAL, IN_FORMATTED)

    call message(' - Reading maneuver data...', LOG_AND_STDOUT)

    !** count number of individual maneuvers (starting with a date tag)
    !-----------------------------------------------------------------

    !** allocate data array for maneuver phases
    allocate(number_phases(MANEUVER_STACK_SIZE))

    do

      read(ich,'(A)',iostat=ios) cbuf

      if(ios /= 0) then
        exit
      else if(index(cbuf,'#') /= 0) then
        cycle
      else if(index(cbuf,'-') /= 0 .and. index(cbuf,':') /= 0) then ! a date tag is available, if a '-' and a ':' are found

        iman = iman + 1

        ! in order to count maneuver phases, save number of phases for last maneuver, before starting to count from anew
        if(iman > 1) then
          !** re-allocate stack if required
          if(iman > size(number_phases)) then  ! re-allocate
            new_size = size(number_phases) + 20
            if(new_size > MANEUVER_STACK_SIZE_MAX) then  ! limit reached
              call setNeptuneError(E_ALLOCATION, FATAL)
              ierr = closeFile(ich)
              return
            else
              deallocate(number_phases)
              allocate(number_phases(new_size))
            end if
          end if
          number_phases(iman - 1) = kphs

        end if

        kphs = 1    ! start counting number of maneuver phases for this maneuver

      else

        kphs = kphs + 1

        if(kphs > MAX_MANEUVER_PHASES) then

          write(ctemp,'(i3)') MAX_MANEUVER_PHASES
          call setNeptuneError(E_MAX_MANEUVERS, FATAL, (/ctemp/))
          ierr = closeFile(ich)
          return

        end if

      end if

    end do

    if(iman == 0) then  !** no maneuvers found

      call setNeptuneError(E_NO_MANEUVERS, REMARK, (/cfile/))
      ierr = closeFile(ich)
      if(isControlled()) then
        call checkOut(csubid)
      end if

      this%flag_no_data = .true.
      return

    end if

    !** save phases for last maneuver, as this has not been counted yet
    number_phases(iman) = kphs

    !** allocated maneuver data structure
    !-------------------------------------
    allocate(this%manv(iman))
    do i = 1, iman

      !** allocate maneuver phases for each maneuver
      this%manv(i)%nphases = number_phases(i)
      allocate(this%manv(i)%phase(this%manv(i)%nphases))

    end do

    !** now finally read data into prepared arrays
    !------------------------------------------------------
    rewind(ich)
    i    = 0
    kphs = 2
    flag_current_maneuver = .false. ! is set to true, as soon as maneuver begins with date tag, set to false as soon as maneuver ends

    do

      read(ich,'(a)',iostat=ios) cbuf
      if(ios /= 0) exit
     if(index(cbuf,'#') /= 0) cycle  ! skip comments

      !** check whether new maneuver begins or only new phase is read
      if(index(cbuf,'-') /= 0 .and. index(cbuf,':') /= 0) then ! a date tag is available, if a '-' and a ':' are found

        i = i + 1
        read(cbuf,*) dateString, duration, (this%manv(i)%phase(1)%acc(j), j=1,3), this%manv(i)%phase(1)%thrust_efficiency

        !** convert we to km/s -- UPDATE: not required, as thrust is also given in kg*m/s**2, which cancels out to kg/s after division...
        !manv(i)%phase(1)%we = manv(i)%phase(1)%we/1.d3

        !** compute start and end date (MJD) for current phase, as well as end date for whole maneuver, if only one phase is available
        call tokenizeDate(dateString, this%manv(i)%start_date)

        this%manv(i)%phase(1)%mjd_start = this%manv(i)%start_date%mjd
        this%manv(i)%phase(1)%mjd_end   = this%manv(i)%start_date%mjd + duration/86400.d0

        if(this%manv(i)%nphases == 1) then
          this%manv(i)%end_date%mjd = this%manv(i)%start_date%mjd + duration/86400.d0
          flag_current_maneuver = .false.
        else
          flag_current_maneuver = .true.
        end if

      else if(flag_current_maneuver) then ! read new phase of current maneuver

        read(cbuf,*) duration, (this%manv(i)%phase(kphs)%acc(j), j=1,3), this%manv(i)%phase(1)%thrust_efficiency

        ! compute dates
        this%manv(i)%phase(kphs)%mjd_start = this%manv(i)%phase(kphs-1)%mjd_end
        this%manv(i)%phase(kphs)%mjd_end   = this%manv(i)%phase(kphs-1)%mjd_end + duration/86400.d0

        ! check whether end of current maneuver has been reached
        if(kphs == this%manv(i)%nphases) then

          ! compute end date of whole maneuver
          this%manv(i)%end_date%mjd = this%manv(i)%phase(kphs)%mjd_end

          kphs = 2  ! reset counter (note that first phase, kphs=1, is already read within 'if' statement)
          flag_current_maneuver = .false. ! end of maneuver

        else

          kphs = kphs + 1

        end if

      end if

    end do

    ierr = closeFile(ich)

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine read_maneuvre_file

!=============================================================================
!
!> @anchor      init_maneuvers_file
!!
!! @brief       Initialization of maneuvers module
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 04.11.2013 (initial design)</li>
!!              </ul>
!!
!! @param[in]   cpath       Path to read maneuver specification file from
!! @param[in]   date_start  propagation span start date (no data prior to this date is read)
!! @param[in]   date_end    propagation span end date (no data beyond this date is read)
!!
!! @details     This routine initializes the maneuvers module. Global parameters are set and
!!              maneuver data is read from a specification file into the specific arrays.
!!              Also the begin and end of the maneuver is prepared in order to provide a
!!              smooth maneuver onset and end.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  subroutine init_maneuvers_file(         &
                              this,      &
                              cpath,     &  ! <-- CHR()   path to read data from
                              date_start,&  ! <-- TYP     propagation start date (no data prior to this date is read)
                              date_end   &  ! <-- TYP     propagation end date (no data beyond this date is read)
                          )

    !** interface
    !-------------------------------------------
    class(Manoeuvres_class)      :: this
    character(len=*), intent(in) :: cpath
    type(time_t),     intent(in) :: date_start
    type(time_t),     intent(in) :: date_end
    !------------------------------------------

    character(len=*), parameter         :: csubid = "init_maneuvers_file"
    integer, dimension(:), allocatable  :: number_phases  ! counting the number of phases for individual maneuvers

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(this%man_initialized) then   ! return if already initialized --> This is not good! A user might call this routine, think to "set" manoeuvers, but actually the call is ignored without feedback
      if(isControlled()) then
        call checkOut(csubid)
      end if
      return
    end if

    !** check input
    if(date_start%mjd > date_end%mjd) then
      call setNeptuneError(E_SGA_INPUT, FATAL)
      return
    end if

    !=====================================================================
    !
    ! Read mnv file
    !
    !---------------------------------------------------------------------

    call this%read_maneuvre_file(cpath,number_phases)

    !=====================================================================
    !
    ! Initialize maneuvre sequence
    !
    !---------------------------------------------------------------------

    call this%init_maneuver_sequence(number_phases)

    this%man_initialized = .true.

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine init_maneuvers_file

!=============================================================================
!
!> @anchor      init_maneuvers_api
!!
!! @brief       Initialization of maneuvers module
!! @author      Christopher Kebschull
!!
!! @date        <ul>
!!                <li> 21.07.2017 (initial design)</li>
!!              </ul>
!!
!! @param[in]   manv_in          an array containing the individual manoeuvres
!! @param[in]   number_phases    array of manoeuvre phases
!!
!! @details     This routine initializes the maneuvers module. Global parameters are set and
!!              maneuver data is read from a specification file into the specific arrays.
!!              Also the begin and end of the maneuver is prepared in order to provide a
!!              smooth maneuver onset and end.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Finish. </li>
!!              </ol>
!!
!!------------------------------------------------------------------------------------------------
  subroutine init_maneuvers_api( this,           &
                                manv_in,         &  ! <-- TYP     List of manoeuvres planned to occur in propagation time
                                number_phases    &
                                )

    !** interface
    !-------------------------------------------
    class(Manoeuvres_class)                                     :: this
    type(maneuver_t), dimension(:), allocatable, intent(in)     :: manv_in           ! maneuver data structure
    integer,          dimension(:), allocatable, intent(in)     :: number_phases     ! counting the number of phases for individual maneuvers
    !------------------------------------------

    character(len=*), parameter :: csubid = "init_maneuvers_api"

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    ! Please not that we allow re-initialization

    !=====================================================================
    !
    ! Transfer maneuvre list
    !
    !---------------------------------------------------------------------
    if (allocated(this%manv)) deallocate(this%manv)

    allocate(this%manv(size(manv_in)))

    this%manv = manv_in

    !=====================================================================
    !
    ! Initialize maneuvre sequence
    !
    !---------------------------------------------------------------------

    call this%init_maneuver_sequence(number_phases)

    this%man_initialized = .true.

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine init_maneuvers_api


  subroutine init_maneuver_sequence(this, number_phases)

    class(Manoeuvres_class)                         :: this
    integer, dimension(:), allocatable, intent(in)  :: number_phases             ! counting the number of phases for individual maneuvers

    character(len=*), parameter     :: csubid = "init_maneuver_sequence"
    integer                         :: i,j,k                                    ! loop counter

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !=====================================================================
    !
    ! Prepare pointer array, in order to know maneuver sequence wrt. time
    !
    !---------------------------------------------------------------------

    if (allocated(this%mnv_sequence)) deallocate(this%mnv_sequence)

    allocate(this%mnv_sequence(size(this%manv)))  ! enough elements to account for transitions

    k = 1
    do i=1,size(this%manv)

      do j=1,this%manv(i)%nphases
        this%mnv_sequence(k)%midx        = i
        this%mnv_sequence(k)%pidx        = j
        this%mnv_sequence(k)%mjd_start   = this%manv(i)%phase(j)%mjd_start
        this%mnv_sequence(k)%mjd_end     = this%manv(i)%phase(j)%mjd_end
        this%mnv_sequence(k)%trans       = NO_TRANSITION

        k = k + 1     ! +1, as between each maneuver phase, also a transition phase will be added later on (see below)

      end do

    end do

    ! allocate(this%mnv_sequence(2*sum(number_phases(1:size(this%manv)))+size(this%manv)))  ! enough elements to account for transitions
    !                                                                             ! between thrust levels (cubic hermite spline, later on)
    ! k = 0

    ! do i=1,size(this%manv)


    !   do j=1,this%manv(i)%nphases

    !     k = k + 2                                                               ! +2, as between each maneuver phase, also a transition phase will be added
    !     this%mnv_sequence(k)%midx        = i
    !     this%mnv_sequence(k)%pidx        = j
    !     this%mnv_sequence(k)%mjd_start   = this%manv(i)%phase(j)%mjd_start + 0.5d0*transDay
    !     this%mnv_sequence(k)%mjd_end     = this%manv(i)%phase(j)%mjd_end   - 0.5d0*transDay
    !     this%mnv_sequence(k)%trans       = NO_TRANSITION

    !     !** transition phase
    !     this%mnv_sequence(k-1)%mjd_start = this%mnv_sequence(k)%mjd_start - transDay
    !     this%mnv_sequence(k-1)%mjd_end   = this%mnv_sequence(k)%mjd_start
    !     this%mnv_sequence(k-1)%midx      = i
    !     this%mnv_sequence(k-1)%pidx      = j

    !     if(j == 1) then

    !       this%mnv_sequence(k-1)%trans = TRANS_START  ! begin of maneuver

    !     else if(j <= this%manv(i)%nphases) then

    !       this%mnv_sequence(k-1)%trans = TRANS_PHASE  ! transition between two phases

    !     end if

    !     if(j == this%manv(i)%nphases) then   ! last phase: add transition interval at the end

    !       this%mnv_sequence(k+1)%mjd_start  = this%manv(i)%phase(j)%mjd_end - 0.5d0*transDay
    !       this%mnv_sequence(k+1)%mjd_end    = this%mnv_sequence(k+1)%mjd_start + transDay
    !       this%mnv_sequence(k+1)%midx       = i
    !       this%mnv_sequence(k+1)%pidx       = j
    !       this%mnv_sequence(k+1)%trans      = TRANS_END

    !     end if

    !   end do

    !   k = k + 1     ! +1, as between each maneuver phase, also a transition phase will be added later on (see below)

    ! end do

    !** final step: sort maneuver sequence
!    do
!
!      swapped = .false.
!
!      do i = 1, k-1
!        if(mnv_sequence(i)%mjd_start > mnv_sequence(i+1)%mjd_start) then
!          !** swap values
!          temp             = mnv_sequence(i)
!          mnv_sequence(i)   = mnv_sequence(i+1)
!          mnv_sequence(i+1) = temp
!          swapped  = .true.
!        end if
!      end do
!
!      if(.not. swapped) exit
!
!    end do

    !================================================================
    !
    ! Check whether there are intersections of phases (not allowed!)
    !
    !----------------------------------------------------------------
!    do i=1,size(mnv_sequence)-1
!        write (*,*) i, mnv_sequence(i+1)%mjd_start, mnv_sequence(i)%mjd_end
!      if(mnv_sequence(i+1)%mjd_start < mnv_sequence(i)%mjd_end) then
!        call setNeptuneError(E_MANEUVER_DEFINITION, FATAL)
!        return
!      end if
!    end do
!
!   do i=1,size(mnv_sequence)
!
!     write(*,'(i,x,2(f15.9,x),3(i,x))') i, mnv_sequence(i)%mjd_start, mnv_sequence(i)%mjd_end, mnv_sequence(i)%midx, mnv_sequence(i)%pidx, mnv_sequence(i)%trans
!
!   end do
!
!   stop

!   do i=1,size(manv)
!
!     write(*,*)
!     write(*,*) "Maneuver ",i
!     write(*,*) "-------------------"
!     write(*,*) " Begin date: ", manv(i)%start_date%mjd
!     write(*,*) " End   date: ", manv(i)%end_date%mjd
!
!     do j=1,size(manv(i)%phase)
!       write(*,*) "<-- Phase ", j, " -->"
!       write(*,*) "  -Begin:  ", manv(i)%phase(j)%mjd_start
!       write(*,*) "  -End:    ", manv(i)%phase(j)%mjd_end
!!       write(*,'(a,3(f8.2,x))') "  -Thrust: ", (manv(i)%phase(j)%acc(kphs), kphs=1,3)
!     end do
!   end do
!
!   write(*,*)

!    do dtemp = mnv_sequence(1)%mjd_start-3.d0*transDay, mnv_sequence(size(mnv_sequence))%mjd_end+3.d0*transday, 0.1d0/86400.d0
!
!      write(44,*) dtemp, mag(get_current_acc(dtemp,get_current_index(dtemp))), get_current_index(dtemp)
!
!    end do

!    stop
    !** finally set initial object mass -- not implemented!!
    !mass_old = getObjectMass()
    !mass_new = mass_old

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine init_maneuver_sequence

!=============================================================================
!
!> @anchor      get_maneuver_acceleration
!!
!! @brief       Computes acceleration due to orbital maneuvers
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 04.11.2013 (initial design)</li>
!!              </ul>
!!
!! @param[in]   r_gcrf            current radius vector in GCRF
!! @param[in]   v_gcrf            current velocity vector in GCRF
!! @param[in]   time_mjd          current time MJD
!! @param[in]   lout              (optional), if .true., then no updates are performed
!!
!! @param[out]  acc_atmosphere    acceleration vector in inertial frame
!!
!! @details     This routine computes the acceleration in the GCRF due to orbital
!!              maneuvers.
!!
!!              \par Overview
!!
!!              <ol>
!!                <li> Finish. </li>
!!              </ol>
!!
!------------------------------------------------------------------------------------------------
  subroutine get_maneuver_acceleration(   this,        &
                                        reduction,   &
                                        r_gcrf,      &   ! <-- DBL(3) current radius vector GCRF / km
                                        v_gcrf,      &   ! <-- DBL(3) current velocity vector GCRF / km/s
                                        time_mjd,    &   ! <-- DBL    current epoch
                                        acc_man      &   ! --> DBL(3) acceleration vector in inertial frame
                                    !    lout         &   ! <-- LOG    (optional), if .true., then no updates are performed, only acceleration for time_mjd given
                                    )                    !            which is important if only output for an intermediate step is required, but the current and
                                                         !            prior states have to remain untouched

    !** interface
    !---------------------------------------------------
    class(Manoeuvres_class)             :: this
    type(Reduction_type)               :: reduction
    real(dp),               intent(in)  :: time_mjd
    real(dp), dimension(3), intent(in)  :: r_gcrf
    real(dp), dimension(3), intent(in)  :: v_gcrf
!    logical,  optional,     intent(in)  :: lout

    real(dp), dimension(3), intent(out) :: acc_man
    !---------------------------------------------------

    character(len=*), parameter :: csubid = "get_maneuver_acceleration"

    real(dp), dimension(3) :: acc_uvw    ! current acceleration for given MJD in UVW frame

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    !** check whether there is maneuver data or not
    if(this%flag_no_data) then

      acc_man = 0.d0

      if(isControlled()) then
        call checkOut(csubid)
      end if

      return

    end if

    !** check whether initialization has been done already...
    if(.not. this%man_initialized) then
      call setNeptuneError(E_MANEUVER_INIT, FATAL)
      return
    end if

    !====================================================
    !
    ! Find out, whether current time corresponds to a
    ! possible maneuver. This is done by computing the
    ! index within the mnv_sequence array. This array holds the
    ! information of each phase of the individual maneuvers.
    ! E.g. mnv_sequence(1) is the start phase (transition) of the first maneuver,
    ! mnv_sequence(2) is the constant thrust phase of the first maneuver,
    ! etc. Now, the index is found by looking at the current date.
    ! The index is given a negative value if it corresponds
    ! to a phase with no thrust - which is actually not mapped by
    ! the mnv_sequence array. Thus, an index of "-1" means that
    ! the integrator is in a phase prior to phase 1, which would be
    ! before the start of the first maneuver. A "-5" would tell that
    ! the time corresponds to a no-thrust phase before phase 5, which itself
    ! would be a start phase of a subsequent maneuver.
    !
    ! <===< ONLY RELEVANT IF MASS UPDATES REQUIRED >===>
    ! The first thing to do, however, is to determine whether there has been a step
    ! forward, or the step size has been reduced and the function is called based
    ! on the same step one more time.
    !
    !----------------------------------------------------
!   if(present(lout)) then
!     lmode = lout
!   else
!     lmode = .false.
!   end if


!    if(lmode) then   ! routine has been called for output purposes - no step updates should be performed!!
!
!      tnew = time_mjd
!      inew = get_current_index(time_mjd)
!
!      if(time_mjd > time_mjd_new) then
!
!       told = time_mjd_new
!       !mold = mass_new
!       iold = idCurrentTime
!
!     else
!
!       told = time_mjd_old
!       !mold = mass_old
!       iold = idLastTime
!
!     end if

!   else

!     if(time_mjd > time_mjd_new) then  ! step forward - update time_mjd_old and time_mjd_new (moving one slot forward)
!                                       ! but also mass_old and the index idLastTime!

!       time_mjd_old = time_mjd_new
!       !mass_old     = mass_new

!       !** update also the mnv_sequence index
!       idLastTime    = idCurrentTime

!     end if

!     time_mjd_new  = time_mjd
!     idCurrentTime = get_current_index(time_mjd)

!     !** and now assign to local variables
!     told = time_mjd_old
!     tnew = time_mjd_new

!     !mold = mass_old

!     iold = idLastTime
!     inew = idCurrentTime

!   end if

!   !** if there was no change in index from subsequent calls, the thrust is zero:
!   if((idCurrentTime <= 0) .and. (idLastTime == idCurrentTime)) then

!     acc_man = 0.d0  ! no mass update
!
!     if(isControlled()) then
!       call checkOut(csubid)
!     end if
!     return

!   end if


    ! 1)  Get current thrust
    acc_uvw = this%get_current_acc(time_mjd, this%get_current_index(time_mjd))

    ! 2)  Convert acceleration to inertial frame, which results in the
    call reduction%uvw2eci(r_gcrf, v_gcrf, acc_uvw, acc_man)

    ! 3)  Update satellite's mass in order to compute current acceleration (stored in module variable mass_new)
    !call updateMass()

    write(*,*) time_mjd, mag(acc_uvw), v_gcrf

    ! 4)  Finally compute acceleration

    !acc_man(1:3) = thrust_gcrf*1.d-3/mass_new   ! 1.d-3 to account for conversion from 'N' to 'kg*km/s**2'

    !** done!
    if(isControlled()) then
      call checkOut(csubid)
    end if

    return

  end subroutine get_maneuver_acceleration

!=============================================================================
!
!> @anchor      getMeanMassFlow
!!
!! @brief       Gets the mean mass flow for a given time interval
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 07.11.2013: initial design</li>
!!              </ul>
!!
!! @detail      This function computes the mean mass flow of an engine for
!!              a given time interval. This is necessary, as the thrust level
!!              may change between subsequent calls of the maneuver acceleration
!!              routine. As the satellite's mass is important for the proper
!!              computation of the acceleration, the first step is to get the
!!              right mass flow
!!
!-----------------------------------------------------------------------------
!  real(dp) function getMeanMassFlow()


    !** make sure that mjd_end > mjd_start
!   if(abs(time_mjd_old - time_mjd_new) < eps9) then  ! same time, so mass flow is for the current step only and does not have to consider an interval

!     mdot = mag(getCurrentThrust(time_mjd_new))/getCurrentExVelocity()
!
!     swap      = mjd_start
!     mjd_start = mjd_end
!     mjd_end   = swap

!   end if

!    getMeanMassFlow = 0.d0

!    return

!  end function getMeanMassFlow


!=============================================================================
!
!> @anchor      get_current_acc
!!
!! @brief       Get the current acceleration for a given MJD in a satellite-based reference frame (UVW)
!! @author      Vitali Braun
!!
!! @param[in]   time_mjd  current MJD
!! @param[in]   idx       current index in mnv_sequence array
!!
!! @date        <ul>
!!                <li> 06.11.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  function get_current_acc(this, time_mjd, idx)

    class(Manoeuvres_class) :: this
    real(dp), intent(in)    :: time_mjd
    integer,  intent(in)    :: idx
    real(dp), dimension(3)  :: get_current_acc

    real(dp) :: ttemp   ! MJD transformed to interval [0:1]
    real(dp), dimension(3) :: y1, y2    ! interpolation points for thrust transitions


    !** if index is negative or zero, thrust = 0
    if(idx <= 0) then

      get_current_acc = 0.d0

    else

      ! if(this%mnv_sequence(idx)%trans == NO_TRANSITION) then  ! constant thrust in this phase

      !   get_current_acc(1:3) = this%manv(this%mnv_sequence(idx)%midx)%phase(this%mnv_sequence(idx)%pidx)%acc(1:3)

      ! else

      !   if(this%mnv_sequence(idx)%trans == TRANS_START) then   ! cubic hermite interpolation for maneuver start

      !     y1 = 0.d0
      !     y2(1:3) = this%manv(this%mnv_sequence(idx)%midx)%phase(this%mnv_sequence(idx)%pidx)%acc(1:3)

      !   else if(this%mnv_sequence(idx)%trans == TRANS_PHASE) then ! cubic hermite interpolation between phases

      !     y2(1:3) = this%manv(this%mnv_sequence(idx)%midx)%phase(this%mnv_sequence(idx)%pidx)%acc(1:3)
      !     y1(1:3) = this%manv(this%mnv_sequence(idx-1)%midx)%phase(this%mnv_sequence(idx-1)%pidx)%acc(1:3)

      !   else if(this%mnv_sequence(idx)%trans == TRANS_END) then ! cubic hermite interpolation for maneuver end

      !     y1(1:3) = this%manv(this%mnv_sequence(idx)%midx)%phase(this%mnv_sequence(idx)%pidx)%acc(1:3)
      !     y2(1:3) = 0.d0

      !   end if

      !   ttemp = (time_mjd - this%mnv_sequence(idx)%mjd_start)/transDay
      !   get_current_acc(1:3) = (2.d0*ttemp**3.d0 - 3.d0*ttemp**2.d0 + 1.d0)*y1(1:3) + (3.d0*ttemp**2.d0 - 2.d0*ttemp**3.d0)*y2(1:3)

      ! end if

      get_current_acc(1:3) = this%manv(this%mnv_sequence(idx)%midx)%phase(this%mnv_sequence(idx)%pidx)%acc(1:3) &
                                * this%manv(this%mnv_sequence(idx)%midx)%phase(this%mnv_sequence(idx)%pidx)%thrust_efficiency

    end if

    return

  end function get_current_acc

!=============================================================================
!
!> @anchor      get_current_index
!!
!! @brief       Get the current index within the mnv_sequence array
!! @author      Vitali Braun
!!
!! @param[in]   mjd   current MJD
!!
!! @date        <ul>
!!                <li> 07.11.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  integer function get_current_index(this, mjd)

    class(Manoeuvres_class) :: this
    real(dp), intent(in)    :: mjd                                              ! current MJD

    integer :: k    ! loop counter

    ! loop through mnv_sequence array to bracket mjd
    do k = 1, size(this%mnv_sequence)

      if(k == size(this%mnv_sequence) .and. mjd >= this%mnv_sequence(k)%mjd_start .and. mjd < this%mnv_sequence(k)%mjd_end) then
        get_current_index = k
        return
      else if(mjd < this%mnv_sequence(k)%mjd_start) then
        exit
      end if

    end do

    if(k == 1) then

      get_current_index = -1  ! prior to first maneuver

    else if(k > size(this%mnv_sequence)) then

      get_current_index =  0  ! after last maneuver phase

    else  ! somewhere in between of all maneuvers...

      get_current_index = k-1
      if((mjd > this%mnv_sequence(k-1)%mjd_end) .and. (mjd <= this%mnv_sequence(k)%mjd_start)) then ! in between of two maneuvers
        get_current_index = -(get_current_index + 1)  ! '+1' between phases, the negative index gives the number of the subsequent phase!
      end if

    end if

    return

  end function get_current_index

!=============================================================================
!
!> @anchor      get_upcoming_manoeuvre_change_epoch
!!
!! @brief       Get the start or end epoch of the next manoeuvre within the mnv_sequence array
!! @author      Christopher Kebschull
!!
!! @param[in]   mjd   current MJD
!!
!! @date        <ul>
!!                <li> 13.07.2020: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  real(dp) function get_upcoming_manoeuvre_change_epoch(this, mjd)

    class(Manoeuvres_class) :: this
    real(dp), intent(in)    :: mjd                                              ! current MJD

    integer :: k    ! loop counter

    get_upcoming_manoeuvre_change_epoch = 0.d0

    ! loop through mnv_sequence array to bracket mjd
    do k = 1, size(this%mnv_sequence)

      if( mjd < this%mnv_sequence(k)%mjd_start) then
        get_upcoming_manoeuvre_change_epoch = this%mnv_sequence(k)%mjd_start
        return
      else if(mjd < this%mnv_sequence(k)%mjd_end) then
        get_upcoming_manoeuvre_change_epoch = this%mnv_sequence(k)%mjd_end
        return
      end if

    end do

  end function


!=============================================================================
!
!> @anchor      get_man_file_name
!!
!! @brief       Get the name of the maneuver specification file
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 04.11.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  character(len=255) function get_man_file_name(this)

    class(Manoeuvres_class) :: this

    get_man_file_name = this%man_file_name
    return

  end function get_man_file_name


!=============================================================================
!
!> @anchor      get_manoeuvre_data
!!
!! @brief       Returns a manoeuvre form the manoeuvres array
!! @author      Christopher Kebschull
!!
!! @date        <ul>
!!                <li> 30.04.2019: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  function get_manoeuvre_data(this,idx) result(manv_record)

    class(Manoeuvres_class) :: this
    integer, intent(in)     :: idx
    type(maneuver_t)        :: manv_record

    character(len=*), parameter :: csubid = 'get_manoeuvre_data'

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    if(idx < 1 .or. idx > size(this%manv)) then

      call setNeptuneError(E_OUT_OF_BOUNDS, WARNING)
      manv_record%nphases = 0
      return

    end if

    manv_record = this%manv(idx)

    if(isControlled()) then
      call checkOut(csubid)
    end if

  end function get_manoeuvre_data


!=============================================================================
!
!> @anchor      get_number_of_manoeuvres
!!
!! @brief       Returns the size of the manoeuvres array
!! @author      Christopher Kebschull
!!
!! @date        <ul>
!!                <li> 30.04.2019: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
  integer function get_number_of_manoeuvres(this)

    class(Manoeuvres_class) :: this

    character(len=*), parameter :: csubid = 'get_number_of_manoeuvres'

    get_number_of_manoeuvres = 0

    if(isControlled()) then
      if(hasToReturn()) return
      call checkIn(csubid)
    end if

    get_number_of_manoeuvres = size(this%manv)

    if(isControlled()) then
      call checkOut(csubid)
    end if

  end function get_number_of_manoeuvres

!=============================================================================
!
!> @anchor      set_man_file_name
!!
!! @brief       Set the name of the maneuver specification file
!! @author      Vitali Braun
!!
!! @param[in]   val   name of maneuver specification file
!!
!! @date        <ul>
!!                <li> 04.11.2013: initial design</li>
!!              </ul>
!!
!-------------------------------------------------------------------------------
  subroutine set_man_file_name(this, val)

    class(Manoeuvres_class)         :: this
    character(len=*), intent(in)    :: val

    this%man_file_name = val(1:min(len(val),len(this%man_file_name)))
    return

  end subroutine set_man_file_name

!=============================================================================
!
!> @anchor      updateMass
!!
!! @brief       Update the mass for the given time step
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li> 08.11.2013: initial design</li>
!!              </ul>
!!
!-----------------------------------------------------------------------------
!  subroutine updateMass()

!   integer  :: k       ! loop counter
!   integer  :: k_start ! starting index
!   integer  :: k_end   ! ending index

!   real(dp) :: dsum    ! sum of product mdot*dt for different phases (basically an integration...)
!   real(dp) :: dt      ! time interval / days
!   real(dp) :: mdot    ! mass flow / kg/s
!   real(dp) :: t1, t2  ! time 1 and 2 for integration / days
!   real(dp) :: thrust_1, thrust_2    ! thrust 1 and 2 for integration / N
!   real(dp) :: we      ! exhaust velocity / km/s

!   dsum = 0.d0

!   write(*,*) "updating mass..."
!   write(*,*) "ID = ", idLastTime, idCurrentTime
!
!   if(idLastTime == 0) then    ! nothing to be done - exit without mass update
!     return
!   end if

!   !** find loop start and end indices from IDs

!   ! 1) k_start
!   if(idLastTime ==  -1) then
!     k_start = 1
!   else if(idLastTime < 0) then
!     k_start = -idLastTime - 1
!   else
!     k_start = idLastTime
!   end if

!   ! 2) k_end
!   if(idCurrentTime == 0) then
!     k_end = size(mnv_sequence)
!   else if(idCurrentTime ==  -1) then
!     k_end = 1
!   else if(idCurrentTime < 0) then
!     k_end = -idCurrentTime - 1
!   else
!     k_end = idCurrentTime
!   end if
!
!   !** loop over relevant phases for update
!   !----------------------------------------------
!   do k = k_start, k_end

!     write(*,*) "k = ", k

!     we   = manv(mnv_sequence(k)%midx)%phase(mnv_sequence(k)%pidx)%we
!     t1   = max(mnv_sequence(k)%mjd_start, time_mjd_old)
!     t2   = min(mnv_sequence(k)%mjd_end, time_mjd_new)
!     dt   = t2 - t1

!     write(*,*) "dt = ", dt*86400.d0
!     !** thrust for phase corresponding to time_mjd_old
!     thrust_1 = mag(getCurrentThrust(t1, k))

!     write(*,*) "thrust1 = ", thrust_1

!     select case(mnv_sequence(k)%trans)
!
!       case(NO_TRANSITION)
!
!         write(*,*) "no_trans"
!         mdot = thrust_1/we

!       case(TRANS_START, TRANS_PHASE, TRANS_END)

!         write(*,*) "start, phase, end"
!         !** compute thrust corresponding to time_mjd_new
!         thrust_2 = mag(getCurrentThrust(t2, k))
!         write(*,*) "thrust2 = ", thrust_2
!         mdot     = 0.5d0*(thrust_2 - thrust_1)/we

!     end select

!     write(*,*) "mdot = ", mdot
!     dsum = dsum + mdot*dt*day_in_sec
!     write(*,*) "dsum = ", dsum

!   end do

!   write(*,*) "mass_new = ", mass_new

!   !** now update mass
!   !-------------------------------
!   mass_new = mass_new - dsum

!   !** also update for satellite model
!   !----------------------------------
!   call setObjectMass(mass_new)

!   write(*,*) "done."
!   return

! end subroutine updateMass

end module maneuvers


