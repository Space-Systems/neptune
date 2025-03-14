!-----------------------------------------------------------------------------------------------
!
!> @brief NEPTUNE testing routine (standalone) for open-mp
!!
!! @author      Christopher Kebschull (CHK)
!!
!! @date        <ul>
!!                <li>ChK:  25.02.2018 (initial implementation)       </li>
!!                <li>ChK:  09.03.2018 (Optimized the OpenMP loop)    </li>
!!              </ul>
!!
!! @details     This is the test routine for the openmp functionality
!!
!!
!! @anchor      openmp-test-sa
!!------------------------------------------------------------------------------------------------
program openmp_test_sa

    use omp_lib
    use slam_time,          only: time_t,mjd2gd,jd2gd,date2string
    use slam_orbit_types,   only: covariance_t, state_t
    use libneptune,         only: propagate,init_neptune
    use neptuneClass,       only: Neptune_class
    use neptune_error_handling, only: getNeptuneErrorMessage
    use slam_rframes,       only: getFrameId
    use slam_strings,       only: toString
    use slam_math,          only: mag
    use slam_astro,         only: getEarthRadius
    use slam_error_handling,only: getLatestError,hasFailed,printTrace,resetError,initErrorHandler
    use slam_io,                  only: STDOUT, LOG_AND_STDOUT, message
    use slam_types,         only: sp,dp

    implicit none
    integer,parameter                         :: num_objects  = 8               ! Number of objects
    integer,parameter                         :: n_threads    = 4               ! Number of available cores
    type(Neptune_class)                       :: single_neptune                 ! Neptune class instances

    type(time_t),dimension(num_objects)       :: startEpoch                     ! Star epoch
    type(time_t),dimension(num_objects)       :: endEpoch                       ! End epoch
    type(state_t),dimension(num_objects)      :: initialState                   ! State available from the object, which i used to initialize the propagator
    type(covariance_t),dimension(num_objects) :: initialCovariance              ! The covariance matrix
    type(state_t),dimension(num_objects)      :: propagatedState                ! State available from the object, which i used to initialize the propagator
    type(covariance_t),dimension(num_objects) :: propagatedCovariance           ! The covariance matrix
    type(time_t),dimension(2)                 :: tmpEpoch                       ! This array holds the start and end epochs as an array

    integer                         :: ineptune,iobject,itime                   ! Neptune counter

    ! For benchmarking
    character(len=12)               :: ctemp
    real(sp)                        :: totalTime
    real(sp), dimension(2)          :: timerArray
    real(dp)                        :: cpu_time_1, cpu_time_2, time, sys_time
    ! Variables for time of processing (CPU/System time)            [seconds]
    integer(kind=8)                 :: sys_time_1, count_rate_1, count_max_1
    integer(kind=8)                 :: sys_time_2, count_rate_2, count_max_2

    call cpu_time(cpu_time_1)
    call system_clock(sys_time_1, count_rate_1, count_max_1)

    call initErrorHandler(control = 'NO', errAction = 'RETURN', traceback = 'NO')

    do iobject=1,num_objects
        startEpoch(iobject)%jd = 2457023.3804976852d0
        !endEpoch%jd = 2457024.3804976852d0
        endEpoch(iobject)%mjd = 57023.8804976852d0

        initialState(iobject)%r(:) = (/ 968.352156207202d0, -176.201355234769d0, 6942.28113053269d0 /)
        initialState(iobject)%v(:) = (/ -0.151841688230462d0, -7.52766998424329d0, -0.171137603666182d0/)
        initialState(iobject)%a(:)   = 0.0d0
        call jd2gd(startEpoch(iobject))
        call mjd2gd(endEpoch(iobject))

        initialState(iobject)%epoch  = startEpoch(iobject)
        initialState(iobject)%frame  = getFrameId('GCRF')
        initialCovariance(iobject)%elem(:,:)   = 0.0d0
        initialCovariance(iobject)%frame = getFrameId('GCRF')
    end do

    call message("First initilization",STDOUT)
    ! Init in the master thread so that the data files are read only once
    single_neptune = Neptune_class()
    call initPropagator(single_neptune,startEpoch(1),endEpoch(1))
    call message("Starting propagation of "//toString(num_objects)//" objects",STDOUT)

do itime=1,5

!$omp parallel num_threads(n_threads) default(none) &
!$omp private(tmpEpoch,iobject) &
!$omp firstprivate(single_neptune) &
!$omp shared(propagatedState,propagatedCovariance,startEpoch,endEpoch, &
!$omp         initialState,initialCovariance)
!$omp do
    do iobject=1,num_objects
        write (*,*) "init propagator ", iobject, num_objects
        ! Note: We need to re-allocate the output arrays within neptune. This is for gfortran, as it does not re-allocate them automatically for every thread.
        ! ifort will do this w/o problems.
        call single_neptune%reallocate()
        call initPropagator(single_neptune,startEpoch(iobject),endEpoch(iobject))
        if (single_neptune%neptuneInitialized) then
            tmpEpoch(1) = startEpoch(iobject)
            tmpEpoch(2) = endEpoch(iobject)
            write (*,*) "starting propagator ", iobject, num_objects
            call propagate(single_neptune,initialState(iobject), &
                            initialCovariance(iobject), &
                            tmpEpoch, &
                            propagatedState(iobject), &
                            propagatedCovariance(iobject), &
                            .true.)
        end if
    end do
!$omp end do
!$omp end parallel

end do

    write(*,*) "propagated states:"
    do ineptune=1,num_objects
        write(*,*) propagatedState(ineptune)%r, propagatedState(ineptune)%v
    end do

    call single_neptune%destroy()

    call cpu_time(cpu_time_2)
    time = cpu_time_2 - cpu_time_1
    call system_clock(sys_time_2, count_rate_2, count_max_2)

    call message(' NEPTUNE finished nominally.', LOG_AND_STDOUT)
    call message(' - Elapsed time for propagation: ', LOG_AND_STDOUT)
    call dtime(timerArray, totalTime) ! There seems to be an issue on ifort with using totalTime ...
    write(ctemp,'(f12.4)') (timerArray(1) + timerArray(2)) ! ... so we just use both of them
    call message('    * total time: '//trim(adjustl(ctemp))//' s', LOG_AND_STDOUT)
    write(ctemp,'(f12.4)') timerArray(1)
    call message('    * user time : '//trim(adjustl(ctemp))//' s', LOG_AND_STDOUT)
    write(ctemp,'(f12.4)') timerArray(2)
    call message('    * system time : '//trim(adjustl(ctemp))//' s', LOG_AND_STDOUT)
    call message('CPU time: '//toString(time)//' sec', STDOUT)
    sys_time = dble(dble(sys_time_2 - sys_time_1)/count_rate_1)
    call message('System time: '//toString(sys_time)//' sec', STDOUT)

  contains

    subroutine initPropagator(neptune,startEpoch,endEpoch)
        use libneptune,             only: init_neptune
        use neptuneClass,           only: Neptune_class
        use slam_time,              only: time_t,jd2gd,mjd2gd,date2string
        use slam_strings,           only: toString
        use slam_error_handling,    only: getLatestError, hasFailed
        use neptune_error_handling, only: getNeptuneErrorMessage
        use slam_types,             only: dp

        implicit none

        type(Neptune_class)                 :: neptune                          ! Neptune class instance
        type(time_t),intent(inout)          :: startEpoch                       ! Star epoch
        type(time_t),intent(inout)          :: endEpoch                         ! End epoch

        !** Auxiliary variable
        integer :: ind                                                          ! Loop variable
        integer :: ind_err                                                      ! Checkup variable for error

        !** Auxiliary characters
        character(len=25)                                   :: cepoch_start
        character(len=25)                                   :: cepoch_end
        character(len=3),dimension(10)                      :: cpert
        character(len=15),dimension(10),parameter           :: perturbationValue = (/ "GEOPOTENTIAL   ", &
                                                                                      "ATMOSPHERE     ", &
                                                                                      "HORIZONTAL_WIND", &
                                                                                      "SUN            ", &
                                                                                      "MOON           ", &
                                                                                      "SRP            ", &
                                                                                      "ALBEDO         ", &
                                                                                      "SOLID_TIDES    ", &
                                                                                      "OCEAN_TIDES    ", &
                                                                                      "MANEUVERS      "/)
        real(dp)                                            :: mass = 156.d0
        real(dp)                                            :: cross_section = 1.d0
        real(dp)                                            :: c_R = 1.3d0
        real(dp)                                            :: c_D = 2.2d0

        ! Perturbations switches, mind the perturbationValue list,
        !  here we set the corresponding values.
        ! Disable everything first, then enable depending on the values from
        !  the database.
        cpert(2:10)  = "ON"
        cpert(1) = toString(6)
        cpert(2) = "ON"
        cpert(3) = "ON"
        cpert(4) = "ON"
        cpert(5) = "ON"
        cpert(6) = "ON"
        cpert(7) = "ON"
        cpert(8) = "ON"
        cpert(9) = "OFF"
        cpert(10) = "OFF"
        ! MANEUVERS are disabled by default.

        do ind=1,10
            ind_err = neptune%setNeptuneVar(perturbationValue(ind),TRIM(cpert(ind)))
        end do

        ind_err = neptune%setNeptuneVar("OPT_HARMONICS","OFF")
        ind_err = neptune%setNeptuneVar("OPT_SHADOW","ON")
        ind_err = neptune%setNeptuneVar("OPT_SRP_CORRECT","ON")

        !** Celestrak space weather
        !ind_err = neptune%setNeptuneVar("FILE_SOLMAG","sw19571001.txt")
        !** ESA's space weather
        ind_err = neptune%setNeptuneVar("FILE_SOLMAG","fap_day.dat")
        !ind_err = neptune%setNeptuneVar("FILE_SOLMAG_MONTHLY","fap_mon.dat")

        !** Options
        ind_err = neptune%setNeptuneVar("OPT_GEO_MODEL","3")    ! being the EIGEN-GL04C model
        ind_err = neptune%setNeptuneVar("OPT_AP_FORECAST","15")

        ind_err = neptune%setNeptuneVar("OPT_SAT_PROPERTIES","1")
        ind_err = neptune%setNeptuneVar("OPT_EOP","ON")
        ind_err = neptune%setNeptuneVar("OPT_INT_LOGFILE","OFF")

        !** Start and end epoch
        cepoch_start = date2string(startEpoch)
        cepoch_end   = date2string(endEpoch)

        ind_err = neptune%setNeptuneVar("EPOCH_START_GD",TRIM(cepoch_start))
        ind_err = neptune%setNeptuneVar("EPOCH_END_GD",TRIM(cepoch_end))

        ind_err = neptune%setNeptuneVar("PAR_MASS",toString(mass))
        ind_err = neptune%setNeptuneVar("PAR_CDRAG",toString(c_D))
        ind_err = neptune%setNeptuneVar("PAR_CROSS_SECTION", toString(cross_section))
        ind_err = neptune%setNeptuneVar("PAR_CREFL",toString(c_R))

        !** Covariance matrix
        ind_err = neptune%setNeptuneVar("COVARIANCE_PROPAGATION", "OFF")
        ind_err = neptune%setNeptuneVar("COVARIANCE_GEOPOTENTIAL","2")
        ind_err = neptune%setNeptuneVar("COVARIANCE_DRAG","OFF")
        ind_err = neptune%setNeptuneVar("COVARIANCE_MOON","OFF")
        ind_err = neptune%setNeptuneVar("COVARIANCE_SUN","OFF")
        ind_err = neptune%setNeptuneVar("COVARIANCE_SRP","OFF")

        ind_err = neptune%setNeptuneVar("PAR_INT_COV_STEP", toString(60.0d0))

        !** Output options
        ind_err = neptune%setNeptuneVar("OUTPUT_FILES","OFF")
        ind_err = neptune%setNeptuneVar("OPT_STORE_DATA", toString(10.d0))
        ind_err = neptune%setNeptuneVar('OPT_PROGRESS', 'OFF')

        !** Relative tolerance
        ind_err = neptune%setNeptuneVar("PAR_INT_RELEPS",toString(1.d-10))

        !** Absolute tolerance
        ind_err = neptune%setNeptuneVar("PAR_INT_ABSEPS",toString(1.d-11))

        ind_err = init_neptune(neptune)

    end subroutine initPropagator

end program openmp_test_sa
