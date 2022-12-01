!----------------------------------------------------------------
!
!> @brief   Perform NEPTUNE validation with Sentinel 1 data
!!
!! @author  Christopher Kebschull
!! @date    <ul> 
!!            <li> 14.07.2020: initial implementation</li>
!!          </ul>
!!
!----------------------------------------------------------------  
program valsent_set

    use slam_io
    use neptuneClass,           only: Neptune_class
    use libneptune,             only: init_neptune, propagate
    use neptuneClass,           only: Neptune_class
    use slam_math,              only: rad2deg
    use slam_orbit_types,       only: kepler_t, state_t, covariance_t
    use slam_reduction_class,   only: Reduction_type
    use slam_astro_conversions, only: rv2coe
    use slam_time,              only: time_t,date2longstring
    use slam_rframes,           only: getFrameId
    use slam_types

    character(len=50)   :: runid
    character(len=255)  :: outfile

    integer :: i,j      ! counter
    integer :: ich, ichr, ichv
    integer :: ndata    ! number of data sets for selected satellite

    logical :: init_flag, flag_backward
    real(dp) :: rms, rms_r, rms_v, max_rms
    real(dp), dimension(3)              :: sent_rUVW, nep_rUVW, sent_vUVW, nep_vUVW, rdiff, vdiff

    type(covariance_t)                  :: dummy1, dummy2
    type(kepler_t)                      :: sentKep
    type(state_t), dimension(100000)    :: sentState
    type(state_t)                       :: state_in, state_out, sent_state
    type(time_t), dimension(2)          :: epoch
    type(time_t), dimension(:), allocatable :: epochs
    type(Reduction_type)                :: reduction_model
    type(Neptune_class)                 :: neptune_instance

    !===================================
    !
    ! Instantiate reduction model and NEPTUNE
    !
    !---------------------------------
    reduction_model = Reduction_type()
    call reduction_model%initEop('data')
    neptune_instance = Neptune_class()

    !===================================
    !
    ! Read input file
    !
    !---------------------------------
    call rdinp(neptune_instance, runid, epoch)

    !============================================================
    !
    ! Read GPS POE data (containing transformation to GCRF)
    !
    !--------------------------------------------------------
    call readPOE(reduction_model, epoch, sentState, ndata)

    !** dump data to output file
    write(outfile,'(a)') 'output/'//trim(adjustl(runid))//'.dmp'
    ich = openFile(outfile, SEQUENTIAL, OUT_FORMATTED_OVERWRITE)

    do i=1,ndata
    write(ich,'(f15.8,3(x,f15.7),3(x,f12.8))') sentState(i)%epoch%mjd, (sentState(i)%r(j), j=1,3), (sentState(i)%v(j), j=1,3)
    end do

    close(ich)

    !** write kepler elements (osculating) for initial orbit
    call rv2coe(sentState(1), sentKep, i)

    !===========================================================
    !
    ! NEPTUNE propagation
    !
    !-----------------------------------------------

    allocate(epochs(ndata))

    if (epoch(1)%mjd < epoch(2)%mjd) then
        flag_backward = .false.
        sent_state = sentState(1)
        epochs(:) = sentState(:ndata)%epoch
    else  ! backwards
        flag_backward = .true.
        sent_state = sentState(ndata)
        epochs(:) = sentState(ndata:1:-1)%epoch
    end if


    i = neptune_instance%setNeptuneVar("INITIAL_STATE", sent_state)   ! setting initial state

    state_in = sent_state  ! initial state is GPS state
    state_in%frame  = getFrameId('GCRF')
    dummy1%frame  = getFrameId('GCRF')

    !** open comparison output file (radius)
    write(outfile,'(a)') 'output/'//trim(adjustl(runid))//'.rad'
    ichr = openFile(outfile, SEQUENTIAL, OUT_FORMATTED_OVERWRITE)
    write(outfile,'(a)') 'output/'//trim(adjustl(runid))//'.vel'
    ichv = openFile(outfile, SEQUENTIAL, OUT_FORMATTED_OVERWRITE)

    init_flag = .true.



    write(*,*) "Propagating from", epochs(1)%mjd, " to ", epochs(ndata)%mjd

    if ((epochs(1)%mjd > epochs(ndata)%mjd) .neqv. flag_backward) then
        write(*,*) "epochs not created correctly!"
        write(*,*) "epoch = ", epoch
        write(*,*) "epochs = ", epochs
        call exit(-1)
    end if

    if (all((epochs(:size(epochs)-1)%mjd > epochs(2:)%mjd) .neqv. flag_backward)) then
        write(*,*) "epochs not created correctly!"
        write(*,*) "epoch = ", epoch
        write(*,*) "epochs = ", epochs
        call exit(-1)
    end if

    write(*,*) "Starting NEPTUNE..."
    call propagate(neptune_instance,state_in, dummy1, epochs, state_out, dummy2, init_flag)
    write(*,*) "Done."

    state_out = neptune_instance%getNeptuneData(1)

    ! write(*,*) "Press enter to run some checks!"
    ! read(*,*)

    max_rms = -huge(0d0)

    do i=1,ndata-1

        ! write(*,*) "--------------"
        ! write(*,*) state_out%epoch%mjd, sent_state%epoch%mjd

        state_out = neptune_instance%getNeptuneData(i+1)
        
        if (flag_backward) then
            sent_state = sentState(ndata-i)
        else
            sent_state = sentState(i+1)
        end if

        ! write(*,*) "i=", i, ", nept=", state_out%epoch%mjd, ", sent=", sent_state%epoch%mjd

        if(abs(state_out%epoch%mjd - sent_state%epoch%mjd) > 1.d-6) then
            write(*,*) "something wrong with the epochs...."
            write(*,*) " - from Neptune: ", state_out%epoch%mjd, i+1
            write(*,*) " - from Sentinel:   ", sent_state%epoch%mjd, i+1
            write(*,*) " diff in secs:   ", (state_out%epoch%mjd - sent_state%epoch%mjd)*86400.d0
            write(*,*) "++ Fix first! ++"
            read(*,*)
            continue
        end if

        !** compute difference in UVW coordinates
        call reduction_model%eci2uvw(sent_state%r, sent_state%v, sent_state%r, sent_rUVW)
        call reduction_model%eci2uvw(sent_state%r, sent_state%v, state_out%r,   nep_rUVW)
        rdiff = nep_rUVW - sent_rUVW

        call reduction_model%eci2uvw(sent_state%r, sent_state%v, sent_state%v, sent_vUVW)
        call reduction_model%eci2uvw(sent_state%r, sent_state%v, state_out%v,   nep_vUVW)
        vdiff = nep_vUVW - sent_vUVW

        !** write differences to output
        rms_r = sqrt(sum(rdiff**2))
        rms_v = sqrt(sum(vdiff**2))
        rms = sqrt(sum(rdiff**2) + sum(vdiff**2))
        write(ichr,'(f15.8,2X,A,4(X,f15.7))') sent_state%epoch%mjd, date2longstring(sent_state%epoch), (rdiff(j), j=1,3)
        write(ichv,'(f15.8,2X,A,4(X,f12.8))') sent_state%epoch%mjd, date2longstring(sent_state%epoch), (vdiff(j), j=1,3)

        write(*,   '(f15.8,2X,A,7(X,f15.7))') sent_state%epoch%mjd, date2longstring(sent_state%epoch), rms, (rdiff(j), j=1,3), (vdiff(j), j=1,3)

        state_in = state_out
        max_rms = max(max_rms, rms)

    end do

    close(ichr)
    close(ichv)

    write(*,*) "Finished nominally."
    write(*,*) "max rms error = ", max_rms

end program valsent_set
