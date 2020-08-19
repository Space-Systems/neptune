!----------------------------------------------------------------
!
!> @brief     Reads POE data from files in SP3-c format
!!
!> @author    Vitali Braun
!> @date      <ul>
!!              <li> 06.02.2014: initial implementation</li>
!!              <li> 14.07.2014: adapted for Sentinel 1b</li>
!!            </ul>
!!
!> @param[in] reduction_model - used for ECEF -> ECI conversion
!!
!> @param[in] epoch - Epoch for which data has to be provided
!!
!> @param[inout] stateGCRF - Holds the POE states
!!
!> @param[inout] ndata - Defines number of states available
!!
!----------------------------------------------------------------  
subroutine readPOE(reduction_model, epoch, stateGCRF, npoe)

    use slam_io
    use slam_math,              only: mag
    use slam_orbit_types,       only: state_t
    use slam_reduction_class,   only: Reduction_type
    use slam_rframes
    use slam_time,              only: gd2mjd, mjd2gd, delta_AT, time_t, mjd2yyddd
    use slam_types

    implicit none

    type(Reduction_type),intent(in)             :: reduction_model
    type(time_t), dimension(2),  intent(in)     :: epoch
    type(state_t), dimension(*), intent(inout)  :: stateGCRF
    integer, intent(out)                        :: npoe


    character(len=255)                          :: cbuf, cbuf2      ! character buffer
    character(len=255), dimension(100)          :: cfile            ! file to read data from

    integer             :: p                    ! counter
    integer             :: ich                  ! channel to read data from
    integer             :: nfiles               ! number of different files to read
    real(dp)            :: dummy, x, y, z, vx, vy, vz
    character(len=35)   :: string_epoch
    logical             :: osv_detected
    integer             :: start_nbr
    integer             :: stop_nbr

    type(state_t) :: state

    write (*,*) "Reading POE"

    osv_detected = .false.

    !=================================================================
    !
    ! Find out data file(s) to read data from based on passed epoch
    !
    !----------------------------------------------------------------

    if(epoch(2)%mjd - epoch(1)%mjd > 10) then
        write(*,*) "Not designed for epochs > 10 days"
        write(*,*) "+++ Program stopped!+++"
        stop
    end if

    ! starting file
    nfiles = 0
    ich = openFile('./poeData/available_files.txt', SEQUENTIAL, IN_FORMATTED)
    do
        call nxtbuf('#', 0, ich, cbuf)
        if (.not. cbuf .eq. '') then
            nfiles = nfiles + 1
            cfile(nfiles) = trim(cbuf)
        else
            exit
        end if
    end do
    close(ich)

    npoe = 0
    do p=1,nfiles
        
        write (*,*) "Reading file ", trim(cfile(p))
        ich = openFile("./poeData/"//trim(cfile(p)), SEQUENTIAL, IN_FORMATTED)
        do
            call nxtbuf('#', 0, ich, cbuf)
            if (.not. cbuf .eq. '') then

                ! Scan for XML-Tag
                if (index(cbuf,'<OSV>')/=0) then
                    osv_detected = .true.
                elseif (index(cbuf,'</OSV>')/=0) then
                    osv_detected = .false.
                    if (npoe < 2) then
                        ! Make sure that duplicate epoch entries are ignored
                        npoe = npoe + 1
                        call reduction_model%earthFixed2inertial(state%r, state%v, state%epoch%mjd ,stateGCRF(npoe)%r, stateGCRF(npoe)%v)
                        stateGCRF(npoe)%epoch = state%epoch
                    else if (state%epoch%mjd > stateGCRF(npoe)%epoch%mjd) then
                        ! Make sure that duplicate epoch entries are ignored
                        npoe = npoe + 1
                        call reduction_model%earthFixed2inertial(state%r, state%v, state%epoch%mjd ,stateGCRF(npoe)%r, stateGCRF(npoe)%v)
                        stateGCRF(npoe)%epoch = state%epoch
                    end if
                endif

                if (osv_detected .eqv. .true.) then
                    if(index(cbuf,'<X unit="m">')/=0 ) then
                        start_nbr = index(cbuf,'<X unit="m">',.true.)+12
                        stop_nbr = index(cbuf,'</X>',.true.) -1
                        read(cbuf(start_nbr:stop_nbr),"(1f20.6)") dummy
                        read(cbuf(start_nbr:stop_nbr),*) x
                        state%r(1) = x/1000.0
                                  
                    elseif(index(cbuf,'<Y unit="m">')/=0 ) then
                        start_nbr = index(cbuf,'<Y unit="m">',.true.)+12
                        stop_nbr = index(cbuf,'</Y>',.true.) -1
                        read(cbuf(start_nbr:stop_nbr),"(1f20.6)") dummy
                        read(cbuf(start_nbr:stop_nbr),*) y 
                        state%r(2) = y/1000.0
                            
                    elseif(index(cbuf,'<Z unit="m">')/=0 ) then
                            
                        start_nbr = index(cbuf,'<Z unit="m">',.true.)+12
                        stop_nbr = index(cbuf,'</Z>',.true.) -1
                        read(cbuf(start_nbr:stop_nbr),"(1f20.6)") dummy
                        read(cbuf(start_nbr:stop_nbr),*) z
                        state%r(3) = z/1000.0 

                    elseif(index(cbuf,'<VX unit="m/s">')/=0 ) then
                            
                        start_nbr = index(cbuf,'<VX unit="m/s">',.true.)+15
                        stop_nbr = index(cbuf,'</VX>',.true.) -1
                        read(cbuf(start_nbr:stop_nbr),"(1f20.6)") dummy
                        read(cbuf(start_nbr:stop_nbr),*) vx
                        state%v(1) = vx/1000.0
                            
                    elseif(index(cbuf,'<VY unit="m/s">')/=0 ) then
                        
                        start_nbr = index(cbuf,'<VY unit="m/s">',.true.)+15
                        stop_nbr = index(cbuf,'</VY>',.true.) -1
                        read(cbuf(start_nbr:stop_nbr),"(1f20.6)") dummy
                        read(cbuf(start_nbr:stop_nbr),*) vy
                        state%v(2) = vy/1000.0

                    elseif(index(cbuf,'<VZ unit="m/s">')/=0 ) then
                        
                        start_nbr = index(cbuf,'<VZ unit="m/s">',.true.)+15
                        stop_nbr = index(cbuf,'</VZ>',.true.)-1
                        read(cbuf(start_nbr:stop_nbr),"(1f20.6)") dummy
                        read(cbuf(start_nbr:stop_nbr),*) vz
                        state%v(3) = vz/1000

                    elseif(index(cbuf,'<UTC>')/=0 ) then
                            
                        start_nbr = index(cbuf,'<UTC>',.true.)+9
                        stop_nbr = index(cbuf,'</UTC>',.true.)-1
                        read(cbuf(start_nbr:stop_nbr),*) string_epoch
                        read(cbuf(start_nbr:stop_nbr),*) string_epoch
                        read(string_epoch(1:4),*)   state%epoch%year  
                        read(string_epoch(6:7),*)   state%epoch%month
                        read(string_epoch(9:10),*)  state%epoch%day
                        read(string_epoch(12:13),*) state%epoch%hour
                        read(string_epoch(15:16),*) state%epoch%minute
                        read(string_epoch(18:19),*) state%epoch%second
                        call gd2mjd(state%epoch)

                    elseif(index(cbuf,'<Quality>')/=0 ) then
                            
                        start_nbr = index(cbuf,'<Quality>',.true.)+9
                        stop_nbr = index(cbuf,'</Quality>',.true.)-1
                        read(cbuf(start_nbr:stop_nbr),*) cbuf2
                    endif
                endif
 
                if (index(cbuf,'</List_of_OSVs>')/=0) then
                    exit
                endif
            else
                exit
            end if
        end do 

        close(ich)

    end do

    return

end subroutine readPOE
