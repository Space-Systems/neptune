!>-----------------------------------------------------------------------------------------------
!!
!! @brief       Reading input file
!!
!! @author      Vitali Braun
!!
!! @date        <ul>
!!                <li>06.02.2014 initial implementation</li>
!!              </ul>
!!
!! @param[out]  runid       run ID is returned for further file operations
!! @param[out]  epoch       begin and end epoch
!!
!> @param[inout] neptune_instance - NEPTUNE class instance to be initialized
!!
!> @param[out]  runid       run ID is returned for further file operations
!!
!> @param[out]  epoch       begin and end epoch
!!
!> @anchor      rdinp
!!
!!------------------------------------------------------------------------------------------------
subroutine rdinp(                   &
                  neptune_instance, & ! <-> Neptune_class
                  runid,            & ! --> CHR       run ID
                  epoch             & ! --> TYP       begin and end epoch
                )

    use slam_io
    use neptuneClass,       only: Neptune_class
    use libneptune,         only: init_neptune
    use slam_orbit_types,   only: state_t
    use slam_time,          only: time_t, tokenizeDate

    type(Neptune_class),intent(inout)           :: neptune_instance
    character(len=*),         intent(inout)     :: runid
    type(time_t), dimension(2), intent(out)     :: epoch

    character(len=255) :: cbuf
    character(len=10)  :: ccdrag
    character(len=10)  :: ccrefl
    character(len=10)  :: ccrsdrag
    character(len=10)  :: ccrssrp
    character(len=30)  :: cepoch_start, cepoch_end, cstart, cend
    character(len=1)   :: cGeoModel
    character(len=10)  :: cmass
    character(len=20)  :: cRelTolerance
    character(len=20)  :: cAbsTolerance
    character(len=3), dimension(10) :: pert

    integer :: ierr
    integer :: i

    !====================================================
    !
    ! Read input file
    !
    !----------------------------------
    open(unit=21, file="input/valsent.inp")

    !** run ID
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) runid

    !** Mass
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cmass

    !** Cross-section drag
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) ccrsdrag

    !** Cross-section SRP
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) ccrssrp

    !** Drag coefficient
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) ccdrag

    !** SRP coefficient
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) ccrefl

    !** Initialization epochs
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cepoch_start

    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cepoch_end

    !** Propagation epochs
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cstart

    call tokenizeDate(trim(cstart), epoch(1))

    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cend

    call tokenizeDate(trim(cend), epoch(2))

    !** Perturbations
    do i=1,size(pert)

    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) pert(i)

    end do

    !** geopotential model to be used
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cGeoModel

    !** relative tolerance
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cRelTolerance

    !** absolute tolerance
    call nxtbuf('#', 0, 21, cbuf)
    read(cbuf,*) cAbsTolerance

    close(21)

    !====================================================
    !
    ! NEPTUNE initialization
    !
    !------------------------------------------

    !** Run ID
    ierr = neptune_instance%setNeptuneVar("RUN_ID", trim(adjustl(runid)))

    !** Perturbations switches
    ierr = neptune_instance%setNeptuneVar("GEOPOTENTIAL",    trim(pert(1)))
    ierr = neptune_instance%setNeptuneVar("ATMOSPHERE",      trim(pert(2)))
    ierr = neptune_instance%setNeptuneVar("HORIZONTAL_WIND", trim(pert(3)))
    ierr = neptune_instance%setNeptuneVar("SUN",             trim(pert(4)))
    ierr = neptune_instance%setNeptuneVar("MOON",            trim(pert(5)))
    ierr = neptune_instance%setNeptuneVar("SRP",             trim(pert(6)))
    ierr = neptune_instance%setNeptuneVar("SOLID_TIDES",     trim(pert(7)))
    ierr = neptune_instance%setNeptuneVar("OCEAN_TIDES",     trim(pert(8)))
    ierr = neptune_instance%setNeptuneVar("MANEUVERS",       trim(pert(9)))
    ierr = neptune_instance%setNeptuneVar("ALBEDO",          trim(pert(10)))

    !** Options
    ierr = neptune_instance%setNeptuneVar("OPT_GEO_MODEL",      trim(cGeoModel))
    ierr = neptune_instance%setNeptuneVar("OPT_AP_FORECAST",    "15")
    ierr = neptune_instance%setNeptuneVar("OPT_SAT_PROPERTIES", "1")
    ierr = neptune_instance%setNeptuneVar("OPT_EOP",            "ON")
    ierr = neptune_instance%setNeptuneVar("OPT_INT_LOGFILE",    "OFF")

    ierr = neptune_instance%setNeptuneVar("FILE_MANEUVERS", "sentinel.mnv")

    !** Parameters
    ierr = neptune_instance%setNeptuneVar("PAR_MASS",     trim(cmass))
    ierr = neptune_instance%setNeptuneVar("PAR_CROSS_SECTION", trim(ccrsdrag))
    ierr = neptune_instance%setNeptuneVar("PAR_CDRAG",    trim(ccdrag))
    ierr = neptune_instance%setNeptuneVar("PAR_CREFL",    trim(ccrefl))

    !** Start and end epoch
    ierr = neptune_instance%setNeptuneVar("EPOCH_START_GD", trim(cepoch_start))
    ierr = neptune_instance%setNeptuneVar("EPOCH_END_GD",   trim(cepoch_end))

    !** Output options
    ierr = neptune_instance%setNeptuneVar("OUTPUT_CSV", "ON")
    ierr = neptune_instance%setNeptuneVar("OUTPUT_ATM", "ON")
    ierr = neptune_instance%setNeptuneVar("OUTPUT_ACR", "ON")
    ierr = neptune_instance%setNeptuneVar("OUTPUT_OSC", "ON")
    !ierr = neptune_instance%setNeptuneVar("OUTPUT_STEP","300")
    ierr = neptune_instance%setNeptuneVar("OPT_STORE_DATA", "10")

    !** relative tolerance
    ierr = neptune_instance%setNeptuneVar("PAR_INT_RELEPS", cRelTolerance)

    !** absolute tolerance
    ierr = neptune_instance%setNeptuneVar("PAR_INT_ABSEPS", cAbsTolerance)

    ierr = init_neptune(neptune_instance)

    return

end subroutine rdinp
