#include "../../libslam/src/oop/enter_exit_macros.inc"
subroutine OPI_Plugin_info(info) bind(c, name="OPI_Plugin_info")
  use OPI
  use ISO_C_BINDING

  implicit none

  type(c_ptr), value :: info
  ! init plugin variables
  call OPI_PluginInfo_init(info, OPI_API_VERSION_MAJOR, OPI_API_VERSION_MINOR, &
    0, 1, 0, OPI_PROPAGATOR_PLUGIN)
  call OPI_PluginInfo_setName(info, "NEPTUNE")
  call OPI_PluginInfo_setAuthor(info, "Vitali Braun @ Institute of Space Systems / TU Braunschweig")
  call OPI_PluginInfo_setDescription(info, "NPI Ephemeris Propagation Tool with Uncertainty Extrapolation")
end subroutine

module mem_copy

  interface
    subroutine memcpy(dest, src, n)bind(C,name='memcpy')
      use iso_c_binding
      type(c_ptr), value        :: dest
      type(c_ptr), value        :: src
      integer(c_size_t),value   :: n
    end subroutine memcpy
  end interface
end module

subroutine OPI_Plugin_init(propagator) bind(c, name="OPI_Plugin_init")

  use OPI
  use ISO_C_BINDING
  use slam_strings, only: toString

  implicit none

  type(c_ptr), value :: propagator
  integer            :: i

  !** perturbations
  call OPI_Module_createProperty(propagator, "geopotential_degree_order", "6")
  call OPI_Module_createProperty(propagator, "atmosphere", "ON")
  call OPI_Module_createProperty(propagator, "horizontal_wind", "OFF")
  call OPI_Module_createProperty(propagator, "solar", "ON")
  call OPI_Module_createProperty(propagator, "moon", "ON")
  call OPI_Module_createProperty(propagator, "srp", "ON")
  call OPI_Module_createProperty(propagator, "albedo", "ON")
  call OPI_Module_createProperty(propagator, "solid_tides", "ON")
  call OPI_Module_createProperty(propagator, "ocean_tides", "OFF")
  call OPI_Module_createProperty(propagator, "manoeuvres", "OFF") !** maneuvers are now "kind of" supported via opi

  !** Neptune errors
  call OPI_Module_createProperty(propagator, "eps_relative", "1.d-10")
  call OPI_Module_createProperty(propagator, "eps_absolute", "1.d-11")

  !** Neptune options
  call OPI_Module_createProperty(propagator, "harmonics", "OFF")
  call OPI_Module_createProperty(propagator, "shadow_model", "conical")
  call OPI_Module_createProperty(propagator, "shadow_boundary_correction", "OFF")
  call OPI_Module_createProperty(propagator, "srp_correct", "OFF")
  call OPI_Module_createProperty(propagator, "solar_activity_file", "sw19571001.txt")
  call OPI_Module_createProperty(propagator, "data_folder", "neptune_data")

  call OPI_Module_createProperty(propagator, "geo_model", "3")                  ! EIGEN-GL04C model
  call OPI_Module_createProperty(propagator, "atm_model", "2")                  ! NRLMSISE-00, 1 is exponential
  call OPI_Module_createProperty(propagator, "const_ap_value", "15")
  call OPI_Module_createProperty(propagator, "const_f107_value", "150.0")


  call OPI_Module_createProperty(propagator, "covariance_propagation_method", "2")
  call OPI_Module_createProperty(propagator, "covariance_propagation_step", "60.d0")
  call OPI_Module_createProperty(propagator, "covariance_prop", "OFF")          ! Enable/Disable covariance propagation
  call OPI_Module_createProperty(propagator, "covariance_prop_geopotential_degree_order", "2")
  call OPI_Module_createProperty(propagator, "covariance_prop_drag", "OFF")
  call OPI_Module_createProperty(propagator, "covariance_prop_sun", "OFF")
  call OPI_Module_createProperty(propagator, "covariance_prop_moon", "OFF")
  call OPI_Module_createProperty(propagator, "covariance_prop_srp", "OFF")
  call OPI_Module_createProperty(propagator, "covariance_prop_correlation_matrix", "OFF")
  call OPI_Module_createProperty(propagator, "covariance_ref_frame", "UVW")

  call OPI_Module_createProperty(propagator, "output_files", "OFF")
  call OPI_Module_createProperty(propagator, "store_data", "0")

  call OPI_Module_createProperty(propagator, "reentry_altitude", "50.d0")      ! km above the surface

  call OPI_Module_createProperty(propagator, "start_epoch_gd", "2009-01-01T00:00:00Z") ! start epoch, TBD
  call OPI_Module_createProperty(propagator, "end_epoch_gd", "2033-01-01T00:00:00Z") ! end epoch, TBD

  call OPI_Module_createProperty(propagator, "sat_prop", "OFF")
  call OPI_Module_createProperty(propagator, "eop", "ON")
  call OPI_Module_createProperty(propagator, "use_pn_lookup_tables", "ON")
  call OPI_Module_createProperty(propagator, "int_logfile", "OFF")

  call OPI_Module_createProperty(propagator, "cheby_polys", "OFF")
  call OPI_Module_createProperty(propagator, "cheby_degree", "36")              ! 36 is max
  call OPI_Module_createProperty(propagator, "keplerian_elements_in", "ON")
  call OPI_Module_createProperty(propagator, "ecef_states_in", "OFF")
  call OPI_Module_createProperty(propagator, "mean_elements_in", "ON")
  call OPI_Module_createProperty(propagator, "keplerian_elements_out", "ON")
  call OPI_Module_createProperty(propagator, "ecef_states_out", "OFF")
  call OPI_Module_createProperty(propagator, "mean_elements_out", "ON")
  call OPI_Module_createProperty(propagator, "return_set_matrix", "OFF")    ! instead of covariance

  do i = 1, 20
    !** manoeuvre --> should be temp, this is a little bit an overkill
    call OPI_Module_createProperty(propagator, "man_mjd_ignition_"//toString(i), "0.0")
    call OPI_Module_createProperty(propagator, "man_duration_"//toString(i), "0.0")
    call OPI_Module_createProperty(propagator, "man_ref_frame_"//toString(i), "UVW")
    call OPI_Module_createProperty(propagator, "man_a1_"//toString(i), "0.0")
    call OPI_Module_createProperty(propagator, "man_a2_"//toString(i), "0.0")
    call OPI_Module_createProperty(propagator, "man_a3_"//toString(i), "0.0")
    call OPI_Module_createProperty(propagator, "man_thrust_uncertainty_"//toString(i), "0.0")
    call OPI_Module_createProperty(propagator, "man_thrust_pointing_uncertainty_"//toString(i), "0.0")
  end do

end subroutine

function OPI_Plugin_propagate(propagator, data, julian_day, dt) result(opi_error_code) bind(c, name="OPI_Plugin_propagate")

    use mem_copy
    use OPI
    use OPI_Types
    use ISO_C_BINDING
    use slam_math
    use slam_io
    use slam_orbit_types
    use slam_interpolation
    use slam_astro_conversions, only: rv2coe,coe2rv,mean2true,true2mean
    use slam_reduction_class,   only: Reduction_type
    use maneuvers,              only: maneuver_t
    use libneptune,             only: init_neptune,propagate
    use neptuneClass,           only: Neptune_class
    use neptune_error_handling
    use slam_rframes,           only: getFrameId
    use slam_units,             only: UNIT_KM, UNIT_KMPS
    use slam_time,              only: jd2gd, jd2mjd, date2longstring, sec_per_day
    use slam_linAlgebra,        only: invert_matrix
    use slam_tehl_class

    implicit none

    !function properties
    integer(c_int)        :: opi_error_code
    type(c_ptr), value    :: data
    type(c_ptr), value    :: propagator
    real(c_double), value :: dt
    real(c_double), value :: julian_day
    integer               :: ierr                                               ! error flag
    logical, save         :: is_initialized = .false.

    !** for the perturbations
    character(len=5),dimension(10)  :: cpert                                    !** to set on/off perturbations
    character(len=50)               :: temp_string                              !** temp string to switch on/off all other options

    character(len=*), parameter     :: csubid = 'OPI_Plugin_propagate'


    !** locals values needed to set up NEPTUNE via the API
    !** some locals for neptune
    character(len=15),dimension(10),parameter :: perturbationValue = (/ "GEOPOTENTIAL   ", &
                                                                        "ATMOSPHERE     ", &
                                                                        "HORIZONTAL_WIND", &
                                                                        "SUN            ", &
                                                                        "MOON           ", &
                                                                        "SRP            ", &
                                                                        "ALBEDO         ", &
                                                                        "SOLID_TIDES    ", &
                                                                        "OCEAN_TIDES    ", &
                                                                        "MANEUVERS      "/)

    !** all those values to enable propagation
    integer(c_int)                                          :: data_size
    type(OPI_Orbit), dimension(:), pointer                  :: object_orbit
    type(OPI_ObjectProperties), dimension(:), pointer       :: object_properties
    type(OPI_Vector3), dimension(:), pointer                :: object_position
    type(OPI_Vector3), dimension(:), pointer                :: object_velocity
    type(OPI_Covariance), dimension(:), pointer             :: object_covariance

    !** states
    type(state_t)               :: initial_state                                ! initial state of the object to be propagated (ECI)
    type(state_t)               :: ecef_state                                   ! state vector (ECEF)
    type(state_t)               :: propagated_state                             ! the propagated state of the object (ECEF)
    type(covariance_t)          :: initial_covariance                           ! Assumed covariance
    type(covariance_t)          :: propagated_covariance                        ! Propagated covariance
    real(dp),dimension(6,6)     :: covariance_matrix_ECI                        ! Covariance matrix in ECI
    real(dp),dimension(6,6)     :: covariance_matrix_UVW                        ! Covariance matrix in UVW
    real(dp),dimension(6,6)     :: jacobi,invJacobi                             ! Covarianice matrix transformation from ECI to UVW and inverted

    !** times
    type(time_t),dimension(2)   :: propagation_epoch                            ! This array holds the start and end epochs as an array

    real(dp)                    :: EA,TA                                        ! For conversion from MA->TA
    real(dp)                    :: dummy                                        ! For conversion from TA->MA
    integer                     :: idummy                                       ! For conversion from TA->MA
    integer                     :: current_error
    integer                     :: neptune_iteration
    real(c_double)              :: temp_real
    real(c_double), save        :: absolute_eps
    real(c_double), save        :: relative_eps

    logical,save                :: use_kepler_elements_in
    logical,save                :: use_kepler_elements_out
    logical,save                :: use_ecef_states_in
    logical,save                :: use_ecef_states_out
    logical,save                :: use_mean_elements_in
    logical,save                :: use_mean_elements_out
    logical,save                :: propagate_covariance
    logical,save                :: propagate_maneuvers

    logical,save                :: create_cheby !** switch to create chebies or not
    integer(c_int),save         :: cheby_degree
    logical,save                :: store_data

    logical,save                :: return_set_matrix    ! switch to return set instead of covariance



    ! trying to access and write into the bytes array
    type(OPI_Char), dimension(:), pointer                   :: bytes_pointer
    real(c_double), target, allocatable, dimension(:,:,:)   :: cheby_polynomial

    integer :: iobject
    integer :: i, j, l

    !** for interpolation
    integer(c_int)                          :: number_of_states
    real(dp), dimension(:,:), allocatable   :: interpol_state
    real(dp), dimension(:),   allocatable   :: interpol_epoch
    real(dp), dimension(:,:), allocatable   :: interpol_ephemeris
    type(state_t)                           :: temp_state
    type(covariance_t)                      :: temp_covariance

    type(Neptune_class), save               :: neptune_instance
    type(Reduction_type), save              :: reduction_model

    logical :: slam_error
    real(c_double), target, allocatable, dimension(:,:,:)   :: ephemeris
    integer :: output_step_size
    integer :: temp_integer

    real(dp) :: temp_fortran_real

    integer :: number_of_maneuvers
    integer, dimension(:), allocatable :: phases_array
    type(maneuver_t), dimension(:), allocatable :: maneuvers
    type(time_t) :: temp_date
    real(dp), dimension(20,7) :: temp_maneuver_array !** contains the temp values for maneuvers: First the number of maneuver, then
                                                    !** temp_maneuver_array(,1): man_mjd_ignition
                                                    !** temp_maneuver_array(,2): man_duration
                                                    !** temp_maneuver_array(,3): man_a1
                                                    !** temp_maneuver_array(,4): man_a2
                                                    !** temp_maneuver_array(,5): man_a3
                                                    !** temp_maneuver_array(,6): man_thrust_uncertainty --> NOT USED YET
                                                    !** temp_maneuver_array(,7): man_thrust_pointing_uncertainty --> NOT USED YET

    !** starting off
    ENTER_STATIC_PROCEDURE("NEPTUNE_Propagatorschnittstelle")

    opi_error_code = 0

    !** init neptune
    if (.not. is_initialized) then

        !===================================
        !
        ! Instantiate reduction model and NEPTUNE
        !
        !---------------------------------
        reduction_model = Reduction_type()
        neptune_instance = Neptune_class()

        !** Set the data folder, which has been defined
        write(temp_string,*) OPI_Module_getPropertyString(propagator,"data_folder")
        call neptune_instance%setDataPath(trim(temp_string))

        !** set error handling to return, in case we need to change the accuracy
        call initErrorHandler(control = "YES", errAction = "RETURN", traceback = "YES")

        is_initialized = .true.
    end if

    !** start with the perturbations. All to zero
    write(cpert(1),*) OPI_Module_getPropertyString(propagator,"geopotential_degree_order")
    write(cpert(2),*) OPI_Module_getPropertyString(propagator,"atmosphere")
    write(cpert(3),*) OPI_Module_getPropertyString(propagator,"horizontal_wind")
    write(cpert(4),*) OPI_Module_getPropertyString(propagator,"solar")
    write(cpert(5),*) OPI_Module_getPropertyString(propagator,"moon")
    write(cpert(6),*) OPI_Module_getPropertyString(propagator,"srp")
    write(cpert(7),*) OPI_Module_getPropertyString(propagator,"albedo")
    write(cpert(8),*) OPI_Module_getPropertyString(propagator,"solid_tides")
    write(cpert(9),*) OPI_Module_getPropertyString(propagator,"ocean_tides")
    write(cpert(10),*) OPI_Module_getPropertyString(propagator,"manoeuvres")

    !** set maneuver switch
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"manoeuvres")
    if (index(temp_string,"OFF") /= 0) then
        propagate_maneuvers = .false.
    else
        propagate_maneuvers = .true.
    end if

    !** move values to neptune
    do i = 1, 10
        ierr = neptune_instance%setNeptuneVar(perturbationValue(i),TRIM(cpert(i)))

        !** check, if something went wrong in libslam/neptune
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

    end do

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"harmonics")
    ierr = neptune_instance%setNeptuneVar("OPT_HARMONICS", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif


    write(temp_string,*) OPI_Module_getPropertyString(propagator,"eps_relative")
    read(temp_string,*) relative_eps
    ierr = neptune_instance%setNeptuneVar("PAR_INT_RELEPS", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"eps_absolute")
    read(temp_string,*) absolute_eps
    ierr = neptune_instance%setNeptuneVar("PAR_INT_ABSEPS", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif


    write(temp_string,*) OPI_Module_getPropertyString(propagator,"shadow_model")
    ierr = neptune_instance%setNeptuneVar("OPT_SHADOW", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        write(*,*) "Error in Neptune during neptune_instance%setNeptuneVar (OPT_SHADOW): "//toString(ierr)
        stop
    end if

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"srp_correct")
    ierr = neptune_instance%setNeptuneVar("OPT_SRP_CORRECT", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"solar_activity_file")
    ierr = neptune_instance%setNeptuneVar("FILE_SOLMAG", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"geo_model")
    ierr = neptune_instance%setNeptuneVar("OPT_GEO_MODEL", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    !** set the atm model to be used
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"atm_model")
    ierr = neptune_instance%setNeptuneVar("OPT_ATMOSPHERE_MODEL", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    !** constant ap value (long term propagation)
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"const_ap_value")
    ierr = neptune_instance%setNeptuneVar("OPT_AP_FORECAST", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    !** constant f10.7 value (long term propagation)
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"const_f107_value")
    ierr = neptune_instance%setNeptuneVar("OPT_SOL_FORECAST", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    !** Enable/Disable covariance propagation
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop")
    if (index(temp_string,"OFF") /= 0) then
        propagate_covariance = .false.
    else
        propagate_covariance = .true.
    end if
    ierr = neptune_instance%setNeptuneVar("COVARIANCE_PROPAGATION", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif


    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop_geopotential_degree_order")
    ierr = neptune_instance%setNeptuneVar("COVARIANCE_GEOPOTENTIAL",trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop_drag")
    ierr = neptune_instance%setNeptuneVar("COVARIANCE_DRAG",trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop_moon")
    ierr = neptune_instance%setNeptuneVar("COVARIANCE_MOON",trim(temp_string))
    !!** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop_sun")
    ierr = neptune_instance%setNeptuneVar("COVARIANCE_SUN",trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop_srp")
    ierr = neptune_instance%setNeptuneVar("COVARIANCE_SRP",trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_prop_correlation_matrix")
    ierr = neptune_instance%setNeptuneVar("CORRELATION_MATRIX",trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_propagation_step")
    ierr = neptune_instance%setNeptuneVar("PAR_INT_COV_STEP", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_propagation_method")
    ierr = neptune_instance%setNeptuneVar("PAR_INT_COV_METHOD", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"output_files")
    ierr = neptune_instance%setNeptuneVar("OUTPUT_FILES", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"reentry_altitude")
    ierr = neptune_instance%setNeptuneVar("PAR_REENTRY", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"start_epoch_gd")
    ierr = neptune_instance%setNeptuneVar("EPOCH_START_GD", trim(adjustl(temp_string)))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"end_epoch_gd")
    ierr = neptune_instance%setNeptuneVar("EPOCH_END_GD", trim(adjustl(temp_string)))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"sat_prop")
    ierr = neptune_instance%setNeptuneVar("OPT_SRP_CORRECT", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"eop")
    ierr = neptune_instance%setNeptuneVar("OPT_EOP", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"use_pn_lookup_tables")
    ierr = neptune_instance%setNeptuneVar("OPT_PN_LOOKUP", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"int_logfile")
    ierr = neptune_instance%setNeptuneVar("OPT_INT_LOGFILE", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"data_folder")
    ierr = neptune_instance%setNeptuneVar("PATH_DATA", trim(temp_string))
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    !** using properties of the satellite
    ierr = neptune_instance%setNeptuneVar("OPT_SAT_PROPERTIES",   "1")
    !** check error
    if (ierr .ne. 0) then
        slam_error = t%check_slam_error()
        if (t%has_to_return()) return
        if (slam_error) then
            call resetError()
        end if
    endif

    !** check if we want to have the polynomials
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"cheby_polys")
    if (index(temp_string,"OFF") /= 0) then
        create_cheby = .false.

        !** We need to store the data, as we are VERY interested in the ephemeris. As a "first try", we use 30 seconds.
        ierr = neptune_instance%setNeptuneVar("OPT_STORE_DATA", "30")
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

    else
        create_cheby = .true.

        !** we need to store data in this case. We use every 30 seconds for now. Needs a check
        ierr = neptune_instance%setNeptuneVar("OPT_STORE_DATA", "30")
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

    end if

    !** check if we should return the set matrix instead of the covariance
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"return_set_matrix")
    if (index(temp_string,"OFF") /= 0) then
        return_set_matrix = .false.

    else
        return_set_matrix = .true.

    end if

    !** Allow saving the states, if wanted
    store_data = .false.
    if (.not. create_cheby) then
        write(temp_string,*) OPI_Module_getPropertyString(propagator,"store_data")
        ierr = neptune_instance%setNeptuneVar("OPT_STORE_DATA", trim(temp_string))
        !!** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        temp_integer = string_to_int(temp_string)

        if ( temp_integer .ne. 0 ) then

        store_data = .true.

        !** store the data
        output_step_size = string_to_int(temp_string)

        end if
    endif

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"cheby_degree")
    read (temp_string,*) cheby_degree

    !** check if we want to use keplerian elements
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"keplerian_elements_in")
    if (index(temp_string,"OFF") /= 0) then
        use_kepler_elements_in = .false.
    else
        use_kepler_elements_in = .true.
    end if

    write(temp_string,*) OPI_Module_getPropertyString(propagator,"keplerian_elements_out")
    if (index(temp_string,"OFF") /= 0) then
        use_kepler_elements_out = .false.
    else
        use_kepler_elements_out = .true.
    end if

    !** check if we want to use ecef states
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"ecef_states_in")
    if (index(temp_string,"OFF") /= 0) then
        use_ecef_states_in = .false.
    else
        use_ecef_states_in = .true.
    end if

    !** check if we want to use ecef states
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"ecef_states_out")
    if (index(temp_string,"OFF") /= 0) then
        use_ecef_states_out = .false.
    else
        use_ecef_states_out = .true.
    end if

    !** check if we want to mean elements as input
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"mean_elements_in")
    if (index(temp_string,"OFF") /= 0) then
        use_mean_elements_in = .false.
    else
        use_mean_elements_in = .true.
    end if

    !** check if we want to use mean elements as output
    write(temp_string,*) OPI_Module_getPropertyString(propagator,"mean_elements_out")
    if (index(temp_string,"OFF") /= 0) then
        use_mean_elements_out = .false.
    else
        use_mean_elements_out = .true.
    end if

    !** have fun with maneuvers
    if (propagate_maneuvers) then
        number_of_maneuvers = 0
        do i = 1, 20

            write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_duration_"//toString(i))
            read (temp_string,*) temp_fortran_real
            ! call message(toString(temp_fortran_real), LOG_AND_STDOUT)
            if (temp_fortran_real .gt. 0.d0) then
                number_of_maneuvers = number_of_maneuvers + 1
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_mjd_ignition_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,1) = temp_fortran_real
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_duration_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,2) = temp_fortran_real
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_a1_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,3) = temp_fortran_real
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_a2_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,4) = temp_fortran_real
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_a3_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,5) = temp_fortran_real
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_thrust_uncertainty_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,6) = temp_fortran_real
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"man_thrust_pointing_uncertainty_"//toString(i))
                read (temp_string,*) temp_fortran_real
                temp_maneuver_array(number_of_maneuvers,7) = temp_fortran_real
            end if
        end do

        !** create the maneuver array
        if (allocated(maneuvers)) deallocate(maneuvers)
        allocate(maneuvers(number_of_maneuvers))
        if (allocated(phases_array)) deallocate(phases_array)
        allocate(phases_array(number_of_maneuvers))
        phases_array(:) = 1
        do i = 1, number_of_maneuvers

            temp_date%mjd = temp_maneuver_array(i,1)
            call mjd2gd(temp_date)
            maneuvers(i)%start_date = temp_date

            temp_date%mjd = temp_maneuver_array(i,1) + temp_maneuver_array(i,2)/86400.d0
            call mjd2gd(temp_date)
            maneuvers(i)%end_date = temp_date

            maneuvers(i)%nphases = 1 !** just one phase now!!!!

            !** allocate that single phase
            allocate(maneuvers(i)%phase(1))

            !** set the alues for that phase
            maneuvers(i)%phase(1)%mjd_start = temp_maneuver_array(i,1)
            maneuvers(i)%phase(1)%mjd_end = temp_maneuver_array(i,1) + temp_maneuver_array(i,2)/86400.d0
            maneuvers(i)%phase(1)%acc(:) = temp_maneuver_array(i,3:5)
            maneuvers(i)%phase(1)%thrust_efficiency = abs(1.d0 - temp_maneuver_array(i,6))

            !** dunno how to use the uncertainties yet

        end do
        call neptune_instance%manoeuvres_model%init_maneuvers(maneuvers, phases_array)

    end if

    data_size = OPI_Population_getSize(data)
    object_orbit => OPI_Population_getOrbit(data)
    object_position => OPI_Population_getPosition(data)
    object_velocity => OPI_Population_getVelocity(data)
    object_covariance => OPI_Population_getCovariance(data)
    object_properties => OPI_Population_getObjectProperties(data)

    if (create_cheby) then
        bytes_pointer => OPI_Population_getBytes(data)
        !** allocate the polynomial
        allocate(cheby_polynomial(1:6,0:cheby_degree,1:data_size))
    end if

    !** loop over the objects and propagate object per object

    !** set the parameters that are needed for every object
    initial_state%epoch%jd = julian_day

    !** update everything in the epoch
    call jd2gd(initial_state%epoch)
    !write (*,*) date2longstring(initial_state%epoch)
    call jd2mjd(initial_state%epoch)
    propagation_epoch(1) = initial_state%epoch
    propagation_epoch(2)%jd = initial_state%epoch%jd + dt

    call jd2gd(propagation_epoch(2))
    call jd2mjd(propagation_epoch(2))
    !write (*,*) date2longstring(propagation_epoch(2))


    initial_state%frame = getFrameId("GCRF")
    initial_state%radius_unit = UNIT_KM
    initial_state%velocity_unit = UNIT_KMPS

    !** covariances are all zero :D
    initial_covariance%elem(:,:) = 0.d0
    initial_covariance%epoch = initial_state%epoch
    initial_covariance%frame = getFrameId("GCRF")
    initial_covariance%unit = UNIT_KM

    do iobject = 1, data_size

        !** set the object properties
        ierr = neptune_instance%setNeptuneVar("PAR_MASS", tostring(object_properties(iobject)%mass))
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        ierr = neptune_instance%setNeptuneVar("PAR_CDRAG", tostring(object_properties(iobject)%drag_coefficient))
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        ierr = neptune_instance%setNeptuneVar("PAR_CROSS_SECTION", tostring(object_properties(iobject)%area_to_mass * object_properties(iobject)%mass))
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        ierr = neptune_instance%setNeptuneVar("PAR_CREFL", tostring(object_properties(iobject)%reflectivity))
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        write(temp_string,*) OPI_Module_getPropertyString(propagator,"eps_relative")
        read(temp_string,*) relative_eps
        ierr = neptune_instance%setNeptuneVar("PAR_INT_RELEPS", trim(temp_string))
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        write(temp_string,*) OPI_Module_getPropertyString(propagator,"eps_absolute")
        read(temp_string,*) absolute_eps
        ierr = neptune_instance%setNeptuneVar("PAR_INT_ABSEPS", trim(temp_string))
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif

        !** init neptune
        ierr = init_neptune(neptune_instance)
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        endif


        if (use_kepler_elements_in) then
            if (use_mean_elements_in) then
                ! Convert from mean to eccentric and true anomaly
                call mean2true( object_orbit(iobject)%eccentricity,     &
                                object_orbit(iobject)%mean_anomaly,     &
                                EA,                                     &
                                TA)
            else
                TA = object_orbit(iobject)%mean_anomaly
            end if

            ! Convert from keplerian elements to cartesian state vector
            call coe2rv(object_orbit(iobject)%semi_major_axis                   &
                            * (1.0d0 - object_orbit(iobject)%eccentricity**2),  &
                        object_orbit(iobject)%eccentricity,                     &
                        object_orbit(iobject)%inclination,                      &
                        object_orbit(iobject)%raan,                             &
                        object_orbit(iobject)%arg_of_perigee,                   &
                        TA,                                                     &
                        TA + object_orbit(iobject)%arg_of_perigee,              &
                        TA + object_orbit(iobject)%arg_of_perigee               &
                                + object_orbit(iobject)%raan,                   &
                        object_orbit(iobject)%raan                              &
                                + object_orbit(iobject)%arg_of_perigee,         &
                        initial_state%r,                                        &
                        initial_state%v)
        else ! states are assumed
            if (use_ecef_states_in) then
                ! convert from ecef to eci
                ecef_state%r(1) = object_position(iobject)%x
                ecef_state%r(2) = object_position(iobject)%y
                ecef_state%r(3) = object_position(iobject)%z
                ecef_state%v(1) = object_velocity(iobject)%x
                ecef_state%v(2) = object_velocity(iobject)%y
                ecef_state%v(3) = object_velocity(iobject)%z
                call reduction_model%earthFixed2inertial(   ecef_state%r,               &
                                                            ecef_state%v,               &
                                                            initial_state%epoch%mjd,    &
                                                            initial_state%r,            &
                                                            initial_state%v)
                !** check error
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if

            else
                ! set the ECI state
                initial_state%r(1) = object_position(iobject)%x
                initial_state%r(2) = object_position(iobject)%y
                initial_state%r(3) = object_position(iobject)%z
                initial_state%v(1) = object_velocity(iobject)%x
                initial_state%v(2) = object_velocity(iobject)%y
                initial_state%v(3) = object_velocity(iobject)%z
            end if
        end if

        if (propagate_covariance) then
            initial_covariance%elem(1,1) = object_covariance(iobject)%k1_k1

            initial_covariance%elem(2,1) = object_covariance(iobject)%k2_k1
            initial_covariance%elem(1,2) = object_covariance(iobject)%k2_k1
            initial_covariance%elem(2,2) = object_covariance(iobject)%k2_k2

            initial_covariance%elem(3,1) = object_covariance(iobject)%k3_k1
            initial_covariance%elem(1,3) = object_covariance(iobject)%k3_k1
            initial_covariance%elem(3,2) = object_covariance(iobject)%k3_k2
            initial_covariance%elem(2,3) = object_covariance(iobject)%k3_k2
            initial_covariance%elem(3,3) = object_covariance(iobject)%k3_k3

            initial_covariance%elem(4,1) = object_covariance(iobject)%k4_k1
            initial_covariance%elem(1,4) = object_covariance(iobject)%k4_k1
            initial_covariance%elem(4,2) = object_covariance(iobject)%k4_k2
            initial_covariance%elem(2,4) = object_covariance(iobject)%k4_k2
            initial_covariance%elem(4,3) = object_covariance(iobject)%k4_k3
            initial_covariance%elem(3,4) = object_covariance(iobject)%k4_k3
            initial_covariance%elem(4,4) = object_covariance(iobject)%k4_k4

            initial_covariance%elem(5,1) = object_covariance(iobject)%k5_k1
            initial_covariance%elem(1,5) = object_covariance(iobject)%k5_k1
            initial_covariance%elem(5,2) = object_covariance(iobject)%k5_k2
            initial_covariance%elem(2,5) = object_covariance(iobject)%k5_k2
            initial_covariance%elem(5,3) = object_covariance(iobject)%k5_k3
            initial_covariance%elem(3,5) = object_covariance(iobject)%k5_k3
            initial_covariance%elem(5,4) = object_covariance(iobject)%k5_k4
            initial_covariance%elem(4,5) = object_covariance(iobject)%k5_k4
            initial_covariance%elem(5,5) = object_covariance(iobject)%k5_k5

            initial_covariance%elem(6,1) = object_covariance(iobject)%k6_k1
            initial_covariance%elem(1,6) = object_covariance(iobject)%k6_k1
            initial_covariance%elem(6,2) = object_covariance(iobject)%k6_k2
            initial_covariance%elem(2,6) = object_covariance(iobject)%k6_k2
            initial_covariance%elem(6,3) = object_covariance(iobject)%k6_k3
            initial_covariance%elem(3,6) = object_covariance(iobject)%k6_k3
            initial_covariance%elem(6,4) = object_covariance(iobject)%k6_k4
            initial_covariance%elem(4,6) = object_covariance(iobject)%k6_k4
            initial_covariance%elem(6,5) = object_covariance(iobject)%k6_k5
            initial_covariance%elem(5,6) = object_covariance(iobject)%k6_k5
            initial_covariance%elem(6,6) = object_covariance(iobject)%k6_k6
            covariance_matrix_UVW = initial_covariance%elem

            !write(*,*) "Initial Covariance (as delivered) in km and km/s"
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,1)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,2)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,3)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,4)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,5)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,6)

            !** check, in what frame the covariance was provided
            write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_ref_frame")
            if ( trim(adjustl(temp_string)) == "UVW" ) then

                ! write(*,*) 'MY MATRIX IS UVW DUDE'

                ! Convert UVW covariance into ECI covariance
                call reduction_model%getJacobianEci2uvw(initial_state%r,initial_state%v,jacobi) ! convert to UVW
                !** check error
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if

                call invert_matrix(jacobi,invJacobi,'LU_DECOMP')
                !** check error
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if

                covariance_matrix_ECI = matmul(matmul(invJacobi,covariance_matrix_UVW),transpose(invJacobi))
                initial_covariance%elem = covariance_matrix_ECI
                !write(*,*) "Initial covariance (in RTN) in km and km/s"
            else
                !write(*,*) "Initial covariance (in GCRF) in km and km/s"
            endif

            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,1)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,2)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,3)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,4)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,5)
            !write(*,'(6(F15.8,1X))') initial_covariance%elem(:,6)
        else
            !write (99,*) "No covariance to propagate"

        end if

        current_error = E_INTEGRATION_ABORT
        neptune_iteration = 0

        !** init neptune
        ierr = init_neptune(neptune_instance)
        !** check error
        if (ierr .ne. 0) then
            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
                call resetError()
            end if
        end if

        current_error = E_INTEGRATION_ABORT
        neptune_iteration = 0

        do while ((current_error .eq. E_INTEGRATION_ABORT) .and. (neptune_iteration .le. 5))

          !** calling neptune
          call propagate(neptune_instance,      &
                        initial_state,          &
                        initial_covariance,     &
                        propagation_epoch,      &
                        propagated_state,       &
                        propagated_covariance,  &
                        .true.) !** set to true, as we always propagate different objects

          !** update error
          current_error = 0

          !** check what happened
          if (hasFailed()) then

            !** check if it was re-entry
            if (getLatestError() .eq. E_MIN_ALTITUDE) then

               !** awesome. we only reentered, therefore:
               object_orbit(iobject)%eol = propagation_epoch(1)%jd  !** <-- this is not very correct. Here actually should stand the decay date, which is not updated in NEPTUNE!

               !** reset error handling
               call resetError()

            elseif (getLatestError() .eq. E_INTEGRATION_ABORT) then

              current_error = getLatestError()

              !** reset error handling
              call resetError()

              neptune_iteration = neptune_iteration + 1

              call t%log_warning("Too many resets in neptune, trying with lower accuracy.")

              temp_real = absolute_eps*10.d0**(neptune_iteration)
              ierr = neptune_instance%setNeptuneVar("PAR_INT_ABSEPS", toString(temp_real))
              if (ierr .ne. 0) then
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if
              end if

              temp_real = relative_eps*10.d0**(neptune_iteration)
              ierr = neptune_instance%setNeptuneVar("PAR_INT_RELEPS",toString(temp_real))
              if (ierr .ne. 0) then
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if
              end if

              !** init neptune
              ierr = init_neptune(neptune_instance)
              if (ierr .ne. 0) then
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if
              end if

            else

              !call initErrorHandler(control = "YES", traceback = "YES")
              slam_error = t%check_slam_error()
              if (t%has_to_return()) return
              if (slam_error) then
                call resetError()
              end if

            endif

          else

            !** no re-entry, eol is 0.d0
            object_orbit(iobject)%eol  = 0.d0

          endif

        enddo

        if (neptune_iteration .gt. 5) then

          call t%log_fatal("Too many iterations of accuracy in Neptune.")
          if (t%has_to_return()) return

        end if

        ! return values to OPI
        if (use_ecef_states_out) then
            call reduction_model%inertial2earthFixed(propagated_state%r,        &
                                                    propagated_state%v,         &
                                                    propagation_epoch(2)%mjd,   &
                                                    ecef_state%r,               &
                                                    ecef_state%v)


            slam_error = t%check_slam_error()
            if (t%has_to_return()) return
            if (slam_error) then
              call resetError()
            end if
            object_position(iobject)%x = ecef_state%r(1)
            object_position(iobject)%y = ecef_state%r(2)
            object_position(iobject)%z = ecef_state%r(3)
            object_velocity(iobject)%x = ecef_state%v(1)
            object_velocity(iobject)%y = ecef_state%v(2)
            object_velocity(iobject)%z = ecef_state%v(3)
        else
            object_position(iobject)%x = propagated_state%r(1)
            object_position(iobject)%y = propagated_state%r(2)
            object_position(iobject)%z = propagated_state%r(3)
            object_velocity(iobject)%x = propagated_state%v(1)
            object_velocity(iobject)%y = propagated_state%v(2)
            object_velocity(iobject)%z = propagated_state%v(3)
        end if

        if (use_kepler_elements_out) then
            call rv2coe(    propagated_state%r,                     &
                            propagated_state%v,                     &
                            dummy,                                  &
                            object_orbit(iobject)%semi_major_axis,  &
                            object_orbit(iobject)%eccentricity,     &
                            object_orbit(iobject)%inclination,      &
                            object_orbit(iobject)%raan,             &
                            object_orbit(iobject)%arg_of_perigee,   &
                            TA,                                     &
                            object_orbit(iobject)%mean_anomaly,     &
                            dummy,                                  &
                            dummy,                                  &
                            dummy,                                  &
                            idummy)

             if (.not. use_mean_elements_out) then
                object_orbit(iobject)%mean_anomaly = TA
             end if
        end if

        if (propagate_covariance) then
            !** check, in what frame the covariance was provided
            write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_ref_frame")
            if ( trim(adjustl(temp_string)) == "UVW" ) then
                ! Convert UVW covariance into ECI covariance
                call reduction_model%getJacobianEci2uvw(propagated_state%r,propagated_state%v,jacobi) ! convert to UVW
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                    call resetError()
                end if
                covariance_matrix_UVW = matmul(matmul(jacobi,propagated_covariance%elem),transpose(jacobi))
                propagated_covariance%elem = covariance_matrix_UVW
            end if
            object_covariance(iobject)%k1_k1 = propagated_covariance%elem(1,1)
            object_covariance(iobject)%k2_k1 = propagated_covariance%elem(2,1)
            object_covariance(iobject)%k2_k2 = propagated_covariance%elem(2,2)
            object_covariance(iobject)%k3_k1 = propagated_covariance%elem(3,1)
            object_covariance(iobject)%k3_k2 = propagated_covariance%elem(3,2)
            object_covariance(iobject)%k3_k3 = propagated_covariance%elem(3,3)
            object_covariance(iobject)%k4_k1 = propagated_covariance%elem(4,1)
            object_covariance(iobject)%k4_k2 = propagated_covariance%elem(4,2)
            object_covariance(iobject)%k4_k3 = propagated_covariance%elem(4,3)
            object_covariance(iobject)%k4_k4 = propagated_covariance%elem(4,4)
            object_covariance(iobject)%k5_k1 = propagated_covariance%elem(5,1)
            object_covariance(iobject)%k5_k2 = propagated_covariance%elem(5,2)
            object_covariance(iobject)%k5_k3 = propagated_covariance%elem(5,3)
            object_covariance(iobject)%k5_k4 = propagated_covariance%elem(5,4)
            object_covariance(iobject)%k5_k5 = propagated_covariance%elem(5,5)
            object_covariance(iobject)%k6_k1 = propagated_covariance%elem(6,1)
            object_covariance(iobject)%k6_k2 = propagated_covariance%elem(6,2)
            object_covariance(iobject)%k6_k3 = propagated_covariance%elem(6,3)
            object_covariance(iobject)%k6_k4 = propagated_covariance%elem(6,4)
            object_covariance(iobject)%k6_k5 = propagated_covariance%elem(6,5)
            object_covariance(iobject)%k6_k6 = propagated_covariance%elem(6,6)

            !write(*,*) "New Covariance (UVW) in km and km/s"
            !write(*,'(6(F15.8,1X))') propagated_covariance%elem(:,1)
            !write(*,'(6(F15.8,1X))') propagated_covariance%elem(:,2)
            !write(*,'(6(F15.8,1X))') propagated_covariance%elem(:,3)
            !write(*,'(6(F15.8,1X))') propagated_covariance%elem(:,4)
            !write(*,'(6(F15.8,1X))') propagated_covariance%elem(:,5)
            !write(*,'(6(F15.8,1X))') propagated_covariance%elem(:,6)
        end if

        if (create_cheby) then
            ! create the actual polynomials
            ! recover the data to interpolate from
            number_of_states = neptune_instance%getNumberOfEphemerides()
            !** from neptune
            if (number_of_states > 100000) then
                call t%log_fatal("Too many ephemerides in NEPTUNE. Change span or step size.")
                if (t%has_to_return()) return
            end if


            ! allocation
            allocate(interpol_state(number_of_states,6))
            allocate(interpol_epoch(number_of_states))
            allocate(interpol_ephemeris(number_of_states,2))

            ! transfer the data to the allocated arrays
            do j = 1, number_of_states

                temp_state = neptune_instance%getNeptuneData(j)
                ! Convert ECI states to ECEF, when requested.
                ! Chebyshev polynomials will then be produced for ECEF instead of ECI states
                if (use_ecef_states_out) then
                    call reduction_model%inertial2earthFixed(temp_state%r,          &
                                                            temp_state%v,           &
                                                            temp_state%epoch%mjd,   &
                                                            interpol_state(j,1:3),  &
                                                            interpol_state(j,4:6))
                    slam_error = t%check_slam_error()
                    if (t%has_to_return()) return
                    if (slam_error) then
                      call resetError()
                    end if

                    interpol_epoch(j) = temp_state%epoch%mjd
                else
                    interpol_state(j,1:3) = temp_state%r(1:3)
                    interpol_state(j,4:6) = temp_state%v(1:3)
                    interpol_epoch(j) = temp_state%epoch%mjd

                end if

            end do

            interpol_ephemeris(:,1) = interpol_epoch(:)

            ! call the interpolation, once for every position
            if (object_orbit(iobject)%eol .eq. 0.d0) then

              do l = 1,6
                interpol_ephemeris(:,2) = interpol_state(:,l)

                call chebyshev_interpolation(cheby_degree, 5, number_of_states, interpol_ephemeris(:,1:2), cheby_polynomial(l,:,iobject))
                slam_error = t%check_slam_error()
                if (t%has_to_return()) return
                if (slam_error) then
                  call resetError()
                end if

              end do

            else

              cheby_polynomial(1:6,:,iobject) = 0.d0

            endif

            ! delete
            deallocate(interpol_state)
            deallocate(interpol_epoch)
            deallocate(interpol_ephemeris)
        end if

        !** we have to move the single states to the byte pointer thing. Exciting.
        if (store_data) then

          !** now get the single states
          number_of_states = neptune_instance%getNumberOfEphemerides()
          !** from neptune
          if (number_of_states > 1000000) then
              call t%log_fatal("Too many ephemerides in NEPTUNE. Change span or step size.")
              if (t%has_to_return()) return
          end if

          !** allocate the ephemeris array
          allocate(ephemeris(1000000, 1:28, data_size)) !** the 1000000 is the maximum allowed at the moment. It is fixed to avoid errors with PICARD. The 28 is: 1 date, 6 state, 21 co-variances
          ephemeris(:,:,:) = 0.d0

          do j = 1, number_of_states

            ! Extract state vector from NEPTUNE API
            temp_state = neptune_instance%getNeptuneData(j)
            ephemeris(j, 1, iobject) = temp_state%epoch%mjd
            ephemeris(j, 2:4, iobject) = temp_state%r(1:3)
            ephemeris(j, 5:7, iobject) = temp_state%v(1:3)
            !write (99,*) temp_state%epoch%mjd, temp_state%r(1:3), temp_state%v(1:3)

            ! Extract covariance from NEPTUNE API
            if (propagate_covariance .and. (.not. return_set_matrix)) then
                temp_covariance = neptune_instance%getNeptuneCovarianceData(j)
                !** check, in what frame the covariance was provided
                write(temp_string,*) OPI_Module_getPropertyString(propagator,"covariance_ref_frame")
                if ( trim(adjustl(temp_string)) == "UVW" ) then
                    ! Convert UVW covariance into ECI covariance
                    call reduction_model%getJacobianEci2uvw(temp_state%r,temp_state%v,jacobi) ! convert to UVW
                    slam_error = t%check_slam_error()
                    if (t%has_to_return()) return
                    if (slam_error) then
                        call resetError()
                    end if
                    covariance_matrix_UVW = matmul(matmul(jacobi,temp_covariance%elem),transpose(jacobi))
                    temp_covariance%elem = covariance_matrix_UVW
                end if
            else if(return_set_matrix .and. propagate_covariance) then
                ! Return SET matrix instead of covariance if requested
                temp_covariance = neptune_instance%getNeptuneSetMatrixData(j)
            else
                temp_covariance%elem = initial_covariance%elem
            end if

            ephemeris(j, 8, iobject) = temp_covariance%elem(1,1)
            ephemeris(j, 9, iobject) = temp_covariance%elem(2,1)
            ephemeris(j,10, iobject) = temp_covariance%elem(2,2)
            ephemeris(j,11, iobject) = temp_covariance%elem(3,1)
            ephemeris(j,12, iobject) = temp_covariance%elem(3,2)
            ephemeris(j,13, iobject) = temp_covariance%elem(3,3)
            ephemeris(j,14, iobject) = temp_covariance%elem(4,1)
            ephemeris(j,15, iobject) = temp_covariance%elem(4,2)
            ephemeris(j,16, iobject) = temp_covariance%elem(4,3)
            ephemeris(j,17, iobject) = temp_covariance%elem(4,4)
            ephemeris(j,18, iobject) = temp_covariance%elem(5,1)
            ephemeris(j,19, iobject) = temp_covariance%elem(5,2)
            ephemeris(j,20, iobject) = temp_covariance%elem(5,3)
            ephemeris(j,21, iobject) = temp_covariance%elem(5,4)
            ephemeris(j,22, iobject) = temp_covariance%elem(5,5)
            ephemeris(j,23, iobject) = temp_covariance%elem(6,1)
            ephemeris(j,24, iobject) = temp_covariance%elem(6,2)
            ephemeris(j,25, iobject) = temp_covariance%elem(6,3)
            ephemeris(j,26, iobject) = temp_covariance%elem(6,4)
            ephemeris(j,27, iobject) = temp_covariance%elem(6,5)
            ephemeris(j,28, iobject) = temp_covariance%elem(6,6)

          end do
          !write (99,*) propagation_epoch(2)%mjd, propagated_state%r(1:3), propagated_state%v(1:3)
        endif

    end do

    if (create_cheby) then

        call memcpy(c_loc(bytes_pointer(1)%bytes), c_loc(cheby_polynomial), int(8*6*(cheby_degree+1)*data_size,8)) !** 21 is for cheby that goes from 0 to 20
        deallocate(cheby_polynomial)

    end if

    if (store_data) then

        !**call message("Calling the bytes",LOG_AND_STDOUT)
        !** get the bytes pointer
        bytes_pointer => OPI_Population_getBytes(data)
        call memcpy(c_loc(bytes_pointer(1)%bytes), c_loc(ephemeris), int(8*28*1000000*data_size,8))
        deallocate(ephemeris)

    end if

    if (allocated(maneuvers)) deallocate(maneuvers)
    if (allocated(phases_array)) deallocate(phases_array)

    !** update the data for opi
    if (opi_error_code .eq. 0) then
        opi_error_code = OPI_Population_update(data, OPI_DATA_ORBIT)
    endif
    if (opi_error_code .eq. 0) then
        opi_error_code = OPI_Population_update(data, OPI_DATA_PROPERTIES)
    endif
    if (opi_error_code .eq. 0) then
        opi_error_code = OPI_Population_update(data, OPI_DATA_POSITION)
    endif
    if (opi_error_code .eq. 0) then
        opi_error_code = OPI_Population_update(data, OPI_DATA_VELOCITY)
    endif
    if (opi_error_code .eq. 0) then
        opi_error_code = OPI_Population_update(data, OPI_DATA_COVARIANCE)
    endif
    if (opi_error_code .eq. 0) then
        opi_error_code = OPI_Population_update(data, OPI_DATA_BYTES)
    endif

    !** EXIT
    EXIT_PROCEDURE()


end function
