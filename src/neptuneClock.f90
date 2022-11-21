!=====================================================================================================
!
!> @brief      NEPTUNE clock manages the integration steps
!!
!> @author     Vitali Braun (VB)
!!
!> @date       <ul>
!!                  <li>VB: 29.04.2017 (initial implementation)</li>
!!                  <li>VB: 24.12.2017 (Adding constructor)</li>
!!              </ul>
!!
!> @copyright   Institute of Space Systems / TU Braunschweig
!!
!> @anchor      neptuneClock
!
! ------------------------------------------------------------------------
module neptuneClock

    use slam_time,      only: time_t
    use slam_types,     only: dp
    use numint,         only: RK4
    use neptuneClass,   only: Neptune_class

    implicit none
    private

    type, public :: Clock_class
        real(dp) :: start_time                                                  !< start time (in seconds)
        real(dp) :: end_time                                                    !< end time (in seconds)
        real(dp) :: step_size                                                   !< step size - usually the output step size (in seconds)
        real(dp) :: cov_step                                                    !< covariance matrix integration step size (seconds)
        real(dp) :: out_counter                                                 !< counter keeping track of output writing (seconds)
        real(dp) :: cov_counter                                                 !< counter keeping track of covariance propagation (seconds)
        real(dp) :: next_step                                                   !< time for the next step the integration aims at (seconds)
        real(dp) :: last_counter_success                                        !< keeps track of the last successful step
        real(dp) :: prop_counter_reset                                          !< keeps track of the time the propagation was reset
        integer  :: cov_half_step_counter                                       !< used to count how many covariance half-steps were made. If the number is even, a full step was made - hence covariance gets updated. On half-step it just saves.
        logical  :: forward                                                     !< .true. if forward propagation, .false. for backwards
        logical  :: flag_output_step                                            !< is .true. during output steps
        logical  :: flag_cov_update                                             !< is .true. during covariance matrix update/integration steps
        logical  :: flag_cov_save                                               !< is .true. when an intermediate step for the covariance needs to be saved (like for RK4)
        logical  :: cov_propagation_flag                                        !< is .true. when covariance matrix is propagated
        real(dp),dimension(:),allocatable  :: step_epochs_sec                   !< intermediate epochs (in seconds)
        logical  :: intermediate_steps_flag                                     !< is .true. when intermediate steps are to be take (requested by user)
        integer  :: intermediate_steps_index                                    !  keeps track of in-array index whilst in intermediate mode
        integer  :: last_index
        real(dp) :: cumulated_time
    contains
        procedure :: init_counter
        procedure :: switch_backward_propagation
        procedure :: switch_forward_propagation
        procedure :: toggle_output_flag
        procedure :: toggle_cov_update_flag
        procedure :: get_next_step
        procedure :: has_finished
        procedure :: has_finished_step
        procedure :: has_to_write_output
        procedure :: has_to_update_covariance
        procedure :: has_to_save_covariance
        procedure :: update_cov_flags
    end type Clock_class

    ! Constructor
    interface Clock_class
        module procedure constructor_time
    end interface

contains

    ! ====================================================================
    !!
    !>  @brief      Constructor that should be used to initialize variables.
    !!
    !!  @author     Christopher Kebschull
    !!  @date       <ul>
    !!                  <li>ChK: 24.12.2017 (initial implementation)</li>
    !!                  <li>ChK: 13.06.2020 (Add intemediate steps support)</li>
    !!              </ul>
    !!  @anchor     constructor_time
    !!
    ! --------------------------------------------------------------------
    type(Clock_class) function constructor_time(start_epoch_sec, end_epoch_sec, step_epochs_sec)
        real(dp),intent(in)             :: start_epoch_sec                      !< start time (in seconds)
        real(dp),intent(in)             :: end_epoch_sec                        !< end time (in seconds)
        real(dp),dimension(:),optional  :: step_epochs_sec                      ! intermediate epochs in seconds (MJD)


        constructor_time%start_time                     = start_epoch_sec
        constructor_time%end_time                       = end_epoch_sec
        constructor_time%step_size                      = 0.d0                  !< step size - usually the output step size (in seconds)
        constructor_time%cov_step                       = 0.d0                  !< covariance matrix integration step size (seconds)
        constructor_time%out_counter                    = 0.d0                  !< counter keeping track of output writing (seconds)
        constructor_time%cov_counter                    = 0.d0                  !< counter keeping track of covariance propagation (seconds)
        constructor_time%next_step                      = 0.d0                  !< time for the next step the integration aims at (seconds)
        constructor_time%last_counter_success           = -1.d0                 !< keeps track of the last successful step
        constructor_time%prop_counter_reset             = huge(1.d0)            !< keeps track of the time the propagation was reset
        constructor_time%cov_half_step_counter          = 0
        constructor_time%forward                        = .true.                !< .true. if forward propagation, .false. for backwards
        constructor_time%flag_output_step               = .false.               !< is .true. during output steps
        constructor_time%flag_cov_update                = .false.               !< is .true. during covariance matrix update/integration steps
        constructor_time%flag_cov_save                  = .false.               !< is .true. when an intermediate step for the covariance needs to be saved (like for RK4)
        constructor_time%intermediate_steps_flag        = .false.
        constructor_time%intermediate_steps_index       = 1
        constructor_time%last_index                     = 1
        constructor_time%cumulated_time                 = 0.d0
        !=====================================================
        !
        ! Consider backward propagation
        !
        !---------------------------------------------
        if(constructor_time%end_time < constructor_time%start_time) then
          call constructor_time%switch_backward_propagation
        else
          call constructor_time%switch_forward_propagation
        end if

        !=====================================================
        !
        ! Consider requested intermediate propagation steps
        !
        !---------------------------------------------
        if (present(step_epochs_sec)) then
          constructor_time%step_epochs_sec = step_epochs_sec
          constructor_time%intermediate_steps_flag = .true.
        end if

    end function constructor_time

    ! ====================================================================
    !!
    !>  @brief      Get the next requested step
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     get_next_step
    !!
    ! --------------------------------------------------------------------
    function get_next_step(this,neptune,current_time) result(step)
        use slam_io,                       only: LOG_AND_STDOUT, message
        use slam_strings,                  only: toString

        class(Clock_class)                  :: this
        type(Neptune_class),intent(inout)   :: neptune
        real(dp)           ,intent(in)      :: current_time

        real(dp)                            :: step
        real(dp), parameter                 :: eps3 = 1.d-3
        real(dp)                            :: cov_step                         !< covariance matrix integration step size / s
        integer                             :: i_index                          !< intermediate step index

        character(len=255)                  :: cmess

        ! update step_size with every request when in intermediate steps mode
        if (this%intermediate_steps_flag) then
            
            ! this%step_size = this%step_epochs_sec(this%intermediate_steps_index)
            
            ! Determine on which time step we are working
            ! -- DO NOT TOUCH --
            ! -- The new indexing-loop implementation now processes cases where 
            ! -- neighbouring array-elements contain the same epoch (due to overlapping TDMs). 
            ! -- Previously, those situations led to wrong in-array indexing and infinite 
            ! -- looping, exhausting the allocated memory. 
            do i_index = this%last_index, size(this%step_epochs_sec)
                this%cumulated_time = this%cumulated_time + this%step_epochs_sec(i_index)
                if ((this%forward       .and. (this%cumulated_time > current_time .or. i_index <= this%last_index)) .or. &
                    (.not. this%forward .and. (this%cumulated_time < current_time .or. i_index <= this%last_index))) then
                    this%last_index = i_index
                    if (i_index < size(this%step_epochs_sec)) then
                        this%last_index = i_index + 1
                    end if
                    exit
                end if
            end do
            ! Extract the step size
            this%step_size = abs(this%step_epochs_sec(this%last_index))
         end if

        ! for the covariance step, it can depend on the integration method. For example, the RK4 method requires
        ! half steps to be saved
        if(neptune%numerical_integrator%get_covariance_integration_method() == RK4) then
            cov_step = 0.5d0*this%cov_step
        else
            cov_step = this%cov_step
        end if

        if(this%forward) then
            if(this%cov_counter < this%out_counter) then
                ! only covariance save/update step
                step                 = this%cov_counter
                this%cov_counter     = this%cov_counter + cov_step
                this%flag_output_step = .false.
                call this%update_cov_flags(neptune)
            else if(abs(this%cov_counter - this%out_counter) < eps3) then
                ! simultaneous output and covariance save/update step
                step                  = min(this%cov_counter, this%out_counter)
                this%out_counter      = this%out_counter + this%step_size
                this%cov_counter      = this%cov_counter + cov_step
                this%flag_output_step = .true.
                call this%update_cov_flags(neptune)
            else
                ! only output step
                step                  = this%out_counter
                this%out_counter      = this%out_counter + this%step_size
                this%flag_output_step = .true.
                if (this%intermediate_steps_flag .and. this%cov_propagation_flag) then
                    ! Always update the set matrix and covariance to the desired time
                    call this%update_cov_flags(neptune)
                else
                    this%flag_cov_update  = .false.
                    this%flag_cov_save    = .false.
                end if
            end if
            if (step > this%end_time) then
                step = this%end_time
            end if
        else
            if(this%cov_counter > this%out_counter) then
                ! only covariance save/update step
                step                  = this%cov_counter
                this%cov_counter      = this%cov_counter - cov_step
                this%flag_output_step = .false.
                call this%update_cov_flags(neptune)
            else if(abs(this%cov_counter - this%out_counter) < eps3) then
                ! simultaneous output and covariance save/update step
                step = min(this%cov_counter, this%out_counter)
                this%out_counter      = this%out_counter - this%step_size
                this%cov_counter      = this%cov_counter - cov_step
                this%flag_output_step = .true.
                call this%update_cov_flags(neptune)
            else
                ! only output step
                step                  = this%out_counter
                this%out_counter      = this%out_counter - this%step_size
                this%flag_output_step = .true.
                if (this%intermediate_steps_flag .and. this%cov_propagation_flag) then
                    ! Always update the set matrix and covariance to the desired time
                    call this%update_cov_flags(neptune)
                else
                    this%flag_cov_update  = .false.
                    this%flag_cov_save    = .false.
                end if
            end if
            if (step < this%end_time) then
                step = this%end_time
            end if
        end if
        ! save the value for later reference (e.g. by has_finished_step)
        this%next_step = step

        ! write(cmess,'(a)') 'this%step_size  = '//toString(this%step_size)//', this%next_step '//toString(this%next_step)
        ! call message(cmess, LOG_AND_STDOUT)        
        ! call message("request epoch: "//toString(this%next_step)//" this%last_index: "//toString(this%last_index), LOG_AND_STDOUT)
        
        return
    end function

    ! ====================================================================
    !
    !>  @brief      Updates the half-step and full step flags for the covariance matrix
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 14.07.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     update_cov_flags
    !
    ! --------------------------------------------------------------------
    subroutine update_cov_flags(this,neptune)
        class(Clock_class)                  :: this
        type(Neptune_class),intent(inout)   :: neptune

        if(neptune%numerical_integrator%get_covariance_integration_method() == RK4) then
            this%cov_half_step_counter = this%cov_half_step_counter + 1
        end if
        ! the logic below makes sure that we either save the half-step or update at full steps
        if(mod(this%cov_half_step_counter,2) == 1) then
            this%flag_cov_save   = .true.
            this%flag_cov_update = .false.
        else
            this%flag_cov_save   = .false.
            this%flag_cov_update = .true.
        end if
        return
    end subroutine

    ! ====================================================================
    !
    !>  @brief      Toggles the output flag (.true. to .false. and .false. to .true.)
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     toggle_output_flag
    !
    ! --------------------------------------------------------------------

    subroutine toggle_output_flag(this)
        class(Clock_class) :: this
        if(this%flag_output_step) then
            this%flag_output_step = .false.
        else
            this%flag_output_step = .true.
        end if
        return
    end subroutine

    ! ====================================================================
    !
    !>  @brief      Toggles the covariance update flag (.true. to .false. and .false. to .true.)
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     toggle_covariance_update_flag
    !
    ! --------------------------------------------------------------------

    subroutine toggle_cov_update_flag(this)
        class(Clock_class) :: this
        if(this%flag_cov_update) then
            this%flag_cov_update = .false.
        else
            this%flag_cov_update = .true.
        end if
        return
    end subroutine


    ! ====================================================================
    !
    !>  @brief      Switch to backward propagation
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     switch_backward_propagation
    !
    ! --------------------------------------------------------------------

    subroutine switch_backward_propagation(this)
        class(Clock_class) :: this
        this%forward = .false.
        return
    end subroutine

    ! ====================================================================
    !
    !>  @brief      Switch to forward propagation
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     switch_backward_propagation
    !
    ! --------------------------------------------------------------------

    subroutine switch_forward_propagation(this)
        class(Clock_class) :: this
        this%forward = .true.
        return
    end subroutine

    ! ====================================================================
    !
    !>  @brief      Initialises the counters
    !!
    !!  @author     Vitali Braun
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor     init_counter
    !
    ! --------------------------------------------------------------------

    subroutine init_counter(this,neptune)
        class(Clock_class)                      :: this
        type(Neptune_class),intent(inout)       :: neptune
        real(dp)                                :: cov_step                     !< to handle the special case of covariance save vs. update

        !** step size / output write counter
        !-------------------------------------------------------
        if(neptune%output%get_output_switch()) then
            if (this%intermediate_steps_flag) then
                ! First step for intermediate steps is defined by the user
                this%step_size = this%step_epochs_sec(1)
            else if(neptune%getStoreDataFlag()) then
                this%step_size = dble(neptune%getStep())
            else
                this%step_size = dble(neptune%get_output_step()) ! in seconds
            end if

            this%flag_output_step = .true.

            if(this%forward) then
                this%out_counter = this%start_time + this%step_size
            else
                this%out_counter = this%start_time - this%step_size
            end if

        else
            if(neptune%getStoreDataFlag()) then
                if (this%intermediate_steps_flag) then
                    ! First step for intermediate steps is defined by the user
                    this%step_size = this%step_epochs_sec(1)
                else
                    this%step_size   = dble(neptune%getStep())
                end if
                if(this%forward) then
                    this%out_counter = this%start_time + this%step_size
                else
                    this%out_counter = this%start_time - this%step_size
                end if
                this%flag_output_step = .true.
            else
                this%step_size   = abs(this%end_time - this%start_time)
                this%out_counter = this%end_time
            end if
        end if

        !** set covariance step counter
        !-----------------------------------------------------
        this%cov_step = neptune%numerical_integrator%getCovarianceIntegrationStep()
        call neptune%numerical_integrator%setCovarianceIntegrationStep(this%cov_step)
        this%cov_propagation_flag = neptune%numerical_integrator%getCovariancePropagationFlag()

        ! for the covariance step, it can depend on the integration method. For example, the RK4 method requires
        ! half steps to be saved
        if(neptune%numerical_integrator%get_covariance_integration_method() == RK4) then
            cov_step = 0.5d0*this%cov_step
        else
            cov_step = this%cov_step
        end if

        if(this%cov_propagation_flag) then
            if(this%forward) then
                this%cov_counter  = this%start_time + cov_step
            else
                this%cov_counter  = this%start_time - cov_step
            end if
        else
            this%cov_counter = this%end_time
        end if

        !** the counter that keeps track of the time of the last reset (initialise with an arbitrary value)
        if(this%forward) then
            this%prop_counter_reset = huge(1.d0)
        else
            this%prop_counter_reset = -1.d0
        end if

        return
    end subroutine

    ! ====================================================================
    !
    !>  @brief      Tells whether the time of the current step has been reached
    !!
    !!  @author     Vitali Braun
    !!
    !!  @param[in]  current_time    Time which is checked against step time (in seconds)
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor    has_finished_step
    !
    ! --------------------------------------------------------------------
    function has_finished_step(this, current_time) result(hfs)
        ! use slam_strings,       only: toString
        ! use slam_io,            only: LOG_AND_STDOUT, message
        Class(Clock_class)            :: this
        real(dp),     intent(in)      :: current_time
        logical                       :: hfs
        character(len=255)            :: cmess
        real(dp)                      :: diff 

        ! diff = abs(this%next_step - current_time)
        ! write(cmess, '(a, D15.3, a)') 'diff_nextStep = ', diff, toString((this%next_step - current_time) < epsilon(1.0d0))
        ! call message(cmess, LOG_AND_STDOUT)
        
        hfs = .false.
        
        if(this%forward) then
            if((this%next_step - current_time) < epsilon(1.0d0)) hfs = .true.
        else
            if((current_time - this%next_step) < epsilon(1.0d0)) hfs = .true.
        end if
        
        ! if (this%intermediate_steps_flag .and. hfs) then 
        !     if (this%intermediate_steps_index < size(this%step_epochs_sec)) then
        !         this%intermediate_steps_index = this%intermediate_steps_index + 1
        !     else 
        !         this%intermediate_steps_index = size(this%step_epochs_sec)
        !     end if
        ! end if

        return
    end function


    ! ====================================================================
    !
    !>  @brief      Tells whether the end epoch has been reached
    !!
    !!  @author     Vitali Braun
    !!
    !!  @param[in]  current_time    Time which is checked against end epoch (in seconds)
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor    has_finished
    !
    ! --------------------------------------------------------------------
    pure function has_finished(this, current_time) result(hf)
        Class(Clock_class), intent(in) :: this
        real(dp),     intent(in) :: current_time
        logical :: hf

        hf = .false.
        if(this%forward) then
            if((this%end_time - current_time) < epsilon(1.0d0)) hf = .true.
        else
            if((current_time - this%end_time) < epsilon(1.0d0)) hf = .true.
        end if
        return
    end function

    ! ====================================================================
    !
    !>  @brief      True if the current step has to write output
    !!
    !!  @author     Vitali Braun
    !!
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor    has_to_write_output
    !
    ! --------------------------------------------------------------------
    pure function has_to_write_output(this) result(h)
        Class(Clock_class), intent(in) :: this
        logical :: h
        h = this%flag_output_step
        return
    end function

    ! ====================================================================
    !
    !>  @brief      True if the current step has to save the current state for a later update of the covariance matrix
    !!
    !!  @author     Vitali Braun
    !!
    !!  @date       <ul>
    !!                  <li>VB: 14.07.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor    has_to_save_covariance
    !
    ! --------------------------------------------------------------------
    pure function has_to_save_covariance(this) result(h)
        Class(Clock_class), intent(in) :: this
        logical :: h
        h = this%flag_cov_save
        return
    end function

    ! ====================================================================
    !
    !>  @brief      True if the current step has to update the covariance matrix
    !!
    !!  @author     Vitali Braun
    !!
    !!  @date       <ul>
    !!                  <li>VB: 29.04.2017 (initial implementation)</li>
    !!              </ul>
    !!  @anchor    has_to_update_covariance
    !
    ! --------------------------------------------------------------------
    pure function has_to_update_covariance(this) result(h)
        Class(Clock_class), intent(in) :: this
        logical :: h
        h = this%flag_cov_update
        return
    end function



end module neptuneClock
