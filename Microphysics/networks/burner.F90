module burner_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use network
  use eos_module
#ifndef SDC
  use actual_burner_module
#else
  use integrator_module
#endif
  use burn_type_module

  logical :: burner_initialized = .false.

contains

  subroutine burner_init() bind(C, name="burner_init")

    implicit none

#ifdef SDC
    call integrator_init()
#else
    call actual_burner_init()
#endif

    burner_initialized = .true.

  end subroutine burner_init



  AMREX_DEVICE subroutine ok_to_burn(state, do_burn)

    !$acc routine seq

    use meth_params_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max

    implicit none

    logical,       intent(inout) :: do_burn
    type (burn_t), intent(in)    :: state

    do_burn = .true.

    if (state % T < react_T_min .or. state % T > react_T_max .or. &
        state % rho < react_rho_min .or. state % rho > react_rho_max) then

       do_burn = .false.

    endif

  end subroutine ok_to_burn



#ifndef SDC
  AMREX_DEVICE subroutine burner(state_in, state_out, dt, time)

    !$acc routine seq

    use amrex_error_module

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    double precision, intent(in) :: dt, time

    logical :: do_burn = .true.

    ! Make sure the network and burner have been initialized.

#if !(defined(ACC) || defined(CUDA))
    if (.NOT. network_initialized) then
       call amrex_error("ERROR in burner: must initialize network first.")
    endif

    if (.NOT. burner_initialized) then
       call amrex_error("ERROR in burner: must initialize burner first.")
    endif
#endif

    ! Initialize the final state by assuming it does not change.
    call copy_burn_t(state_out, state_in)

    ! Do the burning.
    call ok_to_burn(state_in, do_burn)

    if (do_burn) then
       call actual_burner(state_in, state_out, dt, time)
    endif

  end subroutine burner
#endif

end module burner_module
