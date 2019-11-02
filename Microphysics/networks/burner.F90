module burner_module

  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module
  use network
  use eos_module
#ifndef SIMPLIFIED_SDC
  use actual_burner_module
#else
  use integrator_module
#endif
  use burn_type_module

  logical :: burner_initialized = .false.

contains

  subroutine burner_init() bind(C, name="burner_init")

    implicit none

#ifdef SIMPLIFIED_SDC
    call integrator_init()
#else
    call actual_burner_init()
#endif

    burner_initialized = .true.

  end subroutine burner_init



  function ok_to_burn(state) result(ok_to_burnr)

    !$gpu

    use meth_params_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max

    implicit none

    logical       :: ok_to_burnr
    type (burn_t) :: state

    ok_to_burnr = .true.

    if (state % T < react_T_min .or. state % T > react_T_max .or. &
        state % rho < react_rho_min .or. state % rho > react_rho_max) then

       ok_to_burnr = .false.

    endif

  end function ok_to_burn



#ifndef SIMPLIFIED_SDC
  subroutine burner(state_in, state_out, dt, time)

    !$gpu

    use castro_error_module

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    double precision, intent(in) :: dt, time

    ! Make sure the network and burner have been initialized.

#if !(defined(ACC)||defined(AMREX_USE_CUDA))
    if (.NOT. network_initialized) then
       call castro_error("ERROR in burner: must initialize network first.")
    endif

    if (.NOT. burner_initialized) then
       call castro_error("ERROR in burner: must initialize burner first.")
    endif
#endif

    ! Initialize the final state by assuming it does not change.
    call copy_burn_t(state_out, state_in)

    ! Do the burning.

    if (ok_to_burn(state_in)) then
       call actual_burner(state_in, state_out, dt, time)
    endif

  end subroutine burner
#endif

end module burner_module
