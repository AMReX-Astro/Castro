module burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use actual_burner_module
  use burn_type_module

  logical :: burner_initialized = .false.

contains

  subroutine burner_init() bind(C, name="burner_init")

    implicit none

    call actual_burner_init()

    burner_initialized = .true.

  end subroutine burner_init



  function ok_to_burn(state)

    !$acc routine seq

    use meth_params_module, only: react_T_min, react_T_max, react_rho_min, react_rho_max, &
                                  disable_shock_burning

    implicit none

    logical       :: ok_to_burn
    type (burn_t) :: state

    ok_to_burn = .true.

    if (state % T < react_T_min .or. state % T > react_T_max .or. &
        state % rho < react_rho_min .or. state % rho > react_rho_max .or. &
        (disable_shock_burning .eq. 1 .and. state % shock) ) then

       ok_to_burn = .false.

    endif

  end function ok_to_burn



  subroutine burner(state_in, state_out, dt, time)

    !$acc routine seq

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    double precision, intent(in) :: dt, time

    ! Make sure the network and burner have been initialized.

#ifndef ACC
    if (.NOT. network_initialized) then
       call bl_error("ERROR in burner: must initialize network first.")
    endif

    if (.NOT. burner_initialized) then
       call bl_error("ERROR in burner: must initialize burner first.")
    endif
#endif

    ! Initialize the final state by assuming it does not change.

    state_out = state_in

    ! Do the burning.

    if (ok_to_burn(state_in)) then
       call actual_burner(state_in, state_out, dt, time)
    endif

  end subroutine burner

end module burner_module
