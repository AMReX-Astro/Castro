module burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use actual_burner_module
  use burn_type_module
  use meth_params_module, only: react_T_min, react_T_max

contains

  subroutine burner_init() bind(C)

    implicit none

    call actual_burner_init()

  end subroutine burner_init



  subroutine burner(state_in, state_out, dt, time)

    implicit none

    type (burn_t), intent(inout) :: state_in
    type (burn_t), intent(inout) :: state_out
    double precision, intent(in) :: dt, time

    ! Make sure the network has been initialized.

    if (.NOT. network_initialized) then
       call bl_error("ERROR in burner: must initialize network first.")
    endif

    ! Initialize the final state by assuming it does not change.

    state_out = state_in

    ! Do the burning.

    if (state_in % T >= react_T_min .and. state_in % T <= react_T_max) then
       call actual_burner(state_in, state_out, dt, time)
    endif

  end subroutine burner

end module burner_module
