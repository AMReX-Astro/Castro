module burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use actual_burner_module

contains

  subroutine burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t), intent(inout)  :: state_in
    type (eos_t), intent(inout)  :: state_out
    double precision, intent(in) :: dt, time

    integer             :: i
    
    ! Make sure the network has been initialized.
    
    if (.NOT. network_initialized) then
       call bl_error("ERROR in burner: must initialize network first.")
    endif
    
    ! We assume that the valid quantities coming in are (rho, e); do an EOS call
    ! to make sure all other variables are consistent.

    call eos(eos_input_re, state_in)
    
    ! Initialize the final state by assuming it does not change.

    state_out = state_in

    ! Do the burning.
    
    call actual_burner(state_in, state_out, dt, time)

    ! Normalize the mass fractions to unity.

    call normalize_abundances(state_out)
    
    ! Now update the temperature to match the new internal energy.

    call eos(eos_input_re, state_out)

  end subroutine burner

end module burner_module
