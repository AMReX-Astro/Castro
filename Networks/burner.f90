module burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_module
  use specific_burner_module

contains

  subroutine burner(state_in, state_out, dt, time)

    implicit none

    class (eos_type), intent(inout)  :: state_in
    class (eos_type), intent(inout)  :: state_out
    double precision, intent(in) :: dt, time

    integer :: i
    type (eos_t_vector) :: state_vector_in
    type (eos_t_vector) :: state_vector_out
    
    ! Make sure the network has been initialized.
    
    if (.NOT. network_initialized) then
       call bl_error("ERROR in burner: must initialize network first.")
    endif

    ! Initialize the final state by assuming it does not change.

    select type (state_in)
    type is (eos_t_1D)
       select type (state_out)
       type is (eos_t_1D)
          state_out = state_in
       end select

    type is (eos_t_2D)
       select type (state_out)
       type is (eos_t_2D)
          state_out = state_in
       end select

    type is (eos_t_3D)
       select type (state_out)
       type is (eos_t_3D)
          state_out = state_in
       end select

    end select

    ! Get an EOS vector for each case.

    call eos_vector_in(state_vector_in, state_in)
    call eos_vector_in(state_vector_out, state_out)

    ! We assume that the valid quantities coming in are (rho, e); do an EOS call
    ! to make sure all other variables are consistent. This will also
    ! properly

    call eos(eos_input_re, state_vector_in)

    ! Do the burning.
    
    call specific_burner(state_vector_in, state_vector_out, dt, time)

    ! Normalize the mass fractions: they must be individually positive and less than one,
    ! and they must all sum to unity.
    
    do i = 1, state_vector_out % N

       state_vector_out % xn(i,:) = max(smallx, min(ONE, state_vector_out % xn(i,:)))

       state_vector_out % xn(i,:) = state_vector_out % xn(i,:) / sum(state_vector_out % xn(i,:))

    enddo
       
    ! Now update the temperature to match the new internal energy.

    call eos(eos_input_re, state_vector_out)

  end subroutine burner

end module burner_module
