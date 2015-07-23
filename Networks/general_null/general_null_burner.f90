module actual_burner_module

  use bl_types
  use bl_constants_module
  use network
  use eos_type_module
  
contains

  subroutine actual_burner(state_in, state_out, dt, time)

    implicit none

    type (eos_t_vector), intent(in)    :: state_in
    type (eos_t_vector), intent(inout) :: state_out
    double precision, intent(in)       :: dt, time

    ! Do nothing in this burner.
    
  end subroutine actual_burner

end module actual_burner_module
