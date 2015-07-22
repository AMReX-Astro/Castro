module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module
  use specific_eos_module

  implicit none

  public eos_init, eos

contains

  ! EOS initialization routine: read in general EOS parameters, then 
  ! call any specific initialization used by the EOS.

  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module, only: eos_assume_neutral

    implicit none
 
    double precision, optional :: small_temp
    double precision, optional :: small_dens
 
    ! Set up any specific parameters or initialization steps required by the EOS we are using.

    call specific_eos_init

    ! If they exist, save the minimum permitted user temperature and density.
    ! These cannot be less than zero and they also cannot be less than the 
    ! minimum possible EOS quantities.

    if (present(small_temp)) then
       if (small_temp > ZERO) then
          if (small_temp < mintemp) then
             call bl_warn('EOS: small_temp cannot be less than the mintemp allowed by the EOS. Resetting smallt to mintemp.')
             small_temp = mintemp
          endif
          smallt = small_temp
       endif
    endif

    if (present(small_dens)) then
       if (small_dens > ZERO) then
          if (small_dens < mindens) then
             call bl_warn('EOS: small_dens cannot be less than the mindens allowed by the EOS. Resetting smalld to mindens.')
             small_dens = mindens
          endif
          smalld = small_dens
       endif
    endif

    assume_neutral = eos_assume_neutral

  end subroutine eos_init



  subroutine eos(input, state)

    implicit none

    ! Input arguments

    integer,          intent(in   ) :: input
    class (eos_type), intent(inout) :: state

    ! Local variables

    type (eos_t_vector) :: state_vector

    if (.not. initialized) call bl_error('EOS: not initialized')

    ! Convert from the incoming type to the vectorized type we work with in the EOS.

    call eos_vector_in(state_vector, state)

    ! Get abar, zbar, etc.

    call composition(state_vector)

    ! Check to make sure the inputs are valid.

    call check_inputs(input, state_vector)
    
    ! Call the EOS.

    call specific_eos(input, state_vector)

    ! Get dpdX, dedX, dhdX.

    call composition_derivatives(state_vector)

    ! Convert the vectorized state back to the output form.

    call eos_vector_out(state_vector, state)
    
  end subroutine eos

end module eos_module
