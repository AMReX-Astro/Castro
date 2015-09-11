module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use eos_data_module
  use actual_eos_module

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

    call actual_eos_init

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

    initialized = .true.

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

    call actual_eos(input, state_vector)

    ! Get dpdX, dedX, dhdX.

    call composition_derivatives(state_vector)

    ! Convert the vectorized state back to the output form.

    call eos_vector_out(state_vector, state)
    
  end subroutine eos

  

  subroutine check_inputs(input, state)

    implicit none

    integer,             intent(in   ) :: input
    type (eos_t_vector), intent(inout) :: state

    integer :: i, n
    
    do i = 1, state % N

       ! Check the inputs, and do initial setup for iterations.

       do n = 1, nspec
          if (state % xn(i,n) .lt. init_test) call eos_type_error(ierr_init_xn, input)
       enddo

       if ( state % y_e(i) .lt. minye ) then
          print *, 'Y_E = ', state % y_e(i)
          call bl_error('EOS: y_e less than minimum possible electron fraction.')
       endif
       if ( state % y_e(i) .gt. maxye ) then
          print *, 'Y_E = ', state % y_e(i)
          call bl_error('EOS: y_e greater than maximum possible electron fraction.')
       endif       
       
    enddo
       
    ! Our strategy for testing the validity of the inputs is as follows.
    ! First, if the quantities for the given call type haven't been initialized,
    ! then throw an error. Second, if the quantity is rho or T and it has been initialized
    ! but is less than smalld or smallt, reset to smalld or smallt. Third, if the
    ! quantity is something else (e, h, etc.) and is less than zero in a zone,
    ! make sure both T and rho are at least as large as smallt and small d, then
    ! call the EOS in (rho, T) mode just on that zone.
       
    if (input .eq. eos_input_rt) then

       call check_rho(state)
       call check_T(state)
          
    elseif (input .eq. eos_input_rh) then

       call check_rho(state)
       call check_h(state)

    elseif (input .eq. eos_input_tp) then
       
       call check_T(state)
       call check_p(state)

    elseif (input .eq. eos_input_rp) then

       call check_rho(state)
       call check_p(state)

    elseif (input .eq. eos_input_re) then

       call check_rho(state)
       call check_e(state)

    elseif (input .eq. eos_input_ps) then

       call check_p(state)
       call check_s(state)

    elseif (input .eq. eos_input_ph) then

       call check_p(state)
       call check_h(state)

    elseif (input .eq. eos_input_th) then

       call check_t(state)
       call check_h(state)

    endif
    
  end subroutine check_inputs



  subroutine check_rho(state)

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i

    do i = 1, state % N
    
       if (state % rho(i) .lt. init_test) then
          call bl_error('EOS: rho not initialized.')
       endif
    
       if (state % rho(i) .lt. smalld) then         
          if (state % reset) then
             state % rho(i) = smalld
          else
             print *, 'DENS = ', state % rho(i)
             call bl_error('EOS: rho smaller than small_dens and we have not chosen to reset.')
          endif
       endif
    
       if (state % rho(i) .gt. maxdens) then
          print *, 'DENS = ', state % rho(i)
          call bl_error('EOS: dens greater than maximum possible density.')
       endif

    enddo

  end subroutine check_rho


  
  subroutine check_T(state)

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i

    do i = 1, state % N
    
       if (state % T(i) .lt. init_test) then
          call bl_error('EOS: T not initialized.')
       endif
    
       if (state % T(i) .lt. smallt) then         
          if (state % reset) then
             state % T(i) = smallt
          else
             print *, 'DENS = ', state % T(i)
             call bl_error('EOS: T smaller than small_dens and we have not chosen to reset.')
          endif
       endif
    
       if (state % T(i) .gt. maxdens) then
          print *, 'DENS = ', state % T(i)
          call bl_error('EOS: dens greater than maximum possible density.')
       endif

    enddo
       
  end subroutine check_T  

  

  subroutine check_e(state)

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i

    do i = 1, state % N

       if (state % e(i) .lt. init_test) then
          call bl_error('EOS: energy not initialized.')
       endif
    
       if (state % e(i) .lt. ZERO) then
          if (state % reset) then
             state % T(i) = max(smallt, state % T(i))
             state % rho(i) = max(smalld, state % rho(i))
             call eos_reset(state, i)
          else
             call bl_error('EOS: e smaller than zero and we have not chosen to reset.')
          endif
       endif

    enddo
       
  end subroutine check_e
  


  subroutine check_h(state)

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i

    do i = 1, state % N

       if (state % h(i) .lt. init_test) then
          call bl_error('EOS: enthalpy not initialized.')
       endif
    
       if (state % h(i) .lt. ZERO) then
          if (state % reset) then
             state % T(i) = max(smallt, state % T(i))
             state % rho(i) = max(smalld, state % rho(i))             
             call eos_reset(state, i)
          else
             call bl_error('EOS: h smaller than zero and we have not chosen to reset.')
          endif
       endif

    enddo
       
  end subroutine check_h



  subroutine check_s(state)

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i

    do i = 1, state % N

       if (state % s(i) .lt. init_test) then
          call bl_error('EOS: entropy not initialized.')
       endif
    
       if (state % s(i) .lt. ZERO) then
          if (state % reset) then
             state % T(i) = max(smallt, state % T(i))
             state % rho(i) = max(smalld, state % rho(i))             
             call eos_reset(state, i)
          else
             call bl_error('EOS: s smaller than zero and we have not chosen to reset.')
          endif
       endif

    enddo
       
  end subroutine check_s



  subroutine check_p(state)

    implicit none

    type (eos_t_vector), intent(inout) :: state

    integer :: i

    do i = 1, state % N

       if (state % p(i) .lt. init_test) then
          call bl_error('EOS: pressure not initialized.')
       endif
    
       if (state % p(i) .lt. ZERO) then
          if (state % reset) then
             state % T(i) = max(smallt, state % T(i))
             state % rho(i) = max(smalld, state % rho(i))             
             call eos_reset(state, i)
          else
             call bl_error('EOS: p smaller than zero and we have not chosen to reset.')
          endif
       endif

    enddo
       
  end subroutine check_p



  ! Given an EOS vector and an input i,
  ! the code has reset some characteristic of that
  ! state element and we now want to call the EOS just
  ! on that zone to reset its state values.
  
  subroutine eos_reset(state, i)

    use actual_eos_module
    
    implicit none

    type (eos_t_vector), intent(inout) :: state
    integer,             intent(in   ) :: i
    
    type (eos_t) :: reset_state
    type (eos_t_vector) :: reset_state_vector

    if (.not. state % reset) then
       call bl_error('EOS: should not be in eos_reset if reset == .false.')
    endif
    
    call get_eos_t(state, reset_state, i)
    call eos_vector_in(reset_state_vector, reset_state)
    call actual_eos(eos_input_rt, reset_state_vector)
    call eos_vector_out(reset_state_vector, reset_state)
    call put_eos_t(state, reset_state, i)

  end subroutine eos_reset



end module eos_module
