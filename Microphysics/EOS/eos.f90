module eos_module

  use bl_types
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use actual_eos_module

  implicit none

  public eos_init, eos, eos_get_small_temp, eos_get_small_dens

  logical, save :: initialized = .false.  

contains

  subroutine eos_get_small_temp(small_temp_out)

    double precision, intent(out) :: small_temp_out

    small_temp_out = smallt

  end subroutine eos_get_small_temp



  subroutine eos_get_small_dens(small_dens_out)

    double precision, intent(out) :: small_dens_out

    small_dens_out = smalld

  end subroutine eos_get_small_dens



  ! EOS initialization routine: read in general EOS parameters, then 
  ! call any specific initialization used by the EOS.

  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module
    use parallel

    implicit none

    double precision, optional :: small_temp
    double precision, optional :: small_dens

    ! Set up any specific parameters or initialization steps required by the EOS we are using.

    call actual_eos_init

    ! If they exist, save the minimum permitted user temperature and density.
    ! These cannot be less than zero and they also cannot be less than the 
    ! minimum possible EOS quantities.

    ! Note that in this routine we use the Fortran-based parallel_IOProcessor()
    ! command rather than the C++-based version used elsewhere in Castro; this
    ! ensures compatibility with Fortran-based test programs.

    if (present(small_temp)) then
       if (small_temp > ZERO) then
          if (small_temp < mintemp) then
             if (parallel_IOProcessor()) then
                call bl_warn('EOS: small_temp cannot be less than the mintemp allowed by the EOS. Resetting smallt to mintemp.')
             endif
             small_temp = mintemp
          endif
          smallt = small_temp
       endif
    endif

    if (present(small_dens)) then
       if (small_dens > ZERO) then
          if (small_dens < mindens) then
             if (parallel_IOProcessor()) then
                call bl_warn('EOS: small_dens cannot be less than the mindens allowed by the EOS. Resetting smalld to mindens.')
             endif
             small_dens = mindens
          endif
          smalld = small_dens
       endif
    endif

    initialized = .true.

  end subroutine eos_init



  subroutine eos(input, state)

    implicit none

    ! Input arguments

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    logical :: has_been_reset

    ! Local variables

    if (.not. initialized) call bl_error('EOS: not initialized')

    ! Get abar, zbar, etc.

    call composition(state)

    ! Set small variables.

    state % smallt = smallt
    state % smalld = smalld

    ! Check to make sure the inputs are valid.

    call check_inputs(input, state, has_been_reset)

    ! Call the EOS.

    if (.not. has_been_reset) then
       call actual_eos(input, state)
    endif

    ! Get dpdX, dedX, dhdX.

    call composition_derivatives(state)

  end subroutine eos



  subroutine check_inputs(input, state, reset)

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    integer :: n

    reset = .false.

    ! Check the inputs, and do initial setup for iterations.

    do n = 1, nspec
       if (state % xn(n) .lt. init_test) then
          call bl_error('EOS: mass fractions not initialized.')
       endif
    enddo

    if ( state % y_e .lt. minye ) then
       print *, 'Y_E  = ', state % y_e
       print *, 'DENS = ', state % rho
       print *, 'TEMP = ', state % T
       print *, 'X    = ', state % xn
       call bl_error('EOS: y_e less than minimum possible electron fraction.')
    endif
    if ( state % y_e .gt. maxye ) then
       print *, 'Y_E  = ', state % y_e
       print *, 'DENS = ', state % rho
       print *, 'TEMP = ', state % T
       print *, 'X    = ', state % xn
       call bl_error('EOS: y_e greater than maximum possible electron fraction.')
    endif

    ! Our strategy for testing the validity of the inputs is as follows.
    ! First, if the quantities for the given call type haven't been initialized,
    ! then throw an error. Second, if the quantity is rho or T and it has been initialized
    ! but is less than smalld or smallt, reset to smalld or smallt. Third, if the
    ! quantity is something else (e, h, etc.) and is less than zero in a zone,
    ! make sure both T and rho are at least as large as smallt and small d, then
    ! call the EOS in (rho, T) mode just on that zone.

    if (input .eq. eos_input_rt) then

       call check_rho(state, reset)
       call check_T(state, reset)

    elseif (input .eq. eos_input_rh) then

       call check_rho(state, reset)
       call check_h(state, reset)

    elseif (input .eq. eos_input_tp) then

       call check_T(state, reset)
       call check_p(state, reset)

    elseif (input .eq. eos_input_rp) then

       call check_rho(state, reset)
       call check_p(state, reset)

    elseif (input .eq. eos_input_re) then

       call check_rho(state, reset)
       call check_e(state, reset)

    elseif (input .eq. eos_input_ps) then

       call check_p(state, reset)
       call check_s(state, reset)

    elseif (input .eq. eos_input_ph) then

       call check_p(state, reset)
       call check_h(state, reset)

    elseif (input .eq. eos_input_th) then

       call check_t(state, reset)
       call check_h(state, reset)

    endif

  end subroutine check_inputs



  subroutine check_rho(state, reset)

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    if (state % rho .lt. init_test) then
       call bl_error('EOS: rho not initialized.')
    endif

    if (state % rho .lt. smalld .and. state % check_small) then
       if (state % reset) then
          state % rho = smalld
       else
          print *, 'DENS = ', state % rho
          print *, 'TEMP = ', state % T
          print *, 'X    = ', state % xn
          call bl_error('EOS: rho smaller than small_dens and we have not chosen to reset.')
       endif
    endif

    if (state % rho .gt. maxdens) then
       print *, 'DENS = ', state % rho
       print *, 'TEMP = ', state % T
       print *, 'X    = ', state % xn
       call bl_error('EOS: dens greater than maximum possible density.')
    endif

  end subroutine check_rho



  subroutine check_T(state, reset)

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    if (state % T .lt. init_test) then
       call bl_error('EOS: T not initialized.')
    endif

    if (state % T .lt. smallt .and. state % check_small) then
       if (state % reset) then
          state % T = smallt
       else
          print *, 'TEMP = ', state % T
          print *, 'DENS = ', state % rho
          print *, 'X    = ', state % xn
          call bl_error('EOS: T smaller than small_temp and we have not chosen to reset.')
       endif
    endif

    if (state % T .gt. maxdens) then
       print *, 'TEMP = ', state % T
       print *, 'DENS = ', state % rho
       print *, 'X    = ', state % xn
       call bl_error('EOS: T greater than maximum possible temperature.')
    endif

  end subroutine check_T



  subroutine check_e(state, reset)

    use meth_params_module, only: allow_negative_energy

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    if (state % e .lt. init_test) then
       call bl_error('EOS: energy not initialized.')
    endif

    if (state % e .lt. ZERO .and. allow_negative_energy .eq. 0) then
       if (state % reset) then
          state % T = max(smallt, state % T)
          state % rho = max(smalld, state % rho)
          call eos_reset(state, reset)
       else
          call bl_error('EOS: e smaller than zero and we have not chosen to reset.')
       endif
    endif

  end subroutine check_e



  subroutine check_h(state, reset)

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    if (state % h .lt. init_test) then
       call bl_error('EOS: enthalpy not initialized.')
    endif

    if (state % h .lt. ZERO) then
       if (state % reset) then
          state % T = max(smallt, state % T)
          state % rho = max(smalld, state % rho)
          call eos_reset(state, reset)
       else
          call bl_error('EOS: h smaller than zero and we have not chosen to reset.')
       endif
    endif

  end subroutine check_h



  subroutine check_s(state, reset)

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    if (state % s .lt. init_test) then
       call bl_error('EOS: entropy not initialized.')
    endif

    if (state % s .lt. ZERO) then
       if (state % reset) then
          state % T = max(smallt, state % T)
          state % rho = max(smalld, state % rho)
          call eos_reset(state, reset)
       else
          call bl_error('EOS: s smaller than zero and we have not chosen to reset.')
       endif
    endif

  end subroutine check_s



  subroutine check_p(state, reset)

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    if (state % p .lt. init_test) then
       call bl_error('EOS: pressure not initialized.')
    endif

    if (state % p .lt. ZERO) then
       if (state % reset) then
          state % T = max(smallt, state % T)
          state % rho = max(smalld, state % rho)
          call eos_reset(state, reset)
       else
          call bl_error('EOS: p smaller than zero and we have not chosen to reset.')
       endif
    endif

  end subroutine check_p



  ! Given an EOS vector and an input i,
  ! the code has reset some characteristic of that
  ! state element and we now want to call the EOS just
  ! on that zone to reset its state values.

  subroutine eos_reset(state, reset)

    use actual_eos_module

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: reset

    call actual_eos(eos_input_rt, state)

    reset = .true.

  end subroutine eos_reset

end module eos_module
