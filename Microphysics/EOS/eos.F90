module eos_module

  use bl_types, only: dp_t
  use bl_error_module
  use bl_constants_module
  use network, only: nspec, aion, zion
  use eos_type_module
  use actual_eos_module

  implicit none

  public eos_init, eos

  logical, save :: initialized = .false.  

contains

  ! EOS initialization routine: read in general EOS parameters, then 
  ! call any specific initialization used by the EOS.

  subroutine eos_init(small_temp, small_dens)

    use extern_probin_module
    use parallel

    implicit none

    real(dp_t), optional :: small_temp
    real(dp_t), optional :: small_dens

    ! Set up any specific parameters or initialization steps required by the EOS we are using.

    call actual_eos_init

    ! If they exist, save the minimum permitted user temperature and density.
    ! These are only relevant to this module if they are larger than the minimum
    ! possible EOS quantities. We will reset them to be equal to the EOS minimum
    ! if they are smaller than that.

    ! Note that in this routine we use the Fortran-based parallel_IOProcessor()
    ! command rather than the C++-based version used elsewhere in Castro; this
    ! ensures compatibility with Fortran-based test programs.

    if (present(small_temp)) then
       if (small_temp < mintemp) then
          if (parallel_IOProcessor()) then
             call bl_warn('EOS: small_temp cannot be less than the mintemp allowed by the EOS. Resetting small_temp to mintemp.')
          endif
          small_temp = mintemp
       else
          mintemp = small_temp
       endif
    endif

    if (present(small_dens)) then
       if (small_dens < mindens) then
          if (parallel_IOProcessor()) then
             call bl_warn('EOS: small_dens cannot be less than the mindens allowed by the EOS. Resetting small_dens to mindens.')
          endif
          small_dens = mindens
       else
          mindens = small_dens
       endif
    endif

    initialized = .true.

    !$acc update &
    !$acc device(mintemp, maxtemp, mindens, maxdens, minx, maxx, minye, maxye) &
    !$acc device(mine, maxe, minp, maxp, mins, maxs, minh, maxh)

  end subroutine eos_init



  subroutine eos(input, state)

    !$acc routine seq

    implicit none

    ! Input arguments

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    logical :: has_been_reset

    ! Local variables

#ifndef ACC
    if (.not. initialized) call bl_error('EOS: not initialized')
#endif

    ! Get abar, zbar, etc.

    call composition(state)

    ! Force the inputs to be valid.

    has_been_reset = .false.
    call reset_inputs(input, state, has_been_reset)

    ! Call the EOS.

    if (.not. has_been_reset) then
       call actual_eos(input, state)
    endif

    ! Get dpdX, dedX, dhdX.

    call composition_derivatives(state)

  end subroutine eos



  subroutine reset_inputs(input, state, has_been_reset)

    !$acc routine seq

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    ! Reset the input quantities to valid values. For inputs other than rho and T,
    ! this will evolve an EOS call, which will negate the need to do the main EOS call.

    if (input .eq. eos_input_rt) then

       call reset_rho(state, has_been_reset)
       call reset_T(state, has_been_reset)

    elseif (input .eq. eos_input_rh) then

       call reset_rho(state, has_been_reset)
       call reset_h(state, has_been_reset)

    elseif (input .eq. eos_input_tp) then

       call reset_T(state, has_been_reset)
       call reset_p(state, has_been_reset)

    elseif (input .eq. eos_input_rp) then

       call reset_rho(state, has_been_reset)
       call reset_p(state, has_been_reset)

    elseif (input .eq. eos_input_re) then

       call reset_rho(state, has_been_reset)
       call reset_e(state, has_been_reset)

    elseif (input .eq. eos_input_ps) then

       call reset_p(state, has_been_reset)
       call reset_s(state, has_been_reset)

    elseif (input .eq. eos_input_ph) then

       call reset_p(state, has_been_reset)
       call reset_h(state, has_been_reset)

    elseif (input .eq. eos_input_th) then

       call reset_t(state, has_been_reset)
       call reset_h(state, has_been_reset)

    endif

  end subroutine reset_inputs



  ! For density, just ensure that it is within mindens and maxdens.

  subroutine reset_rho(state, has_been_reset)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    state % rho = min(maxdens, max(mindens, state % rho))

  end subroutine reset_rho



  ! For temperature, just ensure that it is within mintemp and maxtemp.

  subroutine reset_T(state, has_been_reset)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    state % T = min(maxtemp, max(mintemp, state % T))

  end subroutine reset_T



  subroutine reset_e(state, has_been_reset)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    if (state % e .lt. mine .or. state % e .gt. maxe) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_e



  subroutine reset_h(state, has_been_reset)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    if (state % h .lt. minh .or. state % h .gt. maxh) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_h



  subroutine reset_s(state, has_been_reset)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    if (state % s .lt. mins .or. state % s .gt. maxs) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_s



  subroutine reset_p(state, has_been_reset)

    !$acc routine seq

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    if (state % p .lt. minp .or. state % p .gt. maxp) then
       call eos_reset(state, has_been_reset)
    endif

  end subroutine reset_p



  ! Given an EOS state, ensure that rho and T are
  ! valid, then call with eos_input_rt.

  subroutine eos_reset(state, has_been_reset)

    !$acc routine seq

    use actual_eos_module

    implicit none

    type (eos_t), intent(inout) :: state
    logical,      intent(inout) :: has_been_reset

    state % T = min(maxtemp, max(mintemp, state % T))
    state % rho = min(maxdens, max(mindens, state % rho))

    call actual_eos(eos_input_rt, state)

    has_been_reset = .true.

  end subroutine eos_reset



#ifndef ACC
  subroutine check_inputs(input, state)

    !$acc routine seq

    implicit none

    integer,      intent(in   ) :: input
    type (eos_t), intent(inout) :: state

    integer :: n

    ! Check the inputs for validity.

    do n = 1, nspec
       if (state % xn(n) .lt. minx) then
          call print_state(state)
          call bl_error('EOS: mass fraction less than minimum possible mass fraction.')
       else if (state % xn(n) .gt. maxx) then
          call print_state(state)
          call bl_error('EOS: mass fraction more than maximum possible mass fraction.')
       endif
    enddo

    if (state % y_e .lt. minye) then
       call print_state(state)
       call bl_error('EOS: y_e less than minimum possible electron fraction.')
    else if (state % y_e .gt. maxye) then
       call print_state(state)
       call bl_error('EOS: y_e greater than maximum possible electron fraction.')
    endif

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

    !$acc routine seq

    implicit none

    type (eos_t), intent(in) :: state

    if (state % rho .lt. mindens) then
       call print_state(state)
       call bl_error('EOS: rho smaller than mindens.')
    else if (state % rho .gt. maxdens) then
       call print_state(state)
       call bl_error('EOS: rho greater than maxdens.')
    endif

  end subroutine check_rho



  subroutine check_T(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(in) :: state

    if (state % T .lt. mintemp) then
       call print_state(state)
       call bl_error('EOS: T smaller than mintemp.')
    else if (state % T .gt. maxdens) then
       call print_state(state)
       call bl_error('EOS: T greater than maxtemp.')
    endif

  end subroutine check_T



  subroutine check_e(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(in) :: state

    if (state % e .lt. mine) then
       call print_state(state)
       call bl_error('EOS: e smaller than mine.')
    else if (state % e .gt. maxe) then
       call print_state(state)
       call bl_error('EOS: e greater than maxe.')
    endif

  end subroutine check_e



  subroutine check_h(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(in) :: state

    if (state % h .lt. minh) then
       call print_state(state)
       call bl_error('EOS: h smaller than minh.')
    else if (state % h .gt. maxh) then
       call print_state(state)
       call bl_error('EOS: h greater than maxh.')
    endif

  end subroutine check_h



  subroutine check_s(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(in) :: state

    if (state % s .lt. mins) then
       call print_state(state)
       call bl_error('EOS: s smaller than mins.')
    else if (state % s .gt. maxs) then
       call print_state(state)
       call bl_error('EOS: s greater than maxs.')
    endif

  end subroutine check_s



  subroutine check_p(state)

    !$acc routine seq

    implicit none

    type (eos_t), intent(in) :: state

    if (state % p .lt. minp) then
       call print_state(state)
       call bl_error('EOS: p smaller than minp.')
    else if (state % p .gt. maxp) then
       call print_state(state)
       call bl_error('EOS: p greater than maxp.')
    endif

  end subroutine check_p
#endif

end module eos_module
