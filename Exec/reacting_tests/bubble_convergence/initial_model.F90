! create a simple hydrostatic initial model that is 4th order accurate
! at cell-centers by using RK4 to integrate:
!
!   dp/dr = rho(p, s) g
!   ds/dr = 0

module initial_model_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec
  use model_parser_module

  implicit none

  type :: model_t
     real(rt) :: xn(nspec)
     real(rt) :: dens_base
     real(rt) :: T_base
  end type model_t

contains

  subroutine generate_initial_model(nx, xmin, xmax, model_params)

    use amrex_constants_module
    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt, eos_input_ps
    use network, only : nspec
    use meth_params_module, only : const_grav

    use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

    implicit none

    integer, intent(in) :: nx
    type(model_t), intent(in) :: model_params
    real(rt) :: xmin, xmax

    integer :: i

    real(rt) :: entropy_fixed

    type (eos_t) :: eos_state

    real(rt) :: k1, k2, k3, k4
    real(rt) :: pnew, p, s, rho, T, xn(nspec)
    real(rt) :: dx

    ! allocate the storage in the model_parser_module
    npts_model = nx

    allocate (model_state(npts_model, nvars_model))
    allocate (model_r(npts_model))

    ! get the base conditions
    eos_state % rho = model_params % dens_base
    eos_state % T = model_params % T_base
    eos_state % xn(:) = model_params % xn(:)

    call eos(eos_input_rt, eos_state)

    entropy_fixed = eos_state % s

    model_state(1, idens_model) = eos_state % rho
    model_state(1, itemp_model) = eos_state % T
    model_state(1, ipres_model) = eos_state % p
    model_state(1, ispec_model:ispec_model-1+nspec) = eos_state % xn(:)

    ! create the grid -- cell centers
    dx = (xmax - xmin)/nx

    do i = 1, nx
       model_r(i) = xmin + (i - HALF)*dx
    end do

    ! do RK 4 integration
    do i = 2, nx

       rho = model_state(i-1, idens_model)
       T = model_state(i-1, itemp_model)
       p = model_state(i-1, ipres_model)
       xn(:) = model_state(i-1, ispec_model:ispec_model-1+nspec)

       ! entropy never changes in this model
       s = entropy_fixed

       k1 = f(p, s, const_grav, rho, T, xn)
       k2 = f(p + HALF*dx*k1, s, const_grav, rho, T, xn)
       k3 = f(p + HALF*dx*k2, s, const_grav, rho, T, xn)
       k4 = f(p + dx*k3, s, const_grav, rho, T, xn)

       pnew = p + SIXTH*dx*(k1 + TWO*k2 + TWO*k3 + k4)


       ! call the EOS to get the remainder of the thermodynamics
       eos_state % T     = model_state(i-1, itemp_model) ! initial guess
       eos_state % rho   = model_state(i-1, idens_model) ! initial guess
       eos_state % xn(:) = model_state(i-1, ispec_model:ispec_model-1+nspec)
       eos_state % p = pnew
       eos_state % s = s

       call eos(eos_input_ps, eos_state)

       ! update the thermodynamics in this zone
       model_state(i, idens_model) = eos_state % rho
       model_state(i, itemp_model) = eos_state % T
       model_state(i, ipres_model) = eos_state % p
       model_state(i, ispec_model:ispec_model-1+nspec) = eos_state % xn(:)

    enddo

  end subroutine generate_initial_model

  function f(p, s, g, rho, T, xn) result (rhs)

    use eos_type_module, only : eos_t, eos_input_ps
    use eos_module, only : eos

    real(rt), intent(in) :: p, s, g, rho, T, xn(nspec)

    real(rt) :: rhs
    type(eos_t) :: eos_state

    eos_state % T = T ! initial guess
    eos_state % rho = rho ! initial guess
    eos_state % xn(:) = xn(:)
    eos_state % p = p
    eos_state % s = s

    call eos(eos_input_ps, eos_state)

    rhs = eos_state % rho * g

  end function f

end module initial_model_module
