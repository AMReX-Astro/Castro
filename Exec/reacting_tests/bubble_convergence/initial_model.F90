! create a simple hydrostatic initial model that is 4th order accurate
! at cell-centers by using RK4 to integrate:
!
!   dp/dr = rho(p, s) g
!   ds/dr = 0

module initial_model_module

  use amrex_fort_module, only : rt => amrex_real
  use network, only : nspec

  implicit none

  ! integer keys for indexing the model_state array
  integer, parameter :: nvars_model = 3 + nspec
  integer, parameter :: idens_model = 1
  integer, parameter :: itemp_model = 2
  integer, parameter :: ipres_model = 3
  integer, parameter :: ispec_model = 4

  ! number of points in the model file
  integer, allocatable :: gen_npts_model

  ! arrays for storing the model data -- we have an extra index here
  ! which is the model number
  real(rt), allocatable, save :: gen_model_state(:,:)
  real(rt), allocatable, save :: gen_model_r(:)

  type :: model_t
     real(rt) :: xn(nspec)
     real(rt) :: dens_base
     real(rt) :: T_base
  end type model_t

#ifdef AMREX_USE_CUDA
  attributes(managed) :: gen_npts_model
  attributes(managed) :: gen_model_state, gen_model_r
#endif

contains

  subroutine init_model_data(nx)

    implicit none

    integer, intent(in) :: nx

    allocate(gen_npts_model)

    ! allocate storage for the model data
    allocate (gen_model_state(nx, nvars_model))
    allocate (gen_model_r(nx))

    gen_npts_model = nx

  end subroutine init_model_data


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

    real (rt) :: entropy_fixed

    type (eos_t) :: eos_state


    ! get the base conditions
    eos_state % rho = model_params % dens_base
    eos_state % T = model_params % T_base
    eos_state % xn(:) = model_params % xn(:)

    call eos(eos_input_rt, eos_state)

    entropy_fixed = eos_state % s

    gen_model_state(1, idens_model) = eos_state % rho
    gen_model_state(1, itemp_model) = eos_state % T
    gen_model_state(1, ipres_model) = eos_state % p
    gen_model_state(1, ispec_model:ispec_model-1+nspec) = eos_state % xn(:)

    ! create the grid -- cell centers
    dx = (xmax - xmin)/nx

    do i = 1, nx
       gen_model_r(i) = xmin + (i - HALF)*dx
    end do

    ! do RK 4 integration
    do n = 2, nx

       p = gen_model_state(i-1, idens_model)

       ! entropy never changes in this model
       s = entropy_fixed

       k1 = f(p, s)

       k2 = f(p + HALF*dx*k1)

       k3 = f(p + HALF*dx*k2)

       k4 = f(p + dx*k3)

       pnew = p + SIXTH*(k1 + TWO*k2 + TWO*k3 + k4)


       ! call the EOS to get the remainder of the thermodynamics
       eos_state % T     = gen_model_state(i-1, itemp_model) ! initial guess
       eos_state % rho   = gen_model_state(i-1, idens_model) ! initial guess
       eos_state % xn(:) = xn(:)
       eos_state % p = pnew
       eos_state % s = s

       call eos(eos_input_ps, eos_state)

       ! update the thermodynamics in this zone
       gen_model_state(i, idens_model) = eos_state % rho
       gen_model_state(i, itemp_model) = eos_state % T
       gen_model_state(i, ipres_model) = eos_state % p
       gen_model_state(i, ispec_model:ispec_model-1+nspec) = eos_state % xn(:)

    enddo

  end subroutine generate_initial_model

  function f(p, s, g, rho, T, xn) result (rhs)

    real(rt), intent(in) :: p, s, g, rho, T, xn(nspec)

    real(rt) :: rhs

    eos_state % T = T ! initial guess
    eos_state % rho = rho ! initial guess
    eos_state % xn(:) = xn(:)
    eos_state % p = p
    eos_state % s = s

    call eos(eos_input_ps, eos_state)

    rhs = eos_state % rho * g

  end function f

end module initial_model_module
