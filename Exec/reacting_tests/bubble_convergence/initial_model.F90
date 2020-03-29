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

  subroutine generate_initial_model(nx, xmin, xmax, model_params, nbuf)

    use amrex_constants_module
    use castro_error_module
    use amrex_fort_module, only : rt => amrex_real

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt, eos_input_ps
    use network, only : nspec
    use meth_params_module, only : const_grav, sdc_order, time_integration_method

    use amrex_paralleldescriptor_module, only: parallel_IOProcessor => amrex_pd_ioprocessor

    implicit none

    integer, intent(in) :: nx, nbuf
    type(model_t), intent(in) :: model_params
    real(rt) :: xmin, xmax

    integer :: i, ibase, itop

    real(rt) :: entropy_fixed, h

    type (eos_t) :: eos_state

    real(rt) :: k1, k2, k3, k4
    real(rt) :: pnew, p, s, rho, T, xn(nspec)
    real(rt) :: dx

    logical :: converged_hse
    real(rt) :: drho, p_want, dens_zone
    integer :: iter
    integer, parameter :: MAX_ITER = 100
    real(rt), parameter :: TOL = 1.e-10_rt

    ! allocate the storage in the model_parser_module
    npts_model = nx + 2*nbuf

    allocate (model_state(npts_model, nvars_model))
    allocate (model_r(npts_model))

    ibase = nbuf + 1
    itop = ibase + nx - 1

    ! get the base conditions
    eos_state % rho = model_params % dens_base
    eos_state % T = model_params % T_base
    eos_state % xn(:) = model_params % xn(:)

    call eos(eos_input_rt, eos_state)

    entropy_fixed = eos_state % s

    ! create the grid -- cell centers
    dx = (xmax - xmin)/nx

    do i = 1, nx + 2*nbuf
       model_r(i) = xmin + (i - HALF - nbuf)*dx
    end do

    ! note, those conditions are the lower boundary.  This means we
    ! need to integrate dx/2 in the first step to get to the first
    ! zone's cell-center.  After that, we integrate from
    ! center-to-center.

    p = eos_state % p


    if (time_integration_method == 3 .and. sdc_order == 4) then

       ! a fourth order accurate method

       ! do RK 4 integration up from the lower boundary
       do i = ibase, itop + nbuf

          ! rho and T here are guesses for the EOS call
          if (i == ibase) then
             rho = eos_state % rho
             T = eos_state % T
          else
             rho = model_state(i-1, idens_model)
             T = model_state(i-1, itemp_model)
          end if

          xn(:) = model_params % xn(:)

          ! step size
          if (i == ibase) then
             h = HALF*dx
          else
             h = dx
          end if

          ! entropy never changes in this model
          s = entropy_fixed

          k1 = f(p, s, const_grav, rho, T, xn)
          k2 = f(p + HALF*h*k1, s, const_grav, rho, T, xn)
          k3 = f(p + HALF*h*k2, s, const_grav, rho, T, xn)
          k4 = f(p + h*k3, s, const_grav, rho, T, xn)

          pnew = p + SIXTH*h*(k1 + TWO*k2 + TWO*k3 + k4)

          ! call the EOS to get the remainder of the thermodynamics
          eos_state % T     = T ! initial guess
          eos_state % rho   = rho ! initial guess
          eos_state % xn(:) = model_params % xn(:)
          eos_state % p = pnew
          eos_state % s = s

          call eos(eos_input_ps, eos_state)

          ! update the thermodynamics in this zone
          model_state(i, idens_model) = eos_state % rho
          model_state(i, itemp_model) = eos_state % T
          model_state(i, ipres_model) = eos_state % p
          model_state(i, ispec_model:ispec_model-1+nspec) = eos_state % xn(:)

          ! reset for the next iteration
          p = pnew

       enddo

       p = model_state(ibase, ipres_model)

       ! now integrate down
       do i = ibase-1, 1, -1

          rho = model_state(i+1, idens_model)
          T = model_state(i+1, itemp_model)
          xn(:) = model_params % xn(:)

          h = -dx

          ! entropy never changes in this model
          s = entropy_fixed

          k1 = f(p, s, const_grav, rho, T, xn)
          k2 = f(p + HALF*h*k1, s, const_grav, rho, T, xn)
          k3 = f(p + HALF*h*k2, s, const_grav, rho, T, xn)
          k4 = f(p + h*k3, s, const_grav, rho, T, xn)

          pnew = p + SIXTH*h*(k1 + TWO*k2 + TWO*k3 + k4)

          ! call the EOS to get the remainder of the thermodynamics
          eos_state % T     = T ! initial guess
          eos_state % rho   = rho ! initial guess
          eos_state % xn(:) = model_params % xn(:)
          eos_state % p = pnew
          eos_state % s = s

          call eos(eos_input_ps, eos_state)

          ! update the thermodynamics in this zone
          model_state(i, idens_model) = eos_state % rho
          model_state(i, itemp_model) = eos_state % T
          model_state(i, ipres_model) = eos_state % p
          model_state(i, ispec_model:ispec_model-1+nspec) = eos_state % xn(:)

          ! reset for the next iteration
          p = pnew

       enddo

    else

       model_state(:, :) = ZERO

       ! a second order accurate scheme
       do i = ibase, itop + nbuf

          if (i == ibase) then
             model_state(i, idens_model) = model_params % dens_base
             model_state(i, ipres_model) = p + HALF*dx* model_params % dens_base * const_grav

             eos_state % rho = model_state(i, idens_model)
             eos_state % T = model_params % dens_base
             eos_state % p = model_state(i, ipres_model)
             eos_state % xn(:) = model_params % xn(:)
             eos_state % s = entropy_fixed

             call eos(eos_input_ps, eos_state)

             model_state(i, idens_model) = eos_state % rho
             model_state(i, itemp_model) = eos_state % T
             model_state(i, ispec_model:ispec_model-1+nspec) = model_params % xn(:)

          else

             dens_zone = model_state(i-1, idens_model)

             converged_hse = .FALSE.

             do iter = 1, MAX_ITER

                p_want = model_state(i-1, ipres_model) + &
                     dx*HALF*(dens_zone + model_state(i-1, idens_model))*const_grav

                ! use the EOS with constant entropy to find corrected state
                eos_state % p = p_want
                eos_state % s = entropy_fixed
                eos_state % xn(:) = model_params % xn(:)

                call eos(eos_input_ps, eos_state)

                drho = eos_state % rho - dens_zone
                dens_zone = eos_state % rho
                if (abs(drho) < TOL*dens_zone) then
                   converged_hse = .TRUE.
                   exit
                end if

             end do

             if (.not. converged_hse) then
                print *, "failed to convergence in initial model generation"
                call castro_error("ERROR")
             end if

             ! initialze zone
             model_state(i, idens_model) = dens_zone
             model_state(i, ipres_model) = p_want
             model_state(i, itemp_model) = eos_state % T
             model_state(i, ispec_model:ispec_model-1+nspec) = model_params % xn(:)

          end if
          print *, i, model_state(i, idens_model), model_state(i, ipres_model)

       end do

       ! now integrate down
       do i = ibase-1, 1, -1

          dens_zone = model_state(i+1, idens_model)

          converged_hse = .FALSE.

          do iter = 1, MAX_ITER

             p_want = model_state(i+1, ipres_model) - &
                  dx*HALF*(dens_zone + model_state(i+1, idens_model))*const_grav

             ! use the EOS with constant entropy to find corrected state
             eos_state % p = p_want
             eos_state % s = entropy_fixed
             eos_state % xn(:) = model_params % xn(:)

             call eos(eos_input_ps, eos_state)

             drho = eos_state % rho - dens_zone
             dens_zone = eos_state % rho
             if (abs(drho) < TOL*dens_zone) then
                converged_hse = .TRUE.
                exit
             end if

          end do

          if (.not. converged_hse) then
             print *, "failed to convergence in initial model generation"
             call castro_error("ERROR")
          end if

          ! initialze zone
          model_state(i, idens_model) = dens_zone
          model_state(i, ipres_model) = p_want
          model_state(i, itemp_model) = eos_state % T
          model_state(i, ispec_model:ispec_model-1+nspec) = model_params % xn(:)

       enddo

    end if

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
