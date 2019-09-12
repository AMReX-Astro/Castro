module model_util_module

  implicit none

contains

  pure function set_species(y) result (xn)

    use network, only : nspec
    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
    use amrex_fort_module, only : rt => amrex_real

    real(rt), intent(in) :: y
    real(rt) :: xn(nspec)

    xn(:) = ZERO
    xn(1) = ONE - fv(y)
    xn(2) = fv(y)

  end function set_species

  pure function fv(y) result (f_v)

    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
    use amrex_fort_module, only : rt => amrex_real

    real(rt), intent(in) :: y
    real(rt) :: f_v

    if (y < 1.9375e0_rt * 4.e8_rt) then
       f_v = ZERO

    else if (y > 2.0625e0_rt * 4.e8_rt) then
       f_v = ONE

    else
       f_v = HALF * (ONE + sin(8.e0_rt * M_PI * (y/4.e8_rt - 2.e0_rt)))
    endif

  end function fv

  function dUdy(y, U) result (dU)

    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
    use amrex_fort_module, only : rt => amrex_real
    use eos_type_module
    use eos_module
    use prescribe_grav_module, only : grav_zone
    use probdata_module, only: gamma1
    use meth_params_module, only : T_guess

    implicit none

    real(rt), intent(in) :: y, U(2)
    real(rt) :: dU(2), gamma, gamma0, dgdy

    type(eos_t) :: eos_state

    ! U(1) = rho
    ! U(2) = p

    eos_state % rho = U(1)
    eos_state % p = U(2)
    eos_state % xn = set_species(y)
    eos_state % T = T_guess

    call eos(eos_input_rp, eos_state)

    gamma0 = eos_state % gam1
    gamma = gamma0 + fv(y) * (gamma1 - gamma0)

    ! dp / dy
    dU(2) = U(1) * grav_zone(y)

    ! drho / dy; this follows from gamma = dlnp/dln rho
    dU(1) = U(1) * dU(2) / (gamma * U(2))

  end function dUdy

  subroutine integrate_model(ny, ymin, ymax, rho0, p0)

    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE, TWO
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : sdc_order, T_guess
    use model_parser_module
    use network, only : nspec
    use eos_type_module, only : eos_t, eos_input_rp
    use eos_module, only : eos

    implicit none

    integer, intent(in) :: ny
    real(rt), intent(in) :: rho0, p0, ymin, ymax

    real(rt) :: ystart, y, dy, k1(2), k2(2), k3(2), k4(2)
    real(rt) :: U_old(2), U_new(2), h
    integer :: j
    type (eos_t) :: eos_state

    ! allocate the storage in the model_parser_module
    npts_model = ny
    allocate (model_state(npts_model, nvars_model))
    allocate (model_r(npts_model))

    ! create the grid -- cell centers
    dy = (ymax - ymin)/ny

    do j = 1, ny
       model_r(j) = ymin + (j - HALF)*dy
    end do

    U_old(1) = rho0
    U_old(2) = p0

    do j = 1, ny
       y = model_r(j)

       ! our integration starts at y - h
       if (j .eq. 0) then
          h = dy * HALF
       else
          h = dy
       endif

       ystart = y - h


       if (sdc_order /= 4) then

          ! do HSE using RK2

          k1(:) = dUdy(ystart, U_old)
          U_new(:) = U_old(:) + h * dUdy(ystart + HALF*h, U_old + HALF*h * k1)

       else

          ! do HSE using RK4

          k1(:) = dUdy(ystart, U_old)
          U_new(:) = U_old(:) + HALF*h * k1(:)

          k2(:) = dUdy(ystart + HALF*h, U_new)
          U_new(:) = U_old(:) + HALF*h * k2(:)

          k3(:) = dUdy(ystart + HALF*h, U_new)
          U_new(:) = U_old(:) + h * k3(:)

          k4(:) = dUdy(ystart + h, U_new)

          U_new = U_old(:) + (1.0_rt/6.0_rt) * h * (k1(:) + TWO*k2(:) + TWO*k3(:) + k4(:))

       end if

       model_state(j, idens_model) = U_new(1)
       model_state(j, ipres_model) = U_new(2)
       model_state(j, ispec_model:ispec_model-1+nspec) = set_species(y)

       eos_state % T = T_guess
       eos_state % rho = model_state(j, idens_model)
       eos_state % xn(:) = model_state(j, ispec_model:ispec_model-1+nspec)
       eos_state % p = model_state(j, ipres_model)

       call eos(eos_input_rp, eos_state)

       model_state(j, itemp_model) = eos_state % T

       U_old(:) = U_new(:)

    end do

  end subroutine integrate_model

end module model_util_module
