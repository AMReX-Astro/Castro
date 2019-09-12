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

  pure function dfvdy(y) result (df_vdy)

    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
    use amrex_fort_module, only : rt => amrex_real

    real(rt), intent(in) :: y
    real(rt) :: df_vdy

    if (y < 1.9375e0_rt * 4.e8_rt) then
       df_vdy = ZERO

    else if (y > 2.0625e0_rt * 4.e8_rt) then
       df_vdy = ZERO

    else
       df_vdy = HALF*8.e0_rt*M_PI * cos(8.e0_rt * M_PI* (y/4.e8_rt - 2.e0_rt))/4.e8_rt
    endif

  end function dfvdy

  function dUdy(y, U) result (dU)

    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
    use amrex_fort_module, only : rt => amrex_real
    use eos_type_module
    use eos_module
    use prescribe_grav_module, only : grav_zone
    use probdata_module, only: gamma1

    real(rt), intent(in) :: y, U(2)
    real(rt) :: dU(2), gamma, gamma0

    type(eos_t) :: eos_state

    ! U(1) = log(rho)
    ! U(2) = log(p)

    eos_state % rho = exp(U(1))
    eos_state % p = exp(U(2))
    eos_state % xn = set_species(y)

    call eos(eos_input_rp, eos_state)

    gamma0 = eos_state % gam1
    gamma = gamma0 + fv(y) * (gamma1 - gamma0)

    dgdy = dfvdy(y) * (gamma1 - gamma0)

    ! dlog p / dy
    dU(2) = exp(U(1)) * grav_zone(y) / exp(U(2))

    ! this follows from gamma = dlnp/dln rho
    dU(1) = dU(2) / gamma

  end function dUdy

  subroutine integrate_model(hi, rho0, p0, ymin, dy, dens, pres)

    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE, TWO
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : sdc_order

    implicit none

    integer, intent(in) :: hi
    real(rt), intent(in) :: rho0, p0, ymin, dy
    real(rt), intent(inout) :: dens(0:hi), pres(0:hi)

    real(rt) :: ystart, y, k1(2), k2(2), k3(2), k4(2)
    real(rt) :: U_old(2), U_new(2), h
    integer :: j

    U_old(1) = log(rho0)
    U_old(2) = log(p0)

    if (sdc_order /= 4) then

       ! do HSE using RK2
       do j = 0, hi
          y = ymin + dy*(dble(j) + HALF)

          ! our integration starts at y - h
          if (j .eq. 0) then
             h = dy * HALF
          else
             h = dy
          endif

          k1(:) = dUdy(y - h, U_old)
          U_new(:) = U_old(:) + h * dUdy(y - HALF*h, U_old + HALF*h * k1)

          dens(j) = exp(U_new(1))
          pres(j) = exp(U_new(2))

          U_old(:) = U_new(:)

       end do

    else

       ! do HSE using RK4
       do j = 0, hi
          y = ymin + dy*(dble(j) + HALF)

          ! our integration starts at y - h
          if (j .eq. 0) then
             h = dy * HALF
          else
             h = dy
          endif

          ystart = y - h

          k1(:) = dUdy(ystart, U_old)
          U_new(:) = U_old(:) + HALF*h * k1(:)

          k2(:) = dUdy(ystart + HALF*h, U_new)
          U_new(:) = U_old(:) + HALF*h * k2(:)

          k3(:) = dUdy(ystart + HALF*h, U_new)
          U_new(:) = U_old(:) + h * k3(:)

          k4(:) = dUdy(ystart + h, U_new)

          U_new = U_old(:) + (1.0_rt/6.0_rt) * h * (k1(:) + TWO*k2(:) + TWO*k3(:) + k4(:))

          dens(j) = exp(U_new(1))
          pres(j) = exp(U_new(2))

          U_old(:) = U_new(:)

       end do

    end if

  end subroutine integrate_model
  
end module model_util_module
