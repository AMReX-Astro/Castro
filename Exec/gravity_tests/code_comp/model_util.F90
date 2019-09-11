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

    ! dlog p / dy
    dU(2) = exp(U(1)) * grav_zone(y) / exp(U(2))

    ! this follows from p = A rho**gamma
    dU(1) = dU(2) / gamma

  end function dUdy

end module model_util_module
