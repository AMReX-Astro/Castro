module model_util_module

  implicit none

contains

  pure function set_species(y) result (xn)

    use network, only : nspec
    use amrex_constants_module, only: HALF, ZERO, M_PI, ONE
    use amrex_fort_module, only : rt => amrex_real

    real(rt), intent(in) :: y
    real(rt) :: xn(nspec)
    real(rt) :: fv

    xn(:) = ZERO

    if (y < 1.9375e0_rt * 4.e8_rt) then
       xn(1) = ONE
       xn(2) = ZERO

    else if (y > 2.0625e0_rt * 4.e8_rt) then
       xn(1) = ZERO
       xn(2) = ONE

    else
       fv = HALF * (ONE + sin(8.e0_rt * M_PI * (y/4.e8_rt - 2.e0_rt)))
       xn(1) = ONE - fv
       xn(2) = fv
    endif

  end function set_species

end module model_util_module
