module binary_module

  ! This module contains a number of routines that are designed to
  ! calculate generic properties of binary orbits.

  use amrex_fort_module, only: rt => amrex_real

  implicit none

contains

  ! Given the mass ratio q of two stars (assumed to be q = M_1 / M_2), 
  ! compute the effective Roche radii of the stars, normalized to unity, 
  ! using the approximate formula of Eggleton (1983). Optionally we can
  ! pass in a distance scale.
  
  subroutine get_roche_radii(mass_ratio, r_1, r_2, a)

    use amrex_constants_module, only: ONE, TWO3RD, THIRD

    implicit none

    real(rt), intent(in   ) :: mass_ratio
    real(rt), intent(inout) :: r_1, r_2
    real(rt), intent(in   ), optional :: a

    real(rt) :: q
    real(rt) :: c1, c2

    real(rt) :: scale

    if (present(a)) then
       scale = a
    else
       scale = ONE
    endif

    c1 = 0.49e0_rt
    c2 = 0.60e0_rt

    q = mass_ratio

    r_1 = scale * c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

    q = ONE / q

    r_2 = scale * c1 * q**(TWO3RD) / (c2 * q**(TWO3RD) + LOG(ONE + q**(THIRD)))

  end subroutine get_roche_radii

end module binary_module
