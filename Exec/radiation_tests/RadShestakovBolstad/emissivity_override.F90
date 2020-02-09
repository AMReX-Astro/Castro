module emissivity_override_module

  ! Allow the user to override the specification for the emissivity.

  implicit none

contains

  subroutine emissivity_override(i, j, k, g, T, kg, dkdT, jg, djdT) &
                                 bind(C, name='emissivity_override')

    use amrex_fort_module, only: rt => amrex_real
    use rad_params_module, only: nugroup, xnu, hplanck, kboltz, clight, pi

    implicit none

    integer,  intent(in   ) :: i, j, k, g
    real(rt), intent(in   ) :: T, kg, dkdT
    real(rt), intent(inout) :: jg, djdT

    real(rt), parameter :: Tf = 1.1604518737857157e6_rt ! 1 (dimensionless)
    real(rt) :: Bg, dBdT, hoverk, cB, nu, num, nup

    !$gpu

    hoverk = hplanck / kboltz
    cB = 8.0_rt * pi * kboltz / clight**3

    nu = nugroup(g)
    num = xnu(g)
    nup = xnu(g+1)

    dBdT = cB * nu**3 * (exp(-hoverk * num / Tf) - exp(-hoverk * nup / Tf))
    Bg = dBdT * T

    jg = Bg * kg
    djdT = dkdT * Bg + dBdT * kg

  end subroutine emissivity_override

end module emissivity_override_module
