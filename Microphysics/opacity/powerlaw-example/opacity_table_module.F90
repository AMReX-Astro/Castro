! This is intended as an example

module opacity_table_module

  implicit none

contains

!     ============================================
!
!     opacity for a given temperature and density
!
!     input:   nu    - radiation frequency (Hz)
!              temp  - temperature  (K)
!              rho   - density    (g cm^-3)
!              rhoYe - rho*Ye
!     output:  kp    - Planck mean opacity  (1/cm)
!              kr    - Rosseland mean opacity  (1/cm)

  subroutine get_opacities(kp, kr, rho, temp, rhoYe, nu, get_Planck_mean, get_Rosseland_mean)

    implicit none

    logical, intent(in) :: get_Planck_mean, get_Rosseland_mean
    double precision, intent(in) :: rho, temp, rhoYe, nu
    double precision, intent(out) :: kp, kr

    double precision, parameter :: const_kappa_p = 0.33d-4
    double precision, parameter :: kappa_p_exp_m = 1.0d0
    double precision, parameter :: kappa_p_exp_n = 0.0d0
    double precision, parameter :: kappa_p_exp_p = 0.0d0
    double precision, parameter :: const_kappa_r = 0.33d0
    double precision, parameter :: kappa_r_exp_m = 1.0d0
    double precision, parameter :: kappa_r_exp_n = 0.0d0
    double precision, parameter :: kappa_r_exp_p = 0.0d0
    double precision, parameter :: tiny = 1.d-50
    double precision, parameter :: tfloor = 0.d0
    double precision, parameter :: kfloor = 0.d0

    double precision :: teff, nup_kpp, nup_kpr

    !$gpu

    nup_kpp = nu**kappa_p_exp_p
    nup_kpr = nu**kappa_r_exp_p

    teff = max(temp, tiny)
    teff = teff + tfloor * exp(-teff / (tfloor + tiny))

    if (get_planck_mean) then
       kp = const_kappa_p * (rho**kappa_p_exp_m) * (teff**(-kappa_p_exp_n)) * nup_kpp
    end if

    if (get_Rosseland_mean) then
       kr = const_kappa_r * (rho**kappa_r_exp_m) * (teff**(-kappa_r_exp_n)) * nup_kpr
       kr = max(kr, kfloor)
    end if

  end subroutine get_opacities


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! dummies

  subroutine get_opacity_emissivity( &
       ab, sc, delta, eta,           &
       rho, rdummy, temp, nu, idummy, comp_ab, comp_sc, ldummy)

    implicit none

    integer, intent(in) :: idummy
    logical, intent(in) :: comp_ab, comp_sc, ldummy
    double precision, intent(in) :: rho, rdummy, temp, nu
    double precision, intent(out) :: ab, sc, delta, eta

    !$gpu

    ab = 0.d0
    sc = 0.d0
    delta = 0.d0
    eta = 0.d0

  end subroutine get_opacity_emissivity


  subroutine prep_opacity(g, inu, er, der)

    implicit none

    integer, intent(in) :: g
    integer, intent(out) :: inu
    double precision, intent(out) :: er, der

    inu = -1
    er = 0.d0
    der = 0.d0

  end subroutine prep_opacity

end module opacity_table_module


subroutine init_opacity_table(iverb)

  use opacity_table_module

  implicit none

  integer :: iverb

end subroutine init_opacity_table

