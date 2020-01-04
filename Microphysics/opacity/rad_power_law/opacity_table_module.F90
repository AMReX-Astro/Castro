! This is an artificial opacity designed for use with radiation
! testing problems.

module opacity_table_module

  use amrex_fort_module, only: rt => amrex_real
  use extern_probin_module, only: const_kappa_p, kappa_p_exp_m, kappa_p_exp_n, kappa_p_exp_p, &
                                  const_kappa_r, kappa_r_exp_m, kappa_r_exp_n, kappa_r_exp_p, &
                                  const_scatter, scatter_exp_m, scatter_exp_n, scatter_exp_p, &
                                  kappa_floor, rad_temp_floor

  implicit none

  real(rt), parameter, private :: tiny = 1.e-50_rt

contains

  ! ============================================
  !
  ! opacity for a given temperature and density
  !
  ! input:   nu    - radiation frequency (Hz)
  !          temp  - temperature  (K)
  !          rho   - density    (g cm^-3)
  !          rhoYe - rho*Ye
  ! output:  kp    - Planck mean opacity  (1/cm)
  !          kr    - Rosseland mean opacity  (1/cm)

  subroutine get_opacities(kp, kr, rho, temp, rhoYe, nu, get_Planck_mean, get_Rosseland_mean)

#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    implicit none

    logical,  intent(in   ) :: get_Planck_mean, get_Rosseland_mean
    real(rt), intent(in   ) :: rho, temp, rhoYe, nu
    real(rt), intent(inout) :: kp, kr

    real(rt) :: teff, nup_kpp, nup_kpr, nup_kps, ks

    !$gpu

    nup_kpp = nu**kappa_p_exp_p
    nup_kpr = nu**kappa_r_exp_p
    nup_kps = nu**scatter_exp_p

    teff = max(temp, tiny)
    teff = teff + rad_temp_floor * exp(-teff / (rad_temp_floor + tiny))

    ks = const_scatter * (rho**scatter_exp_m) * (teff**(-scatter_exp_n)) * nup_kps

    if (get_Planck_mean) then
#ifndef AMREX_USE_GPU
       if (const_kappa_p < 0.0_rt) then
          call castro_error("Must set Planck opacity constant")
       end if
#endif
       kp = const_kappa_p * (rho**kappa_p_exp_m) * (teff**(-kappa_p_exp_n)) * nup_kpp
       kp = max(kp, kappa_floor)
    end if

    if (get_Rosseland_mean) then
#ifndef AMREX_USE_GPU
       if (const_kappa_r < 0.0_rt) then
          call castro_error("Must set Rosseland opacity constant")
       end if
#endif
       kr = const_kappa_r * (rho**kappa_r_exp_m) * (teff**(-kappa_r_exp_n)) * nup_kpr
       kr = max(kr + ks, kappa_floor)
    end if

  end subroutine get_opacities



  ! dummies

  subroutine get_opacity_emissivity( &
       ab, sc, delta, eta,           &
       rho, rdummy, temp, nu, idummy, comp_ab, comp_sc, ldummy)

    implicit none

    integer,  intent(in   ) :: idummy
    logical,  intent(in   ) :: comp_ab, comp_sc, ldummy
    real(rt), intent(in   ) :: rho, rdummy, temp, nu
    real(rt), intent(inout) :: ab, sc, delta, eta

    !$gpu

    ab = 0.d0
    sc = 0.d0
    delta = 0.d0
    eta = 0.d0

  end subroutine get_opacity_emissivity



  subroutine prep_opacity(g, inu, er, der)

    implicit none

    integer,  intent(in   ) :: g
    integer,  intent(inout) :: inu
    real(rt), intent(inout) :: er, der

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
