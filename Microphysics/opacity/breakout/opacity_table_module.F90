module opacity_table_module

  implicit none

contains

!     ============================================
!
!     opacity for a given temperature and density
!
!     input:   nu    - radiation frequency (Hz)
!              temp  - temperature  (K)
!              rhoYe - rho*Ye
!              rho   - density    (g cm^-3)
!     output:  kp    - Planck mean opacity  (1/cm)
!              kr    - Rosseland mean opacity  (1/cm)

  subroutine get_opacities(kp, kr, rho, temp, rhoYe, nu, get_Planck_mean, get_Rosseland_mean)

    implicit none

    logical, intent(in) :: get_Planck_mean, get_Rosseland_mean
    double precision, intent(in) :: rho, temp, rhoYe, nu
    double precision, intent(out) :: kp, kr

    double precision, parameter :: Ksc = 0.4d0 ! Thomson scattering 
    double precision, parameter :: fac = 1.d-4 ! Planck mean is assumed to be fac*Ksc

    !$gpu

    if (get_planck_mean) then
       kp = rhoYe*Ksc * fac
    end if

    if (get_Rosseland_mean) then
       kr = rhoYe*Ksc
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


! dummy
subroutine init_opacity_table(iverb)

  implicit none

  integer :: iverb

end subroutine init_opacity_table


