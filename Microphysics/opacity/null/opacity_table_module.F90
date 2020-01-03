module opacity_table_module

  implicit none

contains

  subroutine prep_opacity(g, inu, er, der)

    implicit none

    integer, intent(in) :: g
    integer, intent(out) :: inu
    double precision, intent(out) :: er, der

    inu = -1
    er = 0.d0
    der = 0.d0

  end subroutine prep_opacity

  subroutine get_opacity_emissivity( &
       ab, sc, delta, eta,           &
       rho, yein, temp, er, inu, comp_ab, comp_sc, comp_eta)

    implicit none

    integer, intent(in) :: inu 
    logical, intent(in) :: comp_ab, comp_sc, comp_eta
    double precision, intent(in) :: rho, yein, temp, er 
    double precision, intent(out) :: ab, sc, delta, eta

    !$gpu

    ab = 0.d0
    sc = 0.d0
    delta = 0.d0
    eta = 0.d0

  end subroutine get_opacity_emissivity

  subroutine get_opacities(kp, kr, rho, temp, rhoYe, nu, get_Planck_mean, get_Rosseland_mean)

    implicit none

    logical, intent(in) :: get_Planck_mean, get_Rosseland_mean
    double precision, intent(in) :: rho, temp, rhoYe, nu
    double precision, intent(out) :: kp, kr

    !$gpu

    kp = 0.d0
    kr = 0.d0

  end subroutine get_opacities

end module opacity_table_module


subroutine init_opacity_table(iverb)

  implicit none

  integer :: iverb

end subroutine init_opacity_table

