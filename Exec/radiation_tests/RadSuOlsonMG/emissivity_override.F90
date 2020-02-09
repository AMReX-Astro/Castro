module emissivity_override_module

  ! Allow the user to override the specification for the emissivity.

  implicit none

contains

  subroutine emissivity_override(i, j, k, g, T, kg, dkdT, jg, djdT) &
                                 bind(C, name='emissivity_override')

    use amrex_fort_module, only: rt => amrex_real
    use rad_params_module, only: arad

    implicit none

    integer,  intent(in   ) :: i, j, k, g
    real(rt), intent(in   ) :: T, kg, dkdT
    real(rt), intent(inout) :: jg, djdT

    real(rt), parameter :: pfc(0:1) = [0.5_rt, 0.5_rt]

    real(rt) :: Bg, dBdT

    !$gpu

    Bg = arad * T**4
    dBdT = 4.e0_rt * arad * T**3

    jg = pfc(g) * Bg * kg
    djdT = pfc(g) * (dkdT * Bg + dBdT * kg)

  end subroutine emissivity_override

end module emissivity_override_module
