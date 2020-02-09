module emissivity_override_module

  ! Allow the user to override the specification for the emissivity.

  implicit none

contains

  subroutine emissivity_override(i, j, k, g, T, kg, dkdT, jg, djdT) &
                                 bind(C, name='emissivity_override')

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: i, j, k, g
    real(rt), intent(in   ) :: T, kg, dkdT
    real(rt), intent(inout) :: jg, djdT

    !$gpu

  end subroutine emissivity_override

end module emissivity_override_module
