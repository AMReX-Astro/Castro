module viscosity_module

  use amrex_fort_module, only : rt => amrex_real
  use network
  use eos_type_module
  use eos_module

  implicit none

contains

  subroutine viscous_coeff(eos_state, visc_coeff)
    
    use extern_probin_module, only: const_viscosity

    type (eos_t), intent(in) :: eos_state
    real (rt), intent(inout) :: visc_coeff

    visc_coeff = const_viscosity
    
  end subroutine viscous_coeff

end module viscosity_module
