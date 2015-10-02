module conductivity_module

  use bl_types
  use network
  use eos_type_module
  use eos_module

  implicit none

contains

  subroutine thermal_conductivity(eos_state, therm_cond)
    
    use extern_probin_module, only: const_conductivity

    type (eos_t), intent(in) :: eos_state
    real (kind=dp_t), intent(inout) :: therm_cond

    therm_cond = const_conductivity
    
  end subroutine thermal_conductivity

end module conductivity_module
