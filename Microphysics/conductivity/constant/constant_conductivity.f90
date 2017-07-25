module actual_conductivity_module

  use amrex_fort_module, only : rt => amrex_real                                
  use eos_type_modulem only : eos_t
  
  implicit none

contains

  subroutine actual_conductivity_init()

    implicit none

  end subroutine actual_conductivity_init

  subroutine actual_conductivity(eos_state, therm_cond)
    
    use extern_probin_module, only: const_conductivity

    implicit none

    type (eos_t), intent(in) :: eos_state
    real (kind=rt), intent(inout) :: therm_cond

    therm_cond = const_conductivity
    
  end subroutine actual_conductivity

end module actual_conductivity_module
