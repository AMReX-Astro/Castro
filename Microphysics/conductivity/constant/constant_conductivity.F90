module actual_conductivity_module

  use amrex_fort_module, only : rt => amrex_real                                
  use eos_type_module, only : eos_t
  
  implicit none

contains

  subroutine actual_conductivity_init()

    implicit none

  end subroutine actual_conductivity_init

  subroutine actual_conductivity(eos_state)
    
    use extern_probin_module, only: const_conductivity

    implicit none

    !$gpu

    type (eos_t), intent(inout) :: eos_state

    eos_state % conductivity = const_conductivity
    
  end subroutine actual_conductivity

end module actual_conductivity_module
