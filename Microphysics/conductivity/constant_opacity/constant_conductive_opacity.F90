module actual_conductivity_module

  use amrex_fort_module, only : rt => amrex_real                                
  use eos_type_module, only : eos_t
  use fundamental_constants_module

  implicit none

contains

  subroutine actual_conductivity_init()

    implicit none

  end subroutine actual_conductivity_init

  subroutine actual_conductivity(eos_state)
    
    use extern_probin_module, only: const_opacity

    type (eos_t), intent(inout) :: eos_state

    !$gpu

    eos_state % conductivity = (16*sigma_SB*(eos_state%T)**3)/(3*const_opacity*eos_state%rho)
    
  end subroutine actual_conductivity

end module actual_conductivity_module
