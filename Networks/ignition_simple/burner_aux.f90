! the burner_aux module contains thermodynamic state variables that are
! needed in the RHS and Jacobian routines

module burner_aux_module

  use bl_types
  use network, only : nspec

  implicit none

  real(kind=dp_t) :: dens_pass
  real(kind=dp_t) :: c_p_pass
  real(kind=dp_t) :: dhdx_pass(nspec)
  real(kind=dp_t) :: X_O16_pass

  !$OMP THREADPRIVATE(dens_pass,c_p_pass,dhdx_pass,X_O16_pass)
  
end module burner_aux_module
