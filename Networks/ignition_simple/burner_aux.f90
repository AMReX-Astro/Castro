! the burner_aux module contains thermodynamic state variables that are
! needed in the RHS and Jacobian routines

module burner_aux_module

  use bl_types
  use network, only : nspec

  implicit none

  real(kind=dp_t), save :: dens_pass
  real(kind=dp_t), save :: c_p_pass
  real(kind=dp_t), save :: dhdx_pass(nspec)
  real(kind=dp_t), save :: X_O16_pass

end module burner_aux_module
