module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), allocatable :: p_l, u_l, v_l, rho_l, p_r, u_r, v_r, rho_r, rhoe_l, rhoe_r, T_l, T_r
  logical, allocatable :: use_Tinit

#ifdef AMREX_USE_CUDA
  attributes(managed) :: p_l, u_l, v_l, rho_l, p_r, u_r, v_r, rho_r, rhoe_l, rhoe_r, frac, T_l, T_r
  attributes(managed) :: use_Tinit
#endif
      
end module probdata_module
