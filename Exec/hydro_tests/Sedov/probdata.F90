module probdata_module

  use amrex_fort_module, only: rt => amrex_real

  real(rt), allocatable :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient
  real(rt), allocatable :: xn_zone(:)
  real(rt), allocatable :: r_init
  integer,  allocatable :: nsub

  real(rt), allocatable :: e_exp ! Explosion energy per zone

#ifdef AMREX_USE_CUDA
  attributes(managed) :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient
  attributes(managed) :: xn_zone, r_init, nsub, e_exp
#endif

end module probdata_module
