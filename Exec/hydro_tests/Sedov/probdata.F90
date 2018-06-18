module probdata_module

  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  real(rt), allocatable, save :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient
  real(rt), allocatable, save :: xn_zone(:)
  real(rt), allocatable, save :: r_init
  integer,  allocatable, save :: nsub

#ifdef CUDA
  attributes(managed) :: p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient
  attributes(managed) :: xn_zone
  attributes(managed) :: r_init
  attributes(managed) :: nsub
#endif

end module probdata_module
