module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), allocatable :: T_l, T_r, dens, cfrac, ofrac, w_T, center_T, smallx, vel, grav_acceleration

  integer,  allocatable :: idir

  integer,  allocatable :: ihe4, ic12, io16
  real(rt), allocatable :: xn(:)

  logical,  allocatable :: fill_ambient_bc

  real(rt), allocatable :: ambient_dens
  real(rt), allocatable :: ambient_temp
  real(rt), allocatable :: ambient_comp(:)
  real(rt), allocatable :: ambient_e_l, ambient_e_r

#ifdef AMREX_USE_CUDA
  attributes(managed) :: T_l, T_r, dens, cfrac, ofrac, w_T, center_T, smallx, vel, grav_acceleration
  attributes(managed) :: idir
  attributes(managed) :: ihe4, ic12, io16
  attributes(managed) :: xn
  attributes(managed) :: fill_ambient_bc
  attributes(managed) :: ambient_dens
  attributes(managed) :: ambient_temp
  attributes(managed) :: ambient_comp
  attributes(managed) :: ambient_e_l, ambient_e_r
#endif

end module probdata_module
