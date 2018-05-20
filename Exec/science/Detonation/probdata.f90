module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), save :: T_l, T_r, dens, cfrac, w_T, center_T, smallx

  integer, save :: idir

  real(rt), save :: center(3)
      
  integer, save :: ihe4, ic12, io16
  real(rt), save, allocatable :: xn(:)

end module probdata_module
