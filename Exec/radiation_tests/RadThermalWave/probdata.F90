module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  real(rt), allocatable, save :: rhocv, T0, Eexp, rexp
  real(rt), allocatable, save :: xmin, xmax, ymin, ymax, zmin, zmax

#ifdef AMREX_USE_CUDA
  attributes(managed) :: rhocv, T0, Eexp, rexp
  attributes(managed) :: xmin, xmax, ymin, ymax, zmin, zmax
#endif

end module probdata_module
