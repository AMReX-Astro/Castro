module probdata_module

  use amrex_fort_module, only : rt => amrex_real

  integer, save :: nx_model
  real(rt), save :: dens_base, temp_base
  real(rt), save :: pert_width
  logical, save :: do_pert

end module probdata_module
