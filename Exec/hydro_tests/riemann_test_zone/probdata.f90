module probdata_module

  use amrex_fort_module, only : rt => c_real
  real(rt)        , save ::  rho_l, u_l, p_l, re_l, gc_l
  real(rt)        , save ::  rho_r, u_r, p_r, re_r, gc_r
  real(rt)        , save ::  cav_s, smallc_s

end module probdata_module
