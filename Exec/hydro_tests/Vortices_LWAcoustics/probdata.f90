module probdata_module

!     variables for initialization; see Lee & Koo AIAA Journal 1995
      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , save :: p_ref, r_0, mach, ratio_c, r_circ 
      real(rt)        , save :: rho_0, c_0, r_c, circ
      real(rt)        , save :: x_c1, y_c1, x_c2, y_c2 
 
end module probdata_module
