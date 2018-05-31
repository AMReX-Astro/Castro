module probdata_module

!     Sod variables
      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , save ::  p_l, u_l, v_l, rho_l, p_r, u_r, v_r, rho_r, rhoe_l, rhoe_r, frac
      real(rt)        , save :: T_l, T_r

      logical, save :: use_Tinit

!     These help specify which specific problem
      integer        , save ::  probtype,idir

      real(rt)        , save :: split(3)
      
end module probdata_module
