module probdata_module

!     Sod variables
      double precision, save ::  p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, frac
      double precision, save :: T_l, T_r

      logical, save :: use_Tinit

!     These help specify which specific problem
      integer        , save ::  probtype,idir

      double precision, save :: prob_center(3)

      
end module probdata_module
