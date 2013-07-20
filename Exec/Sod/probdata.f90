module probdata_module

!     These determine the refinement criteria
      double precision, save :: denerr,  dengrad
      double precision, save :: velerr,  velgrad
      double precision, save :: presserr,pressgrad
      double precision, save :: temperr,tempgrad
      double precision, save :: raderr,radgrad
      integer         , save :: max_denerr_lev   ,max_dengrad_lev
      integer         , save :: max_velerr_lev   ,max_velgrad_lev
      integer         , save :: max_presserr_lev, max_pressgrad_lev
      integer         , save :: max_temperr_lev, max_tempgrad_lev
      integer         , save :: max_raderr_lev, max_radgrad_lev

!     Sod variables
      double precision, save ::  p_l, u_l, rho_l, p_r, u_r, rho_r, rhoe_l, rhoe_r, frac
      double precision, save :: T_l, T_r

!     These help specify which specific problem
      integer        , save ::  probtype,idir

      double precision, save ::  center(3)

      
end module probdata_module
