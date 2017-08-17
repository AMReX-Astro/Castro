module probdata_module

  use network, only: nspec

!     Sod variables
      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , save ::  p_ambient, dens_ambient, exp_energy, temp_ambient, e_ambient
      real(rt)        , save ::  xn_zone(nspec)
      real(rt)        , save ::  r_init
      integer         , save ::  nsub

!     These help specify which specific problem
      integer         , save :: probtype,idir

end module probdata_module
