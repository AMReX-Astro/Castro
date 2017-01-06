module probdata_module

!     Sod variables
      use bl_fort_module, only : rt => c_real
      real(rt)        , save ::  p_ambient, dens_ambient, exp_energy
      real(rt)        , save ::  r_init
      integer         , save ::  nsub

!     These help specify which specific problem
      integer         , save :: probtype,idir

end module probdata_module
