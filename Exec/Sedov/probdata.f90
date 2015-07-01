module probdata_module

!     Sod variables
      double precision, save ::  p_ambient, dens_ambient, exp_energy
      double precision, save ::  r_init
      integer         , save ::  nsub

!     These help specify which specific problem
      integer         , save :: probtype,idir

end module probdata_module
