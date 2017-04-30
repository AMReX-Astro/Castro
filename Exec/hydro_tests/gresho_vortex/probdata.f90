module probdata_module

      use amrex_fort_module, only : rt => amrex_real
      real(rt)        , save ::  p0, rho0, t_r
      real(rt)        , save :: x_r, q_r
      integer, save :: nsub

end module probdata_module
