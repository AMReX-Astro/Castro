module probdata_module

      use bl_fort_module, only : rt => c_real
      real(rt)        , save ::  p0, rho0, t_r
      real(rt)        , save :: x_r, q_r
      integer, save :: nsub

end module probdata_module
