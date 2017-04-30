module probdata_module

      use amrex_fort_module, only : rt => amrex_real
      real(rt)         Pi
      parameter (Pi=3.1415926535897932384e0_rt)

      real(rt)        , save :: rhocv, T0, Eexp, rexp
      
      real(rt)        , save :: xmin,xmax,ymin,ymax,zmin,zmax

end module probdata_module
