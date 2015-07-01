module probdata_module

!     KH variables for initialization; see McNally+ 2012 ApJS 201, 18
      double precision, save ::  rho1, rho2, u1, u2, L, pres
!     These control the amplitude and number of modes in the velocity 
!     perturbation; equivalent to 0.01 and 4 in McNally+ Equation 5
      double precision, save :: vfac, vmode

!     This specifies the direction of the flow; see Prob_?d.f90
      integer        , save ::  idir

end module probdata_module
