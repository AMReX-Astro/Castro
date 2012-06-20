module probdata_module

!     These determine the refinement criteria
      double precision, save :: dengrad
      double precision, save :: velgrad

!     KH variables for initializatoin
      double precision, save ::  u_amb, rho_amb, T_amb, u_pert, rho_pert, T_pert

!     These help specify which specific problem
      integer        , save ::  probtype,idir

!     This is needed for other problems
      double precision, save ::  center(3)
      
end module probdata_module
