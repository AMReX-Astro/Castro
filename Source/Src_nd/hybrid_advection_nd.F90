module hybrid_advection_module

  implicit none

contains

  ! Convert a linear momentum into a "hybrid" momentum that has
  ! an angular momentum component.

  function linear_to_hybrid_momentum(loc, mom_in) result(mom_out)

    implicit none

    double precision :: loc(3), mom_in(3), mom_out(3)

    double precision :: R

    R = sqrt( loc(1)**2 + loc(2)**2 )
    
    mom_out(1) = mom_in(1) * (loc(1) / R) + mom_in(2) * (loc(2) / R)
    mom_out(2) = mom_in(2) * loc(1)       - mom_in(1) * loc(2)
    mom_out(3) = mom_in(3)

  end function linear_to_hybrid_momentum



  ! Convert a "hybrid" momentum into a linear momentum.

  function hybrid_to_linear_momentum(loc, mom_in) result(mom_out)

    implicit none

    double precision :: loc(3), mom_in(3), mom_out(3)

    double precision :: R

    R = sqrt( loc(1)**2 + loc(2)**2 )
    
    mom_out(1) = mom_in(1) * (loc(1) / R)    - mom_in(2) * (loc(2) / R**2)
    mom_out(2) = mom_in(2) * (loc(1) / R**2) + mom_in(1) * (loc(2) / R)
    mom_out(3) = mom_in(3)
    
  end function hybrid_to_linear_momentum

end module hybrid_advection_module

