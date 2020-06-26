module hybrid_advection_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! Avoid the singularity in cylindrical coordinates
  real(rt), parameter :: R_min = epsilon(0.0_rt)

contains

  function linear_to_hybrid(loc, mom_in) result(mom_out)
    ! Convert a linear momentum into the "hybrid" scheme
    ! that has radial and angular components.

    use amrex_constants_module, only: ZERO

    implicit none

    real(rt), intent(in) :: loc(3), mom_in(3)
    real(rt) :: mom_out(3)

    real(rt) :: R

    !$gpu

    R = max(sqrt(loc(1)**2 + loc(2)**2), R_min)

    ! This conversion is Eqs. 25 and 26 in Byerly et al. 2014.
    ! Note that we expect the linear momentum to be consistent
    ! with which frame we're measuring the fluid quantities in.
    ! So we're effectively always using the first form of those
    ! equalities, not the second. If state_in_rotating_frame = 1,
    ! then we're not including the centrifugal term in the angular
    ! momentum anyway, and if state_in_rotating_frame = 0, then
    ! the linear momenta are already expressed in the inertial frame,
    ! so we don't need to explicitly take rotation into account.

    mom_out(1) = mom_in(1) * (loc(1) / R) + mom_in(2) * (loc(2) / R)
    mom_out(2) = mom_in(2) * loc(1) - mom_in(1) * loc(2)
    mom_out(3) = mom_in(3)

  end function linear_to_hybrid

end module hybrid_advection_module
