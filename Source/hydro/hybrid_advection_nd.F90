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


  subroutine set_hybrid_momentum_source(loc, mom, source)

    implicit none

    real(rt), intent(in   ) :: loc(3), source(3)
    real(rt), intent(inout) :: mom(3)

    real(rt) :: R

    !$gpu

    R = max(sqrt(loc(1)**2 + loc(2)**2), R_min)

    ! This is analogous to the conversion of linear momentum to hybrid momentum.

    mom(1) =  source(1) * (loc(1) / R) + source(2) * (loc(2) / R)
    mom(2) =  source(2) * loc(1) - source(1) * loc(2)
    mom(3) =  source(3)

  end subroutine set_hybrid_momentum_source




  subroutine add_hybrid_advection_source(lo, hi, dt, &
                                         update, u_lo, u_hi, &
                                         qx, qx_lo, qx_hi, &
                                         qy, qy_lo, qy_hi, &
                                         qz, qz_lo, qz_hi) bind(C, name="add_hybrid_advection_source")

    use meth_params_module, only: NVAR, NGDNV, GDPRES, UMR
    use prob_params_module, only: center, dx_level
    use castro_util_module, only: position ! function
    use amrinfo_module, only: amr_level

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    integer, intent(in) :: qx_lo(3), qx_hi(3)
    integer, intent(in) :: qy_lo(3), qy_hi(3)
    integer, intent(in) :: qz_lo(3), qz_hi(3)
    real(rt), intent(inout) :: update(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    real(rt), intent(in) :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    real(rt), intent(in) :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    real(rt), intent(in) :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    real(rt), intent(in), value :: dt

    integer  :: i, j, k
    real(rt) :: loc(3), R, dx(3)

    !$gpu

    dx = dx_level(:,amr_level)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = max(sqrt(loc(1)**2 + loc(2)**2), R_min)

             update(i,j,k,UMR) = update(i,j,k,UMR) - ( (loc(1) / R) * (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx(1) + &
                  (loc(2) / R) * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / dx(2) )

          enddo
       enddo
    enddo

  end subroutine add_hybrid_advection_source

end module hybrid_advection_module
