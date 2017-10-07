module sponge_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), save :: sponge_lower_factor, sponge_upper_factor
  real(rt), save :: sponge_lower_radius, sponge_upper_radius
  real(rt), save :: sponge_lower_density, sponge_upper_density
  real(rt), save :: sponge_timescale

  public

  private :: update_factor

contains

  subroutine ca_sponge(lo,hi,state,state_lo,state_hi,source,src_lo,src_hi, &
                       vol,vol_lo,vol_hi,dx,dt,time,mult_factor) &
                       bind(C, name="ca_sponge")

    use prob_params_module,   only: problo, center
    use meth_params_module,   only: URHO, UMX, UMZ, UEDEN, NVAR
    use bl_constants_module,  only: ZERO, HALF, ONE
#ifdef HYBRID_MOMENTUM
    use meth_params_module,   only: UMR, UMP
    use hybrid_advection_module, only: add_hybrid_momentum_source
#endif
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer          :: lo(3),hi(3)
    integer          :: state_lo(3), state_hi(3)
    integer          :: src_lo(3), src_hi(3)
    integer          :: vol_lo(3), vol_hi(3)
    real(rt)         :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    real(rt)         :: source(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
    real(rt)         :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    real(rt)         :: dx(3), dt, time
    real(rt), value  :: mult_factor

    ! Local variables

    real(rt) :: r(3)
    real(rt) :: Sr(3), SrE
    real(rt) :: rho, rhoInv

    integer  :: i, j, k

    real(rt) :: src(NVAR)

    src = ZERO

    do k = lo(3), hi(3)
       r(3) = problo(3) + dble(k + HALF) * dx(3) - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dble(j + HALF) * dx(2) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dble(i + HALF) * dx(1) - center(1)

             rho = state(i,j,k,URHO)
             rhoInv = ONE / rho

             Sr(:) = state(i,j,k,UMX:UMZ) * update_factor(r, rho, dt) * mult_factor / dt

             src(UMX:UMZ) = Sr(:)

             SrE = dot_product(state(i,j,k,UMX:UMZ) * rhoInv, Sr)

             src(UEDEN) = SrE

#ifdef HYBRID_MOMENTUM
             call add_hybrid_momentum_source(r, src(UMR:UMP), src(UMX:UMZ))
#endif

             ! Add terms to the source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_sponge



  real(rt) function update_factor(r, rho, dt)

    use bl_constants_module, only: ZERO, HALF, ONE, M_PI
    use meth_params_module, only: sponge_implicit

    implicit none

    real(rt), intent(in) :: r(3), rho, dt

    real(rt) :: radius, delta_r, delta_rho, alpha, sponge_factor

    ! Radial distance between upper and lower boundaries.

    delta_r = sponge_upper_radius - sponge_lower_radius

    ! Density difference between upper and lower cutoffs.

    delta_rho = sponge_lower_density - sponge_upper_density

    ! alpha is a dimensionless measure of the timestep size; if
    ! sponge_timescale < dt, then the sponge will have a larger effect,
    ! and if sponge_timescale > dt, then the sponge will have a diminished effect.

    if (sponge_timescale > ZERO) then
       alpha = dt / sponge_timescale
    else
       alpha = ZERO
    endif

    ! Apply radial sponge. By default sponge_lower_radius will be zero
    ! so this sponge is applied only if set by the user.

    if (sponge_lower_radius >= ZERO .and. sponge_upper_radius > sponge_lower_radius) then

       radius = sqrt(sum(r**2))

       if (radius < sponge_lower_radius) then
          sponge_factor = sponge_lower_factor
       else if (radius >= sponge_lower_radius .and. radius <= sponge_upper_radius) then
          sponge_factor = sponge_lower_factor + HALF * (sponge_upper_factor - sponge_lower_factor) * &
                                                (ONE - cos(M_PI * (radius - sponge_lower_radius) / delta_r))
       else
          sponge_factor = sponge_upper_factor
       endif

    endif

    ! Apply density sponge. By default sponge_upper_density will be zero
    ! so this sponge is applied only if set by the user.

    ! Note that because we do this second, the density sponge gets priority
    ! over the radial sponge in cases where the two would overlap.

    if (sponge_upper_density > ZERO .and. sponge_lower_density > ZERO) then

       if (rho > sponge_upper_density) then
          sponge_factor = sponge_lower_factor
       else if (rho <= sponge_upper_density .and. rho >= sponge_lower_density) then
          sponge_factor = sponge_lower_factor + HALF * (sponge_upper_factor - sponge_lower_factor) * &
                                                (ONE - cos(M_PI * (rho - sponge_upper_density) / delta_rho))
       else
          sponge_factor = sponge_upper_factor
       endif

    endif

    ! For an explicit update (sponge_implicit /= 1), the source term is given by
    ! -(rho v) * alpha * sponge_factor. We simply add this directly by using the
    ! current value of the momentum.

    ! For an implicit update (sponge_implicit == 1), we choose the (rho v) to be
    ! the momentum after the update. This then leads to an update of the form
    ! (rho v) --> (rho v) * ONE / (ONE + alpha * sponge_factor). To get an equivalent
    ! explicit form of this source term, which we need for the hybrid momentum update,
    ! we can then solve (rho v) + Sr == (rho v) / (ONE + alpha * sponge_factor),
    ! which yields Sr = - (rho v) * (ONE - ONE / (ONE + alpha * sponge_factor)).

    if (sponge_implicit == 1) then
       update_factor = -(ONE - ONE / (ONE + alpha * sponge_factor))
    else
       update_factor = -alpha * sponge_factor
    endif

  end function update_factor

end module sponge_module
