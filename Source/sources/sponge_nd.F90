module sponge_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt), allocatable, save :: sponge_lower_factor, sponge_upper_factor
  real(rt), allocatable, save :: sponge_lower_radius, sponge_upper_radius
  real(rt), allocatable, save :: sponge_lower_density, sponge_upper_density
  real(rt), allocatable, save :: sponge_lower_pressure, sponge_upper_pressure
  real(rt), allocatable, save :: sponge_target_velocity(:)
  real(rt), allocatable, save :: sponge_timescale

#ifdef AMREX_USE_CUDA
  attributes(managed) :: sponge_lower_factor, sponge_upper_factor
  attributes(managed) :: sponge_lower_radius, sponge_upper_radius
  attributes(managed) :: sponge_lower_density, sponge_upper_density
  attributes(managed) :: sponge_lower_pressure, sponge_upper_pressure
  attributes(managed) :: sponge_target_velocity
  attributes(managed) :: sponge_timescale
#endif

  public

  private :: update_factor

contains

  subroutine ca_sponge(lo,hi,state,state_lo,state_hi,source,src_lo,src_hi, &
       vol,vol_lo,vol_hi,dx,dt,time,mult_factor) &
       bind(C, name="ca_sponge")

    use prob_params_module,   only: problo, center
    use meth_params_module,   only: URHO, UMX, UMZ, UEDEN, NVAR
    use amrex_constants_module,  only: ZERO, HALF, ONE
#ifdef HYBRID_MOMENTUM
    use meth_params_module,   only: UMR, UMP
    use hybrid_advection_module, only: set_hybrid_momentum_source
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
    real(rt)         :: dx(3)
    real(rt), value  :: dt, time, mult_factor

    ! Local variables

    real(rt) :: r(3)
    real(rt) :: Sr(3), SrE
    real(rt) :: rho, rhoInv

    integer  :: i, j, k

    real(rt) :: src(NVAR)
    real(rt) :: local_state(NVAR)

    !$gpu

    src(:) = ZERO

    do k = lo(3), hi(3)
       r(3) = problo(3) + dble(k + HALF) * dx(3) - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dble(j + HALF) * dx(2) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dble(i + HALF) * dx(1) - center(1)

             rho = state(i,j,k,URHO)
             rhoInv = ONE / rho

             local_state = state(i,j,k,:)

             Sr(:) = (state(i,j,k,UMX:UMZ) - rho * sponge_target_velocity) * update_factor(r, local_state, dt) * mult_factor / dt

             src(UMX:UMZ) = Sr(:)

             ! Kinetic energy is 1/2 rho u**2, or (rho u)**2 / (2 rho). This means
             ! that d(KE)/dt = u d(rho u)/dt - 1/2 u**2 d(rho)/dt. In this case
             ! the sponge has no contribution to rho, so the kinetic energy source
             ! term, and thus the total energy source term, is u * momentum source.

             SrE = dot_product(state(i,j,k,UMX:UMZ) * rhoInv, Sr)

             src(UEDEN) = SrE

#ifdef HYBRID_MOMENTUM
             call set_hybrid_momentum_source(r, src(UMR:UMP), src(UMX:UMZ))
#endif

             ! Add terms to the source array.

             source(i,j,k,:) = source(i,j,k,:) + src

          enddo
       enddo
    enddo

  end subroutine ca_sponge



  function update_factor(r, state, dt) result(fac)

    use amrex_constants_module, only: ZERO, HALF, ONE, M_PI
    use meth_params_module, only: sponge_implicit, NVAR, URHO, UTEMP, UFS, UFX
    use network, only: nspec, naux
    use eos_type_module, only: eos_t, eos_input_rt
    use eos_module, only: eos

    implicit none

    real(rt), intent(in) :: r(3), state(NVAR), dt

    real(rt) :: radius, rho, rhoInv, p
    real(rt) :: delta_r, delta_rho, delta_p
    real(rt) :: alpha, sponge_factor
    type(eos_t) :: eos_state
    real(rt) :: fac

    !$gpu

    ! Radial distance between upper and lower boundaries.

    delta_r = sponge_upper_radius - sponge_lower_radius

    ! Density difference between upper and lower cutoffs.

    delta_rho = sponge_lower_density - sponge_upper_density

    ! Pressure difference between upper and lower cutoffs.

    delta_p = sponge_lower_pressure - sponge_upper_pressure

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

    sponge_factor = ZERO

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

    ! Apply density sponge. This sponge is applied only if set by the user.

    ! Note that because we do this second, the density sponge gets priority
    ! over the radial sponge in cases where the two would overlap.

    if (sponge_upper_density > ZERO .and. sponge_lower_density > ZERO) then

       rho = state(URHO)

       if (rho > sponge_upper_density) then
          sponge_factor = sponge_lower_factor
       else if (rho <= sponge_upper_density .and. rho >= sponge_lower_density) then
          sponge_factor = sponge_lower_factor + HALF * (sponge_upper_factor - sponge_lower_factor) * &
               (ONE - cos(M_PI * (rho - sponge_upper_density) / delta_rho))
       else
          sponge_factor = sponge_upper_factor
       endif

    endif

    ! Apply pressure sponge. This sponge is applied only if set by the user.

    ! Note that because we do this third, the pressure sponge gets priority
    ! over the radial and density sponges in cases where the two would overlap.

    if (sponge_upper_pressure > ZERO .and. sponge_lower_pressure >= ZERO) then

       rhoInv = ONE / state(URHO)

       eos_state % rho = state(URHO)
       eos_state % T   = state(UTEMP)
       eos_state % xn  = state(UFS:UFS+nspec-1) * rhoInv
       eos_state % aux = state(UFX:UFX+naux-1) * rhoInv

       call eos(eos_input_rt, eos_state)

       p = eos_state % p

       if (p > sponge_upper_pressure) then
          sponge_factor = sponge_lower_factor
       else if (p <= sponge_upper_pressure .and. p >= sponge_lower_pressure) then
          sponge_factor = sponge_lower_factor + HALF * (sponge_upper_factor - sponge_lower_factor) * &
               (ONE - cos(M_PI * (p - sponge_upper_pressure) / delta_p))
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
       fac = -(ONE - ONE / (ONE + alpha * sponge_factor))
    else
       fac = -alpha * sponge_factor
    endif

  end function update_factor

end module sponge_module
