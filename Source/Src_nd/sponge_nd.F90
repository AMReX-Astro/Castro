module sponge_module

  implicit none

  double precision, save :: sponge_lower_radius, sponge_upper_radius
  double precision, save :: sponge_lower_density, sponge_upper_density
  double precision, save :: sponge_timescale  

contains

  subroutine ca_sponge(lo,hi,state,state_lo,state_hi,vol,vol_lo,vol_hi,dx,dt,time,E_added,mom_added) &
                       bind(C, name="ca_sponge")

    use prob_params_module,   only: problo, center
    use meth_params_module,   only: URHO, UMX, UMZ, UEDEN, NVAR, sponge_implicit
    use bl_constants_module,  only: ZERO, HALF, ONE, M_PI
#ifdef HYBRID_MOMENTUM
    use meth_params_module,   only: UMR, UMP
    use hybrid_advection_module, only: linear_to_hybrid
#endif
    
    implicit none
    
    integer          :: lo(3),hi(3)
    integer          :: state_lo(3), state_hi(3)
    integer          :: vol_lo(3), vol_hi(3)
    double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    double precision :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
    double precision :: dx(3), dt, time
    double precision :: E_added, mom_added(3)

    ! Local variables

    double precision :: r(3), radius
    double precision :: ke_old, E_old, mom_old(3)
    double precision :: sponge_factor, alpha
    double precision :: delta_r, delta_rho
    double precision :: rho, rhoInv
    double precision :: update_factor, Sr(3)
    
    integer          :: i, j, k

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

    do k = lo(3), hi(3)
       r(3) = problo(3) + dble(k + HALF) * dx(3) - center(3)

       do j = lo(2), hi(2)
          r(2) = problo(2) + dble(j + HALF) * dx(2) - center(2)

          do i = lo(1), hi(1)
             r(1) = problo(1) + dble(i + HALF) * dx(1) - center(1)

             rho = state(i,j,k,URHO)
             rhoInv = ONE / rho

             ! Starting diagnostic quantities

             E_old   = state(i,j,k,UEDEN)
             ke_old  = HALF * sum(state(i,j,k,UMX:UMZ)**2) * rhoInv
             mom_old = state(i,j,k,UMX:UMZ)

             ! Apply radial sponge. By default sponge_lower_radius will be zero
             ! so this sponge is applied only if set by the user.

             if (sponge_lower_radius > ZERO .and. sponge_upper_radius > ZERO) then

                radius = sqrt(sum(r**2))

                if (radius < sponge_lower_radius) then
                   sponge_factor = ZERO
                else if (radius >= sponge_lower_radius .and. radius <= sponge_upper_radius) then
                   sponge_factor = HALF * (ONE - cos(M_PI * (radius - sponge_lower_radius) / delta_r))
                else
                   sponge_factor = ONE
                endif

             endif

             ! Apply density sponge. By default sponge_upper_density will be zero
             ! so this sponge is applied only if set by the user.

             ! Note that because we do this second, the density sponge gets priority
             ! over the radial sponge in cases where the two would overlap.

             if (sponge_upper_density > ZERO .and. sponge_lower_density > ZERO) then

                if (rho > sponge_upper_density) then
                   sponge_factor = ZERO
                else if (rho <= sponge_upper_density .and. rho >= sponge_lower_density) then
                   sponge_factor = HALF * (ONE - cos(M_PI * (rho - sponge_upper_density) / delta_rho))
                else
                   sponge_factor = ONE
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

             Sr = state(i,j,k,UMX:UMZ) * update_factor

             state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + Sr

#ifdef HYBRID_MOMENTUM
             state(i,j,k,UMR:UMP) = state(i,j,k,UMR:UMP) + linear_to_hybrid(r, Sr)
#endif

             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + HALF * sum(state(i,j,k,UMX:UMZ)**2) * rhoInv - ke_old

             ! Ending diagnostic quantities

             E_added   = E_added   + (state(i,j,k,UEDEN)   - E_old) * vol(i,j,k)
             mom_added = mom_added + (state(i,j,k,UMX:UMZ) - mom_old) * vol(i,j,k)

          enddo
       enddo
    enddo

  end subroutine ca_sponge

end module sponge_module
