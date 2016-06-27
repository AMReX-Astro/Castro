module sponge_module

  implicit none

  double precision, save :: sponge_lower_radius, sponge_upper_radius
  double precision, save :: sponge_lower_density, sponge_upper_density
  double precision, save :: sponge_timescale  

contains

  subroutine ca_sponge(lo,hi,state,state_lo,state_hi,dx,dt,time,E_added,mom_added) &
       bind(C, name="ca_sponge")

    use prob_params_module,   only: problo, center
    use meth_params_module,   only: URHO, UMX, UMZ, UEDEN, NVAR
    use bl_constants_module,  only: ZERO, HALF, ONE, M_PI
    
    implicit none
    
    integer          :: lo(3),hi(3)
    integer          :: state_lo(3), state_hi(3)
    double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
    double precision :: dx(3), dt, time
    double precision :: E_added, mom_added(3)

    ! Local variables

    double precision :: r(3), radius
    double precision :: ke_old, E_old, mom_old(3)
    double precision :: sponge_factor, alpha
    double precision :: delta_r, delta_rho
    double precision :: rho, rhoInv
    
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

             state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) / (ONE + alpha * sponge_factor)

             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + HALF * sum(state(i,j,k,UMX:UMZ)**2) * rhoInv - ke_old

             ! Ending diagnostic quantities

             E_added   = E_added   + state(i,j,k,UEDEN)   - E_old
             mom_added = mom_added + state(i,j,k,UMX:UMZ) - mom_old

          enddo
       enddo
    enddo

  end subroutine ca_sponge

end module sponge_module
