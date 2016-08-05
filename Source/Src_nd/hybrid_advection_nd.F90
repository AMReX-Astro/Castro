module hybrid_advection_module

  implicit none

contains

  ! Takes the initial linear momentum data in a state and converts it
  ! to the hybrid momenta.

  subroutine init_hybrid_momentum(lo, hi, state, s_lo, s_hi) bind(C,name='init_hybrid_momentum')

    use meth_params_module, only: NVAR, UMR, UMP, UMX, UMZ
    use castro_util_module, only: position
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             state(i,j,k,UMR:UMP) = linear_to_hybrid(loc, state(i,j,k,UMX:UMZ))

          enddo
       enddo
    enddo

  end subroutine init_hybrid_momentum



  ! Fill a sources array with the source terms in the hybrid momentum equations.

  subroutine ca_hybrid_hydro_source(lo, hi, state, s_lo, s_hi, ext_src, e_lo, e_hi) bind(C,name='ca_hybrid_hydro_source')

    use meth_params_module, only: NVAR, URHO, UMR, UML
    use prob_params_module, only: center
    use castro_util_module, only: position
    use network, only: nspec, naux
    use eos_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: e_lo(3), e_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: ext_src(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3), R, rhoInv

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             rhoInv = ONE / state(i,j,k,URHO)

             ext_src(i,j,k,UMR) = ext_src(i,j,k,UMR) + (rhoInv / R**3) * state(i,j,k,UML)**2

          enddo
       enddo
    enddo

  end subroutine ca_hybrid_hydro_source



  ! Convert a linear quantity such as momentum into the
  ! "hybrid" scheme that has radial and angular components.

  function linear_to_hybrid(loc, vec_in) result(vec_out)

    implicit none

    double precision, intent(in ) :: loc(3), vec_in(3)
    double precision :: vec_out(3)

    double precision :: R

    ! this conversion is Eqs. 25 and 26 in Byerly et al. 2014.  We
    ! we do not use the corotating frame here
    R = sqrt( loc(1)**2 + loc(2)**2 )

    vec_out(1) = vec_in(1) * (loc(1) / R) + vec_in(2) * (loc(2) / R)
    vec_out(2) = vec_in(2) * loc(1)       - vec_in(1) * loc(2)
    vec_out(3) = vec_in(3)

  end function linear_to_hybrid



  ! Convert a "hybrid" quantity into a linear one.

  function hybrid_to_linear(loc, vec_in) result(vec_out)

    implicit none

    double precision, intent(in ) :: loc(3), vec_in(3)
    double precision :: vec_out(3)

    double precision :: R

    R = sqrt( loc(1)**2 + loc(2)**2 )

    vec_out(1) = vec_in(1) * (loc(1) / R)    - vec_in(2) * (loc(2) / R**2)
    vec_out(2) = vec_in(2) * (loc(1) / R**2) + vec_in(1) * (loc(2) / R)
    vec_out(3) = vec_in(3)

  end function hybrid_to_linear



  ! Update momentum to account for source term

  subroutine add_hybrid_momentum_source(loc, mom, source)

    implicit none

    double precision :: loc(3), mom(3), source(3)

    mom = mom + linear_to_hybrid(loc, source)

  end subroutine add_hybrid_momentum_source



  subroutine compute_hybrid_flux(state, flux, idir, idx)

    use meth_params_module, only: NVAR, NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, UMR, UML, UMP
    use bl_error_module, only: bl_error
    use prob_params_module, only: center
    use castro_util_module, only: position
#ifdef ROTATION
    use amrinfo_module, only: amr_time
    use meth_params_module, only: do_rotation, state_in_rotating_frame
    use rotation_module, only: inertial_to_rotational_velocity
#endif

    implicit none

    double precision :: state(NGDNV)
    double precision :: flux(NVAR)
    integer          :: idir, idx(3)

    double precision :: linear_mom(3), hybrid_mom(3)
    double precision :: loc(3), R

    double precision :: u_adv

#ifdef ROTATION
    double precision :: vel(3)
#endif

    if (idir .eq. 1) then
       loc = position(idx(1),idx(2),idx(3),ccx=.false.) - center
       u_adv = state(GDU)
    else if (idir .eq. 2) then
       loc = position(idx(1),idx(2),idx(3),ccy=.false.) - center
       u_adv = state(GDV)
    else if (idir .eq. 3) then
       loc = position(idx(1),idx(2),idx(3),ccz=.false.) - center
       u_adv = state(GDW)
    else
       call bl_error("Error: unknown direction in compute_hybrid_flux.")
    endif

    ! If we are in the rotating frame but evolving inertial frame variables,
    ! subtract off the coordinate frame velocity here.

    
#ifdef ROTATION
    if (do_rotation .eq. 1 .and. state_in_rotating_frame .ne. 1) then
       vel = state(GDU:GDW)
       call inertial_to_rotational_velocity(idx, amr_time, vel, idir)
       u_adv = vel(idir)
    endif
#endif

    R = sqrt(loc(1)**2 + loc(2)**2)

    linear_mom = state(GDRHO) * state(GDU:GDW)

    hybrid_mom = linear_to_hybrid(loc, linear_mom)

    if (idir .eq. 1) then

       flux(UMR) = hybrid_mom(1) * u_adv
       flux(UML) = hybrid_mom(2) * u_adv - loc(2) * state(GDPRES)
       flux(UMP) = hybrid_mom(3) * u_adv

    else if (idir .eq. 2) then

       flux(UMR) = hybrid_mom(1) * u_adv
       flux(UML) = hybrid_mom(2) * u_adv + loc(1) * state(GDPRES)
       flux(UMP) = hybrid_mom(3) * u_adv

    else if (idir .eq. 3) then

       flux(UMR) = hybrid_mom(1) * u_adv
       flux(UML) = hybrid_mom(2) * u_adv
       flux(UMP) = hybrid_mom(3) * u_adv + state(GDPRES)

    else

       call bl_error("Error: unknown direction in compute_hybrid_flux.")

    endif

  end subroutine compute_hybrid_flux



  subroutine add_hybrid_advection_source(lo, hi, dt, &
                                         update, u_lo, u_hi, &
                                         qx, qx_lo, qx_hi, &
                                         qy, qy_lo, qy_hi, &
                                         qz, qz_lo, qz_hi)

    use meth_params_module, only: NVAR, NGDNV, GDPRES, UMR
    use prob_params_module, only: center, dx_level
    use castro_util_module, only: position
    use amrinfo_module, only: amr_level

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: u_lo(3), u_hi(3)
    integer          :: qx_lo(3), qx_hi(3)
    integer          :: qy_lo(3), qy_hi(3)
    integer          :: qz_lo(3), qz_hi(3)
    double precision :: update(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
    double precision :: qx(qx_lo(1):qx_hi(1),qx_lo(2):qx_hi(2),qx_lo(3):qx_hi(3),NGDNV)
    double precision :: qy(qy_lo(1):qy_hi(1),qy_lo(2):qy_hi(2),qy_lo(3):qy_hi(3),NGDNV)
    double precision :: qz(qz_lo(1):qz_hi(1),qz_lo(2):qz_hi(2),qz_lo(3):qz_hi(3),NGDNV)
    double precision :: dt

    integer          :: i, j, k
    double precision :: loc(3), R, dx(3)

    dx = dx_level(:,amr_level)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             update(i,j,k,UMR) = update(i,j,k,UMR) - ( (loc(1) / R) * (qx(i+1,j,k,GDPRES) - qx(i,j,k,GDPRES)) / dx(1) + &
                                                       (loc(2) / R) * (qy(i,j+1,k,GDPRES) - qy(i,j,k,GDPRES)) / dx(2) )

          enddo
       enddo
    enddo

  end subroutine add_hybrid_advection_source



  ! Update state to account for hybrid advection.

  subroutine hybrid_update(lo, hi, state, state_lo, state_hi) bind(C,name='hybrid_update')

    use bl_constants_module, only: HALF, ONE
    use meth_params_module, only: URHO, UMR, UMP, UMX, UMZ, UEDEN, NVAR
    use castro_util_module, only: position
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: state_lo(3), state_hi(3)
    double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3)
    double precision :: old_ke, new_ke, rhoInv

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             rhoInv = ONE / state(i,j,k,URHO)

             old_ke = HALF * rhoInv * sum(state(i,j,k,UMX:UMZ)**2)

             state(i,j,k,UMX:UMZ) = hybrid_to_linear(loc, state(i,j,k,UMR:UMP))

             new_ke = HALF * rhoInv * sum(state(i,j,k,UMX:UMZ)**2)

             ! We have effectively updated the kinetic energy of the state; to account
             ! for this, add the difference between the old and new KE to the total E.

             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + (new_ke - old_ke)

          enddo
       enddo
    enddo

  end subroutine hybrid_update

end module hybrid_advection_module
