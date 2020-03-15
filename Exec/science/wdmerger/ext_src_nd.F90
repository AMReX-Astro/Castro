
     subroutine ca_ext_src(lo, hi, &
                           old_state, os_lo, os_hi, &
                           new_state, ns_lo, ns_hi, &
                           src, src_lo, src_hi, &
                           problo, dx, time, dt) bind(C, name='ca_ext_src')

       use amrex_fort_module, only: rt => amrex_real
       use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, NSRC
       use prob_params_module, only: center
       use amrex_constants_module, only: ZERO, HALF, ONE, TWO
       use castro_util_module, only: position ! function
       use wdmerger_util_module, only: problem, relaxation_damping_factor, radial_damping_factor, &
                                       t_ff_P, t_ff_S, axis_1, axis_2, axis_3
       use wdmerger_util_module, only: inertial_velocity ! function
#ifdef HYBRID_MOMENTUM
       use hybrid_advection_module, only: linear_to_hybrid ! function
       use meth_params_module, only: UMR, UMP
#endif
#ifdef ROTATION
       use meth_params_module, only: do_rotation, state_in_rotating_frame
       use rotation_module, only: inertial_to_rotational_velocity
#endif

       implicit none

       integer,  intent(in   ) :: lo(3), hi(3)
       integer,  intent(in   ) :: os_lo(3), os_hi(3)
       integer,  intent(in   ) :: ns_lo(3), ns_hi(3)
       integer,  intent(in   ) :: src_lo(3), src_hi(3)
       real(rt), intent(in   ) :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
       real(rt), intent(in   ) :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
       real(rt), intent(inout) :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NSRC)
       real(rt), intent(in   ) :: problo(3), dx(3)
       real(rt), intent(in   ), value :: time, dt

       ! Local variables

       real(rt) :: relaxation_damping_timescale, radial_damping_timescale
       real(rt) :: dynamical_timescale, damping_factor
       real(rt) :: loc(3), R_prp, sinTheta, cosTheta, v_rad, Sr(3)
       integer  :: i, j, k
       real(rt) :: vel(3), mom(3), rhoInv

       !$gpu

       ! First do any relaxation source terms.

       if (problem == 1 .and. relaxation_damping_factor > ZERO) then

          ! The relevant dynamical timescale for determining this source term timescale should be
          ! the smaller of the two WD timescales. Generally this should be the primary, but we'll
          ! be careful just in case.

          dynamical_timescale = min(t_ff_P, t_ff_S)

          ! The relaxation damping factor should be less than unity, so that the damping
          ! timescale is less than the dynamical timescale. This ensures that the stars
          ! are always responding to the damping with quasistatic motion; if the stars
          ! could respond too quickly, they might expand and make contact too early.

          relaxation_damping_timescale = relaxation_damping_factor * dynamical_timescale

          ! Note that we are applying this update implicitly. The implicit and
          ! explicit methods agree in the limit where the damping timescale is
          ! much larger than dt, but the implicit method helps avoid numerical
          ! problems when the damping timescale is shorter than the timestep.
          ! For further information, see Source/sources/sponge_nd.F90.

          damping_factor = -(ONE - ONE / (ONE + dt / relaxation_damping_timescale)) / dt

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   rhoInv = ONE / new_state(i,j,k,URHO)

                   loc = position(i,j,k) - center

                   mom = new_state(i,j,k,UMX:UMZ)

#ifdef ROTATION
                   if (do_rotation == 1 .and. state_in_rotating_frame == 0) then
                      vel = rhoInv * mom
                      call inertial_to_rotational_velocity([i, j, k], time, vel)
                      mom = new_state(i,j,k,URHO) * vel
                   end if
#endif

                   Sr = mom * damping_factor

                   src(i,j,k,UMX:UMZ) = src(i,j,k,UMX:UMZ) + Sr

#ifdef HYBRID_MOMENTUM
                   src(i,j,k,UMR:UMP) = src(i,j,k,UMR:UMP) + linear_to_hybrid(loc, Sr)
#endif

                   ! Do the same thing for the kinetic energy update.

                   src(i,j,k,UEDEN) = src(i,j,k,UEDEN) + dot_product(rhoInv * mom, Sr)

                enddo
             enddo
          enddo

       endif



       ! Now do the radial drift source terms.

       if (problem == 1 .and. radial_damping_factor > ZERO) then

          ! The logic for which dynamical timescale to use and
          ! how to do the implicit coupling follows the reasoning
          ! described above for the relaxation damping.

          dynamical_timescale = max(t_ff_P, t_ff_S)

          radial_damping_timescale = radial_damping_factor * dynamical_timescale

          damping_factor = -(ONE - ONE / (ONE + dt / radial_damping_timescale)) / dt

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   rhoInv = ONE / new_state(i,j,k,URHO)

                   loc      = position(i,j,k) - center
                   R_prp    = sqrt(loc(axis_1)**2 + loc(axis_2)**2)
                   cosTheta = loc(axis_1) / R_prp
                   sinTheta = loc(axis_2) / R_prp

                   mom = new_state(i,j,k,UMX:UMZ)
                   vel = rhoInv * mom
                   vel = inertial_velocity(loc, vel, time)

                   v_rad = cosTheta * vel(axis_1) + sinTheta * vel(axis_2)

                   ! What we want to do is insert a negative radial drift acceleration. If continued
                   ! for long enough, it will eventually drive coalescence of the binary. The
                   ! restriction on how large the acceleration can be is guided by the dynamical
                   ! properties of the system: it needs to be small enough that the WDs can be
                   ! in approximate dynamical equilibrium at all times before significant mass
                   ! transfer begins. So, if we write the force as
                   ! d(v_rad) / dt = -|v_phi| / tau,
                   ! where tau = radial_damping_factor * dynamical_timescale is the timescale
                   ! and |v_phi| is the magnitude of the azimuthal velocity, then
                   ! radial_damping_factor should be much greater than unity.

                   Sr(axis_1) = cosTheta * (new_state(i,j,k,URHO) * abs(v_rad)) * damping_factor
                   Sr(axis_2) = sinTheta * (new_state(i,j,k,URHO) * abs(v_rad)) * damping_factor
                   Sr(axis_3) = ZERO

                   src(i,j,k,UMX:UMZ) = src(i,j,k,UMX:UMZ) + Sr

#ifdef HYBRID_MOMENTUM
                   src(i,j,k,UMR:UMP) = src(i,j,k,UMR:UMP) + linear_to_hybrid(loc, Sr)
#endif

                   ! The kinetic energy source term is v . Sr:

                   src(i,j,k,UEDEN) = src(i,j,k,UEDEN) + dot_product(rhoInv * mom, Sr)

                enddo
             enddo
          enddo

       endif

     end subroutine ca_ext_src
