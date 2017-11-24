
     subroutine ca_ext_src(lo,hi,&
                           old_state,os_lo,os_hi,&
                           new_state,ns_lo,ns_hi,&
                           src,src_lo,src_hi,problo,dx,time,dt)

       use meth_params_module,  only: NVAR, URHO, UMX, UMZ, UEDEN
       use prob_params_module,  only: center
       use bl_constants_module, only: ZERO, HALF, ONE, TWO
       use probdata_module,     only: problem, relaxation_damping_timescale, radial_damping_factor, &
                                      t_ff_P, t_ff_S, axis_1, axis_2, axis_3
       use castro_util_module,  only: position
       use wdmerger_util_module, only: inertial_velocity
#ifdef HYBRID_MOMENTUM
       use hybrid_advection_module, only: linear_to_hybrid
       use meth_params_module, only: UMR, UMP
#endif

       implicit none

       integer          :: lo(3),hi(3)
       integer          :: os_lo(3),os_hi(3)
       integer          :: ns_lo(3),ns_hi(3)
       integer          :: src_lo(3),src_hi(3)
       double precision :: old_state(os_lo(1):os_hi(1),os_lo(2):os_hi(2),os_lo(3):os_hi(3),NVAR)
       double precision :: new_state(ns_lo(1):ns_hi(1),ns_lo(2):ns_hi(2),ns_lo(3):ns_hi(3),NVAR)
       double precision :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),NVAR)
       double precision :: problo(3),dx(3),time,dt

       ! Local variables

       double precision :: radial_damping_timescale
       double precision :: dynamical_timescale, damping_factor
       double precision :: loc(3), R_prp, sinTheta, cosTheta, v_rad, Sr(3)
       integer          :: i, j, k
       double precision :: new_mom(3), old_mom(3), rhoInv

       ! Note that this function exists in a tiling region so we should only 
       ! modify the zones between lo and hi. 

       src(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

       ! The relevant dynamical timescale for determining our source term timescales should be
       ! the larger of the two WD timescales. Generally this should be the secondary, but we'll
       ! be careful just in case.

       dynamical_timescale = max(t_ff_P, t_ff_S)

       ! First do any relaxation source terms.

       if (problem == 3 .and. relaxation_damping_timescale > ZERO) then

          ! Note that we are applying this update implicitly. This helps
          ! avoid numerical problems if the relaxation_damping_timescale
          ! is shorter than the timestep. For further information, see
          ! see Source/sources/sponge_nd.F90.

          damping_factor = -(ONE - ONE / (ONE + dt / relaxation_damping_timescale)) / dt

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   rhoInv = ONE / new_state(i,j,k,URHO)

                   loc = position(i,j,k) - center

                   new_mom = new_state(i,j,k,UMX:UMZ)

                   Sr = new_mom * damping_factor

                   src(i,j,k,UMX:UMZ) = src(i,j,k,UMX:UMZ) + Sr

#ifdef HYBRID_MOMENTUM
                   src(i,j,k,UMR:UMP) = src(i,j,k,UMR:UMP) + linear_to_hybrid(loc, Sr)
#endif

                   ! Do the same thing for the kinetic energy update.

                   src(i,j,k,UEDEN) = src(i,j,k,UEDEN) + dot_product(rhoInv * new_mom, Sr)

                enddo
             enddo
          enddo

       endif



       ! Now do the radial drift source terms.

       if (problem == 3 .and. radial_damping_factor > ZERO) then

          radial_damping_timescale = radial_damping_factor * dynamical_timescale

          damping_factor = ONE / radial_damping_timescale

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   rhoInv = ONE / new_state(i,j,k,URHO)

                   loc      = position(i,j,k) - center
                   R_prp    = sqrt(loc(axis_1)**2 + loc(axis_2)**2)
                   cosTheta = loc(axis_1) / R_prp
                   sinTheta = loc(axis_2) / R_prp

                   old_mom = inertial_velocity(loc, new_state(i,j,k,UMX:UMZ), time)
                   v_rad   = cosTheta * old_mom(UMX + axis_1 - 1) + sinTheta * old_mom(UMX + axis_2 - 1)

                   ! What we want to do is insert a negative radial drift acceleration. If continued
                   ! for long enough, it will eventually drive coalescence of the binary. The
                   ! restriction on how large the acceleration can be is guided by the dynamical
                   ! properties of the system: it needs to be small enough that the WDs can be
                   ! in approximate dynamical equilibrium at all times before significant mass
                   ! transfer begins. So, if we write the force as
                   ! d(v_rad) / dt = -|v_phi| / tau,
                   ! where tau is the timescale and |v_phi| is the magnitude of the azimuthal velocity,
                   ! tau = radial_damping_factor * dynamical_timescale,
                   ! where radial_damping_factor should be much greater than unity.

                   Sr(axis_1) = -cosTheta * abs(v_rad) * damping_factor
                   Sr(axis_2) = -sinTheta * abs(v_rad) * damping_factor
                   Sr(axis_3) = ZERO

                   src(i,j,k,UMX:UMZ) = src(i,j,k,UMX:UMZ) + Sr

#ifdef HYBRID_MOMENTUM
                   src(i,j,k,UMR:UMP) = src(i,j,k,UMR:UMP) + linear_to_hybrid(loc, Sr)
#endif

                   ! The kinetic energy source term is v . Sr:

                   src(i,j,k,UEDEN) = src(i,j,k,UEDEN) + dot_product(rhoInv * old_mom, Sr)

                enddo
             enddo
          enddo

       endif

     end subroutine ca_ext_src
