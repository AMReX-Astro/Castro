   subroutine amrex_probinit(init,name,namlen,problo,probhi) bind(c)

     use problem_io_module, only: initialize_io
     use wdmerger_util_module, only: initialize_problem
     use amrex_fort_module, only: rt => amrex_real

     implicit none

     integer  :: init, namlen
     integer  :: name(namlen)
     real(rt) :: problo(3), probhi(3)

     call initialize_io(name, namlen)
     call initialize_problem(init)

   end subroutine amrex_probinit



   subroutine ca_initdata(lo, hi, &
                          state, state_lo, state_hi, &
                          dx, problo) bind(C, name='ca_initdata')

     use amrex_fort_module, only: rt => amrex_real
     use probdata_module
     use wdmerger_util_module
     use prob_params_module, only: center, dim
     use eos_type_module, only: eos_t
     use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
                                   UEDEN, UEINT, UFS, do_rotation, state_in_rotating_frame
     use network, only: nspec
     use amrex_constants_module
     use model_parser_module, only: idens_model, itemp_model, ipres_model, ispec_model
     use initial_model_module, only: interpolate_3d_from_1d
     use math_module, only: cross_product ! function
     use castro_util_module, only: position ! function
     use rotation_frequency_module, only: get_omega ! function
     use wdmerger_util_module, only: inertial_velocity ! function
     use wdmerger_util_module, only: get_ambient ! function

     implicit none

     integer,  intent(in   ) :: lo(3), hi(3)
     integer,  intent(in   ) :: state_lo(3), state_hi(3)
     real(rt), intent(in   ) :: dx(3), problo(3)
     real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)

     real(rt) :: loc(3), pos(3), omega(3), vel(3), mom(3)
     real(rt) :: rot_loc(3), cap_radius
     real(rt) :: dist_P, dist_S
     real(rt) :: cosTheta, sinTheta, R_prp, mag_vel

     type (eos_t) :: zone_state, ambient_state

     integer :: i, j, k

     real(rt) :: time = ZERO

     !$gpu

     ! Loop through the zones and set the zone state depending on whether we are
     ! inside the primary or secondary (in which case interpolate from the respective model)
     ! or if we are in an ambient zone.

     call get_ambient(ambient_state)

     omega = get_omega(time)

     !$OMP PARALLEL DO PRIVATE(i, j, k, loc, vel, pos, mom) &
     !$OMP PRIVATE(dist_P, dist_S, zone_state) &
     !$OMP PRIVATE(rot_loc, cap_radius)
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              loc = position(i,j,k)

              dist_P = sum((loc - center_P_initial)**2)**HALF
              dist_S = sum((loc - center_S_initial)**2)**HALF

              ! If the zone size is smaller than the stellar radius,
              ! then use interpolation from the 1D model. If the zone
              ! size is larger than the stellar radius, which can happen
              ! on the coarsest levels on large grids, we need to ensure
              ! that some mass gets loaded onto the grid so that we can
              ! refine in the appropriate place. This does not need to be
              ! accurate, as the average down from the fine levels will get
              ! things right on the coarse levels. So we can still use the
              ! interpolation scheme, because it handles this special case
              ! for us by simply using the central zone of the model; we
              ! just need to make sure we catch it.

              if (mass_P > ZERO .and. (dist_P < model_P % radius .or. (model_P % radius <= maxval(dx) .and. dist_P < maxval(dx)))) then
                 pos = loc - center_P_initial
                 call interpolate_3d_from_1d(rho_P, T_P, xn_P, r_P, model_P % npts, pos, model_P % radius, dx, zone_state, nsub)
              else if (mass_S > ZERO .and. (dist_S < model_S % radius .or. (model_S % radius <= maxval(dx) .and. dist_S < maxval(dx)))) then
                 pos = loc - center_S_initial
                 call interpolate_3d_from_1d(rho_S, T_S, xn_S, r_S, model_S % npts, pos, model_S % radius, dx, zone_state, nsub)
              else
                 zone_state = ambient_state
              endif

              state(i,j,k,URHO)  = zone_state % rho
              state(i,j,k,UTEMP) = zone_state % T
              state(i,j,k,UEINT) = zone_state % e * zone_state % rho
              state(i,j,k,UEDEN) = zone_state % e * zone_state % rho
              state(i,j,k,UFS:UFS+nspec-1) = zone_state % rho * zone_state % xn

              ! Set the velocities in each direction equal to the bulk
              ! velocity of the system. By default this is zero so that
              ! the system is at rest in our reference frame.

              state(i,j,k,UMX) = state(i,j,k,URHO) * (bulk_velx + smallu)
              state(i,j,k,UMY) = state(i,j,k,URHO) * (bulk_vely + smallu)
              state(i,j,k,UMZ) = state(i,j,k,URHO) * (bulk_velz + smallu)

              loc = loc - center

              ! Add any additional velocity imparted to the stars, usually
              ! from an eccentric orbit or from a collision calculation.

              dist_P = sum((loc - center_P_initial)**2)**HALF
              dist_S = sum((loc - center_S_initial)**2)**HALF

              if (problem /= 1) then

                 if (dist_P < model_P % radius) then
                    mom = vel_P * state(i,j,k,URHO)
                    state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + mom
                 else if (dist_S < model_S % radius) then
                    mom = vel_S * state(i,j,k,URHO)
                    state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + mom
                 endif

              endif

              ! If we're in the inertial reference frame, use rigid body rotation with velocity omega x r.

              if ( ( (do_rotation .ne. 1) .or. ( (do_rotation .eq. 1) .and. (state_in_rotating_frame .ne. 1) ) ) .and. (problem == 1) ) then

                 rot_loc = loc

                 ! At large enough distances from the center, our rigid body rotation formula gives
                 ! meaningless results, and this is enough to be an issue for the problem sizes of
                 ! interest. We don't want the stars to be plowing through ambient material, though.
                 ! So we need a solution for the ambient material that satisfies both criteria. Our
                 ! solution is to set a cap on the radius used in calculating the rotation velocity:
                 ! the material around the stars will be rotating at the same speed, avoiding unwanted
                 ! numerical effects, but the material near the domain boundaries will still have a
                 ! reasonable velocity. We'll arbitrarily apply the cap at some multiple of the larger
                 ! of the two stellar radii.

                 if (state(i,j,k,URHO) < 1.1d0 * ambient_state % rho) then

                    cap_radius = 1.25d0 * max(model_P % radius + abs(center_P_initial(axis_1)), model_S % radius + abs(center_S_initial(axis_1)))
                    rot_loc(:) = min(cap_radius, abs(rot_loc(:))) * sign(ONE, rot_loc(:))

                 endif

                 vel = cross_product(omega, rot_loc)

                 mom = state(i,j,k,URHO) * vel(:)
                 state(i,j,k,UMX:UMZ) = state(i,j,k,UMX:UMZ) + mom

                 ! In 2D we have to be careful: the third coordinate is an angular
                 ! coordinate, whose unit vector is tangent to the unit circle, so we should
                 ! have the same velocity everywhere along that coordinate to begin with.

                 if (dim .eq. 2) state(i,j,k,UMZ) = abs(state(i,j,k,UMZ))

              endif

              ! If we're doing the merger problem and want to add an initial radial velocity, do that here.
              ! We're only going to impart this velocity to the stars themselves, so that we prevent
              ! artificial infall in the ambient regions.

              if (problem == 1 .and. initial_radial_velocity_factor > ZERO) then

                 if (dist_P < model_P % radius .or. dist_P < model_S % radius) then

                    R_prp    = sqrt(loc(axis_1)**2 + loc(axis_2)**2)
                    cosTheta = loc(axis_1) / R_prp
                    sinTheta = loc(axis_2) / R_prp

                    vel = state(i,j,k,UMX:UMZ)
                    vel = inertial_velocity(loc, vel, time)
                    mag_vel = sqrt( sum( vel**2 ) )

                    state(i,j,k,UMX+axis_1-1) = state(i,j,k,UMX+axis_1-1) - initial_radial_velocity_factor * cosTheta * mag_vel
                    state(i,j,k,UMX+axis_2-1) = state(i,j,k,UMX+axis_2-1) - initial_radial_velocity_factor * sinTheta * mag_vel

                 endif

              endif

              ! Add corresponding kinetic energy from the velocity on the grid.

              state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

           enddo
        enddo
     enddo
     !$OMP END PARALLEL DO

   end subroutine ca_initdata
