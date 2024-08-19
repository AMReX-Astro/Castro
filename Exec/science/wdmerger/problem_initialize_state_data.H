#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <model_parser.H>
#include <ambient.H>
#ifdef ROTATION
#include <Rotation.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{
    // Set the zone state depending on whether we are inside the primary
    // or secondary (in which case interpolate from the respective model)
    // or if we are in an ambient zone.

    auto omega = get_omega();

    GpuArray<Real, 3> loc;
    position(i, j, k, geomdata, loc);

    Real dist_P = std::sqrt(std::pow(loc[0] - problem::center_P_initial[0], 2) +
                            std::pow(loc[1] - problem::center_P_initial[1], 2) +
                            std::pow(loc[2] - problem::center_P_initial[2], 2));

    Real dist_S = std::sqrt(std::pow(loc[0] - problem::center_S_initial[0], 2) +
                            std::pow(loc[1] - problem::center_S_initial[1], 2) +
                            std::pow(loc[2] - problem::center_S_initial[2], 2));

    // If the zone size is smaller than the stellar radius,
    // then use interpolation from the 1D model. If the zone
    // size is larger than the stellar radius, which can happen
    // on the coarsest levels on large grids, we need to ensure
    // that some mass gets loaded onto the grid so that we can
    // refine in the appropriate place. This does not need to be
    // accurate, as the average down from the fine levels will get
    // things right on the coarse levels. So we can still use the
    // interpolation scheme, because it handles this special case
    // for us by simply using the central zone of the model; we
    // just need to make sure we catch it.  Finally, sometimes
    // the stars are so close that the interpolation will overlap.
    // in this case, look at which model has the highest density
    // at that location and use that model.

    const Real* dx = geomdata.CellSize();

    Real max_dx = dx[0];
    if (AMREX_SPACEDIM >= 2) {
        max_dx = amrex::max(max_dx, dx[1]);
    }
    if (AMREX_SPACEDIM == 3) {
        max_dx = amrex::max(max_dx, dx[2]);
    }

    eos_t zone_state;

    bool P_star_test = problem::mass_P > 0.0_rt &&
        (dist_P < problem::radius_P ||
         (problem::radius_P <= max_dx && dist_P < max_dx));

    bool S_star_test = problem::mass_S > 0.0_rt &&
        (dist_S < problem::radius_S ||
         (problem::radius_S <= max_dx && dist_S < max_dx));

    double rho_P{0.0};
    double rho_S{0.0};

    if (P_star_test || S_star_test) {

        Real pos_P[3] = {loc[0] - problem::center_P_initial[0],
                         loc[1] - problem::center_P_initial[1],
                         loc[2] - problem::center_P_initial[2]};

        if (problem::mass_P > 0.0_rt) {
            rho_P = interpolate_3d(pos_P, dx, model::idens, problem::nsub, 0);
        }

        Real pos_S[3] = {loc[0] - problem::center_S_initial[0],
                         loc[1] - problem::center_S_initial[1],
                         loc[2] - problem::center_S_initial[2]};

        if (problem::mass_S > 0.0_rt) {
            rho_S = interpolate_3d(pos_S, dx, model::idens, problem::nsub, 1);
        }

        if (rho_P > rho_S) {
            // use the primary star initialization
            zone_state.rho = rho_P;
            zone_state.T   = interpolate_3d(pos_P, dx, model::itemp, problem::nsub, 0);
            for (int n = 0; n < NumSpec; ++n) {
                zone_state.xn[n] = interpolate_3d(pos_P, dx, model::ispec + n, problem::nsub, 0);
            }

        } else {
            // use the secondary star initialization
            zone_state.rho = rho_S;
            zone_state.T   = interpolate_3d(pos_S, dx, model::itemp, problem::nsub, 1);
            for (int n = 0; n < NumSpec; ++n) {
                zone_state.xn[n] = interpolate_3d(pos_S, dx, model::ispec + n, problem::nsub, 1);
            }
        }

#ifdef AUX_THERMO
        set_aux_comp_from_X(zone_state);
#endif

        eos(eos_input_rt, zone_state);

    } else {

        zone_state.rho = ambient::ambient_state[URHO];
        zone_state.T   = ambient::ambient_state[UTEMP];
        zone_state.e   = ambient::ambient_state[UEINT] / ambient::ambient_state[URHO];
        for (int n = 0; n < NumSpec; ++n) {
            zone_state.xn[n] = ambient::ambient_state[UFS+n] / ambient::ambient_state[URHO];
        }
#ifdef AUX_THERMO
        set_aux_comp_from_X(zone_state);
#endif

    }

    state(i,j,k,URHO)  = zone_state.rho;
    state(i,j,k,UTEMP) = zone_state.T;
    state(i,j,k,UEINT) = zone_state.e * zone_state.rho;
    state(i,j,k,UEDEN) = zone_state.e * zone_state.rho;
    for (int n = 0; n < NumSpec; ++n) {
        state(i,j,k,UFS+n) = zone_state.rho * zone_state.xn[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        state(i,j,k,UFX+n) = zone_state.rho * zone_state.aux[n];
    }
#endif

    // Set the velocities in each direction to zero by default.

    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    for (int n = 0; n < 3; ++n) {
        loc[n] -= problem::center[n];
    }

    // Add any additional velocity imparted to the stars, usually
    // from an eccentric orbit or from a collision calculation.

    if (problem::problem != 1) {

        if (P_star_test || S_star_test) {
            if (rho_P > rho_S && problem::mass_P > 0.0_rt && dist_P < problem::radius_P) {
                state(i,j,k,UMX) += problem::vel_P[0] * state(i,j,k,URHO);
                state(i,j,k,UMY) += problem::vel_P[1] * state(i,j,k,URHO);
                state(i,j,k,UMZ) += problem::vel_P[2] * state(i,j,k,URHO);
            } else if (problem::mass_S > 0.0_rt && dist_S < problem::radius_S) {
                state(i,j,k,UMX) += problem::vel_S[0] * state(i,j,k,URHO);
                state(i,j,k,UMY) += problem::vel_S[1] * state(i,j,k,URHO);
                state(i,j,k,UMZ) += problem::vel_S[2] * state(i,j,k,URHO);
            }
        }
    }

    // If we're in the inertial reference frame, use rigid body rotation with velocity omega x r.

    if (castro::do_rotation != 1 && problem::problem == 1) {

        GpuArray<Real, 3> rot_loc = {loc[0], loc[1], loc[2]};

        // At large enough distances from the center, our rigid body rotation formula gives
        // meaningless results, and this is enough to be an issue for the problem sizes of
        // interest. We don't want the stars to be plowing through ambient material, though.
        // So we need a solution for the ambient material that satisfies both criteria. Our
        // solution is to set a cap on the radius used in calculating the rotation velocity:
        // the material around the stars will be rotating at the same speed, avoiding unwanted
        // numerical effects, but the material near the domain boundaries will still have a
        // reasonable velocity. We'll arbitrarily apply the cap at some multiple of the larger
        // of the two stellar radii.

        if (state(i,j,k,URHO) <= castro::ambient_safety_factor * ambient::ambient_state[URHO]) {

            Real cap_radius = 1.25e0_rt * amrex::max(problem::radius_P + std::abs(problem::center_P_initial[problem::axis_1-1]),
                                                     problem::radius_S + std::abs(problem::center_S_initial[problem::axis_1-1]));

            for (int n = 0; n < 3; ++n) {
                rot_loc[n] = amrex::min(cap_radius, std::abs(rot_loc[n])) * std::copysign(1.0_rt, rot_loc[n]);
            }

        }

        GpuArray<Real, 3> vel;
        cross_product(omega, rot_loc, vel);

        state(i,j,k,UMX) += state(i,j,k,URHO) * vel[0];
        state(i,j,k,UMY) += state(i,j,k,URHO) * vel[1];
        state(i,j,k,UMZ) += state(i,j,k,URHO) * vel[2];

        // In 2D we have to be careful: the third coordinate is an angular
        // coordinate, whose unit vector is tangent to the unit circle, so we should
        // have the same velocity everywhere along that coordinate to begin with.

        if (AMREX_SPACEDIM == 2) {
            state(i,j,k,UMZ) = std::abs(state(i,j,k,UMZ));
        }

    }

    // Add corresponding kinetic energy from the velocity on the grid.

    state(i,j,k,UEDEN) += 0.5_rt * (state(i,j,k,UMX) * state(i,j,k,UMX) +
                                    state(i,j,k,UMY) * state(i,j,k,UMY) +
                                    state(i,j,k,UMZ) * state(i,j,k,UMZ)) / state(i,j,k,URHO);

}

#endif
