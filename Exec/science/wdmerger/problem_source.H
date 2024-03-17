#ifndef problem_source_H
#define problem_source_H

#include <prob_parameters.H>
#ifdef HYBRID_MOMENTUM
#include <hybrid.H>
#endif
#ifdef ROTATION
#include <Rotation.H>
#endif
#include <wdmerger_util.H>
#include <wdmerger_data.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time)
{
    using namespace problem;

    // First do any relaxation source terms.

    if (problem::problem == 1 && relaxation_damping_factor > 0.0_rt) {

        // The relevant dynamical timescale for determining this source term timescale should be
        // the smaller of the two WD timescales. Generally this should be the primary, but we'll
        // be careful just in case.

        const Real dynamical_timescale = amrex::min(t_ff_P, t_ff_S);

        // The relaxation damping factor should be less than unity, so that the damping
        // timescale is less than the dynamical timescale. This ensures that the stars
        // are always responding to the damping with quasistatic motion; if the stars
        // could respond too quickly, they might expand and make contact too early.

        const Real relaxation_damping_timescale = relaxation_damping_factor * dynamical_timescale;

        // Note that we are applying this update implicitly. The implicit and
        // explicit methods agree in the limit where the damping timescale is
        // much larger than dt, but the implicit method helps avoid numerical
        // problems when the damping timescale is shorter than the timestep.
        // For further information, see Source/sources/sponge_nd.F90.

        const auto damping_factor = -(1.0_rt - 1.0_rt / (1.0_rt + dt / relaxation_damping_timescale)) / dt;

        GpuArray<Real, 3> loc;
        GpuArray<Real, 3> mom;
#ifdef HYBRID_MOMENTUM
        GpuArray<Real, 3> hybrid_mom;
#endif
        GpuArray<Real, 3> Sr;

        const auto rhoInv = 1.0_rt / state(i,j,k,URHO);
        position(i, j, k, geomdata, loc);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            loc[dir] -= center[dir];
            mom[dir] = state(i,j,k,UMX+dir);
        }

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            Sr[dir] = mom[dir] * damping_factor;
            src(i,j,k,UMX+dir) += Sr[dir];
        }

#ifdef HYBRID_MOMENTUM
        linear_to_hybrid(loc, Sr, hybrid_mom);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UMR+dir) += hybrid_mom[dir];
        }
#endif

        // Do the same thing for the kinetic energy update.
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UEDEN) += rhoInv * mom[dir] * Sr[dir];
        }
    }

    // This is an inspiral driving force that is inspired by Pakmor et al. 2021
    // (https://academic.oup.com/mnras/article-abstract/503/4/4734/6166776),
    // where a constant drift velocity is applied (see Equation 1 and 2, and note that the
    // denominator in Equation 2 is missing a factor of a). Their approach simplifies to
    // dv/dt = -0.5 * v * (da/dt) / a, so we can achieve a similar effect by choosing as
    // a parameter the ratio [(da/dt) / a], which should be relatively small. To match their
    // method, we choose da/dt to be measured with respect to the initial velocity (relative to
    // the center of mass) so that it is fixed, while the a in the denominator will be the
    // current distance. We do not apply this force until the relaxation phase has completed.

    bool do_radial_damping = false;

    if (problem::problem == 1 && problem::radial_damping_velocity_factor > 0.0_rt && problem::radial_damping_is_done != 1) {
        do_radial_damping = true;
    }

    if (problem::problem == 1 && problem::relaxation_damping_factor > 0.0_rt && problem::relaxation_is_done != 1) {
        do_radial_damping = false;
    }

    if (do_radial_damping) {
        GpuArray<Real, 3> acceleration = {0.0};

        // Calculate the current stellar distance.

        Real a_curr = 0.0_rt;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            a_curr += (problem::com_P[dir] - problem::com_S[dir]) * (problem::com_P[dir] - problem::com_S[dir]);
        }
        a_curr = std::sqrt(a_curr);

        // Determine the initial orbital velocity of the reduced mass system.

        Real v_init = 2.0_rt * M_PI * problem::a / castro::rotational_period;

        if (stellar_mask(i, j, k, geomdata, state(i,j,k,URHO), true) > 0.0_rt) {
            GpuArray<Real, 3> v = {0.0_rt};
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                v[dir] = problem::vel_P[dir];
            }
            if (castro::do_rotation == 1) {
                rotational_to_inertial_velocity(i, j, k, geomdata, time, v);
            }

            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                acceleration[dir] = -0.5_rt * (problem::radial_damping_velocity_factor * v_init / a_curr) * v[dir];
            }
        }
        else if (stellar_mask(i, j, k, geomdata, state(i,j,k,URHO), false) > 0.0_rt) {
            GpuArray<Real, 3> v = {0.0_rt};
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                v[dir] = problem::vel_S[dir];
            }
            if (castro::do_rotation == 1) {
                rotational_to_inertial_velocity(i, j, k, geomdata, time, v);
            }

            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                acceleration[dir] = -0.5_rt * (problem::radial_damping_velocity_factor * v_init / a_curr) * v[dir];
            }
        }

        GpuArray<Real, 3> Sr = {0.0};
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            Sr[dir] = state(i,j,k,URHO) * acceleration[dir];
        }

#ifdef HYBRID_MOMENTUM
        GpuArray<Real, 3> loc;
        position(i, j, k, geomdata, loc);

        GpuArray<Real, 3> hybrid_Sr = {0.0};
        linear_to_hybrid(loc, Sr, hybrid_Sr);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UMR+dir) += hybrid_Sr[dir];
        }
#endif

        Real rhoInv = 1.0_rt / state(i,j,k,URHO);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UMX+dir) += Sr[dir];
            src(i,j,k,UEDEN) += rhoInv * state(i,j,k,UMX+dir) * Sr[dir];
        }
    }
}

#endif
