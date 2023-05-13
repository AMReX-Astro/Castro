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

    // Now do the radial drift source terms.

    if (problem::problem == 1 && radial_damping_factor > 0.0_rt) {

        // The logic for which dynamical timescale to use and
        // how to do the implicit coupling follows the reasoning
        // described above for the relaxation damping.

        const Real dynamical_timescale = amrex::max(t_ff_P, t_ff_S);

        const Real relaxation_damping_timescale = relaxation_damping_factor * dynamical_timescale;

        const auto damping_factor = -(1.0_rt - 1.0_rt / (1.0_rt + dt / relaxation_damping_timescale)) / dt;

        GpuArray<Real, 3> loc;
        GpuArray<Real, 3> vel;
        GpuArray<Real, 3> mom;
#ifdef HYBRID_MOMENTUM
        GpuArray<Real, 3> hybrid_mom;
#endif
        GpuArray<Real, 3> Sr;

        auto rhoInv = 1.0_rt / state(i,j,k,URHO);
        position(i, j, k, geomdata, loc);

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            loc[dir] -= center[dir];
            mom[dir] = state(i,j,k,UMX+dir);
            vel[dir] = rhoInv * mom[dir];
        }
        auto R_prp = std::sqrt(loc[axis_1-1]*loc[axis_1-1] + loc[axis_2-1]*loc[axis_2-1]);
        auto cosTheta = loc[axis_1-1] / R_prp;
        auto sinTheta = loc[axis_2-1] / R_prp;

        GpuArray<Real, 3> vel_i{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, vel_i);
#endif

        auto v_rad = cosTheta * vel_i[axis_1-1] + sinTheta * vel_i[axis_2-1];

        // What we want to do is insert a negative radial drift acceleration. If continued
        // for long enough, it will eventually drive coalescence of the binary. The
        // restriction on how large the acceleration can be is guided by the dynamical
        // properties of the system: it needs to be small enough that the WDs can be
        // in approximate dynamical equilibrium at all times before significant mass
        // transfer begins. So, if we write the force as
        // d(v_rad) / dt = -|v_phi| / tau,
        // where tau = radial_damping_factor * dynamical_timescale is the timescale
        // and |v_phi| is the magnitude of the azimuthal velocity, then
        // radial_damping_factor should be much greater than unity.

        Sr[axis_1-1] = cosTheta * (state(i,j,k,URHO) * std::abs(v_rad)) * damping_factor;
        Sr[axis_2-1] = sinTheta * (state(i,j,k,URHO) * std::abs(v_rad)) * damping_factor;
        Sr[axis_3-1] = 0.0_rt;

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UMX+dir) += Sr[dir];
        }

#ifdef HYBRID_MOMENTUM
        linear_to_hybrid(loc, Sr, hybrid_mom);
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UMR+dir) += hybrid_mom[dir];
        }
#endif

        // The kinetic energy source term is v . Sr:
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            src(i,j,k,UEDEN) += rhoInv * mom[dir] * Sr[dir];
        } 
    }

    // This is a different version of a inspiral driving force that matches the method
    // from Pakmor et al. 2021 (https://academic.oup.com/mnras/article-abstract/503/4/4734/6166776),
    // where a constant drift velocity is applied (see Equation 1 and 2).

    if (problem::problem == 1 && problem::radial_drift_velocity > 0.0_rt && problem::mass_S > 0.0_rt &&
        time < problem::radial_drift_time) {

        GpuArray<Real, 3> acceleration = {0.0};

        Real a_curr = 0.0_rt;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            a_curr += (problem::com_P[dir] - problem::com_S[dir]) * (problem::com_P[dir] - problem::com_S[dir]);
        }
        a_curr = std::sqrt(a_curr);

        Real mass_sum = problem::mass_P + problem::mass_S;
        Real prefactor = 0.5_rt * C::Gconst / (a_curr * mass_sum);

        if (stellar_mask(i, j, k, geomdata, state(i,j,k,URHO), true) > 0.0_rt) {
            GpuArray<Real, 3> v = {0.0_rt};
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                v[dir] = problem::vel_P[dir];
            }
            if (castro::do_rotation == 1) {
                rotational_to_inertial_velocity(i, j, k, geomdata, time, v);
            }

            Real v_sq = 0.0_rt;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                v_sq += v[dir] * v[dir];
            }
            v_sq = amrex::max(v_sq, std::numeric_limits<Real>::epsilon());

            Real M = problem::mass_S;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                acceleration[dir] = -prefactor * (M * M) * (v[dir] / v_sq) * problem::radial_drift_velocity;
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

            Real v_sq = 0.0_rt;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                v_sq += v[dir] * v[dir];
            }
            v_sq = amrex::max(v_sq, std::numeric_limits<Real>::epsilon());

            Real M = problem::mass_P;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                acceleration[dir] = -prefactor * (M * M) * (v[dir] / v_sq) * problem::radial_drift_velocity;
            }
        }

        // if (castro::do_rotation == 1) {
        //     inertial_to_rotational_velocity(i, j, k, geomdata, time, acceleration);
        // }

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
