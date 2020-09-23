#ifndef CASTRO_EXT_SRC_H
#define CASTRO_EXT_SRC_H

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

AMREX_GPU_HOST_DEVICE 
void do_ext_src(int i, int j, int k,
                GeometryData const& geomdata, 
                Array4<const Real> const& state,
                Array4<Real> const& src,
                GpuArray<Real, 3> const& omega,
                const Real dt, const Real time)
{
    // First do any relaxation source terms.

    if (problem == 1 && relaxation_damping_factor > 0.0_rt) {

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
        GpuArray<Real, 3> vel;
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

#ifdef ROTATION
        if (do_rotation == 1 && state_in_rotating_frame == 0) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                vel[dir] = rhoInv * mom[dir];
            }
            inertial_to_rotational_velocity(i, j, k, geomdata, omega, time, vel);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                mom[dir] = state(i,j,k,URHO) * vel[dir];
            }
        }
#endif

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

    if (problem == 1 && radial_damping_factor > 0.0_rt) {

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

        auto vel_i = inertial_velocity(loc, vel, omega);

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
}

#endif
