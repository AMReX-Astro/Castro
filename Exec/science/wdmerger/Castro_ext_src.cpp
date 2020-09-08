#include "Castro.H"
#include "prob_parameters.H"

using namespace amrex;

void
Castro::fill_ext_source (const Real time, const Real dt, 
                         const MultiFab& state_old, 
                         const MultiFab& state_new, MultiFab& ext_src)
{
    // Compute the external sources for all the conservative equations.
    //
    // This is called twice in the evolution:
    //
    // First, for the predictor, it is called with (old, old) states.
    //
    // This is also used in the first pass of the conservative update
    // (adding dt * S there).
    //
    // Next we correct the source terms in the conservative update to
    // time-center them.  Here we call ext_src(old, new), and then
    // in time_center_source_terms we subtract off 1/2 of the first S
    // and add 1/2 of the new S.
    //
    // Therefore, to get a properly time-centered source, generally
    // speaking, you always want to use the "new" state here.  That
    // will be the time n state in the first call and the n+1 in the
    // second call.

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real const> const sold = state_old.array(mfi);
        Array4<Real const> const snew = state_new.array(mfi);
        Array4<Real> const src = ext_src.array(mfi);

        // First do any relaxation source terms.

        if (problem == 1 && relaxation_damping_factor > 0.0_rt) {

            GpuArray<Real, 3> center;
            ca_get_center(center.begin());

            GpuArray<Real, 3> loc;
            GpuArray<Real, 3> vel;
            GpuArray<Real, 3> mom;
            GpuArray<Real, 3> Sr;

            // The relevant dynamical timescale for determining this source term timescale should be
            // the smaller of the two WD timescales. Generally this should be the primary, but we'll
            // be careful just in case.

            auto dynamical_timescale = amrex::min(t_ff_P, t_ff_S);

            // The relaxation damping factor should be less than unity, so that the damping
            // timescale is less than the dynamical timescale. This ensures that the stars
            // are always responding to the damping with quasistatic motion; if the stars
            // could respond too quickly, they might expand and make contact too early.

            auto relaxation_damping_timescale = relaxation_damping_factor * dynamical_timescale;

            // Note that we are applying this update implicitly. The implicit and
            // explicit methods agree in the limit where the damping timescale is
            // much larger than dt, but the implicit method helps avoid numerical
            // problems when the damping timescale is shorter than the timestep.
            // For further information, see Source/sources/sponge_nd.F90.

            auto damping_factor = -(1.0_rt - 1.0_rt / (1.0_rt + dt / relaxation_damping_timescale)) / dt;

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                auto rhoInv = 1.0_rt / new_state(i,j,k,URHO);
                loc = position(i,j,k) - center;
                mom = new_state(i,j,k,UMX:UMZ);

#ifdef ROTATION
                if (do_rotation == 1 && state_in_rotating_frame == 0) {
                    vel = rhoInv * mom
                    inertial_to_rotational_velocity([i, j, k], time, vel)
                    mom = new_state(i,j,k,URHO) * vel
                }
#endif

                auto Sr = mom * damping_factor;

                src(i,j,k,UMX:UMZ) = src(i,j,k,UMX:UMZ) + Sr

#ifdef HYBRID_MOMENTUM
                src(i,j,k,UMR:UMP) = src(i,j,k,UMR:UMP) + linear_to_hybrid(loc, Sr)
#endif

                // Do the same thing for the kinetic energy update.

                src(i,j,k,UEDEN) = src(i,j,k,UEDEN) + dot_product(rhoInv * mom, Sr)
            });
        }

        // Now do the radial drift source terms.

        if (problem == 1 && radial_damping_factor > 0.0_rt) {

            // The logic for which dynamical timescale to use and
            // how to do the implicit coupling follows the reasoning
            // described above for the relaxation damping.

            auto dynamical_timescale = amrex::max(t_ff_P, t_ff_S)

            auto radial_damping_timescale = radial_damping_factor * dynamical_timescale

            auto damping_factor = -(1.0_rt - 1.0_rt / (1.0_rt + dt / radial_damping_timescale)) / dt

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                auto rhoInv = 1.0_rt / new_state(i,j,k,URHO)

                auto loc      = position(i,j,k) - center
                auto R_prp    = sqrt(loc(axis_1)**2 + loc(axis_2)**2)
                auto cosTheta = loc(axis_1) / R_prp
                auto sinTheta = loc(axis_2) / R_prp

                auto mom = new_state(i,j,k,UMX:UMZ)
                auto vel = rhoInv * mom
                vel = inertial_velocity(loc, vel, time)

                auto v_rad = cosTheta * vel(axis_1) + sinTheta * vel(axis_2)

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

                Sr(axis_1) = cosTheta * (new_state(i,j,k,URHO) * abs(v_rad)) * damping_factor;
                Sr(axis_2) = sinTheta * (new_state(i,j,k,URHO) * abs(v_rad)) * damping_factor;
                Sr(axis_3) = 0.0_rt;

                src(i,j,k,UMX:UMZ) = src(i,j,k,UMX:UMZ) + Sr;

#ifdef HYBRID_MOMENTUM
                src(i,j,k,UMR:UMP) = src(i,j,k,UMR:UMP) + linear_to_hybrid(loc, Sr);
#endif

                // The kinetic energy source term is v . Sr:

                src(i,j,k,UEDEN) = src(i,j,k,UEDEN) + dot_product(rhoInv * mom, Sr);
            });
        }
    }
}
