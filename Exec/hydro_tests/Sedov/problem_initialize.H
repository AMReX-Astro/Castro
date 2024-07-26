#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>

AMREX_INLINE
void problem_initialize ()
{
    const Geometry& dgeom = DefaultGeometry();

    const int coord_type = dgeom.Coord();
    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        if (coord_type == 1 || coord_type == 2) {
            problem::center[n] = 0.0_rt;
        }
        else {
            problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
        }
    }

    Real xn_zone[NumSpec] = {0.0};
    xn_zone[0] = 1.0_rt;

    eos_t eos_state;

    // override the pressure with the temperature
    if (problem::temp_ambient > 0.0_rt) {

        eos_state.rho = problem::dens_ambient;
        eos_state.T = problem::temp_ambient;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = xn_zone[n];
        }

        eos(eos_input_rt, eos_state);

        problem::p_ambient = eos_state.p;

    }

    // Calculate ambient state data

    eos_state.rho = problem::dens_ambient;
    eos_state.p   = problem::p_ambient;
    eos_state.T   = 1.e9_rt; // Initial guess for iterations
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = xn_zone[n];
    }

    eos(eos_input_rp, eos_state);

    problem::e_ambient = eos_state.e;
    problem::temp_ambient = eos_state.T;

    // set explosion pressure -- we will convert the point-explosion energy into
    // a corresponding pressure distributed throughout the perturbed volume

    Real vctr = 0.0;

    if (coord_type == 0) {

#if AMREX_SPACEDIM == 1

        amrex::Abort("Sedov problem unsupported in 1D Cartesian geometry.");

#elif AMREX_SPACEDIM == 2

        // Cylindrical problem in Cartesian coordinates

        vctr = M_PI * problem::r_init * problem::r_init;

#else

        // Spherical problem in Cartesian coordinates

        vctr = (4.0_rt / 3.0_rt) * M_PI * problem::r_init * problem::r_init * problem::r_init;

#endif

    }
    else if (coord_type == 1) {

#if AMREX_SPACEDIM == 1

        vctr = M_PI * problem::r_init * problem::r_init;

#elif AMREX_SPACEDIM == 2

        vctr = (4.0_rt / 3.0_rt) * M_PI * problem::r_init * problem::r_init * problem::r_init;

#else
        amrex::Abort("Sedov problem unsupported in 3D axisymmetric geometry.");
#endif

    }
    else if (coord_type == 2) {


#if AMREX_SPACEDIM == 1
        // Must have AMREX_SPACEDIM == 1 for this coord_type.

        vctr = (4.0_rt / 3.0_rt) * M_PI * problem::r_init * problem::r_init * problem::r_init;
#else
        amrex::Abort("Sedov problem unsupported in 2-D or 3-D spherical geometry.");
#endif

    }

    problem::e_exp = problem::exp_energy / vctr / problem::dens_ambient;
}

#endif
