#ifndef problem_initialize_H
#define problem_initialize_H

#include <fundamental_constants.H>
#include <model_parser.H>
#include <ambient.H>

AMREX_INLINE
void problem_initialize ()
{
    // Set the center (where the star will start) to be along the
    // x-axis, halfway to the left edge.

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    problem::center[0] = problo[0] + 0.25_rt * (probhi[0] - problo[0]);
    problem::center[1] = 0.0_rt;
    problem::center[2] = 0.0_rt;

    // Arbitrarily set the composition of the star to 50/50 C/O.

    Real core_comp[NumSpec] = {0.0_rt};
    core_comp[Species::C12-1] = 0.5_rt;
    core_comp[Species::O16-1] = 0.5_rt;

    Real envelope_comp[NumSpec] = {0.0_rt};
    envelope_comp[Species::C12-1] = 0.5_rt;
    envelope_comp[Species::O16-1] = 0.5_rt;

    Real envelope_mass = 0.0_rt;

    Real initial_model_dx = 1.0e6_rt;

    // Generate the WD model.

    problem::wd_mass *= C::M_solar;

    establish_hse(problem::wd_mass, problem::wd_central_density, problem::wd_radius,
                  core_comp, problem::wd_temperature, initial_model_dx,
                  envelope_mass, envelope_comp, 0);

    problem::wd_mass /= C::M_solar;

    amrex::Print() << std::endl;

    amrex::Print() << "Generated initial model for WD of mass " << std::setprecision(3) << problem::wd_mass
                   << " solar masses, central density " << std::setprecision(3) << std::scientific << problem::wd_central_density
                   << " g cm**-3, and radius " << std::setprecision(3) << std::scientific << problem::wd_radius << " cm."
                   << std::endl << std::endl;

    // Safety check: make sure the star is actually inside the computational domain.

    if ((0.5_rt * (probhi[0] - problo[0]) < problem::wd_radius) ||
        (0.5_rt * (probhi[1] - problo[1]) < problem::wd_radius) ||
        (0.5_rt * (probhi[2] - problo[2]) < problem::wd_radius))
    {
        amrex::Error("WD does not fit inside the domain.");
    }

    // Set ambient state composition.

    for (int n = 0; n < NumSpec; ++n) {
        ambient::ambient_state[UFS+n] = ambient::ambient_state[URHO] * envelope_comp[n];
    }
}

#endif
