#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        problem::center[n] = 0.5_rt * (problo[n] + probhi[n]);
    }

    // get the species indices

    problem::ihe4 = network_spec_index("helium-4");
    problem::ic12 = network_spec_index("carbon-12");
    problem::in14 = network_spec_index("nitrogen-14");
    problem::io16 = network_spec_index("oxygen-16");

    if (problem::ihe4 < 0 || problem::ic12 < 0 || problem::io16 < 0) {
        amrex::Error("ERROR: species indices not found");
    }

    // make sure that the carbon fraction falls between 0 and 1

    if (problem::cfrac > 1.0_rt || problem::cfrac < 0.0_rt) {
        amrex::Error("ERROR: cfrac must fall between 0 and 1");
    }

    // make sure that the nitrogen fraction falls between 0 and 1

    if (problem::nfrac > 1.0_rt || problem::nfrac < 0.0_rt) {
        amrex::Error("ERROR: nfrac must fall between 0 and 1");
    }

    // make sure that the oxygen fraction falls between 0 and 1

    if (problem::ofrac > 1.0_rt || problem::ofrac < 0.0_rt) {
        amrex::Error("ERROR: ofrac must fall between 0 and 1");
    }

    // make sure that the C/O fraction sums to no more than 1

    if (problem::cfrac + problem::nfrac + problem::ofrac > 1.0_rt) {
        amrex::Error("ERROR: cfrac + nfrac + ofrac cannot exceed 1.");
    }

    // set the default mass fractions

    for (int n = 0; n < NumSpec; n++) {
        problem::xn[n] = problem::smallx;
    }

    problem::xn[problem::ic12] = amrex::max(problem::cfrac, problem::smallx);
    problem::xn[problem::io16] = amrex::max(problem::ofrac, problem::smallx);

    if (problem::in14 >= 0) {
       problem::xn[problem::in14] = amrex::max(problem::nfrac, problem::smallx);
       problem::xn[problem::ihe4] = 1.0_rt - problem::xn[problem::ic12]
                                           - problem::xn[problem::in14]
                                           - problem::xn[problem::io16]
                                           - (NumSpec - 4) * problem::smallx;
    }
    else {
       problem::xn[problem::ihe4] = 1.0_rt - problem::xn[problem::ic12]
                                           - problem::xn[problem::io16]
                                           - (NumSpec - 3) * problem::smallx;
    }

}

#endif
