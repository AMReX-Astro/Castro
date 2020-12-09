#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>

AMREX_INLINE
void problem_initialize ()
{

    const Geometry& dgeom = DefaultGeometry();

    const Real* problo = dgeom.ProbLo();
    const Real* probhi = dgeom.ProbHi();

    // set center variable

    // read the initial model

    read_model_file(problem::model_name);

    // do the burning

    for (int k = 0; k < model::npts; k++) {

        burn_t burn_state;

        burn_state.rho = model::profile(0).state(k, model::idens);
        burn_state.T = model::profile(0).state(k, model::itemp);
        burn_state.e = 0.0_rt;

        for (int n = 0; n < NumSpec; n++) {
            burn_state.xn[n] = model::profile(0).state(k, model::ispec+n);
        }

        burn_state.success = true;

        burner(burn_state, problem::burn_dt);

        if (! burn_state.success) {
            std::cout << "burning failed for zone " << k << std::endl;
        } else {
            std::cout << "zone, enuc = " << k << " " << burn_state.e / problem::burn_dt << std::endl;
        }
    }

    amrex::Error("done burning");
}
#endif
