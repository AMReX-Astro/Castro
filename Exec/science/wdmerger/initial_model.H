#ifndef initial_model_H
#define initial_model_H

#include <AMReX.H>

#include <network_properties.H>

namespace initial_model {
    const int initial_model_max_npts = 10000;

    struct model {
        // Physical characteristics

        amrex::Real mass;
        amrex::Real envelope_mass;
        amrex::Real central_density;
        amrex::Real central_temp;
        amrex::Real min_density;
        amrex::Real radius;

        // Composition

        amrex::Real core_comp[NumSpec];
        amrex::Real envelope_comp[NumSpec];

        // Model storage

        amrex::Real dx;
        int npts;
        amrex::Real mass_tol, hse_tol;

        amrex::Real r[initial_model_max_npts], rl[initial_model_max_npts], rr[initial_model_max_npts];
        amrex::Real M_enclosed[initial_model_max_npts], g[initial_model_max_npts];
    };
}

extern "C" {

void establish_hse (initial_model::model& model,
                    amrex::Real rho[initial_model::initial_model_max_npts],
                    amrex::Real T[initial_model::initial_model_max_npts],
                    amrex::Real xn[NumSpec][initial_model::initial_model_max_npts],
                    amrex::Real r[initial_model::initial_model_max_npts]);

}

#endif
