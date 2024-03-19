#ifndef problem_initialize_H
#define problem_initialize_H

#include <prob_parameters.H>
#include <eos.H>
#include <network.H>
#include <network_properties.H>
#ifdef AUX_THERMO
#include <actual_network.H>
#endif

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

    bool species_defined = true;

    int ispec1 = network_spec_index(problem::spec1_name);
    if (ispec1 < 0) {
        species_defined = false;
    }

    int ispec2;
    if (!problem::spec2_name.empty()) {
        ispec2 = network_spec_index(problem::spec2_name);
        if (ispec2 < 0) {
            species_defined = false;
        }
    }

    int ispec3;
    if (!problem::spec3_name.empty()) {
        ispec3 = network_spec_index(problem::spec3_name);
        if (ispec3 < 0) {
            species_defined = false;
        }
    }

    int ispec4;
    if (!problem::spec4_name.empty()) {
        ispec4 = network_spec_index(problem::spec4_name);
        if (ispec4 < 0) {
            species_defined = false;
        }
    }

    int ispec5;
    if (!problem::spec5_name.empty()) {
        ispec5 = network_spec_index(problem::spec5_name);
        if (ispec5 < 0) {
            species_defined = false;
        }
    }

    int ispec6;
    if (!problem::spec6_name.empty()) {
        ispec6 = network_spec_index(problem::spec6_name);
        if (ispec6 < 0) {
            species_defined = false;
        }
    }

    if (! species_defined) {
        std::cout << ispec1 << " " << ispec2 << " " << ispec3 << " " << ispec4 << " " << ispec5 << " " << ispec6 << std::endl;
        amrex::Error("ERROR: species not defined");
    }

    Real xn[NumSpec] = {network_rp::small_x};
    // xn[0] = 1.0_rt - (NumSpec-1)*small_x;

    xn[ispec1] = problem::spec1_frac;
    if (!problem::spec2_name.empty()) {
        xn[ispec2] = problem::spec2_frac;
    }
    if (!problem::spec3_name.empty()) {
        xn[ispec3] = problem::spec3_frac;
    }
    if (!problem::spec4_name.empty()) {
        xn[ispec4] = problem::spec4_frac;
    }
    if (!problem::spec5_name.empty()) {
        xn[ispec5] = problem::spec5_frac;
    }
    if (!problem::spec6_name.empty()) {
        xn[ispec6] = problem::spec6_frac;
    }

    // Normalize species
    auto sumX = 0.0_rt;
    for (auto X : xn) {
        sumX += X;
    }
    for (auto & X : xn) {
        X /= sumX;
    }

    eos_t eos_state;

    eos_state.rho = problem::rho0;
    eos_state.T = problem::T0;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = xn[n];
    }

#ifdef AUX_THERMO
    Real aux[NumAux] = {0.0_rt};

    for (int n = 0; n < NumSpec; n++) {
        aux[AuxZero::iye] += xn[n] * zion[n] * aion_inv[n];
        aux[AuxZero::iabar] += xn[n] * aion_inv[n];
        aux[AuxZero::ibea] += xn[n] * network::bion(n+1) * aion_inv[n];
    }
    aux[AuxZero::iabar] = 1.0_rt/aux[AuxZero::iabar];

    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = aux[n];
    }
#endif

    eos(eos_input_rt, eos_state);

    problem::p0 = eos_state.p;
    problem::s0 = eos_state.s;
}

#endif
