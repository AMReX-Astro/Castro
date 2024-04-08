#ifndef problem_initialize_state_data_H
#define problem_initialize_state_data_H

#include <prob_parameters.H>
#include <eos.H>
#include <model_parser.H>
#ifdef NSE_TABLE
#include <nse_table_type.H>
#include <nse_table_check.H>
#include <nse_eos.H>
#endif

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_state_data (int i, int j, int k,
                                    Array4<Real> const& state,
                                    const GeometryData& geomdata)
{

    const Real* dx = geomdata.CellSize();
    const Real* problo = geomdata.ProbLo();

    Real x = problo[0] + dx[0] * (static_cast<Real>(i) + 0.5_rt) - problem::center[0];

    Real y = 0.0;
#if AMREX_SPACEDIM >= 2
    y = problo[1] + dx[1] * (static_cast<Real>(j) + 0.5_rt) - problem::center[1];
#endif

    Real z = 0.0;
#if AMREX_SPACEDIM == 3
    z = problo[2] + dx[2] * (static_cast<Real>(k) + 0.5_rt) - problem::center[2];
#endif

    Real dist = std::sqrt(x * x + y * y + z * z);

    state(i,j,k,URHO) = interpolate(dist, model::idens);
    state(i,j,k,UTEMP) = interpolate(dist, model::itemp);
    Real pres = interpolate(dist, model::ipres);
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = amrex::max(interpolate(dist, model::ispec+n), small_x);
    }


    Real sumX = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
        sumX += state(i,j,k,UFS+n);
    }
    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,UFS+n) / sumX;
    }

    // also get Ye from the model

    state(i,j,k,UFX+AuxZero::iye) = interpolate(dist, model::iaux+AuxZero::iye);

    burn_t burn_state;
    burn_state.rho = state(i,j,k,URHO);
    burn_state.T = state(i,j,k,UTEMP);
    for (int n = 0; n < NumSpec; n++) {
        burn_state.xn[n] = state(i,j,k,UFS+n);
    }

    burn_state.aux[AuxZero::iye] = state(i,j,k,UFX+AuxZero::iye);

#ifdef SIMPLIFIED_SDC
    // if we are doing simplified-SDC + NSE, then the `in_nse()`
    // check will use burn_state.y[], so we need to ensure that
    // those are initialized
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.y[SFS+n] = burn_state.rho * burn_state.xn[n];
    }

    // not used
    burn_state.y[SEINT] = 0.0;
#endif


    // now we use Ye if we are in NSE and use X if we are not.  In
    // either case, we sync up the other composition variables.

    auto nse_check = in_nse(burn_state);

    if (nse_check) {

        nse_table_t nse_state;
        nse_state.T = state(i,j,k,UTEMP);
        nse_state.rho = state(i,j,k,URHO);
        nse_state.Ye = state(i,j,k,UFX+AuxZero::iye);
        nse_interp(nse_state);

        state(i,j,k,UFX+AuxZero::iabar) = nse_state.abar;
        state(i,j,k,UFX+AuxZero::ibea) = nse_state.bea;
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = amrex::max(nse_state.X[n], small_x);
        }

        // we just got new X's, so we need to re-renormalize them
        sumX = 0.0_rt;
        for (int n = 0; n < NumSpec; n++) {
            sumX += state(i,j,k,UFS+n);
        }
        for (int n = 0; n < NumSpec; n++) {
            state(i,j,k,UFS+n) = state(i,j,k,UFS+n) / sumX;
        }

    } else {
        // we are not in NSE, so set the aux quantities from X-- this
        // will fill burn_state.aux[]

        set_aux_comp_from_X(burn_state);

        for (int n = 0; n < NumAux; n++) {
            state(i,j,k,UFX+n) = burn_state.aux[n];
        }
    }

    // now call the EOS -- we need to use an eos_t here since we are
    // going to take pressure as an input

    eos_t eos_state;
    eos_state.rho = state(i,j,k,URHO);
    eos_state.T = state(i,j,k,UTEMP);
    eos_state.p = pres;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = state(i,j,k,UFS+n);
    }
    for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = state(i,j,k,UFX+n);
    }

    if (problem::interpolate_pres == 1) {

        // we need to get T from P, rho, but consistent with NSE (if
        // we are in NSE)
        if (nse_check) {
            // this will also change the composition, since the new T
            // alters the NSE
            nse_T_abar_from_p(eos_state.rho, eos_state.p, eos_state.aux[AuxZero::iye],
                              eos_state.T, eos_state.aux[AuxZero::iabar]);
            state(i,j,k,UTEMP) = eos_state.T;
            state(i,j,k,UFX+AuxZero::iabar) = eos_state.aux[AuxZero::iabar];

            // now call the EOS with the new T and abar to get e
            eos(eos_input_rt, eos_state);

            // finally, get the updated B/A
            nse_table_t nse_state;
            nse_state.T = state(i,j,k,UTEMP);
            nse_state.rho = state(i,j,k,URHO);
            nse_state.Ye = state(i,j,k,UFX+AuxZero::iye);
            nse_interp(nse_state);

            state(i,j,k,UFX+AuxZero::ibea) = nse_state.bea;

        } else {
            eos(eos_input_rp, eos_state);
            state(i,j,k,UTEMP) = eos_state.T;
        }
    } else {
        eos(eos_input_rt, eos_state);
    }

    state(i,j,k,UEINT) = state(i,j,k,URHO) * eos_state.e;
    state(i,j,k,UEDEN) = state(i,j,k,URHO) * eos_state.e;

    // add density scaling to composition

    for (int n = 0; n < NumSpec; n++) {
        state(i,j,k,UFS+n) = state(i,j,k,URHO) * state(i,j,k,UFS+n);
    }

    for (int n = 0; n < NumAux; n++) {
        state(i,j,k,UFX+n) = state(i,j,k,URHO) * state(i,j,k,UFX+n);
    }

    // initial velocities = 0
    state(i,j,k,UMX) = 0.0_rt;
    state(i,j,k,UMY) = 0.0_rt;
    state(i,j,k,UMZ) = 0.0_rt;

    if (problem::perturb_model) {
        // random numbers between -1 and 1
        GpuArray<Real, 27> alpha{};
        GpuArray<Real, 27> beta{};
        GpuArray<Real, 27> gamma{};

        // random numbers between 0 and 2*pi
        GpuArray<Real, 27> phix{};
        GpuArray<Real, 27> phiy{};
        GpuArray<Real, 27> phiz{};

        // compute the norm of k
        GpuArray<Real, 27> normk{};

        for (int kk = 1; kk <= 3; ++kk) {
            for (int jj = 1; jj <= 3; ++jj) {
                for (int ii = 1; ii <= 3; ++ii) {
                    const int n = ii - 1 + 3 * (jj - 1 + 3 * (kk - 1));

                    alpha[n] = std::pow(0.5_rt, ii) * std::pow(0.7_rt, jj) * std::pow(0.3_rt, kk) * std::pow(-1.0_rt, ii);
                    beta[n] = std::pow(0.5_rt, ii) * std::pow(0.3_rt, jj) * std::pow(0.7_rt, kk) * std::pow(-1.0_rt, jj);
                    gamma[n] = std::pow(0.3_rt, ii) * std::pow(0.5_rt, jj) * std::pow(0.7_rt, kk) * std::pow(-1.0_rt, kk);

                    phix[n] = std::pow(0.3_rt, ii) * std::pow(0.7_rt, jj) * std::pow(0.5_rt, kk);
                    phix[n] *= 2.0_rt * M_PI;
                    phiy[n] = std::pow(0.7_rt, ii) * std::pow(0.3_rt, jj) * std::pow(0.5_rt, kk);
                    phiy[n] *= 2.0_rt * M_PI;
                    phiz[n] = std::pow(0.7_rt, ii) * std::pow(0.5_rt, jj) * std::pow(0.3_rt, kk);
                    phiz[n] *= 2.0_rt * M_PI;

                    normk[n] = std::sqrt(Real(ii * ii) + Real(jj * jj) + Real(kk * kk));
                }
            }
        }

        // set perturbational velocity to zero
        Real vpert[3] = {0.0_rt};

        // loop over the 27 combinations of Fourier components
        for (int kk = 1; kk <= 3; ++kk) {
            for (int jj = 1; jj <= 3; ++jj) {
                for (int ii = 1; ii <= 3; ++ii) {
                    // array index
                    const int n = ii - 1 + 3 * (jj - 1 + 3 * (kk - 1));
                    // compute cosines and sines
                    Real cx = std::cos(2.0_rt * M_PI * Real(ii) *
                                        (x + problem::center[0]) /
                                        problem::velpert_scale + phix[n]);
                    Real cy = std::cos(2.0_rt * M_PI * Real(jj) *
                                        (y + problem::center[1]) /
                                        problem::velpert_scale + phiy[n]);
                    Real cz = std::cos(2.0_rt * M_PI * Real(kk) *
                                        (z + problem::center[2]) /
                                        problem::velpert_scale + phiz[n]);

                    Real sx = std::sin(2.0_rt * M_PI * Real(ii) *
                                        (x + problem::center[0]) /
                                        problem::velpert_scale + phix[n]);
                    Real sy = std::sin(2.0_rt * M_PI * Real(jj) *
                                        (y + problem::center[1]) /
                                        problem::velpert_scale + phiy[n]);
                    Real sz = std::sin(2.0_rt * M_PI * Real(kk) *
                                        (z + problem::center[2]) /
                                        problem::velpert_scale + phiz[n]);

                    // compute contribution from perturbation velocity from each mode
                    vpert[0] += (-gamma[n] * Real(jj) * cx * cz * sy +
                                    beta[n] * Real(kk) * cx * cy * sz) / normk[n];

                    vpert[1] += (gamma[n] * Real(ii) * cy * cz * sx -
                                    alpha[n] * Real(kk) * cx * cy * sz) / normk[n];

                    vpert[2] += (-beta[n] * Real(ii) * cy * cz * sx +
                                    alpha[n] * Real(jj) * cx * cz * sy) / normk[n];
                }
            }
        }

        // convolve perturbation by Si28 mass fraction

        // apply the cutoff function to the perturbational velocity
        // add perturbational velocity to background velocity
        for (auto n = 0; n < 3; ++n) {
            vpert[n] *= problem::velpert_amplitude * state(i, j, k, UFS + problem::isi28);

            state(i, j, k, UMX + n) += vpert[n];
        }
    }
}

#endif
