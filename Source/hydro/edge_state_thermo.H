#ifndef CASTRO_EDGE_STATE_THERMO_H
#define CASTRO_EDGE_STATE_THERMO_H

#include <AMReX_Array4.H>
#include <Castro_util.H>
#include <eos.H>

namespace hydro {

// Apply the final thermodynamic repair after the transverse update, including
// its radiation bookkeeping. Keep this in the producing kernel to avoid a
// separate launch and another pass over the edge states.
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void reset_edge_state_thermo(int i, int j, int k,
                             amrex::Array4<amrex::Real> const& qedge,
                             bool use_eos, bool reset_rhoe,
                             amrex::Real small_t, amrex::Real small_p)
{
    if (!use_eos && !reset_rhoe) {
        return;
    }

#ifdef RADIATION
    amrex::Real old_p_state = qedge(i,j,k,QPRES);
#endif

    eos_rep_t eos_state;

    if (reset_rhoe == 1) {
        // if we are still negative, then we need to reset
        if (qedge(i,j,k,QREINT) < 0.0_rt) {

            eos_state.rho = qedge(i,j,k,QRHO);
            eos_state.T = small_t;
            for (int n = 0; n < NumSpec; ++n) {
                eos_state.xn[n] = qedge(i,j,k,QFS+n);
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; ++n) {
                eos_state.aux[n] = qedge(i,j,k,QFX+n);
            }
#endif

            eos(eos_input_rt, eos_state);

            qedge(i,j,k,QREINT) = qedge(i,j,k,QRHO) * eos_state.e;
            qedge(i,j,k,QPRES) = eos_state.p;
        }
    }

    if (use_eos == 1) {
        eos_state.rho = qedge(i,j,k,QRHO);
        eos_state.e   = qedge(i,j,k,QREINT) / qedge(i,j,k,QRHO);
        eos_state.T   = small_t;
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = qedge(i,j,k,QFS+n);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = qedge(i,j,k,QFX+n);
        }
#endif

        eos(eos_input_re, eos_state);

        qedge(i,j,k,QREINT) = eos_state.e * eos_state.rho;
        qedge(i,j,k,QPRES) = amrex::max(eos_state.p, small_p);
    }

#ifdef RADIATION
    // correct the total pressure (gas + radiation) with any
    // change to the gas pressure state
    qedge(i,j,k,QPTOT) = qedge(i,j,k,QPTOT) + (qedge(i,j,k,QPRES) - old_p_state);
#endif
}

}

#endif
