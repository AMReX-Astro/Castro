#include <Castro.H>
#include <Castro_hydro.H>

using namespace amrex;

void
Castro::reset_edge_state_thermo(const Box& bx,
                                Array4<Real> const qedge)
{

    int use_eos = transverse_use_eos;
    int reset_rhoe = transverse_reset_rhoe;

    Real small_t = small_temp;
    Real small_p = small_pres;

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

#ifdef RADIATION
        Real old_p_state = qedge(i,j,k,QPRES);
#endif

        eos_t eos_state;

        if (reset_rhoe == 1) {
            // if we are still negative, then we need to reset
            if (qedge(i,j,k,QREINT) < 0.0_rt) {

                eos_state.rho = qedge(i,j,k,QRHO);
                eos_state.T = small_t;
                for (int n = 0; n < NumSpec; ++n) {
                    eos_state.xn[n] = qedge(i,j,k,QFS+n);
                }
                for (int n = 0; n < NumAux; ++n) {
                    eos_state.aux[n] = qedge(i,j,k,QFX+n);
                }

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
            for (int n = 0; n < NumAux; ++n) {
                eos_state.aux[n] = qedge(i,j,k,QFX+n);
            }

            eos(eos_input_re, eos_state);

            qedge(i,j,k,QREINT) = eos_state.e * eos_state.rho;
            qedge(i,j,k,QPRES) = max(eos_state.p, small_p);
        }

#ifdef RADIATION
        // correct the total pressure (gas + radiation) with any
        // change to the gas pressure state
        qedge(i,j,k,QPTOT) = qedge(i,j,k,QPTOT) + (qedge(i,j,k,QPRES) - old_p_state);
#endif

    });

}
