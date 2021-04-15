#include <Castro.H>
#include <Castro_hydro.H>

using namespace amrex;

void
Castro::reset_edge_state_thermo(const Box& bx,
                                Array4<Real> const& qedge)
{

    int use_eos = transverse_use_eos;
    int reset_rhoe = transverse_reset_rhoe;

    Real small_t = small_temp;
    Real small_p = small_pres;

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

#ifdef RADIATION
        Real old_p_state = qedge(i,j,k,QPRES);
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

    });

}



void
Castro::edge_state_temp_to_pres(const Box& bx,
                                Array4<Real> const& qm,
                                Array4<Real> const& qp)
{

    // use T to define p

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // We just got the extremes corresponding to a particular cell-center, but now
        // we need to assign them to interfaces.

        eos_rep_t eos_state;

        eos_state.rho    = qp(i,j,k,QRHO);
        eos_state.T      = qp(i,j,k,QTEMP);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = qp(i,j,k,QFS+n);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = qp(i,j,k,QFX+n);
        }
#endif

        eos(eos_input_rt, eos_state);

        qp(i,j,k,QPRES) = eos_state.p;
        qp(i,j,k,QREINT) = qp(i,j,k,QRHO) * eos_state.e;
        // should we try to do something about Gamma_! on interface?

        eos_state.rho    = qm(i,j,k,QRHO);
        eos_state.T      = qm(i,j,k,QTEMP);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = qm(i,j,k,QFS+n);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = qm(i,j,k,QFX+n);
        }
#endif

        eos(eos_input_rt, eos_state);

        qm(i,j,k,QPRES) = eos_state.p;
        qm(i,j,k,QREINT) = qm(i,j,k,QRHO) * eos_state.e;
        // should we try to do something about Gamma_! on interface?

    });

}
