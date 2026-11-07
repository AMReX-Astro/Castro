#include <Castro.H>

using namespace amrex;

void
Castro::edge_state_temp_to_pres(const Box& bx,
                                Array4<Real> const& qm,
                                Array4<Real> const& qp)
{

    // use T to define p

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
