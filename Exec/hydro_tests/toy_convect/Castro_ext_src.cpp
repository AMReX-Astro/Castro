#include "Castro.H"
#include "prob_parameters.H"

using namespace amrex;

void
Castro::fill_ext_source (const Real time, const Real dt, 
                         const MultiFab& state_old, 
                         const MultiFab& state_new, MultiFab& ext_src)
{
    // Compute the external sources for all the conservative equations.
    //
    // This is called twice in the evolution:
    //
    // First, for the predictor, it is called with (old, old) states.
    //
    // This is also used in the first pass of the conservative update
    // (adding dt * S there).
    //
    // Next we correct the source terms in the conservative update to
    // time-center them.  Here we call ext_src(old, new), and then
    // in time_center_source_terms we subtract off 1/2 of the first S
    // and add 1/2 of the new S.
    //
    // Therefore, to get a properly time-centered source, generally
    // speaking, you always want to use the "new" state here.  That
    // will be the time n state in the first call and the n+1 in the
    // second call.

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    ext_src.setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(ext_src, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real const> const snew = state_new.array(mfi);
        Array4<Real> const src = ext_src.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            auto rho = snew(i,j,k,URHO);

            auto T6 = snew(i,j,k,UTEMP) / 1.0e6_rt;
            auto T613 = std::pow(T6, 1.0_rt/3.0_rt);

            // CNO abundance
            auto X_CNO = (snew(i,j,k,UFS-1+Species::C12) +
                          snew(i,j,k,UFS-1+Species::N14) +
                          snew(i,j,k,UFS-1+Species::O16)) / rho;

            // H abundance
            auto X_1 = snew(i,j,k,UFS-1+Species::H1) / rho;

            // CNO heating from Kippenhahn & Weigert, Eq. 18.65
            auto g14 = 1.0_rt + 2.7e-3_rt*T613 - 7.78e-3_rt*T613*T613 - 1.49e-4_rt*T6;
            auto eps_CNO = 8.67e27_rt * g14 * X_CNO * X_1 * rho * std::exp(-152.28_rt/T613) / (T613*T613);

            // source terms
            src(i,j,k,UEDEN) = rho * eps_CNO;
            src(i,j,k,UEINT) = rho * eps_CNO;
        });
    }
}
