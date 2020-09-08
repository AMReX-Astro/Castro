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
        Array4<Real> const src = ext_src.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            auto y = prob_lo[1] + (Real(j) + 0.5_rt) * dx[1];

            if (y < 1.125_rt * 4.8e8_rt) {
                auto fheat = std::sin(8.0_rt * M_PI * (y / 4.e8_rt - 1.0_rt));

                src(i,j,k,UEDEN) = heating_factor * fheat;
                src(i,j,k,UEINT) = src(i,j,k,UEDEN);
            }
        });
    }
}
