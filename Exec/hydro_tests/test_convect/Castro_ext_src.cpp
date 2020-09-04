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
    
    const auto y_layer = 1.25e8_rt;
    const auto L_x = 2.5e8_rt;

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

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
            auto x = (Real(i)+0.5e0_rt)*dx[0] + prob_lo[0];
            auto y = (Real(j)+0.5e0_rt)*dx[1] + prob_lo[1];
            auto ey = std::exp(-(y-y_layer)*(y-y_layer)/1.e14_rt);

            auto H = ey * (1.e0_rt + 
                0.00625_rt * std::sin( 2*M_PI*x/L_x) 
              + 0.01875_rt * std::sin((6*M_PI*x/L_x) + M_PI/3.e0_rt) 
              + 0.01250_rt * std::sin((8*M_PI*x/L_x) + M_PI/5.e0_rt));

           // Source terms
           src(i,j,k,UEDEN) = snew(i,j,k,URHO) * H * 2.5e16_rt;
           src(i,j,k,UEINT) = src(i,j,k,UEDEN);
        });
    }
}
