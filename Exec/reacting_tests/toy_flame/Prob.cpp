/* Implementations of functions in Problem.H go here */

#include <Castro.H>

using namespace amrex;

void
Castro::flame_width_properties (Real time, Real& T_max, Real& T_min, Real& grad_T_max)
{
    BL_PROFILE("Castro::flame_width_properties()");

    const auto dx = geom.CellSizeArray();

    auto mf = derive("Temp",time,1);

    BL_ASSERT(mf);

    ReduceOps<ReduceOpMax, ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<Real, Real, Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();

        const auto temp = (*mf)[mfi].array();

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            // Assumes 1D simulation, right now. Also assumes that
            // we have at least one ghost cell in the x dimension.

            Real T = temp(i,j,k);
            Real grad_T = std::abs(temp(i+1,j,k) - temp(i-1,j,k)) / (2.0_rt * dx[0]);

            // Ignore problem zones where we have a negative temperature

            if (T < 0.0_rt) {
                T = 0.0_rt;
                grad_T = 0.0_rt;
            }

            return {T, T, grad_T};
        });
    }

    ReduceTuple hv = reduce_data.value();
    T_max = amrex::max(T_max, amrex::get<0>(hv));
    T_min = amrex::min(T_min, amrex::get<1>(hv));
    grad_T_max = amrex::max(grad_T_max, amrex::get<2>(hv));
}



void
Castro::flame_speed_properties (Real time, Real& rho_fuel_dot)
{
    BL_PROFILE("Castro::flame_speed_properties()");

    const auto dx = geom.CellSizeArray();

    auto mf = derive("rho_omegadot_fuel",time,0);

    BL_ASSERT(mf);

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& box  = mfi.tilebox();

        const auto omegadot = (*mf)[mfi].array();

        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            return {omegadot(i,j,k) * dx[0]};
        });
    }

    ReduceTuple hv = reduce_data.value();
    rho_fuel_dot += amrex::get<0>(hv);
}
