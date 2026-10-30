// Implementations of functions in Problem.H go here

#include <Castro.H>

using namespace amrex;

void
Castro::flame_width_properties (amrex::Real time,
                                amrex::Real& T_max, amrex::Real& T_min,
                                amrex::Real& grad_T_max)
{
    BL_PROFILE("Castro::flame_width_properties()");

    const auto dx = geom.CellSizeArray();

    auto mf = derive("Temp",time,1);

    BL_ASSERT(mf != nullptr);

    ReduceOps<ReduceOpMax, ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(reduce_op);
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

            amrex::Real T = temp(i,j,k);
            amrex::Real grad_T = std::abs(temp(i+1,j,k) - temp(i-1,j,k)) / (2.0_rt * dx[0]);

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
Castro::flame_speed_properties (amrex::Real time, amrex::Real& rho_fuel_dot)
{
    BL_PROFILE("Castro::flame_speed_properties()");

    const auto dx = geom.CellSizeArray();

    std::string name;

    for (const auto & nm : short_spec_names_cxx) {
        if (nm == "He4") {
            name = "rho_omegadot_He4";
            break;
        }

        if (nm == "he4") {
            name = "rho_omegadot_he4";
            break;
        }
    }

  auto mf = derive(name, time, 0);
  BL_ASSERT(mf != nullptr);

  ReduceOps<ReduceOpSum> reduce_op;
  ReduceData<amrex::Real> reduce_data(reduce_op);
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


void
Castro::flame_peak_properties (amrex::Real time,
                               amrex::Real& rho_enuc_max, amrex::Real& peak_x)
{
    BL_PROFILE("Castro::flame_peak_properties()");

    const auto dx = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    auto mf = derive("rho_enuc", time, 0);
    BL_ASSERT(mf != nullptr);

    // find the peak energy generation rate and its index
    // these functions work on the entire MultiFab and do a reduction so
    // the information is available on all processors

    const amrex::Real level_rho_enuc_max = mf->max(0);
    const amrex::IntVect level_max_index = mf->maxIndex(0);

    // now get the physical coordinate of the maximum -- we are assuming
    // the flame is in the x-direction

    const amrex::Real level_peak_x = prob_lo[0] +
        (static_cast<amrex::Real>(level_max_index[0]) + 0.5_rt) * dx[0];

    // find the max over all previous levels

    if (level_rho_enuc_max > rho_enuc_max) {
        rho_enuc_max = level_rho_enuc_max;
        peak_x = level_peak_x;
    }
}
