#include <Castro.H>
#include <Castro_prob_err_F.H>

using namespace amrex;

void
Castro::problem_pre_tagging_hook(TagBoxArray& tags,
                                 int          clearval,
                                 int          tagval,
                                 Real         time)
{

    // Determine if we have any hotspots
    // that are larger than the resolution
    // on this level. Skip it if we've already
    // detected an ignition.

    if (num_zones_ignited > 0 || ignition_level >= 0) return;

    auto ignition_radius_mf = derive("ignition_radius", time, 0);

    const Real* dx = geom.CellSize();

    MultiFab& state = get_new_data(State_Type);

    int num_zones = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:num_zones)
#endif
    for (MFIter mfi(*ignition_radius_mf,true); mfi.isValid(); ++mfi) {

        FArrayBox& ignition_fab = (*ignition_radius_mf)[mfi];

        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        find_ignited_zones(ARLIM_3D(lo), ARLIM_3D(hi),
                           BL_TO_FORTRAN_3D(ignition_fab),
                           BL_TO_FORTRAN_3D(state[mfi]),
                           ZFILL(dx), &num_zones, &level);

    }

    amrex::ParallelDescriptor::ReduceIntSum(num_zones);

    num_zones_ignited = num_zones;

    if (num_zones_ignited > 0)
        ignition_level = level;

    set_num_zones_ignited(&num_zones_ignited, &level);

}

void
Castro::problem_post_tagging_hook(TagBoxArray& tags,
                                  int          clearval,
                                  int          tagval,
                                  Real         time)
{
    // Do nothing here.
}
