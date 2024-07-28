
#include <Castro.H>

#include <diffusion_util.H>

using std::string;

#include <Diffusion.H>

void
Castro::construct_old_diff_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_old_diff_source()");

    amrex::ignore_unused(dt);

    const Real strt_time = ParallelDescriptor::second();

    MultiFab TempDiffTerm(grids, dmap, 1, 0);

    add_temp_diffusion_to_source(source, state_in, TempDiffTerm, time);

    if (verbose > 1)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;
        amrex::Real llevel = level;
#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::construct_old_diff_source() time = " << run_time
                       << " on level " << llevel << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

void
Castro::construct_new_diff_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_new_diff_source()");

    const Real strt_time = ParallelDescriptor::second();

    MultiFab TempDiffTerm(grids, dmap, 1, 0);

    Real mult_factor = 0.5;

    add_temp_diffusion_to_source(source, state_new, TempDiffTerm, time, mult_factor);

    // Time center the source term.

    mult_factor = -0.5;
    Real old_time = time - dt;

    add_temp_diffusion_to_source(source, state_old, TempDiffTerm, old_time, mult_factor);

    if (verbose > 1)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        amrex::Real run_time = ParallelDescriptor::second() - strt_time;
        amrex::Real llevel = level;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        amrex::Print() << "Castro::construct_new_diff_source() time = " << run_time
                       << " on level " << llevel << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

// **********************************************************************************************

void
Castro::add_temp_diffusion_to_source (MultiFab& ext_src, MultiFab& state_in, MultiFab& DiffTerm, Real t, Real mult_factor)
{
    BL_PROFILE("Castro::add_temp_diffusion_to_sources()");

    // Define an explicit temperature update.
    DiffTerm.setVal(0.);
    if (diffuse_temp == 1) {
        getTempDiffusionTerm(t, state_in, DiffTerm);
    }

    if (diffuse_temp == 1) {
       MultiFab::Saxpy(ext_src,mult_factor,DiffTerm,0,UEDEN,1,0);  // NOLINT(readability-suspicious-call-argument)
       MultiFab::Saxpy(ext_src,mult_factor,DiffTerm,0,UEINT,1,0);  // NOLINT(readability-suspicious-call-argument)
    }
}


void
Castro::getTempDiffusionTerm (Real time, MultiFab& state_in, MultiFab& TempDiffTerm)
{
    BL_PROFILE("Castro::getTempDiffusionTerm()");

   // Fill coefficients at this level.
   Vector<std::unique_ptr<MultiFab> > coeffs(AMREX_SPACEDIM);
   for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
       coeffs[dir] = std::make_unique<MultiFab>(getEdgeBoxArray(dir), dmap, 1, 0);
   }

   // Fill temperature at this level.
   MultiFab Temperature(grids, dmap, 1, 1);

   {
       FillPatchIterator fpi(*this, state_in, 1, time, State_Type, 0, NUM_STATE);
       MultiFab& grown_state = fpi.get_mf();

       MultiFab::Copy(Temperature, grown_state, UTEMP, 0, 1, 1);

#ifdef _OPENMP
#pragma omp parallel
#endif
       {
           FArrayBox coeff_cc;

           for (MFIter mfi(grown_state, TilingIfNotGPU()); mfi.isValid(); ++mfi)
           {

               const Box& bx = mfi.tilebox();

               // Create an array for storing cell-centered conductivity data.
               // It needs to have a ghost zone for the next step.

               const Box& obx = amrex::grow(bx, 1);
               coeff_cc.resize(obx, 1);
               Elixir elix_coeff_cc = coeff_cc.elixir();
               Array4<Real> const coeff_arr = coeff_cc.array();

               Array4<Real const> const U_arr = grown_state.array(mfi);

               fill_temp_cond(obx, U_arr, coeff_arr);

               for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

                   const Box& nbx = amrex::surroundingNodes(bx, idir);

                   Array4<Real> const edge_coeff_arr = (*coeffs[idir]).array(mfi);

                   AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                   {

                     if (idir == 0) {
                       edge_coeff_arr(i,j,k) = 0.5_rt * (coeff_arr(i,j,k) + coeff_arr(i-1,j,k));
                     } else if (idir == 1) {
                       edge_coeff_arr(i,j,k) = 0.5_rt * (coeff_arr(i,j,k) + coeff_arr(i,j-1,k));
                     } else {
                       edge_coeff_arr(i,j,k) = 0.5_rt * (coeff_arr(i,j,k) + coeff_arr(i,j,k-1));
                     }
                   });
               }
           }
       }

   }

   MultiFab CrseTemp;

   if (level > 0) {
       // Fill temperature at next coarser level, if it exists.
       const BoxArray& crse_grids = getLevel(level-1).boxArray();
       const DistributionMapping& crse_dmap = getLevel(level-1).DistributionMap();
       CrseTemp.define(crse_grids,crse_dmap,1,1);
       FillPatch(getLevel(level-1),CrseTemp,1,time,State_Type,UTEMP,1);
   }

   diffusion->applyop(level, Temperature, CrseTemp, TempDiffTerm, coeffs);

}
