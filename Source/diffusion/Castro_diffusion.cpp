
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

   if (diffuse_use_amrex_mlmg) {
       // Evaluates ∇ ⋅(k_th ∇T) using AMReX
       diffusion->applyop(level, Temperature, CrseTemp, TempDiffTerm, coeffs);
   } else {
       // Evaluates ∇ ⋅(k_th ∇T) without using AMReX
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(Temperature, TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
           const auto dx = geom.CellSizeArray();
           const auto dxinv = geom.InvCellSizeArray();
           const auto problo = geom.ProbLoArray();
           const int coord = geom.Coord();

           const Box& bx = mfi.tilebox();
           Array4<Real const> const& Temp_array = Temperature.array(mfi);
           Array4<Real> const & TempDiff_array = TempDiffTerm.array(mfi);

           // edged based k_th in different averaging direction.
           const Vector<Array4<Real const>> edge_coeff_arrs {AMREX_D_DECL((*coeffs[0]).array(mfi),
                                                                          (*coeffs[1]).array(mfi),
                                                                          (*coeffs[2]).array(mfi))};

           ParallelFor(bx,
           [=] AMREX_GPU_DEVICE (int i, int j, int k)
           {
               for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {
                   int il = i;
                   int jl = j;
                   int kl = k;

                   int ir = i;
                   int jr = j;
                   int kr = k;

                   if (idir == 0) {
                       il = i - 1;
                       ir = i + 1;
                   } else if (idir == 1) {
                       jl = j - 1;
                       jr = j + 1;
                   } else {
                       kl = k - 1;
                       kr = k + 1;
                   }

                   Real dxinv2 = dxinv[idir]*dxinv[idir];

                   const auto& edge_coeff_arr = edge_coeff_arrs[idir];
                   Real kth_r = edge_coeff_arr(ir,jr,kr);
                   Real kth_l = edge_coeff_arr(i ,j ,k);

#if AMREX_SPACEDIM < 3
                   // Apply geometric terms for curvilinear coordinates

                   if ((coord != 0) && (idir == 0)) {
                       // In curilinear radial direction

                       Real rr = problo[idir] + static_cast<Real>(ir) * dx[idir];
                       Real rl = problo[idir] + static_cast<Real>(i) * dx[idir];
                       Real rc = problo[idir] + (static_cast<Real>(i) + 0.5_rt) * dx[idir];

                       if (coord == 1) {
                           // Cylindrical radial equation looks like: 1/r d(r kth dT/dr)/dr

                           kth_r *= rr;
                           kth_l *= rl;
                           dxinv2 *= 1.0_rt / rc;
                       } else {
                           // Spherical radial equation looks like: 1/r^2 d(r^2 kth dT/dr)/dr

                           kth_r *= rr * rr;
                           kth_l *= rl * rl;
                           dxinv2 *= 1.0_rt / (rc * rc) ;
                       }
                   } else if ((coord == 2) && (idir == 1)) {
                       // In spherical theta direction
                       // Spherical theta equation looks like: 1/r^2 sin(θ) d(sin(θ) kth dT/dθ)/dθ

                       Real rc = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];
                       Real thetar = problo[idir] + static_cast<Real>(jr) * dx[idir];
                       Real thetal = problo[idir] + static_cast<Real>(j) * dx[idir];
                       Real thetac = problo[idir] + (static_cast<Real>(j) + 0.5_rt) * dx[idir];

                       kth_r *= std::sin(thetar);
                       kth_l *= std::sin(thetal);
                       dxinv2 *= 1.0_rt / (rc * rc * std::sin(thetac));
                   }
#endif
                   TempDiff_array(i,j,k) += dxinv2 *
                       (kth_r * (Temp_array(ir,jr,kr) - Temp_array(i ,j ,k )) -
                        kth_l * (Temp_array(i ,j ,k ) - Temp_array(il,jl,kl)));
               }
           });
       }
   }
}
