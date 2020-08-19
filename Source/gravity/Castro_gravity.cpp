#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

#ifdef HYBRID_MOMENTUM
#include "Castro_util.H"
#include "hybrid.H"
#endif

using namespace amrex;

void
Castro::construct_old_gravity(int amr_iteration, int amr_ncycle, Real time)
{
    BL_PROFILE("Castro::construct_old_gravity()");

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& phi_old = get_old_data(PhiGrav_Type);

    // Always set phi to zero initially since some gravity modes
    // don't use it and we want to have valid data.

    if (gravity->get_gravity_type() != "PoissonGrav")
        phi_old.setVal(0.0);

    if (!do_grav) {

        grav_old.setVal(0.0);

        return;

    }

    // Do level solve at beginning of time step in order to compute the
    // difference between the multilevel and the single level solutions.
    // Note that we don't need to do this solve for single-level runs,
    // since the solution at the end of the last timestep won't have changed.

    if (gravity->get_gravity_type() == "PoissonGrav" && parent->finestLevel() > 0)
    {

        // Create a copy of the current (composite) data on this level.

        MultiFab comp_phi;
        Vector<std::unique_ptr<MultiFab> > comp_gphi(BL_SPACEDIM);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() && level < parent->finestLevel() && level <= gravity->get_max_solve_level()) {

            comp_phi.define(phi_old.boxArray(), phi_old.DistributionMap(), phi_old.nComp(), phi_old.nGrow());
            MultiFab::Copy(comp_phi, phi_old, 0, 0, phi_old.nComp(), phi_old.nGrow());

            for (int n = 0; n < BL_SPACEDIM; ++n) {
                comp_gphi[n].reset(new MultiFab(getEdgeBoxArray(n), dmap, 1, 0));
                comp_gphi[n]->copy(*gravity->get_grad_phi_prev(level)[n], 0, 0, 1);
            }

        }

        if (castro::verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "... old-time level Poisson gravity solve at level " << level << std::endl << std::endl;
        }

        int is_new = 0;

        // If we are doing composite solves, then this is a placeholder solve
        // to get the difference between the composite and level solutions. If
        // we are only doing level solves, then this is the main result.

        gravity->solve_for_phi(level,
                               phi_old,
                               amrex::GetVecOfPtrs(gravity->get_grad_phi_prev(level)),
                               is_new);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() && level < parent->finestLevel() && level <= gravity->get_max_solve_level()) {

            // Subtract the level solve from the composite solution.

            gravity->create_comp_minus_level_grad_phi(level,
                                                      comp_phi,
                                                      amrex::GetVecOfPtrs(comp_gphi),
                                                      comp_minus_level_phi,
                                                      comp_minus_level_grad_phi);

            // Copy the composite data back. This way the forcing
            // uses the most accurate data we have.

            MultiFab::Copy(phi_old, comp_phi, 0, 0, phi_old.nComp(), phi_old.nGrow());

            for (int n = 0; n < BL_SPACEDIM; ++n)
                gravity->get_grad_phi_prev(level)[n]->copy(*comp_gphi[n], 0, 0, 1);

        }

        if (gravity->test_results_of_solves() == 1) {

            if (castro::verbose && ParallelDescriptor::IOProcessor()) {
                std::cout << " " << '\n';
                std::cout << "... testing grad_phi_curr after doing single level solve " << '\n';
            }

            gravity->test_level_grad_phi_prev(level);

        }
 
    }

    // Define the old gravity vector.

    gravity->get_old_grav_vector(level, grav_old, time);

}

void
Castro::construct_new_gravity(int amr_iteration, int amr_ncycle, Real time)
{
    BL_PROFILE("Castro::construct_new_gravity()");

    MultiFab& grav_new = get_new_data(Gravity_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    // Always set phi to zero initially since some gravity modes
    // don't use it and we want to have valid data.

    if (gravity->get_gravity_type() != "PoissonGrav")
        phi_new.setVal(0.0);

    if (!do_grav) {

        grav_new.setVal(0.0);

        return;

    }

    // If we're doing Poisson gravity, do the new-time level solve here.

    if (gravity->get_gravity_type() == "PoissonGrav")
    {

        // Use the "old" phi from the current time step as a guess for this solve.

        MultiFab& phi_old = get_old_data(PhiGrav_Type);

        MultiFab::Copy(phi_new, phi_old, 0, 0, 1, phi_new.nGrow());

        // Subtract off the (composite - level) contribution for the purposes
        // of the level solve. We'll add it back later.

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() && level < parent->finestLevel() && level <= gravity->get_max_solve_level())
            phi_new.minus(comp_minus_level_phi, 0, 1, 0);

        if (castro::verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "... new-time level Poisson gravity solve at level " << level << std::endl << std::endl;
        }

        int is_new = 1;

        gravity->solve_for_phi(level,
                               phi_new,
                               amrex::GetVecOfPtrs(gravity->get_grad_phi_curr(level)),
                               is_new);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() == 1 && level < parent->finestLevel() && level <= gravity->get_max_solve_level()) {

            if (gravity->test_results_of_solves() == 1) {

                if (castro::verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << " " << '\n';
                    std::cout << "... testing grad_phi_curr before adding comp_minus_level_grad_phi " << '\n';
                }

                gravity->test_level_grad_phi_curr(level);

            }

            // Add back the (composite - level) contribution. This ensures that
            // if we are not doing a sync solve, then we still get the difference
            // between the composite and level solves added to the force we
            // calculate, so it is slightly more accurate than it would have been.

            phi_new.plus(comp_minus_level_phi, 0, 1, 0);
            for (int n = 0; n < BL_SPACEDIM; ++n)
                gravity->get_grad_phi_curr(level)[n]->plus(*comp_minus_level_grad_phi[n], 0, 1, 0);

            if (gravity->test_results_of_solves() == 1) {

                if (castro::verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << " " << '\n';
                    std::cout << "... testing grad_phi_curr after adding comp_minus_level_grad_phi " << '\n';
                }

                gravity->test_level_grad_phi_curr(level);

            }

        }

    }

    // Define new gravity vector.

    gravity->get_new_grav_vector(level, grav_new, time);

    if (gravity->get_gravity_type() == "PoissonGrav" && level <= gravity->get_max_solve_level()) {

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() == 1 && level < parent->finestLevel()) {

            // Now that we have calculated the force, if we are going to do a sync
            // solve then subtract off the (composite - level) contribution, as it
            // interferes with the sync solve.

            if (gravity->NoSync() == 0) {

                phi_new.minus(comp_minus_level_phi, 0, 1, 0);

                for (int n = 0; n < BL_SPACEDIM; ++n)
                    gravity->get_grad_phi_curr(level)[n]->minus(*comp_minus_level_grad_phi[n], 0, 1, 0);

            }

            // In any event we can now clear this memory, as we no longer need it.

            comp_minus_level_phi.clear();
            comp_minus_level_grad_phi.clear();

        }

    }

}

void Castro::construct_old_gravity_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_old_gravity_source()");

    const Real strt_time = ParallelDescriptor::second();

    const MultiFab& phi_old = get_old_data(PhiGrav_Type);
    const MultiFab& grav_old = get_old_data(Gravity_Type);

    if (!do_grav) return;

    // Gravitational source term for the time-level n data.

#ifdef HYBRID_MOMENTUM
    GeometryData geomdata = geom.data();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());
#endif

    AMREX_ALWAYS_ASSERT(castro::grav_source_type >= 1 && castro::grav_source_type <= 4);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Array4<Real const> const uold = state_in.array(mfi);
        Array4<Real const> const grav = grav_old.array(mfi);
        Array4<Real> const source_arr = source.array(mfi);

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
        {
            // Temporary array for seeing what the new state would be if the update were applied here.

            GpuArray<Real, NUM_STATE> snew;
            for (int n = 0; n < NUM_STATE; ++n) {
                snew[n] = 0.0_rt;
            }

            // Temporary array for holding the update to the state.

            GpuArray<Real, NSRC> src;
            for (int n = 0; n < NSRC; ++n) {
                src[n] = 0.0_rt;
            }

            // Gravitational source options for how to add the work to (rho E):
            // grav_source_type =
            // 1: Original version ("does work")
            // 2: Modification of type 1 that updates the momentum before constructing the energy corrector
            // 3: Puts all gravitational work into KE, not (rho e)
            // 4: Conservative energy formulation

            Real rho    = uold(i,j,k,URHO);
            Real rhoInv = 1.0_rt / rho;

            for (int n = 0; n < NUM_STATE; ++n) {
                snew[n] = uold(i,j,k,n);
            }

            Real old_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoInv;

            GpuArray<Real, 3> Sr;
            for (int n = 0; n < 3; ++n) {
                Sr[n] = rho * grav(i,j,k,n);

                src[UMX+n] = Sr[n];

                snew[UMX+n] += dt * src[UMX+n];
            }

#ifdef HYBRID_MOMENTUM
            GpuArray<Real, 3> loc;
            for (int n = 0; n < 3; ++n) {
                position(i, j, k, geomdata, loc);
                loc[n] -= center[n];
            }

            GpuArray<Real, 3> hybrid_src;

            set_hybrid_momentum_source(loc, Sr, hybrid_src);

            for (int n = 0; n < 3; ++n) {
                 src[UMR+n] = hybrid_src[n];
                 snew[UMR+n] += dt * src[UMR+n];
            }
#endif

            Real SrE;

            if (castro::grav_source_type == 1 || castro::grav_source_type == 2) {

                // Src = rho u dot g, evaluated with all quantities at t^n

                SrE = (uold(i,j,k,UMX) * Sr[0] + uold(i,j,k,UMY) * Sr[1] + uold(i,j,k,UMZ) * Sr[2]) * rhoInv;

            } else if (castro::grav_source_type == 3) {

                Real new_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoInv;
                SrE = new_ke - old_ke;

            } else if (castro::grav_source_type == 4) {

                // The conservative energy formulation does not strictly require
                // any energy source-term here, because it depends only on the
                // fluid motions from the hydrodynamical fluxes which we will only
                // have when we get to the 'corrector' step. Nevertheless we add a
                // predictor energy source term in the way that the other methods
                // do, for consistency. We will fully subtract this predictor value
                // during the corrector step, so that the final result is correct.
                // Here we use the same approach as grav_source_type == 2.

                SrE = (uold(i,j,k,UMX) * Sr[0] + uold(i,j,k,UMY) * Sr[1] + uold(i,j,k,UMZ) * Sr[2]) * rhoInv;

            }

            src[UEDEN] = SrE;

            snew[UEDEN] += dt * SrE;

            // Add to the outgoing source array.

            for (int n = 0; n < NSRC; ++n) {
                source_arr(i,j,k,n) += src[n];
            }

        });

    }

    if (castro::verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_old_gravity_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

}

void Castro::construct_new_gravity_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
    BL_PROFILE("Castro::construct_new_gravity_source()");

    const Real strt_time = ParallelDescriptor::second();

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);

    if (!do_grav) return;

    GpuArray<Real, 3> dx;
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        dx[i] = geom.CellSizeArray()[i];
    }
    for (int i = AMREX_SPACEDIM; i < 3; ++i) {
        dx[i] = 0.0_rt;
    }

#ifdef HYBRID_MOMENTUM
    GeometryData geomdata = geom.data();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());
#endif

    AMREX_ALWAYS_ASSERT(castro::grav_source_type >= 1 && castro::grav_source_type <= 4);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            Array4<Real const> const uold  = state_old.array(mfi);
            Array4<Real const> const unew  = state_new.array(mfi);
            Array4<Real const> const gold  = grav_old.array(mfi);
            Array4<Real const> const gnew  = grav_new.array(mfi);
            Array4<Real const> const vol   = volume.array(mfi);
            Array4<Real const> const flux0 = (*mass_fluxes[0]).array(mfi);
            Array4<Real const> const flux1 = (*mass_fluxes[1]).array(mfi);
            Array4<Real const> const flux2 = (*mass_fluxes[2]).array(mfi);
            Array4<Real> const source_arr  = source.array(mfi);

            amrex::ParallelFor(bx,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
            {
                GpuArray<Real, NSRC> src{};

                Real hdtInv = 0.5_rt / dt;

                // Gravitational source options for how to add the work to (rho E):
                // grav_source_type =
                // 1: Original version ("does work")
                // 2: Modification of type 1 that updates the U before constructing SrEcorr
                // 3: Puts all gravitational work into KE, not (rho e)
                // 4: Conservative gravity approach (discussed in first white dwarf merger paper).

                Real rhoo    = uold(i,j,k,URHO);
                Real rhooinv = 1.0_rt / uold(i,j,k,URHO);

                Real rhon    = unew(i,j,k,URHO);
                Real rhoninv = 1.0_rt / unew(i,j,k,URHO);

                // Temporary array for seeing what the new state would be if the update were applied here.

                GpuArray<Real, NUM_STATE> snew{};
                for (int n = 0; n < NUM_STATE; ++n) {
                    snew[n] = unew(i,j,k,n);
                }

                Real old_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoninv;

                // Define old source terms

                GpuArray<Real, 3> vold;
                for (int n = 0; n < 3; ++n) {
                    vold[n] = uold(i,j,k,UMX+n) * rhooinv;
                }

                GpuArray<Real, 3> Sr_old;
                for (int n = 0; n < 3; ++n) {
                    Sr_old[n] = rhoo * gold(i,j,k,n);
                }

                Real SrE_old = vold[0] * Sr_old[0] + vold[1] * Sr_old[1] + vold[2] * Sr_old[2];

                // Define new source terms

                GpuArray<Real, 3> vnew;
                for (int n = 0; n < 3; ++n) {
                    vnew[n] = snew[UMX+n] * rhoninv;
                }

                GpuArray<Real, 3> Sr_new;
                for (int n = 0; n < 3; ++n) {
                    Sr_new[n] = rhon * gnew(i,j,k,n);
                }

                Real SrE_new = vnew[0] * Sr_new[0] + vnew[1] * Sr_new[1] + vnew[2] * Sr_new[2];

                // Define corrections to source terms

                GpuArray<Real, 3> Srcorr;
                for (int n = 0; n < 3; ++n) {
                    Srcorr[n] = 0.5_rt * (Sr_new[n] - Sr_old[n]);
                }

                // Correct momenta

                for (int n = 0; n < 3; ++n) {
                    src[UMX+n] = Srcorr[n];
                    snew[UMX+n] += dt * src[UMX+n];
                }

#ifdef HYBRID_MOMENTUM
                GpuArray<Real, 3> loc;
                position(i, j, k, geomdata, loc);
                for (int n = 0; n < 3; ++n) {
                    loc[n] -= center[n];
                }

                GpuArray<Real, 3> hybrid_src;

                set_hybrid_momentum_source(loc, Srcorr, hybrid_src);

                for (int n = 0; n < 3; ++n) {
                    src[UMR+n] = hybrid_src[n];
                    snew[UMR+n] += dt * src[UMR+n];
                }
#endif

                // Correct energy

                Real SrEcorr;

                if (castro::grav_source_type == 1) {

                    // If grav_source_type == 1, then we calculated SrEcorr before updating the velocities.

                    SrEcorr = 0.5_rt * (SrE_new - SrE_old);

                } else if (castro::grav_source_type == 2) {

                    // For this source type, we first update the momenta
                    // before we calculate the energy source term.

                    for (int n = 0; n < 3; ++n) {
                        vnew[n] = snew[UMX+n] * rhoninv;
                    }
                    SrE_new = vnew[0] * Sr_new[0] + vnew[1] * Sr_new[1] + vnew[2] * Sr_new[2];

                    SrEcorr = 0.5_rt * (SrE_new - SrE_old);

                } else if (castro::grav_source_type == 3) {

                    // Instead of calculating the energy source term explicitly,
                    // we simply update the kinetic energy.

                    Real new_ke = 0.5_rt * (snew[UMX] * snew[UMX] + snew[UMY] * snew[UMY] + snew[UMZ] * snew[UMZ]) * rhoninv;
                    SrEcorr = new_ke - old_ke;

                } else if (castro::grav_source_type == 4) {

                    // First, subtract the predictor step we applied earlier.

                    SrEcorr = - SrE_old;

                    // For an explanation of this approach, see wdmerger paper I.
                    // The main idea is that we are evaluating the change of the
                    // potential energy at zone edges and applying that in an equal
                    // and opposite sense to the gas energy. The physics is described
                    // in Section 2.4; we are using a version of the formula similar to
                    // Equation 94 in Springel (2010) based on the gradient rather than
                    // the potential because the gradient-version works for all forms
                    // of gravity we use, some of which do not explicitly calculate phi.

                    // Construct the time-averaged edge-centered gravity.

                    GpuArray<Real, 3> g;
                    for (int n = 0; n < 3; ++n) {
                        g[n] = 0.5_rt * (gnew(i,j,k,n) + gold(i,j,k,n));
                    }

                    Real gxl = 0.5_rt * (g[0] + 0.5_rt * (gnew(i-1*dg0,j,k,0) + gold(i-1*dg0,j,k,0)));
                    Real gxr = 0.5_rt * (g[0] + 0.5_rt * (gnew(i+1*dg0,j,k,0) + gold(i+1*dg0,j,k,0)));

                    Real gyl = 0.5_rt * (g[1] + 0.5_rt * (gnew(i,j-1*dg1,k,1) + gold(i,j-1*dg1,k,1)));
                    Real gyr = 0.5_rt * (g[1] + 0.5_rt * (gnew(i,j+1*dg1,k,1) + gold(i,j+1*dg1,k,1)));

                    Real gzl = 0.5_rt * (g[2] + 0.5_rt * (gnew(i,j,k-1*dg2,2) + gold(i,j,k-1*dg2,2)));
                    Real gzr = 0.5_rt * (g[2] + 0.5_rt * (gnew(i,j,k+1*dg2,2) + gold(i,j,k+1*dg2,2)));

                    SrEcorr += hdtInv * (flux0(i      ,j,k) * gxl * dx[0] +
                                         flux0(i+1*dg0,j,k) * gxr * dx[0] +
                                         flux1(i,j      ,k) * gyl * dx[1] +
                                         flux1(i,j+1*dg1,k) * gyr * dx[1] +
                                         flux2(i,j,k      ) * gzl * dx[2] +
                                         flux2(i,j,k+1*dg2) * gzr * dx[2]) / vol(i,j,k);

                }

                src[UEDEN] = SrEcorr;

                snew[UEDEN] += dt * SrEcorr;

                // Add to the outgoing source array.

                for (int n = 0; n < NSRC; ++n) {
                    source_arr(i,j,k,n) += src[n];
                }
            });
        }
    }

    if (castro::verbose > 1)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Castro::construct_new_gravity_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}
