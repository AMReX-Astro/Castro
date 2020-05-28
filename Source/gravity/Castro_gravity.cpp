#include "Castro.H"
#include "Castro_F.H"

#include "Gravity.H"

using namespace amrex;

void Castro::construct_old_gravity(int amr_iteration, int amr_ncycle,
                                   Real time) {
    BL_PROFILE("Castro::construct_old_gravity()");

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& phi_old = get_old_data(PhiGrav_Type);

    // Always set phi to zero initially since some gravity modes
    // don't use it and we want to have valid data.

    if (gravity->get_gravity_type() != "PoissonGrav") phi_old.setVal(0.0);

    if (!do_grav) {

        grav_old.setVal(0.0);

        return;
    }

    // Do level solve at beginning of time step in order to compute the
    // difference between the multilevel and the single level solutions.
    // Note that we don't need to do this solve for single-level runs,
    // since the solution at the end of the last timestep won't have changed.

    if (gravity->get_gravity_type() == "PoissonGrav" &&
        parent->finestLevel() > 0) {

        // Create a copy of the current (composite) data on this level.

        MultiFab comp_phi;
        Vector<std::unique_ptr<MultiFab> > comp_gphi(BL_SPACEDIM);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() &&
            level < parent->finestLevel() &&
            level <= gravity->get_max_solve_level()) {

            comp_phi.define(phi_old.boxArray(), phi_old.DistributionMap(),
                            phi_old.nComp(), phi_old.nGrow());
            MultiFab::Copy(comp_phi, phi_old, 0, 0, phi_old.nComp(),
                           phi_old.nGrow());

            for (int n = 0; n < BL_SPACEDIM; ++n) {
                comp_gphi[n].reset(
                    new MultiFab(getEdgeBoxArray(n), dmap, 1, 0));
                comp_gphi[n]->copy(*gravity->get_grad_phi_prev(level)[n], 0, 0,
                                   1);
            }
        }

        if (castro::verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "... old-time level Poisson gravity solve at level "
                      << level << std::endl
                      << std::endl;
        }

        int is_new = 0;

        // If we are doing composite solves, then this is a placeholder solve
        // to get the difference between the composite and level solutions. If
        // we are only doing level solves, then this is the main result.

        gravity->solve_for_phi(
            level, phi_old,
            amrex::GetVecOfPtrs(gravity->get_grad_phi_prev(level)), is_new);

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() &&
            level < parent->finestLevel() &&
            level <= gravity->get_max_solve_level()) {

            // Subtract the level solve from the composite solution.

            gravity->create_comp_minus_level_grad_phi(
                level, comp_phi, amrex::GetVecOfPtrs(comp_gphi),
                comp_minus_level_phi, comp_minus_level_grad_phi);

            // Copy the composite data back. This way the forcing
            // uses the most accurate data we have.

            MultiFab::Copy(phi_old, comp_phi, 0, 0, phi_old.nComp(),
                           phi_old.nGrow());

            for (int n = 0; n < BL_SPACEDIM; ++n)
                gravity->get_grad_phi_prev(level)[n]->copy(*comp_gphi[n], 0, 0,
                                                           1);
        }

        if (gravity->test_results_of_solves() == 1) {

            if (castro::verbose && ParallelDescriptor::IOProcessor()) {
                std::cout << " " << '\n';
                std::cout << "... testing grad_phi_curr after doing single "
                             "level solve "
                          << '\n';
            }

            gravity->test_level_grad_phi_prev(level);
        }
    }

    // Define the old gravity vector.

    gravity->get_old_grav_vector(level, grav_old, time);
}

void Castro::construct_new_gravity(int amr_iteration, int amr_ncycle,
                                   Real time) {
    BL_PROFILE("Castro::construct_new_gravity()");

    MultiFab& grav_new = get_new_data(Gravity_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    // Always set phi to zero initially since some gravity modes
    // don't use it and we want to have valid data.

    if (gravity->get_gravity_type() != "PoissonGrav") phi_new.setVal(0.0);

    if (!do_grav) {

        grav_new.setVal(0.0);

        return;
    }

    // If we're doing Poisson gravity, do the new-time level solve here.

    if (gravity->get_gravity_type() == "PoissonGrav") {

        // Use the "old" phi from the current time step as a guess for this solve.

        MultiFab& phi_old = get_old_data(PhiGrav_Type);

        MultiFab::Copy(phi_new, phi_old, 0, 0, 1, phi_new.nGrow());

        // Subtract off the (composite - level) contribution for the purposes
        // of the level solve. We'll add it back later.

        if (gravity->NoComposite() != 1 && gravity->DoCompositeCorrection() &&
            level < parent->finestLevel() &&
            level <= gravity->get_max_solve_level())
            phi_new.minus(comp_minus_level_phi, 0, 1, 0);

        if (castro::verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "... new-time level Poisson gravity solve at level "
                      << level << std::endl
                      << std::endl;
        }

        int is_new = 1;

        gravity->solve_for_phi(
            level, phi_new,
            amrex::GetVecOfPtrs(gravity->get_grad_phi_curr(level)), is_new);

        if (gravity->NoComposite() != 1 &&
            gravity->DoCompositeCorrection() == 1 &&
            level < parent->finestLevel() &&
            level <= gravity->get_max_solve_level()) {

            if (gravity->test_results_of_solves() == 1) {

                if (castro::verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << " " << '\n';
                    std::cout << "... testing grad_phi_curr before adding "
                                 "comp_minus_level_grad_phi "
                              << '\n';
                }

                gravity->test_level_grad_phi_curr(level);
            }

            // Add back the (composite - level) contribution. This ensures that
            // if we are not doing a sync solve, then we still get the difference
            // between the composite and level solves added to the force we
            // calculate, so it is slightly more accurate than it would have been.

            phi_new.plus(comp_minus_level_phi, 0, 1, 0);
            for (int n = 0; n < BL_SPACEDIM; ++n)
                gravity->get_grad_phi_curr(level)[n]->plus(
                    *comp_minus_level_grad_phi[n], 0, 1, 0);

            if (gravity->test_results_of_solves() == 1) {

                if (castro::verbose && ParallelDescriptor::IOProcessor()) {
                    std::cout << " " << '\n';
                    std::cout << "... testing grad_phi_curr after adding "
                                 "comp_minus_level_grad_phi "
                              << '\n';
                }

                gravity->test_level_grad_phi_curr(level);
            }
        }
    }

    // Define new gravity vector.

    gravity->get_new_grav_vector(level, grav_new, time);

    if (gravity->get_gravity_type() == "PoissonGrav" &&
        level <= gravity->get_max_solve_level()) {

        if (gravity->NoComposite() != 1 &&
            gravity->DoCompositeCorrection() == 1 &&
            level < parent->finestLevel()) {

            // Now that we have calculated the force, if we are going to do a sync
            // solve then subtract off the (composite - level) contribution, as it
            // interferes with the sync solve.

            if (gravity->NoSync() == 0) {

                phi_new.minus(comp_minus_level_phi, 0, 1, 0);

                for (int n = 0; n < BL_SPACEDIM; ++n)
                    gravity->get_grad_phi_curr(level)[n]->minus(
                        *comp_minus_level_grad_phi[n], 0, 1, 0);
            }

            // In any event we can now clear this memory, as we no longer need it.

            comp_minus_level_phi.clear();
            comp_minus_level_grad_phi.clear();
        }
    }
}

void Castro::construct_old_gravity_source(MultiFab& source, MultiFab& state_in,
                                          Real time, Real dt) {
    BL_PROFILE("Castro::construct_old_gravity_source()");

    const Real strt_time = ParallelDescriptor::second();

    const MultiFab& phi_old = get_old_data(PhiGrav_Type);
    const MultiFab& grav_old = get_old_data(Gravity_Type);

    if (!do_grav) return;

    // Gravitational source term for the time-level n data.

    const Real* dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_in, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
        ca_gsrc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                BL_TO_FORTRAN_ANYD(state_in[mfi]),
                BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                BL_TO_FORTRAN_ANYD(grav_old[mfi]),
                BL_TO_FORTRAN_ANYD(source[mfi]), AMREX_REAL_ANYD(dx), dt, time);
    }

    if (castro::verbose > 1) {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealMax(run_time, IOProc);

            if (ParallelDescriptor::IOProcessor())
                std::cout << "Castro::construct_old_gravity_source() time = "
                          << run_time << "\n"
                          << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}

void Castro::construct_new_gravity_source(MultiFab& source, MultiFab& state_old,
                                          MultiFab& state_new, Real time,
                                          Real dt) {
    BL_PROFILE("Castro::construct_new_gravity_source()");

    const Real strt_time = ParallelDescriptor::second();

    MultiFab& phi_old = get_old_data(PhiGrav_Type);
    MultiFab& phi_new = get_new_data(PhiGrav_Type);

    MultiFab& grav_old = get_old_data(Gravity_Type);
    MultiFab& grav_new = get_new_data(Gravity_Type);

    if (!do_grav) return;

    const Real* dx = geom.CellSize();
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

#pragma gpu box(bx)
            ca_corrgsrc(
                AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
                AMREX_INT_ANYD(domlo), AMREX_INT_ANYD(domhi),
                BL_TO_FORTRAN_ANYD(state_old[mfi]),
                BL_TO_FORTRAN_ANYD(state_new[mfi]),
                BL_TO_FORTRAN_ANYD(phi_old[mfi]),
                BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                BL_TO_FORTRAN_ANYD(grav_old[mfi]),
                BL_TO_FORTRAN_ANYD(grav_new[mfi]),
                BL_TO_FORTRAN_ANYD(volume[mfi]),
                BL_TO_FORTRAN_ANYD((*mass_fluxes[0])[mfi]),
                BL_TO_FORTRAN_ANYD((*mass_fluxes[1])[mfi]),
                BL_TO_FORTRAN_ANYD((*mass_fluxes[2])[mfi]),
                BL_TO_FORTRAN_ANYD(source[mfi]), AMREX_REAL_ANYD(dx), dt, time);
        }
    }

    if (castro::verbose > 1) {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction([=]() mutable {
#endif
            ParallelDescriptor::ReduceRealMax(run_time, IOProc);

            if (ParallelDescriptor::IOProcessor())
                std::cout << "Castro::construct_new_gravity_source() time = "
                          << run_time << "\n"
                          << "\n";
#ifdef BL_LAZY
        });
#endif
    }
}
