#include "Castro.H"
#include "Castro_F.H"
#include "fundamental_constants.H"
#include "Gravity.H"

using namespace amrex;

#ifdef GRAVITY
#ifdef ROTATION
void Castro::scf_relaxation() {

    AMREX_ASSERT(level == 0);

    // First do some sanity checks.

    // Maximum density must be set.

    if (scf_maximum_density <= 0) {
        amrex::Error("castro.scf_maximum_density must be set for SCF relaxation.");
    }

    // Equatorial radius and polar radius must both
    // be positive, and the polar radius cannot be
    // larger than the equatorial radius.

    if (scf_equatorial_radius <= 0) {
        amrex::Error("Equatorial radius must be positive for SCF relaxation.");
    }

    if (scf_polar_radius <= 0) {
        amrex::Error("Polar radius must be positive for SCF relaxation.");
    }

    if (scf_polar_radius > scf_equatorial_radius) {
        amrex::Error("Polar radius must not be larger than equatorial radius for SCF relaxation.");
    }

    // Now do the SCF solve.

    do_hscf_solve();

    // Update the gravitational field. Since we don't need this
    // during the iterations, we only do it once at the end.

    for (int lev = 0; lev <= parent->finestLevel(); ++lev)
    {
        MultiFab& grav_new = getLevel(lev).get_new_data(Gravity_Type);
        gravity->get_new_grav_vector(lev, grav_new, getLevel(lev).state[Gravity_Type].curTime());
    }

    // Perform a final pass to ensure we're consistent with AMR.

    for (int lev = parent->finestLevel()-1; lev >= 0; --lev) {
        getLevel(lev).avgDown();
    }

}

void
Castro::do_hscf_solve()
{

    const int finest_level = parent->finestLevel();
    const int n_levs = finest_level + 1;

    // Do the initial relaxation setup.

    scf_setup_relaxation(scf_equatorial_radius, 1, scf_polar_radius, 3);

    // To do the relaxation loop, we need to know the target maximum
    // enthalpy on the domain. So we'll loop through the grid and
    // calculate what it should be, based on the density, temperature,
    // and composition. We do not specify the latter two -- we assume
    // that these are filled as desired in initdata. So, to get the
    // maximum enthalpy, we'll loop through the grid, using the composition
    // and temperature in each zone, combined with the requested maximum
    // density, call the EOS to calculate the enthalpy for that combination,
    // and then take the maximum of that. (Note that the current density on
    // the grid is treated as an initial guess and does not need to reach
    // the requested maximum.) The SCF procedure is only intended for a setup
    // with uniform temperature and composition, so in that case this procedure
    // should yield the same result as just picking a zone at random and using
    // its temperature and composition. The reason for doing it this way is so
    // that the user doesn't have to request temperature and composition in their
    // inputs file (which is messy in general) -- they can tell us what composition
    // they want by filling it in the initialization.

    Real target_h_max = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev) {

        MultiFab& state_new = getLevel(lev).get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(max:target_h_max)
#endif
        for (MFIter mfi(state_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();

            scf_calculate_target_h_max(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                       BL_TO_FORTRAN_ANYD(state_new[mfi]),
                                       scf_maximum_density,
                                       &target_h_max);

        }

    }

    ParallelDescriptor::ReduceRealMax(target_h_max);

    Vector< std::unique_ptr<MultiFab> > psi(n_levs);
    Vector< std::unique_ptr<MultiFab> > enthalpy(n_levs);
    Vector< std::unique_ptr<MultiFab> > state_vec(n_levs);
    Vector< std::unique_ptr<MultiFab> > phi(n_levs);
    Vector< std::unique_ptr<MultiFab> > phi_rot(n_levs);

    // Iterate until the system is relaxed by filling the level data
    // and then doing a multilevel gravity solve.

    int j = 1;

    while (j <= scf_max_iterations) {

        Real time = getLevel(0).state[State_Type].curTime();

        // Construct a local MultiFab for the rotational psi.
        // This does not change over the loop iterations, but
        // this (and the data to follow) must be constructed
        // inside the loop on each iteration because of regrids.

        for (int lev = 0; lev <= finest_level; ++lev) {

            psi[lev].reset(new MultiFab(getLevel(lev).grids, getLevel(lev).dmap, 1, 0));

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*psi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                fill_rotational_psi(bx, (*psi[lev]).array(mfi), time);

            }

        }

        // Construct a local MultiFab for the enthalpy.

        for (int lev = 0; lev <= finest_level; ++lev) {
            enthalpy[lev].reset(new MultiFab(getLevel(lev).grids, getLevel(lev).dmap, 1, 0));
        }

        // Construct a local MultiFab for the state data. We do
        // this because we have to mask out the state several times
        // in the below calculation, and it's easiest to have a scratch
        // copy of the data to work with.

        for (int lev = 0; lev <= finest_level; ++lev) {
            state_vec[lev].reset(new MultiFab(getLevel(lev).grids, getLevel(lev).dmap, NUM_STATE, 0));
            phi[lev].reset(new MultiFab(getLevel(lev).grids, getLevel(lev).dmap, 1, 0));
            phi_rot[lev].reset(new MultiFab(getLevel(lev).grids, getLevel(lev).dmap, 1, 0));
        }

        // Copy in the state data. Mask it out on coarse levels.

        for (int lev = 0; lev <= finest_level; ++lev) {

            MultiFab::Copy((*state_vec[lev]), getLevel(lev).get_new_data(State_Type), 0, 0, NUM_STATE, 0);
            MultiFab::Copy((*phi[lev]), getLevel(lev).get_new_data(PhiGrav_Type), 0, 0, 1, 0);
            MultiFab::Copy((*phi_rot[lev]), getLevel(lev).get_new_data(PhiRot_Type), 0, 0, 1, 0);

            if (lev < finest_level) {
                const MultiFab& mask = getLevel(lev+1).build_fine_mask();

                for (int n = 0; n < NUM_STATE; ++n) {
                    MultiFab::Multiply((*state_vec[lev]), mask, 0, n, 1, 0);
                }

                MultiFab::Multiply((*phi[lev]), mask, 0, 0, 1, 0);
                MultiFab::Multiply((*phi_rot[lev]), mask, 0, 0, 1, 0);
            }

        }

        // First step is to find the rotational frequency.

        Real phi_A = 0.0;
        Real psi_A = 0.0;
        Real phi_B = 0.0;
        Real psi_B = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:phi_A, psi_A, phi_B, psi_B)
#endif
            for (MFIter mfi((*phi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                scf_update_for_omegasq(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                       BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
                                       BL_TO_FORTRAN_ANYD((*psi[lev])[mfi]),
                                       AMREX_ZFILL(dx),
                                       &phi_A, &psi_A, &phi_B, &psi_B);

            }

        }

        ParallelDescriptor::ReduceRealSum(phi_A);
        ParallelDescriptor::ReduceRealSum(psi_A);
        ParallelDescriptor::ReduceRealSum(phi_B);
        ParallelDescriptor::ReduceRealSum(psi_B);

        // Now update the square of the rotation frequency, following Hachisu (Equation 16).
        // Deal carefully with the special case where phi_A and phi_B are equal -- we assume
        // this is the non-rotating case, and we want to handle that gracefully. Since the
        // rotation module knows to treat rotational_period = 0 as effectively disabling
        // rotation, we'll use that approach.

        if (std::abs(scf_polar_radius - scf_equatorial_radius) / std::abs(scf_equatorial_radius) < 1.e-6) {

            rotational_period = 0.0;

        } else {

            Real omegasq = -(phi_A - phi_B) / (psi_A - psi_B);

            if (omegasq < 0.0 && ParallelDescriptor::IOProcessor()) {
                std::cerr << "Omega squared = " << omegasq << " is negative in the relaxation step; aborting." << std::endl;
                amrex::Error();
            }

            Real omega = sqrt(omegasq);

            // Rotational period is 2 pi / omega.

            rotational_period = 2.0 * M_PI / omega;

        }

        // Now save the updated rotational frequency in the Fortran module.

        set_rot_period(&rotational_period);

        // With the updated period, we can construct the updated rotational
        // potential, which will be used in the remaining steps below.

        for (int lev = 0; lev <= finest_level; ++lev) {

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*phi_rot[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                fill_rotational_potential(bx, (*phi_rot[lev]).array(mfi), time);

            }

        }



        // Second step is to evaluate the Bernoulli constant.

        Real bernoulli = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:bernoulli)
#endif
            for (MFIter mfi((*phi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                scf_get_bernoulli_const(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                        BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
                                        BL_TO_FORTRAN_ANYD((*phi_rot[lev])[mfi]),
                                        AMREX_ZFILL(dx), &bernoulli);

            }

        }

        ParallelDescriptor::ReduceRealSum(bernoulli);



        // Third step is to construct the enthalpy field and
        // find the maximum enthalpy for the star.

        for (int lev = 0; lev <= finest_level; ++lev) {

            enthalpy[lev]->setVal(0.0);

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi((*phi[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                scf_construct_enthalpy(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                       BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
                                       BL_TO_FORTRAN_ANYD((*phi_rot[lev])[mfi]),
                                       BL_TO_FORTRAN_ANYD((*enthalpy[lev])[mfi]),
                                       AMREX_ZFILL(dx), bernoulli);

            }

        }

        Real actual_h_max = 0.0;
        Real actual_rho_max = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            actual_h_max = std::max(actual_h_max, enthalpy[lev]->max(0));
            actual_rho_max = std::max(actual_rho_max, state_vec[lev]->max(URHO));

        }

        Real Linf_norm = 0.0;

        // Finally, update the density using the enthalpy field.

        for (int lev = 0; lev <= finest_level; ++lev) {

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(max:Linf_norm)
#endif
            for (MFIter mfi((*state_vec[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                scf_update_density(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                   BL_TO_FORTRAN_ANYD((*state_vec[lev])[mfi]),
                                   BL_TO_FORTRAN_ANYD((*enthalpy[lev])[mfi]),
                                   AMREX_ZFILL(dx), actual_rho_max,
                                   actual_h_max, target_h_max,
                                   &Linf_norm);

            }

        }

        ParallelDescriptor::ReduceRealMax(Linf_norm);

        // Copy state data back to its source, and synchronize it on coarser levels.

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(getLevel(lev).get_new_data(State_Type), (*state_vec[lev]), 0, 0, NUM_STATE, 0);
            MultiFab::Copy(getLevel(lev).get_new_data(PhiGrav_Type), (*phi[lev]), 0, 0, 1, 0);
            MultiFab::Copy(getLevel(lev).get_new_data(PhiRot_Type), (*phi_rot[lev]), 0, 0, 1, 0);
        }

        for (int lev = finest_level-1; lev >= 0; --lev) {
            getLevel(lev).avgDown();
        }

        // Since we've changed the density distribution on the grid, regrid.

        bool do_io = false;
        parent->RegridOnly(time, do_io);

        // Update the gravitational field -- only after we've completed cleaning up the state above.

        gravity->multilevel_solve_for_new_phi(0, finest_level);

        // Update diagnostic quantities.

        Real kin_eng = 0.0;
        Real pot_eng = 0.0;
        Real int_eng = 0.0;
        Real mass = 0.0;

        for (int lev = 0; lev <= finest_level; ++lev) {

            const Real* dx = parent->Geom(lev).CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:kin_eng,pot_eng,int_eng,mass)
#endif
            for (MFIter mfi((*state_vec[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.tilebox();

                scf_diagnostics(AMREX_ARLIM_ANYD(bx.loVect()), AMREX_ARLIM_ANYD(bx.hiVect()),
                                BL_TO_FORTRAN_ANYD((*state_vec[lev])[mfi]),
                                BL_TO_FORTRAN_ANYD((*phi[lev])[mfi]),
                                BL_TO_FORTRAN_ANYD((*phi_rot[lev])[mfi]),
                                AMREX_ZFILL(dx),
                                &kin_eng, &pot_eng, &int_eng, &mass);

            }

        }

        ParallelDescriptor::ReduceRealSum(kin_eng);
        ParallelDescriptor::ReduceRealSum(pot_eng);
        ParallelDescriptor::ReduceRealSum(int_eng);
        ParallelDescriptor::ReduceRealSum(mass);

        Real virial_error = std::abs(2.0 * kin_eng + pot_eng + 3.0 * int_eng) / std::abs(pot_eng);

        // Now check to see if we're converged.

        int is_relaxed = 0;

        if (Linf_norm < scf_relax_tol) {
            is_relaxed = 1;
        }

        if (ParallelDescriptor::IOProcessor()) {

            // Grab the value for the solar mass.

            std::cout << std::endl << std::endl;
            std::cout << "   Relaxation iterations completed: " << j << std::endl;
            std::cout << "   L-infinity norm of residual (relative to old state): " << Linf_norm << std::endl;
            std::cout << "   Rotational period (s): " << rotational_period << std::endl;
            std::cout << "   Kinetic energy: " << kin_eng << std::endl;
            std::cout << "   Potential energy: " << pot_eng << std::endl;
            std::cout << "   Internal energy: " << int_eng << std::endl;
            std::cout << "   Virial error: " << virial_error << std::endl;
            std::cout << "   Mass: " << mass / C::M_solar << " solar masses" << std::endl;

            if (is_relaxed == 1) {
                std::cout << "  Relaxation completed!" << std::endl;
            }

            std::cout << std::endl << std::endl;

        }

        if (is_relaxed == 1) break;

        j++;

    }

}
#endif
#endif
