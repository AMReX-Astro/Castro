#include <Castro.H>
#include <Castro_F.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

using namespace amrex;

#ifdef TRUE_SDC
void
Castro::cons_to_prim(const Real time)
{

    BL_PROFILE("Castro::cons_to_prim()");
    
#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& qbx = mfi.growntilebox(NUM_GROW);

        Array4<Real> const q_arr = q.array(mfi);
        Array4<Real> const qaux_arr = qaux.array(mfi);
        Array4<Real> const Sborder_arr = Sborder.array(mfi);
#ifdef RADIATION
        Array4<Real> const Erborder_arr = Erborder.array(mfi);
        Array4<Real> const lamborder_arr = lamborder.array(mfi);
#endif

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.

        ctoprim(qbx,
                time,
                Sborder_arr,
#ifdef RADIATION
                Erborder_arr,
                lamborder_arr,
#endif
                q_arr,
                qaux_arr);

    }

}
#endif

// Convert a MultiFab with conservative state data u to a primitive MultiFab q.
void
Castro::cons_to_prim(MultiFab& u, MultiFab& q_in, MultiFab& qaux_in, Real time)
{

    BL_PROFILE("Castro::cons_to_prim()");

    BL_ASSERT(u.nComp() == NUM_STATE);
    BL_ASSERT(q_in.nComp() == NQ);
    BL_ASSERT(u.nGrow() >= q_in.nGrow());

    int ng = q_in.nGrow();

#ifdef RADIATION
    AmrLevel::FillPatch(*this, Erborder, NUM_GROW, time, Rad_Type, 0, Radiation::nGroups);

    MultiFab lamborder(grids, dmap, Radiation::nGroups, NUM_GROW);
    if (radiation->pure_hydro) {
      lamborder.setVal(0.0, NUM_GROW);
    }
    else {
      radiation->compute_limiter(level, grids, Sborder, Erborder, lamborder);
    }
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(u, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.growntilebox(ng);

        auto u_arr = u.array(mfi);
#ifdef RADIATION
        auto Erborder_arr = Erborder.array(mfi);
        auto lamborder_arr = lamborder.array(mfi);
#endif
        auto q_in_arr = q_in.array(mfi);
        auto qaux_in_arr = qaux_in.array(mfi);

        ctoprim(bx,
                time,
                u_arr,
#ifdef RADIATION
                Erborder_arr,
                lamborder_arr,
#endif
                q_in_arr,
                qaux_in_arr);

    }

}

#ifdef TRUE_SDC
void
Castro::cons_to_prim_fourth(const Real time)
{

    BL_PROFILE("Castro::cons_to_prim_fourth()");
    
    // convert the conservative state cell averages to primitive cell
    // averages with 4th order accuracy

    auto domain_lo = geom.Domain().loVect3d();
    auto domain_hi = geom.Domain().hiVect3d();

    MultiFab& S_new = get_new_data(State_Type);

    // we don't support radiation here
#ifdef RADIATION
    amrex::Abort("radiation not supported to fourth order");
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox U_cc;

        for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

            const Box& qbx = mfi.growntilebox(NUM_GROW);
            const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

            // note: these conversions are using a growntilebox, so it
            // will include ghost cells

            // convert U_avg to U_cc -- this will use a Laplacian
            // operation and will result in U_cc defined only on
            // NUM_GROW-1 ghost cells at the end.

            U_cc.resize(qbx, NUM_STATE);
            Elixir elix_u_cc = U_cc.elixir();
            auto const U_cc_arr = U_cc.array();

            make_cell_center(qbxm1, Sborder.array(mfi), U_cc_arr, domain_lo, domain_hi);

            // enforce the minimum density on the new cell-centered state
            do_enforce_minimum_density(qbxm1, U_cc.array(), verbose);

            // and ensure that the internal energy is positive
            reset_internal_energy(qbxm1, U_cc.array());

            // convert U_avg to q_bar -- this will be done on all NUM_GROW
            // ghost cells.
            auto Sborder_arr = Sborder.array(mfi);
            auto q_bar_arr = q_bar.array(mfi);
            auto qaux_bar_arr = qaux_bar.array(mfi);

            ctoprim(qbx,
                    time, 
                    Sborder_arr,
                    q_bar_arr,
                    qaux_bar_arr);

            // this is what we should construct the flattening coefficient
            // from

            // convert U_cc to q_cc (we'll store this temporarily in q,
            // qaux).  This will remain valid only on the NUM_GROW-1 ghost
            // cells.
            auto q_arr = q.array(mfi);
            auto qaux_arr = qaux.array(mfi);

            ctoprim(qbxm1,
                    time,
                    U_cc_arr,
                    q_arr,
                    qaux_arr);
        }
    }

    // check for NaNs
#ifndef AMREX_USE_GPU
    check_for_nan(q);
    check_for_nan(q_bar);
#endif

#ifdef DIFFUSION
    // we need the cell-center temperature for the diffusion stencil,
    // so save it here, by copying from q (which is cell-center at the
    // moment).
    MultiFab::Copy(T_cc, q, QTEMP, 0, 1, NUM_GROW-1);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // now convert q, qaux into 4th order accurate averages
      // this will create q, qaux in NUM_GROW-1 ghost cells, but that's
      // we need here

      make_fourth_average(qbxm1, q.array(mfi), q_bar.array(mfi), domain_lo, domain_hi);

      // not sure if we need to convert qaux this way, or if we can
      // just evaluate it (we may not need qaux at all actually)

      make_fourth_average(qbxm1, qaux.array(mfi), qaux_bar.array(mfi), domain_lo, domain_hi);

    }

#endif // RADIATION
}
#endif

void
Castro::check_for_cfl_violation(const MultiFab& State, const Real dt)
{

    BL_PROFILE("Castro::check_for_cfl_violation()");

    auto dx = geom.CellSizeArray();

    Real dtdx = dt / dx[0];

    Real dtdy = 0.0_rt;
    if (AMREX_SPACEDIM >= 2) {
      dtdy = dt / dx[1];
    }

    Real dtdz = 0.0_rt;
    if (AMREX_SPACEDIM == 3) {
      dtdz = dt / dx[2];
    }

    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

        auto U = State.array(mfi);

        reduce_op.eval(bx, reduce_data,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            // Compute running max of Courant number over grids

            Real rho = U(i,j,k,URHO);
            Real rhoInv = 1.0 / rho;

            Real u = U(i,j,k,UMX) * rhoInv;
            Real v = U(i,j,k,UMY) * rhoInv;
            Real w = U(i,j,k,UMZ) * rhoInv;

            eos_rep_t eos_state;

            eos_state.rho = U(i,j,k,URHO);
            eos_state.T = U(i,j,k,UTEMP);
            eos_state.e = U(i,j,k,UEINT) * rhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = U(i,j,k,UFS+n) * rhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = U(i,j,k,UFX+n) * rhoInv;
            }
#endif

            eos(eos_input_re, eos_state);

            Real cs = eos_state.cs;

            Real courx = (cs + std::abs(u)) * dtdx;
            Real coury = (cs + std::abs(v)) * dtdy;
            Real courz = (cs + std::abs(w)) * dtdz;

            if (castro::time_integration_method == 0) {

#ifndef AMREX_USE_GPU
                if (verbose == 1) {

                    if (courx > 1.0_rt) {
                        std::cout << std::endl;
                        std::cout << "Warning:: CFL violation in check_for_cfl_violation" << std::endl;
                        std::cout << ">>> ... (u+c) * dt / dx > 1 " << courx << std::endl;
                        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
                        std::cout << ">>> ... u = " << u << " c = " << cs << std::endl;
                        std::cout << ">>> ... density = " << rho << std::endl;
                    }

                    if (coury > 1.0_rt) {
                        std::cout << std::endl;
                        std::cout << "Warning:: CFL violation in check_for_cfl_violation" << std::endl;
                        std::cout << ">>> ... (v+c) * dt / dx > 1 " << coury << std::endl;
                        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
                        std::cout << ">>> ... v = " << v << " c = " << cs << std::endl;
                        std::cout << ">>> ... density = " << rho << std::endl;
                    }

                    if (courz > 1.0_rt) {
                        std::cout << std::endl;
                        std::cout << "Warning:: CFL violation in check_for_cfl_violation" << std::endl;
                        std::cout << ">>> ... (w+c) * dt / dx > 1 " << courz << std::endl;
                        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
                        std::cout << ">>> ... w = " << w << " c = " << cs << std::endl;
                        std::cout << ">>> ... density = " << rho << std::endl;
                    }

                }
#endif

                // CTU integration constraint

                return {amrex::max(courx, coury, courz)};

            }
            else {

                // method-of-lines constraint
                Real courtmp = courx;
                if (AMREX_SPACEDIM >= 2) {
                    courtmp += coury;
                }
                if (AMREX_SPACEDIM == 3) {
                    courtmp += courz;
                }

#ifndef AMREX_USE_GPU
                if (verbose == 1) {

                    // note: it might not be 1 for all RK integrators
                    if (courtmp > 1.0_rt) {
                        std::cout << std::endl;
                        std::cout << "Warning:: CFL violation in check_for_cfl_violation" << std::endl;
                        std::cout << ">>> ... at cell i = " << i << " j = " << j << " k = " << k << std::endl;
                        std::cout << ">>> ... u = " << u << " v = " << v
                                  << " w = " << w << " c = " << cs << std::endl;
                        std::cout << ">>> ... density = " << rho << std::endl;
                    }

                }
#endif

                return {courtmp};

            }

        });

    }

    ReduceTuple hv = reduce_data.value();
    Real courno = amrex::get<0>(hv);

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
        amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

        cfl_violation = 1;
    }

}
