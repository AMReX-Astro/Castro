#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

using namespace amrex;

void
Castro::do_sdc_update(int m_start, int m_end, Real dt)
{

    BL_PROFILE("Castro::do_sdc_update()");

    // NOTE: dt here is the full dt not the dt between time nodes

    // this routine needs to do the update from time node m to m+1
    //
    // We come in with:
    //   A_new[m_start] : this is the advective update at node m_start
    //   A_old[:] : this is the advective source for all nodes at the old iterate
    //
    //   R_old[:] : this is the reaction source for all nodes at the old iterate

    // If we do advection only, then the update is explicit.  If we do
    // reactions, then the update is implicit within a zone.

    // for 4th order reactive SDC, we need to first compute the source, C
    // and do a ghost cell fill on it

    const int* domain_lo = geom.Domain().loVect();
    const int* domain_hi = geom.Domain().hiVect();

    // the timestep from m to m+1
    Real dt_m = (dt_sdc[m_end] - dt_sdc[m_start]) * dt;

#ifdef REACTIONS
    // SDC_Source_Type is only defined for 4th order
    MultiFab tmp;
    MultiFab& C_source = (sdc_order == 4) ? get_new_data(SDC_Source_Type) : tmp;

    if (sdc_order == 4)
    {

        // for 4th order reacting flow, we need to create the "source" C
        // as averages and then convert it to cell centers.  The cell-center
        // version needs to have 2 ghost cells
        for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();
            amrex::Array4<amrex::Real> const& C_source_arr=C_source.array(mfi);

            amrex::Array4<const amrex::Real> const& A_new_arr=(A_new[m_start])->array(mfi);
            amrex::Array4<const amrex::Real> const& A_old_0_arr=(A_old[0])->array(mfi);
            amrex::Array4<const amrex::Real> const& A_old_1_arr=(A_old[1])->array(mfi);
			amrex::Array4<const amrex::Real> const& A_old_2_arr=(A_old[2])->array(mfi);
            amrex::Array4<const amrex::Real> const& R_old_0_arr=(R_old[0])->array(mfi);
            amrex::Array4<const amrex::Real> const& R_old_1_arr=(R_old[1])->array(mfi);
			amrex::Array4<const amrex::Real> const& R_old_2_arr=(R_old[2])->array(mfi);
            if (sdc_quadrature == 0)
            {

                ca_sdc_compute_C4_lobatto(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
										  A_old_2_arr, R_old_0_arr, R_old_1_arr, R_old_2_arr,
										  C_source_arr, m_start);

            }
            else
            {

				amrex::Array4<const amrex::Real> const& A_old_3_arr=(A_old[3])->array(mfi);
				amrex::Array4<const amrex::Real> const& R_old_3_arr=(R_old[3])->array(mfi);

				ca_sdc_compute_C4_radau(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
										A_old_3_arr, R_old_0_arr, R_old_1_arr, R_old_2_arr, R_old_3_arr,
										C_source_arr, m_start);

            }
        }

        // need to construct the time for this stage -- but it is not really
        // at a single instance in time.  For single level this does not matter,
        Real time = state[SDC_Source_Type].curTime();
        AmrLevel::FillPatch(*this, C_source, C_source.nGrow(), time,
                            SDC_Source_Type, 0, NUM_STATE);


        // we'll also construct an initial guess for the nonlinear solve,
        // and store this in the Sburn MultiFab.  We'll use S_new as the
        // staging place so we can do a FillPatch
        MultiFab& S_new = get_new_data(State_Type);

        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();

            amrex::Array4<const amrex::Real> const& k_new_m_start_arr=
                (k_new[m_start])->array(mfi);
            amrex::Array4<const amrex::Real> const& k_new_m_end_arr=(k_new[m_end])->array(
                                                                        mfi);
            amrex::Array4<const amrex::Real> const& A_old_arr=(A_old[m_start])->array(mfi);
            amrex::Array4<const amrex::Real> const& R_old_arr=(R_old[m_start])->array(mfi);
            amrex::Array4<amrex::Real> const& S_new_arr=S_new.array(mfi);

            ca_sdc_compute_initial_guess(bx, k_new_m_start_arr, k_new_m_end_arr,
                                         A_old_arr, R_old_arr, S_new_arr,
                                         dt_m, sdc_iteration);


        }

        const Real cur_time = state[State_Type].curTime();
        expand_state(Sburn, cur_time, 2);

    }
#endif

    // main update loop -- we are updating k_new[m_start] to
    // k_new[m_end]

    FArrayBox U_center;
    FArrayBox C_center;
    FArrayBox U_new_center;
    FArrayBox R_new;

    FArrayBox C2;

    for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi)
    {

        const Box& bx = mfi.tilebox();
        const Box& bx1 = mfi.growntilebox(1);

#ifdef REACTIONS
        // advection + reactions
        if (sdc_order == 2)
        {

            // second order SDC reaction update -- we don't care about
            // the difference between cell-centers and averages

            // first compute the source term, C -- this differs depending
            // on whether we are Lobatto or Radau
            C2.resize(bx, NUM_STATE);
            Elixir elix_C2 = C2.elixir();
            amrex::Array4<amrex::Real> const& C2_arr=C2.array();

            amrex::Array4<const amrex::Real> const& A_new_arr=(A_new[m_start])->array(mfi);
            amrex::Array4<const amrex::Real> const& A_old_0_arr=(A_old[0])->array(mfi);
            amrex::Array4<const amrex::Real> const& A_old_1_arr=(A_old[1])->array(mfi);
            amrex::Array4<const amrex::Real> const& R_old_0_arr=(R_old[0])->array(mfi);
            amrex::Array4<const amrex::Real> const& R_old_1_arr=(R_old[1])->array(mfi);

            if (sdc_quadrature == 0)
            {

                ca_sdc_compute_C2_lobatto(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
                                          R_old_0_arr, R_old_1_arr, C2_arr, m_start);

            }
            else
            {

                amrex::Array4<const amrex::Real> const& A_old_2_arr=(A_old[2])->array(mfi);
                amrex::Array4<const amrex::Real> const& R_old_2_arr=(R_old[2])->array(mfi);
                ca_sdc_compute_C2_radau(bx, dt_m, dt, A_new_arr, A_old_0_arr, A_old_1_arr,
                                        A_old_2_arr,
                                        R_old_0_arr, R_old_1_arr, R_old_2_arr, C2_arr, m_start);

            }

            ca_sdc_update_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                             BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                             BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                             BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                             BL_TO_FORTRAN_3D((*R_old[m_start])[mfi]),
                             BL_TO_FORTRAN_3D(C2),
                             &sdc_iteration,
                             &m_start);
        }
        else
        {

            // fourth order SDC reaction update -- we need to respect the
            // difference between cell-centers and averages

			amrex::Array4<const amrex::Real> const& k_new_m_start_arr=
				(k_new[m_start])->array(mfi);
			amrex::Array4<amrex::Real> const& k_new_m_end_arr=(k_new[m_end])->array(mfi);
			amrex::Array4<amrex::Real> const& C_source_arr=C_source.array(mfi);
            amrex::Array4<const amrex::Real> const& R_new_arr=R_new.array();
			
            // convert the starting U to cell-centered on a fab-by-fab basis
            // -- including one ghost cell
            U_center.resize(bx1, NUM_STATE);
            ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                                BL_TO_FORTRAN_FAB(Sborder[mfi]),
                                BL_TO_FORTRAN_FAB(U_center),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

            // sometimes the Laplacian can make the species go negative near discontinuities
            ca_normalize_species(AMREX_INT_ANYD(bx1.loVect()), AMREX_INT_ANYD(bx1.hiVect()),
                                 BL_TO_FORTRAN_ANYD(U_center));

            // convert the C source to cell-centers
            C_center.resize(bx1, NUM_STATE);
            ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                                BL_TO_FORTRAN_FAB(C_source[mfi]),
                                BL_TO_FORTRAN_FAB(C_center),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

            // solve for the updated cell-center U using our cell-centered C -- we
            // need to do this with one ghost cell
            U_new_center.resize(bx1, NUM_STATE);

            // initialize U_new with our guess for the new state, stored as
            // an average in Sburn
            ca_make_cell_center(BL_TO_FORTRAN_BOX(bx1),
                                BL_TO_FORTRAN_FAB(Sburn[mfi]),
                                BL_TO_FORTRAN_FAB(U_new_center),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

            ca_sdc_update_centers_o4(BL_TO_FORTRAN_BOX(bx1), &dt_m,
                                     BL_TO_FORTRAN_3D(U_center),
                                     BL_TO_FORTRAN_3D(U_new_center),
                                     BL_TO_FORTRAN_3D(C_center),
                                     &sdc_iteration);

            // compute R_i and in 1 ghost cell and then convert to <R> in
            // place (only for the interior)
            R_new.resize(bx1, NUM_STATE);
            ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx1),
                                   BL_TO_FORTRAN_3D(U_new_center),
                                   BL_TO_FORTRAN_3D(R_new));

            ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_FAB(R_new),
                                    AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

            // now do the conservative update using this <R> to get <U>
            // We'll also need to pass in <C>
            ca_sdc_conservative_update(bx, dt_m, k_new_m_start_arr, k_new_m_end_arr,
									   C_source_arr,R_new_arr);

        }
#else
        amrex::Array4<const amrex::Real> const& k_new_m_start_arr=
            (k_new[m_start])->array(mfi);
        amrex::Array4<amrex::Real> const& k_new_m_end_arr=(k_new[m_end])->array(mfi);
        amrex::Array4<const amrex::Real> const& A_new_arr=(A_new[m_start])->array(mfi);
        amrex::Array4<const amrex::Real> const& A_old_0_arr=(A_old[0])->array(mfi);
        amrex::Array4<const amrex::Real> const& A_old_1_arr=(A_old[1])->array(mfi);
        // pure advection
        if (sdc_order == 2)
        {

            if (sdc_quadrature == 0)
            {
                ca_sdc_update_advection_o2_lobatto(bx, dt_m, dt, k_new_m_start_arr,
                                                   k_new_m_end_arr,
                                                   A_new_arr, A_old_0_arr, A_old_1_arr,
                                                   m_start);

            }
            else
            {
                amrex::Array4<const amrex::Real> const& A_old_2_arr=(A_old[2])->array(mfi);
                ca_sdc_update_advection_o2_radau(bx, dt_m, dt, k_new_m_start_arr,
                                                 k_new_m_end_arr,
                                                 A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
                                                 m_start);

            }

        }
        else
        {
            amrex::Array4<const amrex::Real> const& A_old_2_arr=(A_old[2])->array(mfi);
            if (sdc_quadrature == 0)
            {
                ca_sdc_update_advection_o4_lobatto(bx, dt_m, dt, k_new_m_start_arr,
                                                   k_new_m_end_arr,
                                                   A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
                                                   m_start);

            }
            else
            {
                amrex::Array4<const amrex::Real> const& A_old_3_arr=(A_old[3])->array(mfi);
                ca_sdc_update_advection_o4_radau(bx, dt_m, dt, k_new_m_start_arr,
                                                 k_new_m_end_arr,
                                                 A_new_arr, A_old_0_arr, A_old_1_arr, A_old_2_arr,
                                                 A_old_3_arr, m_start);

            }

        }
#endif

    }
}


#ifdef REACTIONS
void
Castro::construct_old_react_source(amrex::MultiFab& U_state,
                                   amrex::MultiFab& R_source,
                                   const bool input_is_average)
{

    BL_PROFILE("Castro::construct_old_react_source()");

    const int* domain_lo = geom.Domain().loVect();
    const int* domain_hi = geom.Domain().hiVect();

    if (sdc_order == 4 && input_is_average)
    {
        // we have cell-averages
        // Note: we cannot tile these operations

        FArrayBox U_center;
        FArrayBox R_center;

        for (MFIter mfi(U_state); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();
            const Box& obx = mfi.growntilebox(1);

            // Convert to centers
            U_center.resize(obx, NUM_STATE);
            ca_make_cell_center(BL_TO_FORTRAN_BOX(obx),
                                BL_TO_FORTRAN_FAB(U_state[mfi]),
                                BL_TO_FORTRAN_FAB(U_center),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

            // burn, including one ghost cell
            R_center.resize(obx, NUM_STATE);
            ca_instantaneous_react(BL_TO_FORTRAN_BOX(obx),
                                   BL_TO_FORTRAN_3D(U_center),
                                   BL_TO_FORTRAN_3D(R_center));

            // at this point, we have the reaction term on centers,
            // including a ghost cell.  Save this into Sburn so we can use
            // it later for the plotfile filling
            Sburn[mfi].copy(R_center, obx, 0, obx, 0, NUM_STATE);

            // convert R to averages (in place)
            ca_make_fourth_in_place(BL_TO_FORTRAN_BOX(bx),
                                    BL_TO_FORTRAN_FAB(R_center),
                                    AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

            // copy this to the center
            R_source[mfi].copy(R_center, bx, 0, bx, 0, NUM_STATE);
        }

    }
    else
    {
        // we are cell-centers

        for (MFIter mfi(U_state); mfi.isValid(); ++mfi)
        {

            const Box& bx = mfi.tilebox();

            // construct the reactive source term
            ca_instantaneous_react(BL_TO_FORTRAN_BOX(bx),
                                   BL_TO_FORTRAN_3D(U_state[mfi]),
                                   BL_TO_FORTRAN_3D(R_source[mfi]));


        }
    }
}
#endif
