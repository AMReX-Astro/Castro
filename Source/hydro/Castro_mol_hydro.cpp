#include <Castro.H>
#include <Castro_F.H>
#include <Castro_util.H>

#ifdef DIFFUSION
#include <diffusion_util.H>
#endif

#include <fourth_center_average.H>

using namespace amrex;

///
/// Construct the hydrodynamic source at the specified time
/// (essentially the flux divergence).  This source is suitable for
/// method of lines or SDC integration.  The output, as a MultiFab
/// is stored in A_update, which comes through the argument list.
///
void
Castro::construct_mol_hydro_source(Real time, Real dt, MultiFab& A_update)
{

#ifdef RADIATION
  amrex::Abort("Error: radiation not supported for the MOL hydro source term");
#else

  BL_PROFILE("Castro::construct_mol_hydro_source()");


  const Real strt_time = ParallelDescriptor::second();

  if (verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << "... construct advection term, SDC iteration: " << sdc_iteration << "; current node: " << current_sdc_node << std::endl;
  }


  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  int coord = geom.Coord();

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

  GpuArray<int, 3> domain_lo = geom.Domain().loVect3d();
  GpuArray<int, 3> domain_hi = geom.Domain().hiVect3d();

#ifdef HYBRID_MOMENTUM
  GeometryData geomdata = geom.data();
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

    // Declare local storage now. This should be done outside the MFIter loop,
    // and then we will resize the Fabs in each MFIter loop iteration. Then,
    // we apply an Elixir to ensure that their memory is saved until it is no
    // longer needed (only relevant for the asynchronous case, usually on GPUs).

    FArrayBox flatn;
    FArrayBox cond;
    FArrayBox dq;
    FArrayBox src_q;
    FArrayBox shk;
    FArrayBox qm, qp;
    FArrayBox div;
    FArrayBox q_int;
    FArrayBox q_avg;
    FArrayBox q_fc;
    FArrayBox f_avg;
    FArrayBox flux[AMREX_SPACEDIM];
    FArrayBox qe[AMREX_SPACEDIM];
#if AMREX_SPACEDIM <= 2
    FArrayBox pradial;
#endif
    FArrayBox avis;

    // The fourth order stuff cannot do tiling because of the Laplacian corrections
    for (MFIter mfi(S_new, (sdc_order == 4) ? no_tile_size : hydro_tile_size); mfi.isValid(); ++mfi)
      {
        const Box& bx  = mfi.tilebox();

        const Box& obx = amrex::grow(bx, 1);
        const Box& obx2 = amrex::grow(bx, 2);

        FArrayBox &statein  = Sborder[mfi];
        Array4<Real const> const uin_arr = statein.array();

        FArrayBox &stateout = S_new[mfi];

        FArrayBox &source_in  = sources_for_hydro[mfi];
        auto source_in_arr = source_in.array();

        // the output of this will be stored in the correct stage MF
        FArrayBox &source_out = A_update[mfi];
        auto source_out_arr = source_out.array();

        Real stage_weight = 1.0;

        if (time_integration_method == SpectralDeferredCorrections) {
          stage_weight = node_weights[current_sdc_node];
        }

        // get the flattening coefficient
        flatn.resize(obx, 1);
        Elixir elix_flatn = flatn.elixir();

        Array4<Real const> const q_arr = q.array(mfi);
        Array4<Real> const flatn_arr = flatn.array();

        if (first_order_hydro == 1) {
          amrex::ParallelFor(obx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
            flatn_arr(i,j,k) = 0.0;
          });
        } else if (use_flattening == 1) {
          uflatten(obx, q_arr, flatn_arr, QPRES);
        } else {
          amrex::ParallelFor(obx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
            flatn_arr(i,j,k) = 1.0;
          });
        }

        // get the interface states and shock variable

        shk.resize(obx, 1);
        Elixir elix_shk = shk.elixir();

        Array4<Real> const shk_arr = shk.array();

        // Multidimensional shock detection
        // Used for the hybrid Riemann solver

#ifdef SHOCK_VAR
        bool compute_shock = true;
#else
        bool compute_shock = false;
#endif

        if (hybrid_riemann == 1 || compute_shock) {
          shock(obx, q_arr, shk_arr);
        }
        else {
          amrex::ParallelFor(obx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
            shk_arr(i,j,k) = 0.0;
          });
        }

        const Box& xbx = amrex::surroundingNodes(bx, 0);
        const Box& gxbx = amrex::grow(xbx, 1);
#if AMREX_SPACEDIM >= 2
        const Box& ybx = amrex::surroundingNodes(bx, 1);
        const Box& gybx = amrex::grow(ybx, 1);
#endif
#if AMREX_SPACEDIM == 3
        const Box& zbx = amrex::surroundingNodes(bx, 2);
        const Box& gzbx = amrex::grow(zbx, 1);
#endif

        auto qaux_arr = qaux.array(mfi);

        flux[0].resize(xbx, NUM_STATE);
        Elixir elix_flux_x = flux[0].elixir();

        qe[0].resize(gxbx, NGDNV);
        Elixir elix_qe_x = qe[0].elixir();

#if AMREX_SPACEDIM >= 2
        flux[1].resize(ybx, NUM_STATE);
        Elixir elix_flux_y = flux[1].elixir();

        qe[1].resize(gybx, NGDNV);
        Elixir elix_qe_y = qe[1].elixir();
#endif

#if AMREX_SPACEDIM == 3
        flux[2].resize(zbx, NUM_STATE);
        Elixir elix_flux_z = flux[2].elixir();

        qe[2].resize(gzbx, NGDNV);
        Elixir elix_qe_z = qe[2].elixir();
#endif

        avis.resize(obx, 1);
        auto avis_arr = avis.array();
        Elixir elix_avis = avis.elixir();

#ifndef AMREX_USE_GPU
        if (sdc_order == 4) {

          // -----------------------------------------------------------------
          // fourth order method
          // -----------------------------------------------------------------

          Box ibx[AMREX_SPACEDIM];
          ibx[0] = amrex::grow(amrex::surroundingNodes(bx, 0), IntVect(AMREX_D_DECL(0,1,1)));
#if AMREX_SPACEDIM >= 2
          ibx[1] = amrex::grow(amrex::surroundingNodes(bx, 1), IntVect(AMREX_D_DECL(1,0,1)));
#endif
#if AMREX_SPACEDIM == 3
          ibx[2] = amrex::grow(amrex::surroundingNodes(bx, 2), IntVect(AMREX_D_DECL(1,1,0)));
#endif
          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const Box& nbx = amrex::surroundingNodes(bx, idir);
            const Box& nbx1 = amrex::grow(nbx, 1);

            qm.resize(obx2, NQ);
            Elixir elix_qm = qm.elixir();
            auto qm_arr = qm.array();

            qp.resize(obx2, NQ);
            Elixir elix_qp = qp.elixir();
            auto qp_arr = qp.array();

            q_int.resize(nbx1, 1);
            Elixir elix_qint = q_int.elixir();
            auto q_int_arr = q_int.array();

            q_avg.resize(ibx[idir], NQ);
            Elixir elix_qavg = q_avg.elixir();
            auto q_avg_arr = q_avg.array();

            q_fc.resize(nbx, NQ);
            Elixir elix_qfc = q_fc.elixir();
            auto q_fc_arr = q_fc.array();

            f_avg.resize(ibx[idir], NUM_STATE);
            Elixir elix_favg = f_avg.elixir();
            auto f_avg_arr = f_avg.array();

            for (int n = 0; n < NQ; n++) {

              // construct the interface states in the idir direction
              // operate on nbx1
              fourth_interfaces(nbx1,
                                idir, n,
                                q_arr, q_int_arr);

              // compute the limited interface states
              // operate on obx -- this loop is over cell-centers
              states(obx,
                     idir, n,
                     q_arr, q_int_arr, flatn_arr,
                     qm_arr, qp_arr);
            }

            // get the face-averaged state and flux, <q> and F(<q>),
            // in the idir direction by solving the Riemann problem
            // operate on ibx[idir]

            cmpflx_plus_godunov(ibx[idir],
                                qm_arr, qp_arr,
                                f_avg_arr, q_avg_arr,
                                qaux_arr,
                                shk_arr,
                                idir, true);

            if (do_hydro == 0) {
              amrex::ParallelFor(nbx, NUM_STATE,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
              {
                f_avg_arr(i,j,k,n) = 0.0;
              });

            }

#ifdef DIFFUSION
            // add diffusive flux to F(<q>) if needed
            // Note: this can act even if do_hydro = 0
            // operate on ibx[idir]
            if (diffuse_temp == 1) {
              bool is_avg = true;
              fourth_add_diffusive_flux(ibx[idir],
                                        q_arr, QTEMP,
                                        q_avg_arr, f_avg_arr,
                                        idir, is_avg);
            }
#endif

            // we now have the face-average interface states and
            // fluxes evaluated with these

#if AMREX_SPACEDIM == 1
            // for 1-d, we are done here, so just copy f_avg to flux[0],
            // since there is no face averaging
            Array4<Real> const flux_arr = (flux[0]).array();

            amrex::ParallelFor(nbx, NUM_STATE,
            [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
            {
              flux_arr(i,j,k,n) = f_avg_arr(i,j,k,n);
            });

#endif

#if AMREX_SPACEDIM >= 2
            // construct the face-center interface states q_fc
            const int* lo_bc = phys_bc.lo();
            const int* hi_bc = phys_bc.hi();

            GpuArray<bool, AMREX_SPACEDIM> lo_periodic;
            GpuArray<bool, AMREX_SPACEDIM> hi_periodic;
            for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
              lo_periodic[idir] = lo_bc[idir] == Interior;
              hi_periodic[idir] = hi_bc[idir] == Interior;
            }

            amrex::ParallelFor(nbx, NQ,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
              bool test = (n == QGC) || (n == QTEMP);

              if (test) return;

              Real lap = trans_laplacian(i, j, k, n, idir, q_avg_arr,
                                         lo_periodic, hi_periodic, domain_lo, domain_hi);

              q_fc_arr(i,j,k,n) = q_avg_arr(i,j,k,n) - (1.0_rt/24.0_rt) * lap;

            });

            // compute the face-center fluxes F(q_fc)
            Array4<Real> const f_arr = (flux[idir]).array();

            compute_flux_from_q(nbx, q_fc_arr, f_arr, idir);

            if (do_hydro == 0) {

              amrex::ParallelFor(nbx, NUM_STATE,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
              {
                f_arr(i,j,k,n) = 0.0;
              });

            }


#ifdef DIFFUSION
            if (diffuse_temp == 1) {
              // add the diffusive flux to F(q_fc) if needed
              bool is_avg = false;
              fourth_add_diffusive_flux(nbx,
                                        T_cc.array(mfi), 0,
                                        q_fc_arr, flux[idir].array(),
                                        idir, is_avg);
            }
#endif


            // compute the final fluxes
            Array4<Real> const flux_arr = (flux[idir]).array();

            amrex::ParallelFor(nbx, NUM_STATE,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {

              Real lap = trans_laplacian(i, j, k, n, idir, f_avg_arr,
                                         lo_periodic, hi_periodic, domain_lo, domain_hi);

              flux_arr(i,j,k,n) += (1.0_rt/24.0_rt) * lap;

            });

#endif

            // add artifical viscosity
            if (do_hydro == 1) {

              // avisc_coefficient is the coefficent we use.  The
              // McCorquodale & Colella paper suggest alpha = 0.3, but
              // our other hydro solvers use a coefficient on the
              // divergence that defaults to 0.1, so we normalize to
              // that value, to allow for adjustments
              const Real alpha = 0.3;

              Real avisc_coeff = alpha * (difmag / 0.1);

              auto q_bar_arr = q_bar.array(mfi);
              auto qaux_bar_arr = qaux_bar.array(mfi);

              fourth_avisc(nbx,
                           q_bar_arr, qaux_bar_arr,
                           avis_arr, idir);

              Array4<Real const> const avis_arr = avis.array();

              amrex::ParallelFor(nbx, NUM_STATE,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k, int n)
              {
                  if (n == UTEMP) {
                    flux_arr(i,j,k,n) = 0.0;
#ifdef SHOCK_VAR
                  } else if (n == USHK) {
                    flux_arr(i,j,k,n) == 0.0;
#endif
                  } else {

                    if (idir == 0) {
                      flux_arr(i,j,k,n) = flux_arr(i,j,k,n) +
                        avisc_coeff * avis_arr(i,j,k,0) * (uin_arr(i,j,k,n) - uin_arr(i-1,j,k,n));

                    } else if (idir == 1) {
                      flux_arr(i,j,k,n) = flux_arr(i,j,k,n) +
                        avisc_coeff * avis_arr(i,j,k,0) * (uin_arr(i,j,k,n) - uin_arr(i,j-1,k,n));

                    } else {
                      flux_arr(i,j,k,n) = flux_arr(i,j,k,n) +
                        avisc_coeff * avis_arr(i,j,k,0) * (uin_arr(i,j,k,n) - uin_arr(i,j,k-1,n));

                    }
                  }

                });

            }

            // store the Godunov state
            auto qe_arr = (qe[idir]).array();
            store_godunov_state(nbx, q_avg_arr, qe_arr);

          }


        } else {
#endif   // AMREX_USE_GPU

          // -----------------------------------------------------------------
          // second order method
          // -----------------------------------------------------------------

          // get div{U} -- we'll use this for artificial viscosity
          div.resize(obx, 1);
          Elixir elix_div = div.elixir();

          auto div_arr = div.array();

          if (do_hydro) {
            divu(obx, q_arr, div_arr);
          }

          const Box& tbx = amrex::grow(bx, 2);

          qm.resize(tbx, NQ);
          Elixir elix_qm = qm.elixir();
          Array4<Real> const qm_arr = qm.array();

          qp.resize(tbx, NQ);
          Elixir elix_qp = qp.elixir();
          Array4<Real> const qp_arr = qp.array();

          // compute the fluxes and add artificial viscosity

          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            if (do_hydro) {

              if (ppm_type == 0) {

                dq.resize(obx, NQ);
                Elixir elix_dq = dq.elixir();
                auto dq_arr = dq.array();

                // for well-balancing, we need to primitive variable
                // source terms
                const Box& qbx = amrex::grow(bx, NUM_GROW);
                src_q.resize(qbx, NQSRC);
                Elixir elix_src_q = src_q.elixir();
                Array4<Real> const src_q_arr = src_q.array();

                Array4<Real> const src_arr = sources_for_hydro.array(mfi);

                src_to_prim(qbx, q_arr, src_arr, src_q_arr);

                mol_plm_reconstruct(obx, idir,
                                    q_arr, flatn_arr, src_q_arr,
                                    dq_arr,
                                    qm_arr, qp_arr);

              } else {

                mol_ppm_reconstruct(obx, idir,
                                    q_arr, flatn_arr,
                                    qm_arr, qp_arr);
              }

              // ppm_temp_fix = 1
              if (ppm_temp_fix == 1) {
                  edge_state_temp_to_pres(obx, qm.array(), qp.array());
              }

              const Box& nbx = amrex::surroundingNodes(bx, idir);

              auto qe_arr = (qe[idir]).array();
              auto flux_arr = (flux[idir]).array();

              cmpflx_plus_godunov
                (nbx,
                 qm_arr, qp_arr,
                 flux_arr, qe_arr,
                 qaux_arr, shk_arr,
                 idir, false);

              // set UTEMP and USHK fluxes to zero
              Array4<Real const> const uin_arr = Sborder.array(mfi);

              amrex::ParallelFor(nbx,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
              {
                flux_arr(i,j,k,UTEMP) = 0.e0;
#ifdef SHOCK_VAR
                flux_arr(i,j,k,USHK) = 0.e0;
#endif
              });


              // apply artificial viscosity
              apply_av(nbx, idir, div_arr, uin_arr, flux_arr);

            } else {
              // we are not doing hydro, so simply zero out the fluxes
              const Box& nbx = amrex::surroundingNodes(bx, idir);
              const Box& gbx = amrex::grow(nbx, 1);

              Array4<Real> const flux_arr = (flux[idir]).array();
              Array4<Real> const qe_arr = (qe[idir]).array();
              const int nstate = NUM_STATE;

              AMREX_HOST_DEVICE_FOR_4D(gbx, nstate, i, j, k, n,
                                       {
                                         flux_arr(i,j,k,n) = 0.e0;
                                       });

              const int ncomp = NGDNV;

              AMREX_HOST_DEVICE_FOR_4D(gbx, ncomp, i, j, k, n,
                                       {
                                         qe_arr(i,j,k,n) = 0.e0;
                                       });
            } // end do_hydro

            // add a diffusive flux
            cond.resize(obx, 1);
            Elixir elix_cond = cond.elixir();
            auto cond_arr = cond.array();

#ifdef DIFFUSION
            fill_temp_cond(obx, Sborder.array(mfi), cond_arr);

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            mol_diffusive_flux(nbx, idir,
                               uin_arr, cond_arr,
                               flux[idir].array());

#endif
          } // end idir loop

#ifndef AMREX_USE_GPU
        } // end of 4th vs 2nd order MOL update
#endif

        if (do_hydro) {

          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            Array4<Real> const flux_arr = (flux[idir]).array();

            // apply the density flux limiter
            if (limit_fluxes_on_small_dens == 1) {
                limit_hydro_fluxes_on_small_dens
                  (nbx, idir,
                   Sborder.array(mfi),
                   q.array(mfi),
                   volume.array(mfi),
                   flux[idir].array(),
                   area[idir].array(mfi),
                   dt);
            }

            // ensure that the species fluxes are normalized
            normalize_species_fluxes(nbx, flux_arr);

          }

        }

        // do the conservative update -- and store the shock variable
        mol_consup(bx,
#ifdef SHOCK_VAR
                   shk_arr,
#endif
                   source_in_arr,
                   source_out_arr,
                   dt,
                   flux[0].array(),
#if AMREX_SPACEDIM >= 2
                   flux[1].array(),
#endif
#if AMREX_SPACEDIM == 3
                   flux[2].array(),
#endif
                   area[0].array(mfi),
#if AMREX_SPACEDIM >= 2
                   area[1].array(mfi),
#endif
#if AMREX_SPACEDIM == 3
                   area[2].array(mfi),
#endif
#if AMREX_SPACEDIM <= 2
                   qe[0].array(),
#endif
                   volume.array(mfi));


#ifdef HYBRID_MOMENTUM
        auto dx_arr = geom.CellSizeArray();

        amrex::ParallelFor(bx,
        [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
        {

          GpuArray<Real, 3> loc;

          position(i, j, k, geomdata, loc);

          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            loc[dir] -= problem::center[dir];

          Real R = amrex::max(std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]), R_min);
          Real RInv = 1.0_rt / R;

          source_out_arr(i,j,k,UMR) -= ((loc[0] * RInv) * (qx_arr(i+1,j,k,GDPRES) - qx_arr(i,j,k,GDPRES)) / dx_arr[0] +
                                        (loc[1] * RInv) * (qy_arr(i,j+1,k,GDPRES) - qy_arr(i,j,k,GDPRES)) / dx_arr[1]);

        });
#endif

        // scale the fluxes
#if AMREX_SPACEDIM <= 2
        if (!Geom().IsCartesian()) {
          pradial.resize(xbx, 1);
        }
        Elixir elix_pradial = pradial.elixir();

        Array4<Real> pradial_fab = pradial.array();
        Array4<Real> const qex_arr = qe[0].array();
#endif

        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          Array4<Real> const flux_arr = (flux[idir]).array();
          Array4<Real const> const area_arr = (area[idir]).array(mfi);

          scale_flux(nbx,
#if AMREX_SPACEDIM == 1
                     qex_arr,
#endif
                     flux_arr, area_arr, dt);


          if (idir == 0) {
            // get the scaled radial pressure -- we need to treat this specially
            Array4<Real> const qex_fab = qe[idir].array();
            const int prescomp = GDPRES;

#if AMREX_SPACEDIM == 1
            if (!Geom().IsCartesian()) {
              amrex::ParallelFor(nbx,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
              {
                pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
              });
            }
#endif

#if AMREX_SPACEDIM == 2
            if (!mom_flux_has_p(0, 0, coord)) {
              amrex::ParallelFor(nbx,
              [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
              {
                pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
              });
            }
#endif
          }
        }


        // Store the fluxes from this advance -- we weight them by the
        // integrator weight for this stage

        // For SDC, we store node 0 the only time we enter here (the
        // first iteration) and we store the other nodes only on the
        // last iteration.
        if (time_integration_method == SpectralDeferredCorrections &&
             (current_sdc_node == 0 || sdc_iteration == sdc_order+sdc_extra-1)) {

          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            Array4<Real> const flux_fab = (flux[idir]).array();
            Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);
            const int numcomp = NUM_STATE;

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(idir), numcomp, i, j, k, n,
            {
                fluxes_fab(i,j,k,n) += stage_weight * flux_fab(i,j,k,n);
            });

          }

#if AMREX_SPACEDIM <= 2
          if (!Geom().IsCartesian()) {

            Array4<Real> P_radial_fab = P_radial.array(mfi);
            const Real scale = stage_weight;

            AMREX_HOST_DEVICE_FOR_4D(mfi.nodaltilebox(0), 1, i, j, k, n,
            {
                P_radial_fab(i,j,k,0) += scale * pradial_fab(i,j,k,0);
            });

          }
#endif
        }

      } // MFIter loop

  }  // end of omp parallel region


  BL_PROFILE_VAR_STOP(CA_UMDRV);

  // Flush Fortran output

  if (verbose)
    flush_output();


  if (print_update_diagnostics) {
      evaluate_and_print_source_change(A_update, dt, "hydro source");
    }

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
          std::cout << "Castro::construct_mol_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
        });
#endif
    }

#endif // radiation
}
