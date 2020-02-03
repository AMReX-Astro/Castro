#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

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


  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

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
	FArrayBox &stateout = S_new[mfi];

	FArrayBox &source_in  = sources_for_hydro[mfi];

	// the output of this will be stored in the correct stage MF
	FArrayBox &source_out = A_update[mfi];

        Real stage_weight = 1.0;

        if (time_integration_method == SpectralDeferredCorrections) {
          stage_weight = node_weights[current_sdc_node];
        }

        // get the flattening coefficient
        flatn.resize(obx, 1);
        Elixir elix_flatn = flatn.elixir();

        Array4<Real> const flatn_arr = flatn.array();

        if (first_order_hydro == 1) {
          AMREX_PARALLEL_FOR_3D(obx, i, j, k, { flatn_arr(i,j,k) = 0.0; });
        } else if (use_flattening == 1) {
#pragma gpu box(obx)
          ca_uflatten
            (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(flatn), QPRES+1);
        } else {
          AMREX_PARALLEL_FOR_3D(obx, i, j, k, { flatn_arr(i,j,k) = 1.0; });
        }

        // get the interface states and shock variable

        shk.resize(obx, 1);
        Elixir elix_shk = shk.elixir();

        // Multidimensional shock detection
        // Used for the hybrid Riemann solver

#ifdef SHOCK_VAR
        bool compute_shock = true;
#else
        bool compute_shock = false;
#endif

        if (hybrid_riemann == 1 || compute_shock) {
#pragma gpu box(obx)
          ca_shock(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(shk),
                   AMREX_REAL_ANYD(dx));
        }
        else {
          shk.setVal(0.0);
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
        Elixir elix_avis = avis.elixir();

#ifndef AMREX_USE_CUDA
        if (sdc_order == 4) {

          // -----------------------------------------------------------------
          // fourth order method
          // -----------------------------------------------------------------

          const int* lo = bx.loVect();
          const int* hi = bx.hiVect();

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

            qp.resize(obx2, NQ);
            Elixir elix_qp = qp.elixir();

            q_int.resize(nbx1, 1);
            Elixir elix_qint = q_int.elixir();

            q_avg.resize(ibx[idir], NQ);
            Elixir elix_qavg = q_avg.elixir();

            q_fc.resize(nbx, NQ);
            Elixir elix_qfc = q_fc.elixir();

            f_avg.resize(ibx[idir], NUM_STATE);
            Elixir elix_favg = f_avg.elixir();

            int idir_f = idir + 1;

            for (int n = 0; n < NQ; n++) {

              int ncomp_f = n + 1;

              // construct the interface states in the idir direction
              // operate on nbx1
              ca_fourth_interfaces(AMREX_INT_ANYD(nbx1.loVect()), AMREX_INT_ANYD(nbx1.hiVect()),
                                   idir_f, ncomp_f,
                                   BL_TO_FORTRAN_ANYD(q[mfi]),
                                   BL_TO_FORTRAN_ANYD(q_int),
                                   ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

              // compute the limited interface states
              // operate on obx -- this loop is over cell-centers
              ca_states(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                        idir_f, ncomp_f,
                        BL_TO_FORTRAN_ANYD(q[mfi]),
                        BL_TO_FORTRAN_ANYD(q_int),
                        BL_TO_FORTRAN_ANYD(flatn),
                        BL_TO_FORTRAN_ANYD(qm),
                        BL_TO_FORTRAN_ANYD(qp),
                        ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));
            }

            // get the face-averaged state and flux, <q> and F(<q>),
            // in the idir direction by solving the Riemann problem
            // operate on ibx[idir]
            riemann_state(AMREX_INT_ANYD(ibx[idir].loVect()), AMREX_INT_ANYD(ibx[idir].hiVect()),
                          BL_TO_FORTRAN_ANYD(qm),
                          BL_TO_FORTRAN_ANYD(qp),
                          BL_TO_FORTRAN_ANYD(q_avg),
                          BL_TO_FORTRAN_ANYD(qaux[mfi]),
                          idir_f, 0,
                          ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

            compute_flux_q(AMREX_INT_ANYD(ibx[idir].loVect()), AMREX_INT_ANYD(ibx[idir].hiVect()),
                           BL_TO_FORTRAN_ANYD(q_avg),
                           BL_TO_FORTRAN_ANYD(f_avg),
                           idir_f, 0);


            if (do_hydro == 0) {
              Array4<Real> const f_avg_arr = f_avg.array();

              AMREX_PARALLEL_FOR_4D(nbx, NUM_STATE, i, j, k, n, {
                f_avg_arr(i,j,k,n) = 0.0;});

            }

#ifdef DIFFUSION
            // add diffusive flux to F(<q>) if needed
            // Note: this can act even if do_hydro = 0
            // operate on ibx[idir]
            if (diffuse_temp == 1) {
              int is_avg = 1;
              add_diffusive_flux(AMREX_INT_ANYD(ibx[idir].loVect()), AMREX_INT_ANYD(ibx[idir].hiVect()),
                                 BL_TO_FORTRAN_ANYD(q[mfi]), NQ, QTEMP+1,
                                 BL_TO_FORTRAN_ANYD(q_avg),
                                 BL_TO_FORTRAN_ANYD(f_avg),
                                 ZFILL(dx), idir_f, is_avg);
            }
#endif

            // we now have the face-average interface states and
            // fluxes evaluated with these

#if AMREX_SPACEDIM == 1
            // for 1-d, we are done here, so just copy f_avg to flux[0],
            // since there is no face averaging
            Array4<Real> const flux_arr = (flux[0]).array();
            Array4<const Real> const f_avg_arr = f_avg.array();

            AMREX_PARALLEL_FOR_4D(nbx, NUM_STATE, i, j, k, n, {
                flux_arr(i,j,k,n) = f_avg_arr(i,j,k,n);});

#endif

#if AMREX_SPACEDIM >= 2
            // construct the face-center interface states q_fc
            Array4<Real> const q_fc_arr = q_fc.array();
            Array4<const Real> const q_avg_arr = q_avg.array();

            AMREX_PARALLEL_FOR_4D(nbx, NQ, i, j, k, n, {
                bool test = (n == QGAME) || (n == QGC) || (n == QTEMP);

                if (test) continue;

                Real lap = 0.0;
                int ncomp_f = n + 1;

                trans_laplacian(i, j, k, ncomp_f,
                                idir_f,
                                BL_TO_FORTRAN_ANYD(q_avg), NQ,
                                &lap,
                                ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

                q_fc_arr(i,j,k,n) = q_avg_arr(i,j,k,n) - 1.0/24.0 * lap;

              });

            // compute the face-center fluxes F(q_fc)
            compute_flux_q(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                           BL_TO_FORTRAN_ANYD(q_fc),
                           BL_TO_FORTRAN_ANYD(flux[idir]),
                           idir_f, 0);

            if (do_hydro == 0) {
              Array4<Real> const f_arr = (flux[idir]).array();

              AMREX_PARALLEL_FOR_4D(nbx, NUM_STATE, i, j, k, n, {
                f_arr(i,j,k,n) = 0.0;});

            }


#ifdef DIFFUSION
            if (diffuse_temp == 1) {
              // add the diffusive flux to F(q_fc) if needed
              int is_avg = 0;
              add_diffusive_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                                 BL_TO_FORTRAN_ANYD(T_cc[mfi]), 1, 1,
                                 BL_TO_FORTRAN_ANYD(q_fc),
                                 BL_TO_FORTRAN_ANYD(flux[idir]),
                                 ZFILL(dx), idir_f, is_avg);
            }
#endif


            // compute the final fluxes
            Array4<Real> const flux_arr = (flux[idir]).array();

            AMREX_PARALLEL_FOR_4D(nbx, NUM_STATE, i, j, k, n, {

                Real lap = 0.0;
                int ncomp_f = n + 1;

                trans_laplacian(i, j, k, ncomp_f,
                                idir_f,
                                BL_TO_FORTRAN_ANYD(f_avg), NUM_STATE,
                                &lap,
                                ARLIM_3D(domain_lo), ARLIM_3D(domain_hi));

                flux_arr(i,j,k,n) = flux_arr(i,j,k,n) + 1.0/24.0 * lap;

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

              avisc(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                    BL_TO_FORTRAN_ANYD(q_bar[mfi]),
                    BL_TO_FORTRAN_ANYD(qaux_bar[mfi]),
                    ZFILL(dx),
                    BL_TO_FORTRAN_ANYD(avis),
                    idir_f);


              Array4<Real const> const uin_arr = statein.array();
              Array4<Real const> const avis_arr = avis.array();

              AMREX_PARALLEL_FOR_4D(nbx, NUM_STATE, i, j, k, n, {
                  if (n == Temp) {
                    flux_arr(i,j,k,n) = 0.0;
#ifdef SHOCK_VAR
                  } else if (n == Shock) {
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
            ca_store_godunov_state(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                                   BL_TO_FORTRAN_ANYD(q_avg),
                                   BL_TO_FORTRAN_ANYD(qe[idir]));

          }


        } else {
#endif   // AMREX_USE_CUDA

          // -----------------------------------------------------------------
          // second order method
          // -----------------------------------------------------------------

          // get div{U} -- we'll use this for artificial viscosity
          div.resize(obx, 1);
          Elixir elix_div = div.elixir();

          if (do_hydro) {

#pragma gpu box(obx)
            divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 AMREX_REAL_ANYD(dx),
                 BL_TO_FORTRAN_ANYD(div));

          }

          const Box& tbx = amrex::grow(bx, 2);

          qm.resize(tbx, NQ);
          Elixir elix_qm = qm.elixir();

          qp.resize(tbx, NQ);
          Elixir elix_qp = qp.elixir();

          // compute the fluxes and add artificial viscosity

          q_int.resize(obx, NQ);
          Elixir elix_q_int = q_int.elixir();

          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const int idir_f = idir + 1;

            if (do_hydro) {

              if (ppm_type == 0) {

                dq.resize(obx, NQ);
                Elixir elix_dq = dq.elixir();

#pragma gpu box(obx)
                ca_mol_plm_reconstruct
                  (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                   idir_f,
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(flatn),
                   BL_TO_FORTRAN_ANYD(dq),
                   BL_TO_FORTRAN_ANYD(qm),
                   BL_TO_FORTRAN_ANYD(qp),
                   AMREX_REAL_ANYD(dx),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

              } else {

#pragma gpu box(obx)
                ca_mol_ppm_reconstruct
                  (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                   idir_f,
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(flatn),
                   BL_TO_FORTRAN_ANYD(qm),
                   BL_TO_FORTRAN_ANYD(qp),
                   AMREX_REAL_ANYD(dx));
              }

              const Box& nbx = amrex::surroundingNodes(bx, idir);

#pragma gpu box(nbx)
              cmpflx_plus_godunov
                (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                 BL_TO_FORTRAN_ANYD(qm),
                 BL_TO_FORTRAN_ANYD(qp),
                 BL_TO_FORTRAN_ANYD(flux[idir]),
                 BL_TO_FORTRAN_ANYD(q_int),
                 BL_TO_FORTRAN_ANYD(qe[idir]),
                 BL_TO_FORTRAN_ANYD(qaux[mfi]),
                 BL_TO_FORTRAN_ANYD(shk),
                 idir_f, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

              // set UTEMP and USHK fluxes to zero
              Array4<Real> const flux_arr = (flux[idir]).array();
              const int temp_comp = Temp;
#ifdef SHOCK_VAR
              const int shk_comp = Shock;
#endif

              AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                                    {
                                      flux_arr(i,j,k,temp_comp) = 0.e0;
#ifdef SHOCK_VAR
                                      flux_arr(i,j,k,shk_comp) = 0.e0;
#endif
                                    });


              // apply artificial viscosity
#pragma gpu box(nbx)
              apply_av
                (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                 idir_f, AMREX_REAL_ANYD(dx),
                 BL_TO_FORTRAN_ANYD(div),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(flux[idir]));

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

#ifdef DIFFUSION
            ca_fill_temp_cond
              (AMREX_ARLIM_ANYD(obx.loVect()), AMREX_ARLIM_ANYD(obx.hiVect()),
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(cond));

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            ca_mol_diffusive_flux
              (AMREX_ARLIM_ANYD(nbx.loVect()), AMREX_ARLIM_ANYD(nbx.hiVect()),
               idir_f,
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(cond),
               BL_TO_FORTRAN_ANYD(flux[idir]),
               AMREX_ZFILL(dx));

#endif
          } // end idir loop

#ifndef AMREX_USE_CUDA
        } // end of 4th vs 2nd order MOL update
#endif

        if (do_hydro) {

          for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            int idir_f = idir + 1;


            // apply the density flux limiter
            if (limit_fluxes_on_small_dens == 1) {
#pragma gpu box(nbx)
              limit_hydro_fluxes_on_small_dens
                (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                 idir_f,
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(volume[mfi]),
                 BL_TO_FORTRAN_ANYD(flux[idir]),
                 BL_TO_FORTRAN_ANYD(area[idir][mfi]),
                 dt, AMREX_REAL_ANYD(dx));
            }

            // ensure that the species fluxes are normalized
#pragma gpu box(nbx)
            normalize_species_fluxes
              (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
               BL_TO_FORTRAN_ANYD(flux[idir]));

          }

        }

        // do the conservative update -- and store the shock variable
#pragma gpu box(bx)
        ca_mol_consup
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           BL_TO_FORTRAN_ANYD(shk),
           BL_TO_FORTRAN_ANYD(statein),
           BL_TO_FORTRAN_ANYD(source_in),
           BL_TO_FORTRAN_ANYD(source_out),
           AMREX_REAL_ANYD(dx), dt,
           BL_TO_FORTRAN_ANYD(flux[0]),
#if AMREX_SPACEDIM >= 2
           BL_TO_FORTRAN_ANYD(flux[1]),
#endif
#if AMREX_SPACEDIM == 3
           BL_TO_FORTRAN_ANYD(flux[2]),
#endif
           BL_TO_FORTRAN_ANYD(area[0][mfi]),
#if AMREX_SPACEDIM >= 2
           BL_TO_FORTRAN_ANYD(area[1][mfi]),
#endif
#if AMREX_SPACEDIM == 3
           BL_TO_FORTRAN_ANYD(area[2][mfi]),
#endif
           BL_TO_FORTRAN_ANYD(qe[0]),
#if AMREX_SPACEDIM >= 2
           BL_TO_FORTRAN_ANYD(qe[1]),
#endif
#if AMREX_SPACEDIM == 3
           BL_TO_FORTRAN_ANYD(qe[2]),
#endif
           BL_TO_FORTRAN_ANYD(volume[mfi]));


        // scale the fluxes
#if AMREX_SPACEDIM <= 2
        if (!Geom().IsCartesian()) {
          pradial.resize(xbx, 1);
        }
        Elixir elix_pradial = pradial.elixir();

        Array4<Real> pradial_fab = pradial.array();
#endif

        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

#pragma gpu box(nbx)
          scale_flux(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
#if AMREX_SPACEDIM == 1
                     BL_TO_FORTRAN_ANYD(qe[idir]),
#endif
                     BL_TO_FORTRAN_ANYD(flux[idir]),
                     BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);


          if (idir == 0) {
            // get the scaled radial pressure -- we need to treat this specially
            Array4<Real> const qex_fab = qe[idir].array();
            const int prescomp = GDPRES;

#if AMREX_SPACEDIM == 1
            if (!Geom().IsCartesian()) {
              AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
                                    {
                                      pradial_fab(i,j,k) = qex_fab(i,j,k,prescomp) * dt;
                                    });
            }
#endif

#if AMREX_SPACEDIM == 2
            if (!mom_flux_has_p[0][0]) {
              AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
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


  if (print_update_diagnostics)
    {

      bool local = true;
      Vector<Real> hydro_update = evaluate_source_change(A_update, dt, local);

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceRealSum(hydro_update.dataPtr(), hydro_update.size(), ParallelDescriptor::IOProcessorNumber());

	  if (ParallelDescriptor::IOProcessor())
	    std::cout << std::endl << "  Contributions to the state from the hydro source:" << std::endl;

	  print_source_change(hydro_update);

#ifdef BL_LAZY
	});
#endif
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
