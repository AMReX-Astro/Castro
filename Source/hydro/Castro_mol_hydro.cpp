#include "Castro.H"
#include "Castro_F.H"
#include "Castro_hydro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

void
Castro::construct_mol_hydro_source(Real time, Real dt, MultiFab& A_update)
{

  BL_PROFILE("Castro::construct_mol_hydro_source()");

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a

  // update to the state, is stored in the multifab A_update, which is
  // passed in

  const Real strt_time = ParallelDescriptor::second();

  if (verbose && ParallelDescriptor::IOProcessor()) {
    if (time_integration_method == MethodOfLines) {
      std::cout << "... hydro MOL stage " << mol_iteration << std::endl;
    } else if (time_integration_method == SpectralDeferredCorrections) {
      std::cout << "... SDC iteration: " << sdc_iteration << "; current node: " << current_sdc_node << std::endl;
    }
  }

  // we'll add each stage's contribution to -div{F(U)} as we compute them
  // (I don't think we need hydro_source anymore)
  if ((time_integration_method == MethodOfLines && mol_iteration == 0) ||
      (time_integration_method == SpectralDeferredCorrections && current_sdc_node == 0)) {
    hydro_source.setVal(0.0);
  }


  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
  MultiFab& Er_new = get_new_data(Rad_Type);

  if (!Radiation::rad_hydro_combined) {
    amrex::Abort("Castro::construct_mol_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
  }

  int nstep_fsp = -1;
#endif

  BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

  const int* domain_lo = geom.Domain().loVect();
  const int* domain_hi = geom.Domain().hiVect();

#ifndef AMREX_USE_CUDA

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp)
#endif
#endif
  {

    FArrayBox flux[AMREX_SPACEDIM];
#if (AMREX_SPACEDIM <= 2)
    FArrayBox pradial(Box::TheUnitBox(),1);
#endif
#ifdef RADIATION
    FArrayBox rad_flux[AMREX_SPACEDIM];
#endif

#ifdef RADIATION
    int priv_nstep_fsp = -1;
#endif
    // The fourth order stuff cannot do tiling because of the Laplacian corrections
    for (MFIter mfi(S_new, (fourth_order) ? no_tile_size : hydro_tile_size); mfi.isValid(); ++mfi)
      {
	const Box& bx  = mfi.tilebox();

	const int* lo = bx.loVect();
	const int* hi = bx.hiVect();

	FArrayBox &statein  = Sborder[mfi];
	FArrayBox &stateout = S_new[mfi];

	FArrayBox &source_in  = sources_for_hydro[mfi];

	// the output of this will be stored in the correct stage MF
	FArrayBox &source_out = A_update[mfi];
	FArrayBox &source_hydro_only = hydro_source[mfi];

#ifdef RADIATION
	FArrayBox &Er = Erborder[mfi];
	FArrayBox &lam = lamborder[mfi];
	FArrayBox &Erout = Er_new[mfi];
#endif

	FArrayBox& vol = volume[mfi];

        Real stage_weight = 1.0;

        if (time_integration_method == MethodOfLines) {
          stage_weight = b_mol[mol_iteration];
        }

	// Allocate fabs for fluxes
	for (int i = 0; i < AMREX_SPACEDIM ; i++)  {
	  const Box& bxtmp = amrex::surroundingNodes(bx,i);
	  flux[i].resize(bxtmp,NUM_STATE);
#ifdef RADIATION
	  rad_flux[i].resize(bxtmp,Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  pradial.resize(amrex::surroundingNodes(bx,0),1);
	}
#endif
        if (fourth_order) {
          ca_fourth_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &stage_weight,
             BL_TO_FORTRAN_ANYD(statein), 
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(q_bar[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(qaux_bar[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
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
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);

        } else {
          ca_mol_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &stage_weight,
             BL_TO_FORTRAN_ANYD(statein), 
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
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
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);
        }

	// Store the fluxes from this advance -- we weight them by the
	// integrator weight for this stage
	for (int i = 0; i < AMREX_SPACEDIM ; i++) {
	  (*fluxes[i])[mfi].saxpy(stage_weight, flux[i], 
                                  mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
#ifdef RADIATION
	  (*rad_fluxes[i])[mfi].saxpy(stage_weight, rad_flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  P_radial[mfi].saxpy(stage_weight, pradial,
                              mfi.nodaltilebox(0), mfi.nodaltilebox(0), 0, 0, 1);
	}
#endif
      } // MFIter loop

#ifdef RADIATION
    nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif
  }  // end of omp parallel region

#else
  // CUDA version
  // TODO: add radiation

  MultiFab& k_stage = *k_mol[mol_iteration];

#ifndef RADIATION

#ifdef _OPENMP
#pragma omp parallel
#endif
  {

  FArrayBox flatn;
  FArrayBox div;
  FArrayBox qm, qp;
  FArrayBox shk;
  FArrayBox flux[AMREX_SPACEDIM];
  FArrayBox qe[AMREX_SPACEDIM];
  FArrayBox qi[AMREX_SPACEDIM];

  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();
      const Box& obx = amrex::grow(bx, 1);
      const Box& tbx = amrex::grow(bx, 2);

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

      flatn.resize(obx, 1);
      Elixir elix_flatn = flatn.elixir();

      div.resize(obx, 1);
      Elixir elix_div = div.elixir();

      qm.resize(tbx, AMREX_SPACEDIM*NQ);
      Elixir elix_qm = qm.elixir();

      qp.resize(tbx, AMREX_SPACEDIM*NQ);
      Elixir elix_qp = qp.elixir();

      shk.resize(obx, 1);
      Elixir elix_shk = shk.elixir();

      flux[0].resize(gxbx, NUM_STATE);
      Elixir elix_flux_x = flux[0].elixir();

      qe[0].resize(gxbx, NGDNV);
      Elixir elix_qe_x = qe[0].elixir();

      qi[0].resize(gxbx, NQ);
      Elixir elix_qi_x = qi[0].elixir();

#if AMREX_SPACEDIM >= 2
      flux[1].resize(gybx, NUM_STATE);
      Elixir elix_flux_y = flux[1].elixir();

      qe[1].resize(gybx, NGDNV);
      Elixir elix_qe_y = qe[1].elixir();

      qi[1].resize(gybx, NQ);
      Elixir elix_qi_y = qi[1].elixir();
#endif

#if AMREX_SPACEDIM == 3
      flux[2].resize(gzbx, NUM_STATE);
      Elixir elix_flux_z = flux[2].elixir();

      qe[2].resize(gzbx, NGDNV);
      Elixir elix_qe_z = qe[2].elixir();

      qi[2].resize(gzbx, NQ);
      Elixir elix_qi_z = qi[2].elixir();
#endif

      // Compute divergence of velocity field.
#pragma gpu
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div));

      // Compute flattening coefficient for slope calculations.
#pragma gpu
      ca_uflatten
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           BL_TO_FORTRAN_ANYD(flatn), QPRES+1);

      // Do PPM reconstruction to the zone edges.
      int put_on_edges = 1;

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

        int idir_f = idir + 1;

#pragma gpu
        ca_ppm_reconstruct
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()), put_on_edges, idir_f,
           BL_TO_FORTRAN_ANYD(q[mfi]), NQ, 1, NQ,
           BL_TO_FORTRAN_ANYD(flatn),
           BL_TO_FORTRAN_ANYD(qm),
           BL_TO_FORTRAN_ANYD(qp), NQ, 1, NQ);

      }

      // Compute the shk variable
#pragma gpu
      ca_shock
        (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
         BL_TO_FORTRAN_ANYD(q[mfi]),
         BL_TO_FORTRAN_ANYD(shk),
         AMREX_REAL_ANYD(dx));

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& nbx = amrex::surroundingNodes(bx, idir);

          int idir_f = idir + 1;

#pragma gpu
          cmpflx_plus_godunov
              (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
               BL_TO_FORTRAN_ANYD(qm),
               BL_TO_FORTRAN_ANYD(qp), AMREX_SPACEDIM, idir_f,
               BL_TO_FORTRAN_ANYD(flux[idir]),
               BL_TO_FORTRAN_ANYD(qi[idir]),
               BL_TO_FORTRAN_ANYD(qe[idir]),
               BL_TO_FORTRAN_ANYD(qaux[mfi]),
               BL_TO_FORTRAN_ANYD(shk),
               idir_f, AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

#pragma gpu
          apply_av
              (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
               idir_f, AMREX_REAL_ANYD(dx),
               BL_TO_FORTRAN_ANYD(div),
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(flux[idir]));

#pragma gpu
          normalize_species_fluxes
              (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
               BL_TO_FORTRAN_ANYD(flux[idir]));

#pragma gpu
        scale_flux
            (AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
#if AMREX_SPACEDIM == 1
             BL_TO_FORTRAN_ANYD(qe[idir]),
#endif
             BL_TO_FORTRAN_ANYD(flux[idir]),
             BL_TO_FORTRAN_ANYD(area[idir][mfi]), dt);

          // Store the fluxes from this advance -- we weight them by the
          // integrator weight for this stage

          Array4<Real> const flux_fab = (flux[idir]).array();
          Array4<Real> fluxes_fab = (*fluxes[idir]).array(mfi);
          const int numcomp = NUM_STATE;
          const Real scale = b_mol[mol_iteration];

          AMREX_HOST_DEVICE_FOR_4D(nbx, numcomp, i, j, k, n,
          {
              fluxes_fab(i,j,k,n) = flux_fab(i,j,k,n);
          });

      }

#pragma gpu
      ca_construct_hydro_update_cuda
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           AMREX_REAL_ANYD(dx), dt,
           BL_TO_FORTRAN_ANYD(qe[0]),
           BL_TO_FORTRAN_ANYD(qe[1]),
           BL_TO_FORTRAN_ANYD(qe[2]),
           BL_TO_FORTRAN_ANYD(flux[0]),
           BL_TO_FORTRAN_ANYD(flux[1]),
           BL_TO_FORTRAN_ANYD(flux[2]),
           BL_TO_FORTRAN_ANYD(area[0][mfi]),
           BL_TO_FORTRAN_ANYD(area[1][mfi]),
           BL_TO_FORTRAN_ANYD(area[2][mfi]),
           BL_TO_FORTRAN_ANYD(volume[mfi]),
           BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
           BL_TO_FORTRAN_ANYD(k_stage[mfi]));

  } // MFIter loop

  } // OpenMP loop

#endif // RADIATION

#endif // CUDA check

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


}
