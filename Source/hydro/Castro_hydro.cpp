#include "Castro.H"
#include "Castro_F.H"

#ifdef RADIATION
#include "Radiation.H"
#endif

using namespace amrex;

#ifndef AMREX_USE_CUDA
void
Castro::construct_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_hydro_source()");

  const Real strt_time = ParallelDescriptor::second();

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using the CTU framework for unsplit hydrodynamics

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Entering hydro advance" << std::endl << std::endl;

    hydro_source.setVal(0.0);

    int finest_level = parent->finestLevel();

    const Real *dx = geom.CellSize();

    const int* domain_lo = geom.Domain().loVect();
    const int* domain_hi = geom.Domain().hiVect();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef RADIATION
    MultiFab& Er_new = get_new_data(Rad_Type);

    if (!Radiation::rad_hydro_combined) {
      amrex::Abort("Castro::construct_hydro_source -- we don't implement a mode where we have radiation, but it is not coupled to hydro");
    }

    int nstep_fsp = -1;
#endif

    // note: the radiation consup currently does not fill these
    Real mass_lost       = 0.;
    Real xmom_lost       = 0.;
    Real ymom_lost       = 0.;
    Real zmom_lost       = 0.;
    Real eden_lost       = 0.;
    Real xang_lost       = 0.;
    Real yang_lost       = 0.;
    Real zang_lost       = 0.;

    BL_PROFILE_VAR("Castro::advance_hydro_ca_umdrv()", CA_UMDRV);

#ifdef _OPENMP
#ifdef RADIATION
#pragma omp parallel reduction(max:nstep_fsp)
#else
#pragma omp parallel reduction(+:mass_lost,xmom_lost,ymom_lost,zmom_lost) \
		     reduction(+:eden_lost,xang_lost,yang_lost,zang_lost)
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

      int priv_nstep_fsp = -1;

      int is_finest_level = (level == finest_level) ? 1 : 0;

      for (MFIter mfi(S_new,hydro_tile_size); mfi.isValid(); ++mfi)
      {
	  const Box& bx    = mfi.tilebox();

	  const int* lo = bx.loVect();
	  const int* hi = bx.hiVect();

	  FArrayBox &statein  = Sborder[mfi];
	  FArrayBox &stateout = S_new[mfi];

	  FArrayBox &source_out = hydro_source[mfi];

#ifdef RADIATION
          FArrayBox &Er = Erborder[mfi];
          FArrayBox &lam = lamborder[mfi];
          FArrayBox &Erout = Er_new[mfi];
#endif

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

	  ca_ctu_update
	    (ARLIM_3D(lo), ARLIM_3D(hi), &is_finest_level, &time,
	     ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
	     BL_TO_FORTRAN_ANYD(statein), 
	     BL_TO_FORTRAN_ANYD(stateout),
#ifdef RADIATION
	     BL_TO_FORTRAN_ANYD(Er), 
	     BL_TO_FORTRAN_ANYD(Erout),
#endif
	     BL_TO_FORTRAN_ANYD(q[mfi]),
	     BL_TO_FORTRAN_ANYD(qaux[mfi]),
	     BL_TO_FORTRAN_ANYD(src_q[mfi]),
	     BL_TO_FORTRAN_ANYD(source_out),
	     ZFILL(dx), &dt,
	     D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
		    BL_TO_FORTRAN_ANYD(flux[1]),
		    BL_TO_FORTRAN_ANYD(flux[2])),
#ifdef RADIATION
	     D_DECL(BL_TO_FORTRAN_ANYD(rad_flux[0]),
		    BL_TO_FORTRAN_ANYD(rad_flux[1]),
		    BL_TO_FORTRAN_ANYD(rad_flux[2])),
#endif
	     D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
		    BL_TO_FORTRAN_ANYD(area[1][mfi]),
		    BL_TO_FORTRAN_ANYD(area[2][mfi])),
#if (AMREX_SPACEDIM < 3)
	     BL_TO_FORTRAN_ANYD(pradial),
	     BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
	     BL_TO_FORTRAN_ANYD(volume[mfi]),
	     verbose,
#ifdef RADIATION
	     &priv_nstep_fsp,
#endif
	     mass_lost, xmom_lost, ymom_lost, zmom_lost,
	     eden_lost, xang_lost, yang_lost, zang_lost);

	  // Store the fluxes from this advance.
	  // For normal integration we want to add the fluxes from this advance
	  // since we may be subcycling the timestep. But for SDC integration
	  // we want to copy the fluxes since we expect that there will not be
	  // subcycling and we only want the last iteration's fluxes.

	  for (int i = 0; i < AMREX_SPACEDIM ; i++) {
#ifndef SDC
	    (*fluxes    [i])[mfi].plus(    flux[i],mfi.nodaltilebox(i),0,0,NUM_STATE);
#ifdef RADIATION
	    (*rad_fluxes[i])[mfi].plus(rad_flux[i],mfi.nodaltilebox(i),0,0,Radiation::nGroups);
#endif
#else
	    (*fluxes    [i])[mfi].copy(    flux[i],mfi.nodaltilebox(i),0,mfi.nodaltilebox(i),0,NUM_STATE);
#ifdef RADIATION
	    (*rad_fluxes[i])[mfi].copy(rad_flux[i],mfi.nodaltilebox(i),0,mfi.nodaltilebox(i),0,Radiation::nGroups);
#endif	    
#endif
            (*mass_fluxes[i])[mfi].copy(flux[i],mfi.nodaltilebox(i),Density,mfi.nodaltilebox(i),0,1);
	  }

#if (AMREX_SPACEDIM <= 2)
	  if (!Geometry::IsCartesian()) {
#ifndef SDC
	    P_radial[mfi].plus(pradial,mfi.nodaltilebox(0),0,0,1);
#else
	    P_radial[mfi].copy(pradial,mfi.nodaltilebox(0),0,mfi.nodaltilebox(0),0,1);
#endif
	  }
#endif
      } // MFIter loop

#ifdef RADIATION
      nstep_fsp = std::max(nstep_fsp, priv_nstep_fsp);
#endif
    }  // end of omp parallel region

    BL_PROFILE_VAR_STOP(CA_UMDRV);

#ifdef RADIATION
    if (radiation->verbose>=1) {
#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceIntMax(nstep_fsp, ParallelDescriptor::IOProcessorNumber());
	  if (ParallelDescriptor::IOProcessor() && nstep_fsp > 0) {
	    std::cout << "Radiation f-space advection on level " << level
		      << " takes as many as " << nstep_fsp;
	    if (nstep_fsp == 1) {
	      std::cout<< " substep.\n";
	    }
	    else {
	      std::cout<< " substeps.\n";
	    }
	  }
#ifdef BL_LAZY
	});
#endif
    }
#else
    // Flush Fortran output

    if (verbose)
      flush_output();

    if (track_grid_losses)
    {
	material_lost_through_boundary_temp[0] += mass_lost;
	material_lost_through_boundary_temp[1] += xmom_lost;
	material_lost_through_boundary_temp[2] += ymom_lost;
	material_lost_through_boundary_temp[3] += zmom_lost;
	material_lost_through_boundary_temp[4] += eden_lost;
	material_lost_through_boundary_temp[5] += xang_lost;
	material_lost_through_boundary_temp[6] += yang_lost;
	material_lost_through_boundary_temp[7] += zang_lost;
    }

    if (print_update_diagnostics)
    {

	bool local = true;
	Vector<Real> hydro_update = evaluate_source_change(hydro_source, dt, local);

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
#endif

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... Leaving hydro advance" << std::endl << std::endl;

    if (verbose > 0)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

	if (ParallelDescriptor::IOProcessor())
	  std::cout << "Castro::construct_hydro_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
	});
#endif
    }

}
#endif



void
Castro::construct_mol_hydro_source(Real time, Real dt)
{

  BL_PROFILE("Castro::construct_mol_hydro_source()");

  // this constructs the hydrodynamic source (essentially the flux
  // divergence) using method of lines integration.  The output, as a
  // update to the state, is stored in the k_mol array of multifabs.

  const Real strt_time = ParallelDescriptor::second();

  if (verbose && ParallelDescriptor::IOProcessor())
    std::cout << "... hydro MOL stage " << mol_iteration << std::endl;


  // we'll add each stage's contribution to -div{F(U)} as we compute them
  if (mol_iteration == 0) {
    hydro_source.setVal(0.0);
  }


  const Real *dx = geom.CellSize();

  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& k_stage = *k_mol[mol_iteration];

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

    int priv_nstep_fsp = -1;
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
	FArrayBox &source_out = k_stage[mfi];
	FArrayBox &source_hydro_only = hydro_source[mfi];

#ifdef RADIATION
	FArrayBox &Er = Erborder[mfi];
	FArrayBox &lam = lamborder[mfi];
	FArrayBox &Erout = Er_new[mfi];
#endif

	// All cate fabs for fluxes
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
             &(b_mol[mol_iteration]),
             BL_TO_FORTRAN_ANYD(statein), 
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(q_bar[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
             D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2])),
             D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi])),
#if (AMREX_SPACEDIM < 3)
             BL_TO_FORTRAN_ANYD(pradial),
             BL_TO_FORTRAN_ANYD(dLogArea[0][mfi]),
#endif
             BL_TO_FORTRAN_ANYD(volume[mfi]),
             verbose);

        } else {
          ca_mol_single_stage
            (ARLIM_3D(lo), ARLIM_3D(hi), &time, ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
             &(b_mol[mol_iteration]),
             BL_TO_FORTRAN_ANYD(statein), 
             BL_TO_FORTRAN_ANYD(stateout),
             BL_TO_FORTRAN_ANYD(q[mfi]),
             BL_TO_FORTRAN_ANYD(qaux[mfi]),
             BL_TO_FORTRAN_ANYD(source_in),
             BL_TO_FORTRAN_ANYD(source_out),
             BL_TO_FORTRAN_ANYD(source_hydro_only),
             ZFILL(dx), &dt,
             D_DECL(BL_TO_FORTRAN_ANYD(flux[0]),
                    BL_TO_FORTRAN_ANYD(flux[1]),
                    BL_TO_FORTRAN_ANYD(flux[2])),
             D_DECL(BL_TO_FORTRAN_ANYD(area[0][mfi]),
                    BL_TO_FORTRAN_ANYD(area[1][mfi]),
                    BL_TO_FORTRAN_ANYD(area[2][mfi])),
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
	  (*fluxes    [i])[mfi].saxpy(b_mol[mol_iteration], flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, NUM_STATE);
#ifdef RADIATION
	  (*rad_fluxes[i])[mfi].saxpy(b_mol[mol_iteration], rad_flux[i], 
				      mfi.nodaltilebox(i), mfi.nodaltilebox(i), 0, 0, Radiation::nGroups);
#endif
	}

#if (AMREX_SPACEDIM <= 2)
	if (!Geometry::IsCartesian()) {
	  P_radial[mfi].saxpy(b_mol[mol_iteration], pradial,
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

  MultiFab flatn;
  flatn.define(grids, dmap, 1, 1);

  MultiFab div;
  div.define(grids, dmap, 1, 1);

  MultiFab qm;
  qm.define(grids, dmap, AMREX_SPACEDIM*NQ, 2);

  MultiFab qp;
  qp.define(grids, dmap, AMREX_SPACEDIM*NQ, 2);

  MultiFab shk;
  shk.define(grids, dmap, 1, 1);

  // can this be replaced with a single FAB in the loop?
  MultiFab qi_tmp;
  qi_tmp.define(grids, dmap, NQ, 1);

  MultiFab flux[AMREX_SPACEDIM];
  MultiFab qe[AMREX_SPACEDIM];

  for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      flux[i].define(getEdgeBoxArray(i), dmap, NUM_STATE, 0);
      qe[i].define(getEdgeBoxArray(i), dmap, NGDNV, 0);
  }


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& obx = mfi.growntilebox(1);

      // Compute divergence of velocity field.
#pragma gpu
      divu(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           AMREX_REAL_ANYD(dx),
           BL_TO_FORTRAN_ANYD(div[mfi]));

      // Compute flattening coefficient for slope calculations.
#pragma gpu
      ca_uflatten
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
           BL_TO_FORTRAN_ANYD(q[mfi]),
           BL_TO_FORTRAN_ANYD(flatn[mfi]), QPRES+1);

      // Do PPM reconstruction to the zone edges.
      int put_on_edges = 1;

#pragma gpu
      ca_ppm_reconstruct
          (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()), put_on_edges,
           BL_TO_FORTRAN_ANYD(q[mfi]), NQ, 1, NQ,
           BL_TO_FORTRAN_ANYD(flatn[mfi]),
           BL_TO_FORTRAN_ANYD(qm[mfi]),
           BL_TO_FORTRAN_ANYD(qp[mfi]), NQ, 1, NQ);

      // Compute the shk variable
#pragma gpu
      ca_shock
        (AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
         BL_TO_FORTRAN_ANYD(q[mfi]),
         BL_TO_FORTRAN_ANYD(shk[mfi]),
         AMREX_REAL_ANYD(dx));

  } // MFIter loop



#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

          const Box& ebx = mfi.nodaltilebox(idir);

          int idir_f = idir + 1;

#pragma gpu
          ca_construct_flux_cuda
              (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
               AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
               AMREX_REAL_ANYD(dx), dt,
               idir_f,
               BL_TO_FORTRAN_ANYD(Sborder[mfi]),
               BL_TO_FORTRAN_ANYD(div[mfi]),
               BL_TO_FORTRAN_ANYD(qaux[mfi]),
               BL_TO_FORTRAN_ANYD(shk[mfi]),
               BL_TO_FORTRAN_ANYD(qm[mfi]),
               BL_TO_FORTRAN_ANYD(qp[mfi]),
               BL_TO_FORTRAN_ANYD(qi_tmp[mfi]),
               BL_TO_FORTRAN_ANYD(flux[idir][mfi]),
               BL_TO_FORTRAN_ANYD(area[idir][mfi]));

#pragma gpu
          ca_store_godunov_state
            (AMREX_INT_ANYD(ebx.loVect()), AMREX_INT_ANYD(ebx.hiVect()),
             BL_TO_FORTRAN_ANYD(qi_tmp[mfi]),
             BL_TO_FORTRAN_ANYD(qe[idir][mfi]));

          // Store the fluxes from this advance -- we weight them by the
          // integrator weight for this stage
#ifdef AMREX_USE_CUDA
          Cuda::Device::synchronize();  // because saxpy below is run on cpu
#endif
          (*fluxes[idir])[mfi].saxpy(b_mol[mol_iteration], flux[idir][mfi], ebx, ebx, 0, 0, NUM_STATE);

      }

  } // MFIter loop


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& bx = mfi.tilebox();

#pragma gpu
      ca_construct_hydro_update_cuda
          (AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
           AMREX_REAL_ANYD(dx), dt,
           BL_TO_FORTRAN_ANYD(qe[0][mfi]),
           BL_TO_FORTRAN_ANYD(qe[1][mfi]),
           BL_TO_FORTRAN_ANYD(qe[2][mfi]),
           BL_TO_FORTRAN_ANYD(flux[0][mfi]),
           BL_TO_FORTRAN_ANYD(flux[1][mfi]),
           BL_TO_FORTRAN_ANYD(flux[2][mfi]),
           BL_TO_FORTRAN_ANYD(area[0][mfi]),
           BL_TO_FORTRAN_ANYD(area[1][mfi]),
           BL_TO_FORTRAN_ANYD(area[2][mfi]),
           BL_TO_FORTRAN_ANYD(volume[mfi]),
           BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
           BL_TO_FORTRAN_ANYD(k_stage[mfi]));

  } // MFIter loop

#endif

  BL_PROFILE_VAR_STOP(CA_UMDRV);

  // Flush Fortran output

  if (verbose)
    flush_output();


  if (print_update_diagnostics)
    {

      bool local = true;
      Vector<Real> hydro_update = evaluate_source_change(k_stage, dt, local);

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



void
Castro::cons_to_prim(const Real time)
{

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

        // Convert the conservative state to the primitive variable state.
        // This fills both q and qaux.

#pragma gpu
        ca_ctoprim(AMREX_INT_ANYD(qbx.loVect()), AMREX_INT_ANYD(qbx.hiVect()),
                   BL_TO_FORTRAN_ANYD(Sborder[mfi]),
#ifdef RADIATION
                   BL_TO_FORTRAN_ANYD(Erborder[mfi]),
                   BL_TO_FORTRAN_ANYD(lamborder[mfi]),
#endif
                   BL_TO_FORTRAN_ANYD(q[mfi]),
                   BL_TO_FORTRAN_ANYD(qaux[mfi]));

        // Convert the source terms expressed as sources to the conserved state to those
        // expressed as sources for the primitive state.
#ifndef AMREX_USE_CUDA
        if (do_ctu) {
          ca_srctoprim(BL_TO_FORTRAN_BOX(qbx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       BL_TO_FORTRAN_ANYD(sources_for_hydro[mfi]),
                       BL_TO_FORTRAN_ANYD(src_q[mfi]));
        }
#endif

#ifndef RADIATION

        // Add in the reactions source term; only done in SDC.

#ifdef SDC
#ifdef REACTIONS
        MultiFab& SDC_react_source = get_new_data(SDC_React_Type);

        if (do_react)
	    src_q[mfi].plus(SDC_react_source[mfi],qbx,qbx,0,0,QVAR);
#endif
#endif
#endif
      
    }

}


#ifndef AMREX_USE_CUDA
void
Castro::cons_to_prim_fourth(const Real time)
{
  // convert the conservative state cell averages to primitive cell
  // averages with 4th order accuracy

    MultiFab& S_new = get_new_data(State_Type);

    // we don't support radiation here
#ifdef RADIATION
    amrex::Abort("radiation not supported to fourth order");
#else
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbx = mfi.growntilebox(NUM_GROW);
      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // note: these conversions are using a growntilebox, so it
      // will include ghost cells

      // convert U_avg to U_cc -- this will use a Laplacian
      // operation and will result in U_cc defined only on
      // NUM_GROW-1 ghost cells at the end.
      FArrayBox U_cc;
      U_cc.resize(qbx, NUM_STATE);

      ca_make_cell_center(BL_TO_FORTRAN_BOX(qbxm1),
                          BL_TO_FORTRAN_FAB(Sborder[mfi]),
                          BL_TO_FORTRAN_FAB(U_cc));

      // convert U_avg to q_bar -- this will be done on all NUM_GROW
      // ghost cells.
      FArrayBox qaux_bar;
      qaux_bar.resize(qbx, NQAUX);

      ca_ctoprim(BL_TO_FORTRAN_BOX(qbx),
                 BL_TO_FORTRAN_ANYD(Sborder[mfi]),
                 BL_TO_FORTRAN_ANYD(q_bar[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux_bar));

      // this is what we should construct the flattening coefficient
      // from

      // convert U_cc to q_cc (we'll store this temporarily in q,
      // qaux).  This will remain valid only on the NUM_GROW-1 ghost
      // cells.
      ca_ctoprim(BL_TO_FORTRAN_BOX(qbxm1),
                 BL_TO_FORTRAN_ANYD(U_cc),
                 BL_TO_FORTRAN_ANYD(q[mfi]),
                 BL_TO_FORTRAN_ANYD(qaux[mfi]));
    }


    // check for NaNs
    check_for_nan(q);
    check_for_nan(q_bar);


#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

      const Box& qbxm1 = mfi.growntilebox(NUM_GROW-1);

      // now convert q, qaux into 4th order accurate averages
      // this will create q, qaux in NUM_GROW-1 ghost cells, but that's
      // we need here

      ca_make_fourth_average(BL_TO_FORTRAN_BOX(qbxm1),
                             BL_TO_FORTRAN_FAB(q[mfi]),
                             BL_TO_FORTRAN_FAB(q_bar[mfi]));

      // not sure if we need to convert qaux this way, or if we can
      // just evaluate it (we may not need qaux at all actually)

    }

    check_for_nan(q_bar);
#endif // RADIATION
}
#endif



void
Castro::check_for_cfl_violation(const Real dt)
{

    Real courno = -1.0e+200;

    const Real *dx = geom.CellSize();

    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel reduction(max:courno)
#endif
    for (MFIter mfi(S_new, hydro_tile_size); mfi.isValid(); ++mfi) {

        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_compute_cfl(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(q[mfi]),
                       BL_TO_FORTRAN_ANYD(qaux[mfi]),
                       dt, AMREX_REAL_ANYD(dx), AMREX_MFITER_REDUCE_MAX(&courno), print_fortran_warnings);

    }

    ParallelDescriptor::ReduceRealMax(courno);

    if (courno > 1.0) {
        amrex::Print() << "WARNING -- EFFECTIVE CFL AT LEVEL " << level << " IS " << courno << std::endl << std::endl;

        cfl_violation = 1;
    }

}
