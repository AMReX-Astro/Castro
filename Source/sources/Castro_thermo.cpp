#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::construct_old_thermo_source(MultiFab& source, MultiFab& state_in, Real time, Real dt)
{

#ifndef MHD
  if (!(time_integration_method == SpectralDeferredCorrections)) return;
#endif

  const Real strt_time = ParallelDescriptor::second();

  MultiFab thermo_src(grids, dmap, source.nComp(), 0);

  thermo_src.setVal(0.0);

  fill_thermo_source(state_in, thermo_src);

  Real mult_factor = 1.0;

  MultiFab::Saxpy(source, mult_factor, thermo_src, 0, 0, source.nComp(), 0);

  if (verbose > 1)
  {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      if (ParallelDescriptor::IOProcessor())
          std::cout << "Castro::construct_old_thermo_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
      });
#endif
    }

}



void
Castro::construct_new_thermo_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{

#ifndef MHD
  if (!(time_integration_method == SpectralDeferredCorrections)) return;

  amrex::Abort("you should not get here!");
#else

  const Real strt_time = ParallelDescriptor::second();

  MultiFab thermo_src(grids, dmap, source.nComp(), 0);

  thermo_src.setVal(0.0);

  //Substract off the old-time value first
  Real old_time = time - dt;

  fill_thermo_source(state_old, thermo_src);

  Real mult_factor = -0.5;

  MultiFab::Saxpy(source, mult_factor, thermo_src, 0, 0, source.nComp(), 0);

  //Time center with the new data

  thermo_src.setVal(0.0);

  mult_factor = 0.5;

  //state_new needs ghost cell
  FillPatchIterator fpi(*this, state_new, 1, time, State_Type, 0, NUM_STATE);
  MultiFab& grown_state = fpi.get_mf();

  fill_thermo_source(grown_state, thermo_src);

  MultiFab::Saxpy(source, mult_factor, thermo_src, 0, 0, source.nComp(), 0);

  if (verbose > 1)
  {
      const int IOProc   = ParallelDescriptor::IOProcessorNumber();
      Real      run_time = ParallelDescriptor::second() - strt_time;

#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
      ParallelDescriptor::ReduceRealMax(run_time,IOProc);

      if (ParallelDescriptor::IOProcessor())
          std::cout << "Castro::construct_new_thermo_source() time = " << run_time << "\n" << "\n";
#ifdef BL_LAZY
      });
#endif
    }


#endif //MHD
}



void
Castro::fill_thermo_source (MultiFab& state_in, MultiFab& thermo_src)
{

  // Compute thermodynamic sources for the internal energy equation.
  // At the moment, this is only the -p div{U} term in the internal
  // energy equation, and only for method-of-lines integration,
  // including the new SDC method (the `-` is because it is on the RHS
  // of the equation)

  const auto dx = geom.CellSizeArray();
  const auto prob_lo = geom.ProbLoArray();

  auto coord = geom.Coord();


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(thermo_src, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    Array4<Real const> const U = state_in.array(mfi);
    Array4<Real> const src = thermo_src.array(mfi);


    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

      // radius for non-Cartesian
      Real rp = prob_lo[0] + (static_cast<Real>(i) + 1.5_rt)*dx[0];
      Real rm = prob_lo[0] + (static_cast<Real>(i) - 0.5_rt)*dx[0];
      Real r = 0.5_rt*(rm + rp);

      // compute -div{U}
      if (coord == 0) {
        src(i,j,k,UEINT) = -0.5_rt*(U(i+1,j,k,UMX)/U(i+1,j,k,URHO)  -
                                    U(i-1,j,k,UMX)/U(i-1,j,k,URHO))/dx[0];

      } else if (coord == 1) {
        // axisymmetric
        src(i,j,k,UEINT) = -0.5_rt*(rp*U(i+1,j,k,UMX)/U(i+1,j,k,URHO) -
                                    rm*U(i-1,j,k,UMX)/U(i-1,j,k,URHO))/(r*dx[0]);

      } else if (coord == 2) {
        // spherical
        src(i,j,k,UEINT) = -0.5_rt*(rp*rp*U(i+1,j,k,UMX)/U(i+1,j,k,URHO) -
                                    rm*rm*U(i-1,j,k,UMX)/U(i-1,j,k,URHO))/(r*r*dx[0]);
      }

#if BL_SPACEDIM >= 2
      src(i,j,k,UEINT) += -0.5_rt*(U(i,j+1,k,UMY)/U(i,j+1,k,URHO) -
                                   U(i,j-1,k,UMY)/U(i,j-1,k,URHO))/dx[1];
#endif
#if BL_SPACEDIM == 3
      src(i,j,k,UEINT) += -0.5_rt*(U(i,j,k+1,UMZ)/U(i,j,k+1,URHO) -
                                   U(i,j,k-1,UMZ)/U(i,j,k-1,URHO))/dx[2];
#endif

      // we now need the pressure -- we will assume that the
      // temperature is consistent with the input state
      eos_rep_t eos_state;
      eos_state.rho = U(i,j,k,URHO);
      eos_state.T = U(i,j,k,UTEMP);
      for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = U(i,j,k,UFS+n)/U(i,j,k,URHO);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
        eos_state.aux[n] = U(i,j,k,UFX+n)/U(i,j,k,URHO);
      }
#endif

      eos(eos_input_rt, eos_state);

      // final source term, -p div{U}
      src(i,j,k,UEINT) = eos_state.p * src(i,j,k,UEINT);

    });
  }
}
