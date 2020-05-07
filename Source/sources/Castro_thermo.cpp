#include "Castro.H"
#include "Castro_F.H"

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

  fill_thermo_source(time, dt, state_in, state_in, thermo_src);

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

  fill_thermo_source(old_time, dt, state_old, state_old, thermo_src);

  Real mult_factor = -0.5;

  MultiFab::Saxpy(source, mult_factor, thermo_src, 0, 0, source.nComp(), 0);

  //Time center with the new data
  
  thermo_src.setVal(0.0);

  mult_factor = 0.5;
  
  //state_new needs ghost cell  
  FillPatchIterator fpi(*this, state_new, 1, time, State_Type, 0, NUM_STATE);
  MultiFab& grown_state = fpi.get_mf();
  
  fill_thermo_source(time, dt, state_old, grown_state, thermo_src);
  
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
Castro::fill_thermo_source (Real time, Real dt,
                            MultiFab& state_old, MultiFab& state_new,
                            MultiFab& thermo_src)
{

  // Compute thermodynamic sources for the internal energy equation.
  // At the moment, this is only the -p div{U} term in the internal
  // energy equation, and only for method-of-lines integration,
  // including the new SDC method (the `-` is because it is on the RHS
  // of the equation)

  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

  auto coord = geom.Coord(); 


#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(thermo_src, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

    Array4<Real const> const old_state = state_old.array(mfi);
    Array4<Real const> const new_state = state_new.array(mfi);
    Array4<Real> const src = thermo_src.array(mfi);


    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

      // radius for non-Cartesian
      Real rp = prob_lo[0] + (static_cast<Real>(i) + 1.5_rt)*dx[0];
      Real rm = prob_lo[0] + (static_cast<Real>(i) - 0.5_rt)*dx[0];
      Real r = 0.5_rt*(rm + rp);

      // compute -div{U}
      if (coord == 0) {
        src(i,j,k,UEINT) = -0.25_rt*((old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) +
                                      new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) -
                                     (old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) +
                                      new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/dx[0];

      } else if (coord == 1) {
        // axisymmetric
        src(i,j,k,UEINT) = -0.25_rt*(rp*(old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) +
                                         new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) -
                                     rm*(old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) +
                                         new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/(r*dx[0]);

      } else if (coord == 2) {
        // spherical
        src(i,j,k,UEINT) = -0.25_rt*(rp*rp*(old_state(i+1,j,k,UMX)/old_state(i+1,j,k,URHO) +
                                            new_state(i+1,j,k,UMX)/new_state(i+1,j,k,URHO)) -
                                     rm*rm*(old_state(i-1,j,k,UMX)/old_state(i-1,j,k,URHO) +
                                            new_state(i-1,j,k,UMX)/new_state(i-1,j,k,URHO)))/(r*r*dx[0]);
      }

#if BL_SPACEDIM >= 2
      src(i,j,k,UEINT) += -0.25_rt*((old_state(i,j+1,k,UMY)/old_state(i,j+1,k,URHO) +
                                     new_state(i,j+1,k,UMY)/new_state(i,j+1,k,URHO)) -
                                    (old_state(i,j-1,k,UMY)/old_state(i,j-1,k,URHO) +
                                     new_state(i,j-1,k,UMY)/new_state(i,j-1,k,URHO)))/dx[1];
#endif
#if BL_SPACEDIM == 3
      src(i,j,k,UEINT) += -0.25_rt*((old_state(i,j,k+1,UMZ)/old_state(i,j,k+1,URHO) +
                                     new_state(i,j,k+1,UMZ)/new_state(i,j,k+1,URHO)) -
                                    (old_state(i,j,k-1,UMZ)/old_state(i,j,k-1,URHO) +
                                     new_state(i,j,k-1,UMZ)/new_state(i,j,k-1,URHO)))/dx[2];
#endif

      // we now need the pressure -- we will assume that the
      // temperature is consistent with the input state
      eos_t eos_state_old;
      eos_state_old.rho = old_state(i,j,k,URHO);
      eos_state_old.T = old_state(i,j,k,UTEMP);
      for (int n = 0; n < NumSpec; n++) {
        eos_state_old.xn[n] = old_state(i,j,k,UFS+n)/old_state(i,j,k,URHO);
      }
      for (int n = 0; n < NumAux; n++) {
        eos_state_old.aux[n] = old_state(i,j,k,UFX+n)/old_state(i,j,k,URHO);
      }

      eos(eos_input_rt, eos_state_old);

      eos_t eos_state_new;
      eos_state_new.rho = new_state(i,j,k,URHO);
      eos_state_new.T = new_state(i,j,k,UTEMP);
      for (int n = 0; n < NumSpec; n++) {
        eos_state_new.xn[n] = new_state(i,j,k,UFS+n)/new_state(i,j,k,URHO);
      }
      for (int n = 0; n < NumAux; n++) {
        eos_state_new.aux[n] = new_state(i,j,k,UFX+n)/new_state(i,j,k,URHO);
      }

      eos(eos_input_rt, eos_state_new);

      // final source term, -p div{U}
      src(i,j,k,UEINT) = 0.5_rt*(eos_state_old.p + eos_state_new.p)*src(i,j,k,UEINT);

    });
  }
}
