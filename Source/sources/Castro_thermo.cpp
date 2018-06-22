#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::construct_old_thermo_source(MultiFab& source, MultiFab& state, Real time, Real dt)
{
  // we only include p divU in method of lines integration
  if (do_ctu) return;

  MultiFab thermo_src(grids, dmap, NUM_STATE, 0);

  thermo_src.setVal(0.0);

  fill_thermo_source(time, dt, state, state, thermo_src);

  Real mult_factor = 1.0;

  MultiFab::Saxpy(source, mult_factor, thermo_src, 0, 0, NUM_STATE, 0);
}



void
Castro::construct_new_thermo_source(MultiFab& source, MultiFab& state_old, MultiFab& state_new, Real time, Real dt)
{
  // we only include p divU in method of lines integration
  if (do_ctu) return;

  amrex::Abort("you should not get here!");
}



void
Castro::fill_thermo_source (Real time, Real dt,
                            MultiFab& state_old, MultiFab& state_new,
                            MultiFab& thermo_src)
{
  const Real* dx = geom.CellSize();
  const Real* prob_lo = geom.ProbLo();

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (MFIter mfi(thermo_src, true); mfi.isValid(); ++mfi)
    {

      const Box& bx = mfi.tilebox();

#pragma gpu
      ca_thermo_src(AMREX_ARLIM_ARG(bx.loVect()), AMREX_ARLIM_ARG(bx.hiVect()),
                    BL_TO_FORTRAN_3D(state_old[mfi]),
                    BL_TO_FORTRAN_3D(state_new[mfi]),
                    BL_TO_FORTRAN_3D(thermo_src[mfi]),
                    ZFILL(prob_lo),ZFILL(dx),time,dt);
    }
}
