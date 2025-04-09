
#include <AMReX_ParmParse.H>
#include <AMReX_LO_BCTYPES.H>
#include <RadBndry.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

int         RadBndry::first = 1;
Vector<int>  RadBndry::bcflag(2*AMREX_SPACEDIM);
Vector<Real> RadBndry::bcval(2*AMREX_SPACEDIM);
Real        RadBndry::time = 0.0;
int         RadBndry::correction = 0;

RadBndry::RadBndry(const BoxArray& _grids, const DistributionMapping& _dmap,
                   const Geometry& _geom) :
    NGBndry(_grids,_dmap,1,_geom)
{
  if (first)
    init();
  first = 0;

  const BoxArray& grids = boxes();
  const Box& domain = geom.Domain();
  //  const Real* dx = geom.CellSize();
  //  const Real* xlo = geom.ProbLo();

  // It is desirable that the type array be set up after static init.
  // This is part of the reason this step is not in NGBndry.
  for (OrientationIter fi; fi; ++fi) {
    Orientation face = fi();
    if (bcflag[face] == 2) {
      bctypearray[face].resize(grids.size());
      for (FabSetIter bi(bndry[face]); bi.isValid(); ++bi) {
        int igrid = bi.index();
        if (domain[face] == boxes()[igrid][face] &&
            !geom.isPeriodic(face.coordDir())) {
          const Box& face_box = bndry[face][bi].box();
          bctypearray[face][igrid].reset(new BaseFab<int>(face_box));
          // We don't care about the bndry values here, only the type array.
#if 0
          FORT_RADBNDRY2(BL_TO_FORTRAN(bndry[face][bi]),
                         bctypearray[face][igrid]->dataPtr(),
                         AMREX_ARLIM(domain.loVect()), AMREX_ARLIM(domain.hiVect()), dx, xlo, time);
#endif
        }
      }
    }
  }
}

RadBndry::RadBndry(const BoxArray& _grids, const DistributionMapping& _dmap,
                   const Geometry& _geom, Real bv) :
    NGBndry(_grids,_dmap,1,_geom)
{
  if (first)
    init(bv);
  first = 0;

  setBndryValues(bv);
}

void RadBndry::init()
{
  // obsolete implementation of the Marshak boundary condition requires
  // us to scale incoming fluxes by c here, since the solver classes
  // do not themselves know about c.  This c compensates
  // for the c in the diffusion coefficient at the boundaries:

//  c = Radiation::clight;

  ParmParse pp("radiation");

  Vector<int> lo_bcflag(AMREX_SPACEDIM, 0), hi_bcflag(AMREX_SPACEDIM, 0);
  pp.queryarr("lo_bcflag",lo_bcflag,0,AMREX_SPACEDIM);
  pp.queryarr("hi_bcflag",hi_bcflag,0,AMREX_SPACEDIM);

  Vector<Real> lo_bcval(AMREX_SPACEDIM, 0.0), hi_bcval(AMREX_SPACEDIM, 0.0);
  pp.queryarr("lo_bcval",lo_bcval,0,AMREX_SPACEDIM);
  pp.queryarr("hi_bcval",hi_bcval,0,AMREX_SPACEDIM);

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    int dir = face.coordDir();
    if (face.isLow()) {
      bcflag[face] = lo_bcflag[dir];
      bcval[face]  = lo_bcval[dir];
    }
    else {
      bcflag[face] = hi_bcflag[dir];
      bcval[face]  = hi_bcval[dir];
    }
  }
}

void RadBndry::init(Real bv)
{
  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    if (face.isLow()) {
      bcflag[face] = 0;
      bcval[face]  = bv;
    }
    else {
      bcflag[face] = 0;
      bcval[face]  = bv;
    }
  }
}

void RadBndry::setBndryConds(const BCRec& bc,
                             const Geometry& geom, IntVect& ratio)
{

//  NOTE: ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL DIMENSIONS
//        *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
  const BoxArray& grids = boxes();
  int ngrds = grids.size();
  const Real* dx = geom.CellSize();
  const Box& domain = geom.Domain();

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    Vector<Real> &bloc = bcloc[face];
    Vector<RadBoundCond> &bctag = bcond[face];

    int dir = face.coordDir();
    Real delta = dx[dir]*ratio[dir];
    int p_bc = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

    for (int i = 0; i < ngrds; i++) {
      const Box& grd = grids[i];

      if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
/*
        // All physical bc values are located on face
        if (p_bc == amrex::BCType::ext_dir) {
          bctag[i] = AMREX_LO_DIRICHLET;
          bloc[i] = 0.;
        }
        else if (p_bc == EXTRAP || p_bc == amrex::BCType::hoextrap || p_bc == amrex::BCType::reflect_even) {
          bctag[i] = AMREX_LO_NEUMANN;
          bloc[i] = 0.;
        }
        else if (p_bc == amrex::BCType::reflect_odd) {
          bctag[i] = AMREX_LO_REFLECT_ODD;
          bloc[i] = 0.;
        }
*/
        if (p_bc == AMREX_LO_DIRICHLET   || p_bc == AMREX_LO_NEUMANN ||
            p_bc == AMREX_LO_REFLECT_ODD) {
          bctag[i] = p_bc;
          bloc[i] = 0.;
        }
        else if (p_bc == AMREX_LO_MARSHAK || p_bc == AMREX_LO_SANCHEZ_POMRANING) {
          bctag[i] = p_bc;
          //gives asymmetric, second-order version of Marshak b.c.
          // (worked for bbmg, works with nonsymmetric hypre solvers):
          bloc[i] = 0.;
          //gives symmetric version of Marshak b.c.
          //(hypre symmetric solvers ignore bloc and do this automatically):
          //bloc[i] = -0.5 * dx[dir];
        }
        else {
          std::cerr << "RadBndry---Not a recognized boundary condition" << std::endl;
          exit(1);
        }
      }
      else {
        // internal bndry
        bctag[i] = AMREX_LO_DIRICHLET;
        bloc[i] = 0.5*delta;
      }
    }
  }
}

void RadBndry::setBndryFluxConds(const BCRec& bc, const BC_Mode phys_bc_mode)
{
  const BoxArray& grids = boxes();
//  int ngrds = grids.size();
  const Real* dx = geom.CellSize();
  const Real* xlo = geom.ProbLo();
  const Box& domain = geom.Domain();

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
//    Vector<Real> &bloc = bcloc[face];
//    Vector<RadBoundCond> &bctag = bcond[face];

    int dir = face.coordDir();
    int p_bc = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

    int p_bcflag = (phys_bc_mode == Inhomogeneous_BC && !correction)
      ? bcflag[face] : 0;
    Real value = (phys_bc_mode == Inhomogeneous_BC && !correction)
      ? bcval[face] : 0.0;

    for (FabSetIter bi(bndry[face]); bi.isValid(); ++bi) {
      int i = bi.index();
      const Box& grd = grids[i];

      if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
        if (bcflag[face] <= 1) {
          if (p_bc == AMREX_LO_MARSHAK   || p_bc == AMREX_LO_SANCHEZ_POMRANING ||
              p_bc == AMREX_LO_DIRICHLET || p_bc == AMREX_LO_NEUMANN) {
              setValue(face, i, value);
          }
        }
        else {
          if(ParallelDescriptor::IOProcessor()) {
            std::cout << "RadBndry ERROR 2: feature not implemented" << std::endl;
          }
          exit(2);

#if 0
          FArrayBox& bnd_fab = bndry[face][bi];
          BaseFab<int>& tfab = bctypearray[face][i];
          FORT_RADBNDRY2(BL_TO_FORTRAN(bnd_fab),
                         tfab.dataPtr(),
                         AMREX_ARLIM(domain.loVect()), AMREX_ARLIM(domain.hiVect()), dx, xlo, time);
          if (p_bcflag == 0) {
            // Homogeneous case.  We called RADBNDRY2 only to set tfab right.
            setValue(face, i, value);
          }
#endif
        }
      }
    }
  }
}

// *************************************************************************

void RadBndry::setHomogValues(const BCRec& bc, IntVect& ratio)
{

  setBndryConds(bc, geom, ratio);

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    for (FabSetIter bi(bndry[face]); bi.isValid(); ++bi) {
      FArrayBox& bnd_fab = bndry[face][bi];
      bnd_fab.setVal<RunOn::Host>(0.0);
    }
  }
}
