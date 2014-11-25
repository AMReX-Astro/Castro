
#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include "RadBndry.H"

#include <Using.H>

#include "LHH.H"

#undef BL_USE_ARLIM

#include "RAD_F.H"

int         RadBndry::first = 1;
Array<int>  RadBndry::bcflag(2*BL_SPACEDIM);
Array<Real> RadBndry::bcval(2*BL_SPACEDIM);
Real        RadBndry::time = 0.0;
int         RadBndry::correction = 0;

RadBndry::RadBndry(const BoxArray& _grids, const Geometry& _geom) :
  NGBndry(_grids,1,_geom)
{
  if (first)
    init();
  first = 0;

  const BoxArray& grids = boxes();
  const Box& domain = geom.Domain();
  const Real* dx = geom.CellSize();
  const Real* xlo = geom.ProbLo();

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
	  bctypearray[face].set(igrid, new BaseFab<int>(face_box));
          // We don't care about the bndry values here, only the type array.
#if 0
          FORT_RADBNDRY2(bndry[face][bi].dataPtr(), dimlist(face_box),
                         bctypearray[face][igrid].dataPtr(),
                         dimlist(domain), dx, xlo, time);
#endif
        }
      }
    }
  }
}

RadBndry::RadBndry(const BoxArray& _grids, const Geometry& _geom, Real bv) :
  NGBndry(_grids,1,_geom)
{
  if (first)
    init(bv);
  first = 0;

  setBndryValues(bv);
}

RadBndry::~RadBndry()
{
  for (OrientationIter fi; fi; ++fi) {
    Orientation face = fi();
    if (bcflag[face] == 2) {
      int len = grids.size();
      for (int igrid = 0; igrid < len; igrid++) {
	if (bctypearray[face].defined(igrid)) {
	  delete bctypearray[face].remove(igrid);
	}
      }
    }
  }
}

void RadBndry::init()
{
  // obsolete implementation of the Marshak boundary condition requires
  // us to scale incoming fluxes by c here, since the solver classes
  // do not themselves know about c.  This c compensates
  // for the c in the diffusion coefficient at the boundaries:

//  c = Radiation::clight;

  ParmParse pp("radiation");

  Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
  pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
  pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);

  Array<int> lo_bcflag(BL_SPACEDIM, 0), hi_bcflag(BL_SPACEDIM, 0);
  pp.queryarr("lo_bcflag",lo_bcflag,0,BL_SPACEDIM);
  pp.queryarr("hi_bcflag",hi_bcflag,0,BL_SPACEDIM);

  Array<Real> lo_bcval(BL_SPACEDIM, 0.0), hi_bcval(BL_SPACEDIM, 0.0);
  pp.queryarr("lo_bcval",lo_bcval,0,BL_SPACEDIM);
  pp.queryarr("hi_bcval",hi_bcval,0,BL_SPACEDIM);

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    int dir = face.coordDir();
    int bctype;
    if (face.isLow()) {
      bctype       = lo_bc[dir];
      bcflag[face] = lo_bcflag[dir];
      bcval[face]  = lo_bcval[dir];
    }
    else {
      bctype       = hi_bc[dir];
      bcflag[face] = hi_bcflag[dir];
      bcval[face]  = hi_bcval[dir];
    }
  }
}

void RadBndry::init(Real bv)
{
  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    int dir = face.coordDir();
    int bctype;
    if (face.isLow()) {
      bctype       = LO_DIRICHLET;
      bcflag[face] = 0;
      bcval[face]  = bv;
    }
    else {
      bctype       = LO_DIRICHLET;
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
    Array<Real> &bloc = bcloc[face];
    Array<RadBoundCond> &bctag = bcond[face];

    int dir = face.coordDir();
    Real delta = dx[dir]*ratio[dir];
    int p_bc = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

    for (int i = 0; i < ngrds; i++) {
      const Box& grd = grids[i];

      if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
/*
	// All physical bc values are located on face
	if (p_bc == EXT_DIR) {
	  bctag[i] = LO_DIRICHLET;
	  bloc[i] = 0.;
	}
	else if (p_bc == EXTRAP || p_bc == HOEXTRAP || p_bc == REFLECT_EVEN) {
	  bctag[i] = LO_NEUMANN;
	  bloc[i] = 0.;
	}
	else if (p_bc == REFLECT_ODD) {
	  bctag[i] = LO_REFLECT_ODD;
	  bloc[i] = 0.;
	}
*/
	if (p_bc == LO_DIRICHLET   || p_bc == LO_NEUMANN ||
	    p_bc == LO_REFLECT_ODD) {
	  bctag[i] = p_bc;
	  bloc[i] = 0.;
	}
	else if (p_bc == LO_MARSHAK || p_bc == LO_SANCHEZ_POMRANING) {
	  bctag[i] = p_bc;
	  //gives asymmetric, second-order version of Marshak b.c.
          // (worked for bbmg, works with nonsymmetric hypre solvers):
	  bloc[i] = 0.;
	  //gives symmetric version of Marshak b.c.
          //(hypre symmetric solvers ignore bloc and do this automatically):
	  //bloc[i] = -0.5 * dx[dir];
	}
	else {
	  cerr << "RadBndry---Not a recognized boundary condition" << endl;
	  exit(1);
	}
      }
      else {
	// internal bndry
	bctag[i] = LO_DIRICHLET;
	bloc[i] = 0.5*delta;
      }
    }
  }
}

void RadBndry::setBndryFluxConds(const BCRec& bc, const BC_Mode phys_bc_mode)
{
  const BoxArray& grids = boxes();
  int ngrds = grids.size();
  const Real* dx = geom.CellSize();
  const Real* xlo = geom.ProbLo();
  const Box& domain = geom.Domain();

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    Array<Real> &bloc = bcloc[face];
    Array<RadBoundCond> &bctag = bcond[face];

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
	  if (p_bc == LO_MARSHAK   || p_bc == LO_SANCHEZ_POMRANING || 
	      p_bc == LO_DIRICHLET || p_bc == LO_NEUMANN) {	      
	    if (p_bcflag == 0) {
	      setValue(face, i, value);
	    }
	    else {
	      Fab& bnd_fab = bndry[face][i];
	      const Box& bnd_box = bnd_fab.box();
	      int iface = face.isLow() ? 0 : 1;
	      FORT_RADBNDRY(bnd_fab.dataPtr(), dimlist(bnd_box),
			    dimlist(domain), dx, xlo, time, dir, iface);
	    }
	  }
	}
	else {
	  if(ParallelDescriptor::IOProcessor()) {
	    cout << "RadBndry ERROR 2: feature not implemented" << endl;
	  }
	  exit(2);

	  Fab& bnd_fab = bndry[face][bi];
	  const Box& bnd_box = bnd_fab.box();
	  BaseFab<int>& tfab = bctypearray[face][i];
#if 0
	  FORT_RADBNDRY2(bnd_fab.dataPtr(), dimlist(bnd_box),
			 tfab.dataPtr(), dimlist(domain), dx, xlo, time);
#endif
	  if (p_bcflag == 0) {
	    // Homogeneous case.  We called RADBNDRY2 only to set tfab right.
	    setValue(face, i, value);
	  }
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
      bnd_fab.setVal(0.0);
    }
  }
}
