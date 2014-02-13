// multigroup version of RadBndry.cpp

#include <ParmParse.H>
#include <LO_BCTYPES.H>
#include "MGRadBndry.H"

#include <Using.H>

#include "LHH.H"

#undef BL_USE_ARLIM

#include "RAD_F.H"

int         MGRadBndry::ngroups = 1;
int         MGRadBndry::first = 1;
Array<int>  MGRadBndry::bcflag(2*BL_SPACEDIM);
Array< Array<Real> > MGRadBndry::bcval(2*BL_SPACEDIM);
Real        MGRadBndry::time = 0.0;
int         MGRadBndry::correction = 0;

MGRadBndry::MGRadBndry(const BoxArray& _grids,
                       const int _ngroups,
                       const Geometry& _geom) :
  NGBndry(_grids,_ngroups,_geom)
{
  if (first)
    init(_ngroups);
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
	  const Box& face_box = bndry[face][igrid].box();
	  bctypearray[face].set(igrid, new BaseFab<int>(face_box));
          // We don't care about the bndry values here, only the type array.
#if 0
          FORT_RADBNDRY2(bndry[face][igrid].dataPtr(), dimlist(face_box),
                         bctypearray[face][igrid].dataPtr(),
                         dimlist(domain), dx, xlo, time);
#endif
        }
      }
    }
  }
}

MGRadBndry::~MGRadBndry()
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

void MGRadBndry::init(const int _ngroups)
{
  // obsolete implementation of the Marshak boundary condition requires
  // us to scale incoming fluxes by c here, since the solver classes
  // do not themselves know about c.  This c compensates
  // for the c in the diffusion coefficient at the boundaries:

  ngroups = _ngroups;

  ParmParse pp("radiation");

  Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
  pp.getarr("lo_bc",lo_bc,0,BL_SPACEDIM);
  pp.getarr("hi_bc",hi_bc,0,BL_SPACEDIM);

  Array<int> lo_bcflag(BL_SPACEDIM, 0), hi_bcflag(BL_SPACEDIM, 0);
  pp.queryarr("lo_bcflag",lo_bcflag,0,BL_SPACEDIM);
  pp.queryarr("hi_bcflag",hi_bcflag,0,BL_SPACEDIM);

  Array< Array<Real> > lo_bcval(BL_SPACEDIM), hi_bcval(BL_SPACEDIM);
  lo_bcval[0].resize(ngroups, 0.0);
  hi_bcval[0].resize(ngroups, 0.0);
  pp.queryarr("lo_bcval0", lo_bcval[0], 0, ngroups);
  pp.queryarr("hi_bcval0", hi_bcval[0], 0, ngroups);

#if (BL_SPACEDIM >= 2)
  lo_bcval[1].resize(ngroups, 0.0);
  hi_bcval[1].resize(ngroups, 0.0);
  pp.queryarr("lo_bcval1", lo_bcval[1], 0, ngroups);
  pp.queryarr("hi_bcval1", hi_bcval[1], 0, ngroups);
#endif

#if (BL_SPACEDIM == 3)
  lo_bcval[2].resize(ngroups, 0.0);
  hi_bcval[2].resize(ngroups, 0.0);
  pp.queryarr("lo_bcval2", lo_bcval[2], 0, ngroups);
  pp.queryarr("hi_bcval2", hi_bcval[2], 0, ngroups);
#endif

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    int dir = face.coordDir();
    int bctype;
    bcval[face].resize(ngroups);
    if (face.isLow()) {
      bctype       = lo_bc[dir];
      bcflag[face] = lo_bcflag[dir];
      for(int igroup = 0; igroup < ngroups; igroup++) {
        bcval[face][igroup] = lo_bcval[dir][igroup];
      }
    }
    else {
      bctype       = hi_bc[dir];
      bcflag[face] = hi_bcflag[dir];
      for(int igroup = 0; igroup < ngroups; igroup++) {
        bcval[face][igroup] = hi_bcval[dir][igroup];
      }
    }
  }
}

void MGRadBndry::setBndryConds(const BCRec& bc,
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
	  cerr << "MGRadBndry---Not a recognized boundary condition" << endl;
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

void MGRadBndry::setBndryFluxConds(const BCRec& bc, const BC_Mode phys_bc_mode)
{
  const BoxArray& grids = boxes();
  int ngrds = grids.size();
  const Real* dx = geom.CellSize();
  const Real* xlo = geom.ProbLo();
  const Box& domain = geom.Domain();

  // variables for multigroup
  Array<Real> value_nu;
  value_nu.resize(ngroups);

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    Array<Real> &bloc = bcloc[face];
    Array<RadBoundCond> &bctag = bcond[face];

    int dir = face.coordDir();
    int p_bc = (face.isLow() ? bc.lo(dir) : bc.hi(dir));

    int p_bcflag = (phys_bc_mode == Inhomogeneous_BC && !correction)
      ? bcflag[face] : 0;
    if(phys_bc_mode == Inhomogeneous_BC && !correction) {
      for(int igroup = 0; igroup < ngroups; igroup++) {
	value_nu[igroup] = bcval[face][igroup];
      }
    }
    else {
      for(int igroup = 0; igroup < ngroups; igroup++) {
	value_nu[igroup] = 0.0;
      }
    }

    for (FabSetIter bi(bndry[face]); bi.isValid(); ++bi) {
      int i = bi.index();
      const Box& grd = grids[i];

      if (domain[face] == grd[face] && !geom.isPeriodic(dir)) {
	if (bcflag[face] <= 1) {
	  if (p_bc == LO_MARSHAK   || p_bc == LO_SANCHEZ_POMRANING || 
	      p_bc == LO_DIRICHLET || p_bc == LO_NEUMANN) {
	    if (p_bcflag == 0) {
	      for(int igroup = 0; igroup < ngroups; igroup++) {
		bndry[face][i].setVal(value_nu[igroup], igroup);
	      }
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
	    cout << "MGRadBndry ERROR 2: feature not implemented" << endl;
	  }
	  exit(2);

	  Fab& bnd_fab = bndry[face][i];
	  const Box& bnd_box = bnd_fab.box();
	  BaseFab<int>& tfab = bctypearray[face][i];
#if 0
	  FORT_RADBNDRY2(bnd_fab.dataPtr(), dimlist(bnd_box),
			 tfab.dataPtr(), dimlist(domain), dx, xlo, time);
#endif
	  if (p_bcflag == 0) {
	    // Homogeneous case.  We called RADBNDRY2 only to set tfab right.
	    //setValue(face, i, value);
	  }
	}
      }
    }
  }
}

// *************************************************************************

void MGRadBndry::setHomogValues(const BCRec& bc, IntVect& ratio)
{

  setBndryConds(bc, geom, ratio);

  for (OrientationIter fi; fi; ++fi) {
    Orientation face(fi());
    for (FabSetIter bi(bndry[face]); bi.isValid(); ++bi) {
      int grd = bi.index();
      FArrayBox& bnd_fab = bndry[face][grd];
      bnd_fab.setVal(0.);
    }
  }
}
