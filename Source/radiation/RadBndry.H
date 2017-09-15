#ifndef _RadBndry_H_
#define _RadBndry_H_

#include "NGBndry.H"

class RadBndry : public NGBndry
{
public:
  RadBndry(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
	   const amrex::Geometry& _geom);
  RadBndry(const amrex::BoxArray& _grids, const amrex::DistributionMapping& _dmap,
	   const amrex::Geometry& _geom, amrex::Real bv);

  ~RadBndry();

  virtual void setBndryConds(const amrex::BCRec& phys_bc,
			     const amrex::Geometry& geom, amrex::IntVect& ratio);

  virtual void setBndryFluxConds(const amrex::BCRec& bc,
		       const BC_Mode phys_bc_mode = Inhomogeneous_BC);

  virtual int mixedBndry(const amrex::Orientation& _face) const {
    return (bcflag[_face] > 1) ? 1 : 0;
  }

  void setHomogValues(const amrex::BCRec& bc, amrex::IntVect& ratio);

  static void setTime(amrex::Real Time) {
    time = Time;
  }
  // Warning: hidden state change, use carefully!
  static void setCorrection() {
    correction = 1;
  }
  static void unsetCorrection() {
    correction = 0;
  }

protected:
  static void init();
  static void init(amrex::Real bv);

  static int first;            // only set up bcval once
  static amrex::Array<int> bcflag;    // 0 = use constant bcval; 1 = use function
  static amrex::Array<amrex::Real> bcval;    // constant bcval across each face
  static amrex::Real time;            // time value passed to bc function
  static int correction;       // ignore bcval when this is set
};

class RadBndryBld : public NGBndryBld
{
public:
  NGBndry* operator()(const amrex::BoxArray& _grids,
		      const amrex::DistributionMapping& _dmap,
		      int _ncomp,
		      const amrex::Geometry& _geom) const {
      return new RadBndry(_grids, _dmap, _geom);
  }
};

#endif
