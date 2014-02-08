#ifndef _MGRadBndry_H_
#define _MGRadBndry_H_

#include "NGBndry.H"

class MGRadBndry : public NGBndry
{
public:
  MGRadBndry(const BoxArray& _grids,
             const int _ngroups,
             const Geometry& _geom);

  ~MGRadBndry();

  virtual void setBndryConds(const BCRec& phys_bc,
			     const Geometry& geom, IntVect& ratio);

  virtual void setBndryFluxConds(const BCRec& bc,
		       const BC_Mode phys_bc_mode = Inhomogeneous_BC);

  virtual int mixedBndry(const Orientation& _face) const {
    return (bcflag[_face] > 1) ? 1 : 0;
  }

  void setHomogValues(const BCRec& bc, IntVect& ratio);

  static void setTime(Real Time) {
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
  static void init(const int _ngroups);
  static int ngroups;

  static int first;            // only set up bcval once
  static Array<int> bcflag;    // 0 = use constant bcval; 1 = use function
  //  static Array<Real> bcval;    // constant bcval across each face
  static Array< Array<Real> > bcval;
  static Real time;            // time value passed to bc function
  static Real c;               // speed of light
  static int correction;       // ignore bcval when this is set
};

class MGRadBndryBld : public NGBndryBld
{
public:
  NGBndry* operator()(const BoxArray& _grids, int _ngroups,
		      const Geometry& _geom) const {
    return new MGRadBndry(_grids, _ngroups, _geom);
  }
};

#endif
