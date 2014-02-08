#ifndef _HypreExtMultiABec_H_
#define _HypreExtMultiABec_H_

#include "HypreMultiABec.H"


class HypreExtMultiABec : public HypreMultiABec {

 public:

  HypreExtMultiABec(int _crse_level, int _fine_level, int _solver_flag)
    : HypreMultiABec(_crse_level, _fine_level, _solver_flag),
      a2coefs(fine_level+1),
      ccoefs( fine_level+1),
      d1coefs(fine_level+1),
      d2coefs(fine_level+1),
      alpha2(0.0), gamma(0.0), delta1(0.0), delta2(0.0),
      boundary_Xi(fine_level+1)
  {
  }

  ~HypreExtMultiABec();

  Real& a2Multiplier() {
    return alpha2;
  }
  Real& cMultiplier() {
    return gamma;
  }
  Real& d1Multiplier() {
    return delta1;
  }
  Real& d2Multiplier() {
    return delta2;
  }

  void a2Coefficients(int level, const MultiFab &a2, int dir);
  void cCoefficients( int level, const MultiFab &c,  int dir);
  void d1Coefficients(int level, const MultiFab &d1, int dir);
  void d2Coefficients(int level, const MultiFab &d2, int dir);

  const MultiFab& a2Coefficients(int level, int dir) {
    return a2coefs[level][dir];
  }
  const MultiFab& cCoefficients( int level, int dir) {
    return ccoefs[level][dir];
  }
  const MultiFab& d1Coefficients(int level, int dir) {
    return d1coefs[level][dir];
  }
  const MultiFab& d2Coefficients(int level, int dir) {
    return d2coefs[level][dir];
  }

  void boundaryXi(int level, const MultiFab &Xi, int dir);

  void loadLevelVectors(int level,
                        MultiFab& dest,
                        int icomp,
                        MultiFab& rhs, // will not be altered
                        BC_Mode inhom);
  void loadLevelVectorB(int level,
                        MultiFab& rhs, // will not be altered
                        BC_Mode inhom);

  void loadMatrix();  // once all level coeffs and scalars have been set

  void boundaryDterm(int level,
		     MultiFab* Dterm,
		     MultiFab& Er,
		     int icomp);

 protected:

  PArray< Tuple<MultiFab, BL_SPACEDIM> > a2coefs; // face-based
  PArray< Tuple<MultiFab, BL_SPACEDIM> > ccoefs;  // face-based, 2 component
  PArray< Tuple<MultiFab, BL_SPACEDIM> > d1coefs; // cell-based but directional
  PArray< Tuple<MultiFab, BL_SPACEDIM> > d2coefs; // face-based
  Real alpha2, gamma, delta1, delta2; // multipliers for the above

  // the following used at Neumann and Marshak boundaries:
  PArray< Tuple<PArray<Fab>, 2 * BL_SPACEDIM> > boundary_Xi;
};

#endif
