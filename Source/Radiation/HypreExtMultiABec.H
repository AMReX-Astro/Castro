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
      alpha2(0.0), gamma(0.0), delta1(0.0), delta2(0.0)
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

  void a2Coefficients(int level, const amrex::MultiFab &a2, int dir);
  void cCoefficients( int level, const amrex::MultiFab &c,  int dir);
  void d1Coefficients(int level, const amrex::MultiFab &d1, int dir);
  void d2Coefficients(int level, const amrex::MultiFab &d2, int dir);

  const amrex::MultiFab& a2Coefficients(int level, int dir) {
      return (*a2coefs[level])[dir];
  }
  const amrex::MultiFab& cCoefficients( int level, int dir) {
      return (*ccoefs[level])[dir];
  }
  const amrex::MultiFab& d1Coefficients(int level, int dir) {
      return (*d1coefs[level])[dir];
  }
  const amrex::MultiFab& d2Coefficients(int level, int dir) {
      return (*d2coefs[level])[dir];
  }

  void loadLevelVectors(int level,
                        amrex::MultiFab& dest,
                        int icomp,
                        amrex::MultiFab& rhs, // will not be altered
                        BC_Mode inhom);
  void loadLevelVectorB(int level,
                        amrex::MultiFab& rhs, // will not be altered
                        BC_Mode inhom);

  void loadMatrix();  // once all level coeffs and scalars have been set

  void boundaryDterm(int level,
		     amrex::MultiFab* Dterm,
		     amrex::MultiFab& Er,
		     int icomp);

 protected:

  amrex::Array<std::unique_ptr<amrex::Tuple<amrex::MultiFab, BL_SPACEDIM> > > a2coefs; // face-based
  amrex::Array<std::unique_ptr<amrex::Tuple<amrex::MultiFab, BL_SPACEDIM> > > ccoefs;  // face-based, 2 component
  amrex::Array<std::unique_ptr<amrex::Tuple<amrex::MultiFab, BL_SPACEDIM> > > d1coefs; // cell-based but directional
  amrex::Array<std::unique_ptr<amrex::Tuple<amrex::MultiFab, BL_SPACEDIM> > > d2coefs; // face-based
  Real alpha2, gamma, delta1, delta2; // multipliers for the above
};

#endif
