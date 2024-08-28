#ifndef CASTRO_HYPREEXTMULTIABEC_H
#define CASTRO_HYPREEXTMULTIABEC_H

#include <HypreMultiABec.H>


///
/// @class HypreExtMultiABec
/// @brief
///
class HypreExtMultiABec : public HypreMultiABec {

 public:


///
/// @param _crse_level
/// @param _fine_level
/// @param _solver_flag
///
  HypreExtMultiABec(int _crse_level, int _fine_level, int _solver_flag)
    : HypreMultiABec(_crse_level, _fine_level, _solver_flag),
      a2coefs(fine_level+1),
      ccoefs( fine_level+1),
      d1coefs(fine_level+1),
      d2coefs(fine_level+1)
  {}

  ~HypreExtMultiABec();

  HypreExtMultiABec(const HypreExtMultiABec& src) = delete;
  HypreExtMultiABec(const HypreExtMultiABec&& src) = delete;

  HypreExtMultiABec& operator= (const HypreExtMultiABec& src) = delete;
  HypreExtMultiABec& operator= (const HypreExtMultiABec&& src) = delete;

  amrex::Real& a2Multiplier() {
    return alpha2;
  }
  amrex::Real& cMultiplier() {
    return gamma;
  }
  amrex::Real& d1Multiplier() {
    return delta1;
  }
  amrex::Real& d2Multiplier() {
    return delta2;
  }


///
/// @param level
/// @param &a2
/// @param dir
///
  void a2Coefficients(int level, const amrex::MultiFab &a2, int dir);

///
/// @param level
/// @param &c
/// @param dir
///
  void cCoefficients( int level, const amrex::MultiFab &c,  int dir);

///
/// @param level
/// @param &d1
/// @param dir
///
  void d1Coefficients(int level, const amrex::MultiFab &d1, int dir);

///
/// @param level
/// @param &d2
/// @param dir
///
  void d2Coefficients(int level, const amrex::MultiFab &d2, int dir);


///
/// @param level
/// @param dir
///
  const amrex::MultiFab& a2Coefficients(int level, int dir) {
      return (*a2coefs[level])[dir];
  }

///
/// @param level
/// @param dir
///
  const amrex::MultiFab& cCoefficients( int level, int dir) {
      return (*ccoefs[level])[dir];
  }

///
/// @param level
/// @param dir
///
  const amrex::MultiFab& d1Coefficients(int level, int dir) {
      return (*d1coefs[level])[dir];
  }

///
/// @param level
/// @param dir
///
  const amrex::MultiFab& d2Coefficients(int level, int dir) {
      return (*d2coefs[level])[dir];
  }

  void loadLevelVectors(int level,
                        amrex::MultiFab& dest,
                        int icomp,
                        amrex::MultiFab& rhs, ///< will not be altered
                        BC_Mode inhom);
  void loadLevelVectorB(int level,
                        amrex::MultiFab& rhs, ///< will not be altered
                        BC_Mode inhom);

  void loadMatrix();  ///< once all level coeffs and scalars have been set


///
/// @param level
/// @param Dterm
/// @param Er
/// @param icomp
///
  void boundaryDterm(int level,
                     amrex::MultiFab* Dterm,
                     amrex::MultiFab& Er,
                     int icomp);

 protected:

  amrex::Vector<std::unique_ptr<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> > > a2coefs; ///< face-based
  amrex::Vector<std::unique_ptr<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> > > ccoefs;  ///< face-based, 2 component
  amrex::Vector<std::unique_ptr<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> > > d1coefs; ///< cell-based but directional
  amrex::Vector<std::unique_ptr<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> > > d2coefs; ///< face-based
  amrex::Real alpha2{}, gamma{}, delta1{}, delta2{}; ///< multipliers for the above
};

#endif
