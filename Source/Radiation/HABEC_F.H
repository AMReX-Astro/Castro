#ifndef _HABEC_F_H_
#define _HABEC_F_H_

#include "LHH.H"

extern "C" {
  void hacoef(Real* mat, Real* acoefs, intS, intS,
	      const Real& alpha);

  void hbcoef(Real* mat, Real* bcoefs, intS, intS,
	      const Real& beta, const Real* dx, const int& n);

  void hbmat(Real* mat, intS,
	     const int& cdir, const RadBoundCond& bct, const Real& bcl,
	     const int* mask, intS,
	     Real* bcoefs, intS,
	     const Real& beta, const Real* dx);

  void hbmat2(Real* mat, intS,
	      const int& cdir, const int* tf, const Real& bcl, intS,
	      const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx);

  void hbmat3(Real* mat, intS,
	      const int& cdir, const int& bctype,
	      const int* tf, const Real& bcl, intS,
	      const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      const Real* spa, intS);

  void hbvec(Real* vec, intS,
	     const int& cdir, const RadBoundCond& bct, const int& bho,
	     const Real& bcl,
	     const Real* bcval, intS, const int* mask, intS,
	     Real* bcoefs, intS,
	     const Real& beta, const Real* dx);

  void hbvec2(Real* vec, intS,
	      const int& cdir, const int* tf,
	      const int& bho, const Real& bcl,
	      const Real* bcval, intS, const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx);

  void hbvec3(Real* vec, intS,
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl,
	      const Real* bcval, intS, const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx, const Real* r);

  void hbflx(Real* flux, intS,
	     Real* soln, intS, intS,
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const Real& bcl,
	     const Real* bcval, intS, const int* mask, intS,
	     Real* bcoefs, intS,
	     const Real& beta, const Real* dx, const int& inhom);

  void hbflx2(Real* flux, intS,
	      Real* soln, intS, intS,
	      const int& cdir, const int* tf,
	      const int& bho, const Real& bcl,
	      const Real* bcval, intS, const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx, const int& inhom);

  void hbflx3(Real* flux, intS,
	      Real* soln, intS, intS,
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl,
	      const Real* bcval, intS, const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      const int& inhom,
	      const Real* spa, intS);

  void hdterm(Real* dterm, intS,
	      Real* soln, intS, intS,
	      const int& cdir, const RadBoundCond& bct,
	      const Real& bcl,
	      const Real* bcval, intS, const int* mask, intS,
	      Real* dcoefs, intS,
	      const Real* dx);

  void hdterm3(Real* dterm, intS,
	       Real* soln, intS, intS,
	       const int& cdir, const int& bctype, const int* tf,
	       const Real& bcl,
	       const Real* bcval, intS, const int* mask, intS,
	       Real* dcoefs, intS,
	       const Real* dx);

  void hmac(Real* mat, Real* acoefs, intS, intS,
	    const Real& alpha);

  void hmbc(Real* mat, Real* bcoefs, intS, intS,
	    const Real& beta, const Real* dx, const int& n);

  void hma2c(Real* mat, Real* a2coefs, intS, intS,
	     const Real& alpha2, const int& n);

  void hmcc(Real* mat, Real* ccoefs, intS, intS,
	    const Real& gamma, const Real* dx, const int& n);

  void hmd1c(Real* mat, Real* d1coefs, intS, intS,
	     const Real& delta1, const Real* dx, const int& n);

  void hmd2c(Real* mat, Real* d2coefs, intS, intS,
	     const Real& delta2, const Real* dx, const int& n);

  void hmmat(Real* mat, intS,
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const Real& bcl,
	     const int* mask, intS,
	     Real* bcoefs, intS,
	     const Real& beta, const Real* dx);

  void hmmat2(Real* mat, intS,
	      const int& cdir, const int* tf,
	      const int& bho, const Real& bcl, intS,
	      const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx);

  void hmmat3(Real* mat, intS,
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl, intS,
	      const int* mask, intS,
	      Real* bcoefs, intS,
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      Real* spa, intS);

  void add_ccoef_flux(const int& n,
		      Real* Erborder, intS,
		      Real* ccoef, intS,
		      const Real& gamma,
		      const Real* dx,
		      Real* Flux, intS);

  void set_abec_flux(intS, int*,
		     Real*, intS,
		     Real*, intS,
		     const Real*,
		     const Real*,
		     Real*, intS);
};

#endif
