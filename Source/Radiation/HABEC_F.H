#ifndef _HABEC_F_H_
#define _HABEC_F_H_

#include "LHH.H"
#include <ArrayLim.H>

extern "C" {
  void hacoef(Real* mat, Real* acoefs, ARLIM_P(alo), ARLIM_P(ahi), 
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const Real& alpha);

  void hbcoef(Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi), 
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const Real& beta, const Real* dx, const int& n);

  void hbmat(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	     const int& cdir, const RadBoundCond& bct, const Real& bcl,
	     const int* mask, ARLIM_P(mlo),ARLIM_P(mhi),
	     Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	     const Real& beta, const Real* dx);

  // void hbmat2(Real* mat, intS,
  // 	      const int& cdir, const int* tf, const Real& bcl, intS,
  // 	      const int* mask, intS,
  // 	      Real* bcoefs, intS,
  // 	      const Real& beta, const Real* dx);

  void hbmat3(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype,
	      const int* tf, const Real& bcl, ARLIM_P(flo), ARLIM_P(fhi),
	      const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	      Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      const Real* spa, ARLIM_P(slo), ARLIM_P(shi));

  void hbvec(Real* vec, ARLIM_P(reglo), ARLIM_P(reghi),
	     const int& cdir, const RadBoundCond& bct, const int& bho,
	     const Real& bcl,
	     const Real* bcval, ARLIM_P(flo), ARLIM_P(fhi), 
	     const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	     Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	     const Real& beta, const Real* dx);

  // void hbvec2(Real* vec, intS,
  // 	      const int& cdir, const int* tf,
  // 	      const int& bho, const Real& bcl,
  // 	      const Real* bcval, intS, const int* mask, intS,
  // 	      Real* bcoefs, intS,
  // 	      const Real& beta, const Real* dx);

  void hbvec3(Real* vec, ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl,
	      const Real* bcval, ARLIM_P(flo), ARLIM_P(fhi), 
	      const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	      Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	      const Real& beta, const Real* dx, const Real* r);

  void hbflx(Real* flux, ARLIM_P(flo), ARLIM_P(fhi),
	     Real* soln, ARLIM_P(slo), ARLIM_P(shi),
	     ARLIM_P(reglo), ARLIM_P(reghi), 
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const Real& bcl,
	     const Real* bcval, ARLIM_P(fslo), ARLIM_P(fshi), 
	     const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	     Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	     const Real& beta, const Real* dx, const int& inhom);

  // void hbflx2(Real* flux, intS,
  // 	      Real* soln, intS, intS,
  // 	      const int& cdir, const int* tf,
  // 	      const int& bho, const Real& bcl,
  // 	      const Real* bcval, intS, const int* mask, intS,
  // 	      Real* bcoefs, intS,
  // 	      const Real& beta, const Real* dx, const int& inhom);

  void hbflx3(Real* flux, ARLIM_P(flo), ARLIM_P(fhi),
	      Real* soln, ARLIM_P(slo), ARLIM_P(shi),
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl,
	      const Real* bcval, ARLIM_P(fslo), ARLIM_P(fshi),
	      const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	      Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      const int& inhom,
	      const Real* spa, ARLIM_P(splo), ARLIM_P(sphi));

  void hdterm(Real* dterm, ARLIM_P(dtlo), ARLIM_P(dthi),
	      Real* soln, ARLIM_P(slo), ARLIM_P(shi), 
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const RadBoundCond& bct,
	      const Real& bcl,
	      const Real* bcval, ARLIM_P(bclo), ARLIM_P(bchi), 
	      const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	      Real* dcoefs, ARLIM_P(dclo), ARLIM_P(dchi),
	      const Real* dx);

  void hdterm3(Real* dterm, ARLIM_P(dtlo), ARLIM_P(dthi),
	       Real* soln, ARLIM_P(slo), ARLIM_P(shi),
	       ARLIM_P(reglo), ARLIM_P(reghi),
	       const int& cdir, const int& bctype, const int* tf,
	       const Real& bcl,
	       const Real* bcval, ARLIM_P(bclo), ARLIM_P(bchi), 
	       const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	       Real* dcoefs, ARLIM_P(dclo), ARLIM_P(dchi),
	       const Real* dx);

  void hmac(Real* mat, Real* acoefs, ARLIM_P(alo), ARLIM_P(ahi),
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    const Real& alpha);

  void hmbc(Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    const Real& beta, const Real* dx, const int& n);

  void hma2c(Real* mat, Real* a2coefs, ARLIM_P(alo), ARLIM_P(ahi),
	     ARLIM_P(reglo), ARLIM_P(reghi), 
	     const Real& alpha2, const int& n);

  void hmcc(Real* mat, Real* ccoefs, ARLIM_P(clo), ARLIM_P(chi), 
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    const Real& gamma, const Real* dx, const int& n);

  void hmd1c(Real* mat, Real* d1coefs, ARLIM_P(d1lo), ARLIM_P(d1hi),
	     ARLIM_P(reglo), ARLIM_P(reghi),
	     const Real& delta1, const Real* dx, const int& n);

  void hmd2c(Real* mat, Real* d2coefs, ARLIM_P(d2lo), ARLIM_P(d2hi),
	     ARLIM_P(reglo), ARLIM_P(reghi),
	     const Real& delta2, const Real* dx, const int& n);

  void hmmat(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const Real& bcl,
	     const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	     Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	     const Real& beta, const Real* dx);

  // void hmmat2(Real* mat, intS,
  // 	      const int& cdir, const int* tf,
  // 	      const int& bho, const Real& bcl, intS,
  // 	      const int* mask, intS,
  // 	      Real* bcoefs, intS,
  // 	      const Real& beta, const Real* dx);

  void hmmat3(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl, ARLIM_P(fslo), ARLIM_P(fshi),
	      const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	      Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      Real* spa, ARLIM_P(splo), ARLIM_P(sphi));

  // not used?
  void add_ccoef_flux(const int& n,
		      Real* Erborder, ARLIM_P(erlo), ARLIM_P(erhi),
		      Real* ccoef, ARLIM_P(cclo), ARLIM_P(cchi),
		      const Real& gamma,
		      const Real* dx,
		      Real* Flux, ARLIM_P(flo), ARLIM_P(fhi));

  void set_abec_flux(ARLIM_P(reglo), ARLIM_P(reghi), 
		     int*,
		     Real*, ARLIM_P(dlo), ARLIM_P(dhi),
		     Real*, ARLIM_P(dclo), ARLIM_P(dchi),
		     const Real*,
		     const Real*,
		     Real*, ARLIM_P(flo), ARLIM_P(fhi));
};

#endif
