#ifndef _HABEC_F_H_
#define _HABEC_F_H_

#include "LHH.H"
#include <ArrayLim.H>
#include <BLFort.H>

extern "C" {
  void hacoef(Real* mat, 
	      BL_FORT_FAB_ARG(acoefs),
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const Real& alpha);

  void hbcoef(Real* mat, 
	      BL_FORT_FAB_ARG(bcoefs),
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const Real& beta, const Real* dx, const int& n);

  void hbmat(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	     const int& cdir, const RadBoundCond& bct, const Real& bcl,
	     const BL_FORT_IFAB_ARG(mask),
	     BL_FORT_FAB_ARG(bcoefs),
	     const Real& beta, const Real* dx);

  void hbmat3(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype,
	      const int* tf, const Real& bcl, ARLIM_P(flo), ARLIM_P(fhi),
	      const BL_FORT_IFAB_ARG(mask),
	      BL_FORT_FAB_ARG(bcoefs),
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      const Real* spa, ARLIM_P(slo), ARLIM_P(shi));

  void hbvec(Real* vec, ARLIM_P(reglo), ARLIM_P(reghi),
	     const int& cdir, const RadBoundCond& bct, const int& bho,
	     const Real& bcl,
	     const BL_FORT_FAB_ARG(bcval),
	     const BL_FORT_IFAB_ARG(mask),
	     BL_FORT_FAB_ARG(bcoefs),
	     const Real& beta, const Real* dx);

  void hbvec3(Real* vec, ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl,
	      const BL_FORT_FAB_ARG(bcval),
	      const BL_FORT_IFAB_ARG(mask),
	      BL_FORT_FAB_ARG(bcoefs),
	      const Real& beta, const Real* dx, const Real* r);

  void hbflx(BL_FORT_FAB_ARG(flux),
	     BL_FORT_FAB_ARG(soln),
	     ARLIM_P(reglo), ARLIM_P(reghi), 
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const Real& bcl,
	     const Real* bcval, ARLIM_P(fslo), ARLIM_P(fshi), 
	     const BL_FORT_IFAB_ARG(mask),
	     BL_FORT_FAB_ARG(bcoefs),
	     const Real& beta, const Real* dx, const int& inhom);

  void hbflx3(BL_FORT_FAB_ARG(flux),
	      BL_FORT_FAB_ARG(soln),
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const Real& bcl,
	      const BL_FORT_FAB_ARG(bcval),
	      const BL_FORT_IFAB_ARG(mask),
	      BL_FORT_FAB_ARG(bcoefs),
	      const Real& beta, const Real* dx,
	      const Real& flux_factor, const Real* r,
	      const int& inhom,
	      const Real* spa, ARLIM_P(splo), ARLIM_P(sphi));

  void hdterm(BL_FORT_FAB_ARG(dterm),
	      BL_FORT_FAB_ARG(soln),
	      ARLIM_P(reglo), ARLIM_P(reghi),
	      const int& cdir, const RadBoundCond& bct,
	      const Real& bcl,
	      const BL_FORT_FAB_ARG(bcval),
	      const BL_FORT_IFAB_ARG(mask),
	      BL_FORT_FAB_ARG(dcoefs),
	      const Real* dx);

  void hdterm3(BL_FORT_FAB_ARG(dterm),
	       BL_FORT_FAB_ARG(soln),
	       ARLIM_P(reglo), ARLIM_P(reghi),
	       const int& cdir, const int& bctype, const int* tf,
	       const Real& bcl,
	       const BL_FORT_FAB_ARG(bcval),
	       const BL_FORT_IFAB_ARG(mask),
	       BL_FORT_FAB_ARG(dcoefs),
	       const Real* dx);

  void hmac(Real* mat, Real* acoefs, ARLIM_P(alo), ARLIM_P(ahi),
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    const Real& alpha);

  void hmbc(Real* mat, Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    const Real& beta, const Real* dx, const int& n);

  void hma2c(Real* mat, 
	     BL_FORT_FAB_ARG(a2coefs), 
	     ARLIM_P(reglo), ARLIM_P(reghi), 
	     const Real& alpha2, const int& n);

  void hmcc(Real* mat, 
	    BL_FORT_FAB_ARG(ccoefs),
	    ARLIM_P(reglo), ARLIM_P(reghi),
	    const Real& gamma, const Real* dx, const int& n);

  void hmd1c(Real* mat,
	     BL_FORT_FAB_ARG(d1coefs),
	     ARLIM_P(reglo), ARLIM_P(reghi),
	     const Real& delta1, const Real* dx, const int& n);

  void hmd2c(Real* mat,
	     BL_FORT_FAB_ARG(d2coefs),
	     ARLIM_P(reglo), ARLIM_P(reghi),
	     const Real& delta2, const Real* dx, const int& n);

  void hmmat(Real* mat, ARLIM_P(reglo), ARLIM_P(reghi),
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const Real& bcl,
	     const int* mask, ARLIM_P(mlo), ARLIM_P(mhi),
	     Real* bcoefs, ARLIM_P(blo), ARLIM_P(bhi),
	     const Real& beta, const Real* dx);

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
		     BL_FORT_FAB_ARG(c), 
		     BL_FORT_FAB_ARG(d), 
		     const Real*,
		     const Real*,
		     BL_FORT_FAB_ARG(flux));
};

#endif
