#ifndef _HABEC_F_H_
#define _HABEC_F_H_


extern "C" {
  void hacoef(amrex::Real* mat, amrex::Real* acoefs, intS, intS,
	      const amrex::Real& alpha);

  void hbcoef(amrex::Real* mat, amrex::Real* bcoefs, intS, intS,
	      const amrex::Real& beta, const amrex::Real* dx, const int& n);

  void hbmat(amrex::Real* mat, intS,
	     const int& cdir, const RadBoundCond& bct, const amrex::Real& bcl,
	     const int* mask, intS,
	     amrex::Real* bcoefs, intS,
	     const amrex::Real& beta, const amrex::Real* dx);

  void hbmat2(amrex::Real* mat, intS,
	      const int& cdir, const int* tf, const amrex::Real& bcl, intS,
	      const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx);

  void hbmat3(amrex::Real* mat, intS,
	      const int& cdir, const int& bctype,
	      const int* tf, const amrex::Real& bcl, intS,
	      const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx,
	      const amrex::Real& flux_factor, const amrex::Real* r,
	      const amrex::Real* spa, intS);

  void hbvec(amrex::Real* vec, intS,
	     const int& cdir, const RadBoundCond& bct, const int& bho,
	     const amrex::Real& bcl,
	     const amrex::Real* bcval, intS, const int* mask, intS,
	     amrex::Real* bcoefs, intS,
	     const amrex::Real& beta, const amrex::Real* dx);

  void hbvec2(amrex::Real* vec, intS,
	      const int& cdir, const int* tf,
	      const int& bho, const amrex::Real& bcl,
	      const amrex::Real* bcval, intS, const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx);

  void hbvec3(amrex::Real* vec, intS,
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const amrex::Real& bcl,
	      const amrex::Real* bcval, intS, const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx, const amrex::Real* r);

  void hbflx(amrex::Real* flux, intS,
	     amrex::Real* soln, intS, intS,
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const amrex::Real& bcl,
	     const amrex::Real* bcval, intS, const int* mask, intS,
	     amrex::Real* bcoefs, intS,
	     const amrex::Real& beta, const amrex::Real* dx, const int& inhom);

  void hbflx2(amrex::Real* flux, intS,
	      amrex::Real* soln, intS, intS,
	      const int& cdir, const int* tf,
	      const int& bho, const amrex::Real& bcl,
	      const amrex::Real* bcval, intS, const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx, const int& inhom);

  void hbflx3(amrex::Real* flux, intS,
	      amrex::Real* soln, intS, intS,
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const amrex::Real& bcl,
	      const amrex::Real* bcval, intS, const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx,
	      const amrex::Real& flux_factor, const amrex::Real* r,
	      const int& inhom,
	      const amrex::Real* spa, intS);

  void hdterm(amrex::Real* dterm, intS,
	      amrex::Real* soln, intS, intS,
	      const int& cdir, const RadBoundCond& bct,
	      const amrex::Real& bcl,
	      const amrex::Real* bcval, intS, const int* mask, intS,
	      amrex::Real* dcoefs, intS,
	      const amrex::Real* dx);

  void hdterm3(amrex::Real* dterm, intS,
	       amrex::Real* soln, intS, intS,
	       const int& cdir, const int& bctype, const int* tf,
	       const amrex::Real& bcl,
	       const amrex::Real* bcval, intS, const int* mask, intS,
	       amrex::Real* dcoefs, intS,
	       const amrex::Real* dx);

  void hmac(amrex::Real* mat, amrex::Real* acoefs, intS, intS,
	    const amrex::Real& alpha);

  void hmbc(amrex::Real* mat, amrex::Real* bcoefs, intS, intS,
	    const amrex::Real& beta, const amrex::Real* dx, const int& n);

  void hma2c(amrex::Real* mat, amrex::Real* a2coefs, intS, intS,
	     const amrex::Real& alpha2, const int& n);

  void hmcc(amrex::Real* mat, amrex::Real* ccoefs, intS, intS,
	    const amrex::Real& gamma, const amrex::Real* dx, const int& n);

  void hmd1c(amrex::Real* mat, amrex::Real* d1coefs, intS, intS,
	     const amrex::Real& delta1, const amrex::Real* dx, const int& n);

  void hmd2c(amrex::Real* mat, amrex::Real* d2coefs, intS, intS,
	     const amrex::Real& delta2, const amrex::Real* dx, const int& n);

  void hmmat(amrex::Real* mat, intS,
	     const int& cdir, const RadBoundCond& bct,
	     const int& bho, const amrex::Real& bcl,
	     const int* mask, intS,
	     amrex::Real* bcoefs, intS,
	     const amrex::Real& beta, const amrex::Real* dx);

  void hmmat2(amrex::Real* mat, intS,
	      const int& cdir, const int* tf,
	      const int& bho, const amrex::Real& bcl, intS,
	      const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx);

  void hmmat3(amrex::Real* mat, intS,
	      const int& cdir, const int& bctype, const int* tf,
	      const int& bho, const amrex::Real& bcl, intS,
	      const int* mask, intS,
	      amrex::Real* bcoefs, intS,
	      const amrex::Real& beta, const amrex::Real* dx,
	      const amrex::Real& flux_factor, const amrex::Real* r,
	      amrex::Real* spa, intS);

  void add_ccoef_flux(const int& n,
		      amrex::Real* Erborder, intS,
		      amrex::Real* ccoef, intS,
		      const amrex::Real& gamma,
		      const amrex::Real* dx,
		      amrex::Real* Flux, intS);

  void set_abec_flux(intS, int*,
		     amrex::Real*, intS,
		     amrex::Real*, intS,
		     const amrex::Real*,
		     const amrex::Real*,
		     amrex::Real*, intS);
};

#endif
