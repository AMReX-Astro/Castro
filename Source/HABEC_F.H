#ifndef _HABEC_F_H_
#define _HABEC_F_H_

#ifdef BL_FORT_USE_UNDERSCORE
#  define   FORT_HACOEF  hacoef_
#  define   FORT_HBCOEF  hbcoef_
#  define   FORT_HBMAT   hbmat_
#  define   FORT_HBMAT2  hbmat2_
#  define   FORT_HBMAT3  hbmat3_
#  define   FORT_HBVEC   hbvec_
#  define   FORT_HBVEC2  hbvec2_
#  define   FORT_HBVEC3  hbvec3_
#  define   FORT_HBFLX   hbflx_
#  define   FORT_HBFLX2  hbflx2_
#  define   FORT_HBFLX3  hbflx3_
#  define   FORT_HDTERM   hdterm_
#  define   FORT_HDTERM3   hdterm3_
#  define   FORT_HMAC    hmac_
#  define   FORT_HMBC    hmbc_
#  define   FORT_HMA2C   hma2c_
#  define   FORT_HMCC    hmcc_
#  define   FORT_HMD1C   hmd1c_
#  define   FORT_HMD2C   hmd2c_
#  define   FORT_HMMAT   hmmat_
#  define   FORT_HMMAT2  hmmat2_
#  define   FORT_HMMAT3  hmmat3_
#  define   FORT_ADD_CCOEF_FLUX   add_ccoef_flux_
#  define   FORT_SET_ABEC_FLUX    setabecflux_
#elif defined(BL_FORT_USE_UPPERCASE)
#  define   FORT_HACOEF  HACOEF
#  define   FORT_HBCOEF  HBCOEF
#  define   FORT_HBMAT   HBMAT
#  define   FORT_HBMAT2  HBMAT2
#  define   FORT_HBMAT3  HBMAT3
#  define   FORT_HBVEC   HBVEC
#  define   FORT_HBVEC2  HBVEC2
#  define   FORT_HBVEC3  HBVEC3
#  define   FORT_HBFLX   HBFLX
#  define   FORT_HBFLX2  HBFLX2
#  define   FORT_HBFLX3  HBFLX3
#  define   FORT_HDTERM   HDTERM
#  define   FORT_HDTERM3   HDTERM3
#  define   FORT_HMAC    HMAC
#  define   FORT_HMBC    HMBC
#  define   FORT_HMA2C   HMA2C
#  define   FORT_HMCC    HMCC
#  define   FORT_HMD1C   HMD1C
#  define   FORT_HMD2C   HMD2C
#  define   FORT_HMMAT   HMMAT
#  define   FORT_HMMAT2  HMMAT2
#  define   FORT_HMMAT3  HMMAT3
#  define   FORT_ADD_CCOEF_FLUX   ADD_CCOEF_FLUX
#  define   FORT_SET_ABEC_FLUX    SETABECFLUX
#  endif

extern "C" {
  void FORT_HACOEF(Real* mat, Real* acoefs, intS, intS,
		   const Real& alpha);

  void FORT_HBCOEF(Real* mat, Real* bcoefs, intS, intS,
		   const Real& beta, const Real* dx, const int& n);

  void FORT_HBMAT(Real* mat, intS,
		  const int& cdir, const RadBoundCond& bct, const Real& bcl,
		  const int* mask, intS,
		  Real* bcoefs, intS,
		  const Real& beta, const Real* dx);

  void FORT_HBMAT2(Real* mat, intS,
		   const int& cdir, const int* tf, const Real& bcl, intS,
		   const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx);

  void FORT_HBMAT3(Real* mat, intS,
		   const int& cdir, const int& bctype,
                   const int* tf, const Real& bcl, intS,
		   const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx,
                   const Real& flux_factor, const Real* r,
		   const Real* spa, intS);

  void FORT_HBVEC(Real* vec, intS,
		  const int& cdir, const RadBoundCond& bct, const int& bho,
                  const Real& bcl,
		  const Real* bcval, intS, const int* mask, intS,
		  Real* bcoefs, intS,
		  const Real& beta, const Real* dx);

  void FORT_HBVEC2(Real* vec, intS,
		   const int& cdir, const int* tf,
                   const int& bho, const Real& bcl,
		   const Real* bcval, intS, const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx);

  void FORT_HBVEC3(Real* vec, intS,
		   const int& cdir, const int& bctype, const int* tf,
                   const int& bho, const Real& bcl,
		   const Real* bcval, intS, const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx, const Real* r);

  void FORT_HBFLX(Real* flux, intS,
		  Real* soln, intS, intS,
		  const int& cdir, const RadBoundCond& bct,
                  const int& bho, const Real& bcl,
		  const Real* bcval, intS, const int* mask, intS,
		  Real* bcoefs, intS,
		  const Real& beta, const Real* dx, const int& inhom);

  void FORT_HBFLX2(Real* flux, intS,
		   Real* soln, intS, intS,
		   const int& cdir, const int* tf,
                   const int& bho, const Real& bcl,
		   const Real* bcval, intS, const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx, const int& inhom);

  void FORT_HBFLX3(Real* flux, intS,
		   Real* soln, intS, intS,
		   const int& cdir, const int& bctype, const int* tf,
                   const int& bho, const Real& bcl,
		   const Real* bcval, intS, const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx,
                   const Real& flux_factor, const Real* r,
                   const int& inhom,
		   const Real* spa, intS);

  void FORT_HDTERM(Real* dterm, intS,
		   Real* soln, intS, intS,
		   const int& cdir, const RadBoundCond& bct,
		   const Real& bcl,
		   const Real* bcval, intS, const int* mask, intS,
		   Real* dcoefs, intS,
		   const Real* dx);

  void FORT_HDTERM3(Real* dterm, intS,
		    Real* soln, intS, intS,
		    const int& cdir, const int& bctype, const int* tf,
		    const Real& bcl,
		    const Real* bcval, intS, const int* mask, intS,
		    Real* dcoefs, intS,
		    const Real* dx);

  void FORT_HMAC(Real* mat, Real* acoefs, intS, intS,
		 const Real& alpha);

  void FORT_HMBC(Real* mat, Real* bcoefs, intS, intS,
		 const Real& beta, const Real* dx, const int& n);

  void FORT_HMA2C(Real* mat, Real* a2coefs, intS, intS,
                  const Real& alpha2, const int& n);

  void FORT_HMCC(Real* mat, Real* ccoefs, intS, intS,
		 const Real& gamma, const Real* dx, const int& n);

  void FORT_HMD1C(Real* mat, Real* d1coefs, intS, intS,
                  const Real& delta1, const Real* dx, const int& n);

  void FORT_HMD2C(Real* mat, Real* d2coefs, intS, intS,
                  const Real& delta2, const Real* dx, const int& n);

  void FORT_HMMAT(Real* mat, intS,
		  const int& cdir, const RadBoundCond& bct,
		  const int& bho, const Real& bcl,
                  const int* mask, intS,
		  Real* bcoefs, intS,
		  const Real& beta, const Real* dx);

  void FORT_HMMAT2(Real* mat, intS,
		   const int& cdir, const int* tf,
                   const int& bho, const Real& bcl, intS,
		   const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx);

  void FORT_HMMAT3(Real* mat, intS,
		   const int& cdir, const int& bctype, const int* tf,
                   const int& bho, const Real& bcl, intS,
		   const int* mask, intS,
		   Real* bcoefs, intS,
		   const Real& beta, const Real* dx,
                   const Real& flux_factor, const Real* r,
		   Real* spa, intS);

  void FORT_ADD_CCOEF_FLUX(const int& n,
                           Real* Erborder, intS,
                           Real* ccoef, intS,
                           const Real& gamma,
                           const Real* dx,
                           Real* Flux, intS);

  void FORT_SET_ABEC_FLUX(intS, int*,
			  Real*, intS,
			  Real*, intS,
			  const Real*,
			  const Real*,
			  Real*, intS);
};

#endif
