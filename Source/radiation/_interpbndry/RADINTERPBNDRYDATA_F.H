#ifndef _RADINTERPBNDRY_F_H_
#define _RADINTERPBNDRY_F_H_

/*************************************************************/

#if defined(RAD_INTERP)

#if defined(BL_LANG_FORT)

#  define FORT_BDINTERPXLO  bdintrpxlo2
#  define FORT_BDINTERPXHI  bdintrpxhi2
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo2
#    define FORT_BDINTERPYHI  bdintrpyhi2
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo2
#    define FORT_BDINTERPZHI  bdintrpzhi2
#  endif

#else

#if defined(BL_FORT_USE_LOWERCASE)
#  define FORT_BDINTERPXLO  bdintrpxlo2
#  define FORT_BDINTERPXHI  bdintrpxhi2
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo2
#    define FORT_BDINTERPYHI  bdintrpyhi2
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo2
#    define FORT_BDINTERPZHI  bdintrpzhi2
#  endif
#elif defined(BL_FORT_USE_UPPERCASE)
#  define FORT_BDINTERPXLO  BDINTRPXLO2
#  define FORT_BDINTERPXHI  BDINTRPXHI2
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  BDINTRPYLO2
#    define FORT_BDINTERPYHI  BDINTRPYHI2
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  BDINTRPZLO2
#    define FORT_BDINTERPZHI  BDINTRPZHI2
#  endif
#elif defined(BL_FORT_USE_UNDERSCORE)
#  define FORT_BDINTERPXLO  bdintrpxlo2_
#  define FORT_BDINTERPXHI  bdintrpxhi2_
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo2_
#    define FORT_BDINTERPYHI  bdintrpyhi2_
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo2_
#    define FORT_BDINTERPZHI  bdintrpzhi2_
#  endif
#endif

#include <AMReX_ArrayLim.H>

typedef void BDInterpFunc(amrex::Real* bndr, ARLIM_P(blo), ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          ARLIM_P(cblo), ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  ARLIM_P(mlo), ARLIM_P(mhi),
                          const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                          amrex::Real* derives);
#ifdef __cplusplus
extern "C"
{
#endif

  BDInterpFunc FORT_BDINTERPXLO;
  BDInterpFunc FORT_BDINTERPXHI;

#if (BL_SPACEDIM > 1)    
    BDInterpFunc FORT_BDINTERPYLO;
    BDInterpFunc FORT_BDINTERPYHI;
#endif    

#if (BL_SPACEDIM > 2)    
    BDInterpFunc FORT_BDINTERPZLO;
    BDInterpFunc FORT_BDINTERPZHI;
#endif
#ifdef __cplusplus
}
#endif

#endif

/*************************************************************/

#else /*RAD_INTERP*/

/*************************************************************/

#if defined(BL_LANG_FORT)

#  define FORT_BDINTERPXLO  bdintrpxlo
#  define FORT_BDINTERPXHI  bdintrpxhi
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo
#    define FORT_BDINTERPYHI  bdintrpyhi
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo
#    define FORT_BDINTERPZHI  bdintrpzhi
#  endif

#else

#if defined(BL_FORT_USE_LOWERCASE)
#  define FORT_BDINTERPXLO  bdintrpxlo
#  define FORT_BDINTERPXHI  bdintrpxhi
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo
#    define FORT_BDINTERPYHI  bdintrpyhi
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo
#    define FORT_BDINTERPZHI  bdintrpzhi
#  endif
#elif defined(BL_FORT_USE_UPPERCASE)
#  define FORT_BDINTERPXLO  BDINTRPXLO
#  define FORT_BDINTERPXHI  BDINTRPXHI
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  BDINTRPYLO
#    define FORT_BDINTERPYHI  BDINTRPYHI
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  BDINTRPZLO
#    define FORT_BDINTERPZHI  BDINTRPZHI
#  endif
#elif defined(BL_FORT_USE_UNDERSCORE)
#  define FORT_BDINTERPXLO  bdintrpxlo_
#  define FORT_BDINTERPXHI  bdintrpxhi_
#  if (BL_SPACEDIM > 1)
#    define FORT_BDINTERPYLO  bdintrpylo_
#    define FORT_BDINTERPYHI  bdintrpyhi_
#  endif
#  if (BL_SPACEDIM > 2)
#    define FORT_BDINTERPZLO  bdintrpzlo_
#    define FORT_BDINTERPZHI  bdintrpzhi_
#  endif
#endif

#include <AMReX_ArrayLim.H>

typedef void BDInterpFunc(amrex::Real* bndr, ARLIM_P(blo), ARLIM_P(bhi),
                          const int* lo, const int* hi,
                          ARLIM_P(cblo), ARLIM_P(cbhi),
                          const int* nvar, const int* ratio,
                          const int* not_covered,
                          const int* mask,  ARLIM_P(mlo), ARLIM_P(mhi),
                          const amrex::Real* crse, ARLIM_P(clo), ARLIM_P(chi),
                          amrex::Real* derives);
#ifdef __cplusplus
extern "C"
{
#endif

  BDInterpFunc FORT_BDINTERPXLO;
  BDInterpFunc FORT_BDINTERPXHI;

#if (BL_SPACEDIM > 1)    
    BDInterpFunc FORT_BDINTERPYLO;
    BDInterpFunc FORT_BDINTERPYHI;
#endif    

#if (BL_SPACEDIM > 2)    
    BDInterpFunc FORT_BDINTERPZLO;
    BDInterpFunc FORT_BDINTERPZHI;
#endif    
#ifdef __cplusplus
}
#endif

#endif

/*************************************************************/

#endif /*RAD_INTERP*/

/*************************************************************/

#endif
