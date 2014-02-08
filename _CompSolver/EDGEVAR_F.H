#ifndef _EDGEVAR_F_H_
#define _EDGEVAR_F_H_

#if defined(BL_LANG_FORT)

#    define FORT_EDGE_AVG_DOWN edge_avg_down
#    define FORT_EDGE_INTERP   edge_interp
#    define FORT_SET_GRADIENT set_gradient
#    define FORT_AVG_DIV avg_div

#else

#include "LHH.H"
#  if defined(BL_FORT_USE_LOWERCASE)
#    define FORT_EDGE_AVG_DOWN edge_avg_down
#    define FORT_EDGE_INTERP   edge_interp
#    define FORT_SET_GRADIENT set_gradient
#    define FORT_AVG_DIV avg_div
#  elif defined (BL_FORT_USE_UPPERCASE)
#    define FORT_EDGE_AVG_DOWN EDGE_AVG_DOWN
#    define FORT_EDGE_INTERP   EDGE_INTERP
#    define FORT_SET_GRADIENT SET_GRADIENT
#    define FORT_AVG_DIV AVG_DIV
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_EDGE_AVG_DOWN edge_avg_down_
#    define FORT_EDGE_INTERP   edge_interp_
#    define FORT_SET_GRADIENT set_gradient_
#    define FORT_AVG_DIV avg_div_
#  endif

   extern "C" {

      void FORT_SET_GRADIENT( const Real*, intS,
			      const Real *,
			      const int *,
#if   (BL_SPACEDIM == 3)
			      const Real*, intS,
			      const Real*, intS,
#elif (BL_SPACEDIM == 2)
			      const Real*, intS,
#endif
			      const Real*, intS) ;

      void FORT_EDGE_AVG_DOWN( const Real*, intS, 
#if   (BL_SPACEDIM == 3)
                               const Real*, intS, 
                               const Real*, intS,
#elif (BL_SPACEDIM == 2)
                               const Real*, intS, 
#endif
                               const int *,
#if   (BL_SPACEDIM == 3)
                               const Real*, intS,
                               const Real*, intS, 
#elif (BL_SPACEDIM == 2)
                               const Real*, intS, 
#endif
                               const Real*, intS) ;

      void FORT_EDGE_INTERP( const Real*, intS, 
#if   (BL_SPACEDIM == 3)
			     const Real*, intS, 
			     const Real*, intS,
#elif (BL_SPACEDIM == 2)
			     const Real*, intS, 
#endif
			     const int *,
#if   (BL_SPACEDIM == 3)
			     const Real*, intS,
			     const Real*, intS, 
#elif (BL_SPACEDIM == 2)
			     const Real*, intS, 
#endif
			     const Real*, intS) ;

      void FORT_AVG_DIV( const Real*, intS,
#if   (BL_SPACEDIM == 3)
                         const Real*, intS,
			 const Real*, intS,
#elif (BL_SPACEDIM == 2)
                         const Real*, intS,
#endif
                         const Real *,
                         const int *,
                         const Real*, intS ) ;

   } ;

#endif

#endif
