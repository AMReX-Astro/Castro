#ifndef _COMPSOLVER_F_H_
#define _COMPSOLVER_F_H_

#if defined(BL_LANG_FORT)

#   define FORT_EDGE_AVG_DOWN   edgeavgdown
#   define FORT_EDGE_HARM_AVG   edgeharmavg
#   define FORT_AVG_DOWN        avg_down
#   define FORT_SET_ABEC_FLUX   setabecflux
#   define FORT_CONS_INTERP     consintrp
#   define FORT_PC_INTERP       pcintrp
#   define FORT_PC_RESTRICT     pcrestrict
#   define FORT_PL_INTERP       plinterp
#   define FORT_PL_RESTRICT     plrestrict

#else

#  include "LHH.H"
#  if defined(BL_FORT_USE_LOWERCASE)
#    define FORT_EDGE_AVG_DOWN   edgeavgdown
#    define FORT_EDGE_HARM_AVG   edgeharmavg
#    define FORT_AVG_DOWN        avg_down
#    define FORT_SET_ABEC_FLUX   setabecflux
#    define FORT_CONS_INTERP     consintrp
#    define FORT_PC_INTERP       pcintrp
#    define FORT_PC_RESTRICT     pcrestrict
#    define FORT_PL_INTERP       plinterp
#    define FORT_PL_RESTRICT     plrestrict
#  elif defined(BL_FORT_USE_UPPERCASE)
#    define FORT_EDGE_AVG_DOWN  EDGEAVGDOWN
#    define FORT_EDGE_HARM_AVG  EDGEHARMAVG
#    define FORT_AVG_DOWN       AVG_DOWN
#    define FORT_SET_ABEC_FLUX  SETABECFLUX
#    define FORT_CONS_INTERP    CONSINTRP
#    define FORT_PC_INTERP      PCINTRP
#    define FORT_PC_RESTRICT    PCRESTRICT
#    define FORT_PL_INTERP      PLINTERP
#    define FORT_PL_RESTRICT    PLRESTRICT
#  elif defined(BL_FORT_USE_UNDERSCORE)
#    define FORT_EDGE_AVG_DOWN  edgeavgdown_
#    define FORT_EDGE_HARM_AVG  edgeharmavg_
#    define FORT_AVG_DOWN       avg_down_
#    define FORT_SET_ABEC_FLUX  setabecflux_
#    define FORT_CONS_INTERP    consintrp_
#    define FORT_PC_INTERP      pcintrp_
#    define FORT_PC_RESTRICT    pcrestrict_
#    define FORT_PL_INTERP      plinterp_
#    define FORT_PL_RESTRICT    plrestrict_
#  endif

   extern "C" {

     void FORT_EDGE_AVG_DOWN( int*,
			      const Real*, intS, 
			      const int *,
			      const Real*, intS);

     void FORT_EDGE_HARM_AVG( int*,
			      const Real*, intS, 
			      const int *,
			      const Real*, intS);

     void FORT_AVG_DOWN(const Real*, intS,
                        const int* nvar,
			const Real*, intS,
			const int* lo, const int* hi,
			const int* ratio);

     void FORT_SET_ABEC_FLUX(
			     int*,
			     Real*, intS,
			     Real*, intS,
			     const Real*,
			     const Real*,
			     Real*, intS
			     );

    void FORT_CONS_INTERP (
			   Real*, intS,
			   Real*, intS,
			   const int*, const int*,
			   const int*,
			   const int*, const int*
			   );

    void FORT_PC_INTERP (
			 Real*, intS,
			 Real*, intS,
			 const int*, const int*,
			 const int*
			 );

    void FORT_PC_RESTRICT (
			   const Real*, intS,
			   Real*, intS,
			   const int*, const int*,
			   const int*
			   );

    void FORT_PL_INTERP (
			 const Real*, intS,
                         Real*, intS,
                         const int*, const int*,
                         const int*
                         );

    void FORT_PL_RESTRICT (
			   const Real*, intS,
			   Real*, intS,
			   const int*, const int*,
			   const int*
			   );
   } ;

#endif

#endif
