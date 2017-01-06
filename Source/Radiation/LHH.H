#ifndef _LHH_H_
#define _LHH_H_

typedef IntVect   Iv;
typedef FArrayBox Fab;

#if (BL_SPACEDIM == 1)
#  define intS  const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.bigEnd(0)
#elif (BL_SPACEDIM == 2)
#  define intS  const int&, const int&, const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.smallEnd(1), a.bigEnd(0), a.bigEnd(1)
#elif (BL_SPACEDIM == 3)
#  define intS  const int&, const int&, const int&, const int&, \
                const int&, const int&
#  define dimlist(a) a.smallEnd(0), a.smallEnd(1), a.smallEnd(2), \
                     a.bigEnd(0),   a.bigEnd(1),   a.bigEnd(2)
#endif


#endif
