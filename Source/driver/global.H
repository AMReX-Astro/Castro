#ifndef GLOBAL_H
#define GLOBAL_H

#include <AMReX_Amr.H>
#ifdef RADIATION
#include <Radiation.H>
#endif

namespace global
{
    AMREX_INLINE Amr* the_amr_ptr;
#ifdef RADIATION
    AMREX_INLINE Radiation* the_radiation_ptr;
#endif
}

#endif
