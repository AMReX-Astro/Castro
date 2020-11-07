#ifndef ambient_H
#define ambient_H

#include <AMReX_REAL.H>
#include <state_indices.H>

namespace ambient {
    extern AMREX_GPU_MANAGED amrex::Real ambient_state[NUM_STATE];
}

#endif
