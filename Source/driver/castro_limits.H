#ifndef CASTRO_LIMITS_H
#define CASTRO_LIMITS_H

constexpr int MAX_LEV{15};

// NUM_GROW is the number of ghost cells needed for the hyperbolic
// portions -- note that this includes the flattening, which
// generally requires 4 ghost cells
#ifdef MHD
constexpr int NUM_GROW{6};
#else
constexpr int NUM_GROW{4};
#endif

// NUM_GROW_SRC is for quantities that will be reconstructed, but
// don't need the full stencil required for flattening
#ifdef MHD
constexpr int NUM_GROW_SRC{6};
#elif TRUE_SDC
constexpr int NUM_GROW_SRC{NUM_GROW};
#else
constexpr int NUM_GROW_SRC{4};
#endif

#endif
