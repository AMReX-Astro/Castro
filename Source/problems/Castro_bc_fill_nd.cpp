
#include <AMReX_BLFort.H>
#include <Castro_bc_fill_nd.H>
#include <Castro_bc_fill_nd_F.H>

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    void ca_hypfill(Real* adv, const int* adv_lo, const int* adv_hi,
                    const int* domlo, const int* domhi, const Real* dx, const Real* xlo,
                    const Real* time, const int* bc)
    {
        int lo[3];
        int hi[3];

        for (int i = 0; i < 3; ++i) {
            lo[i] = adv_lo[i];
            hi[i] = adv_hi[i];
        }

        hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, dx, xlo, time, bc);
    }

#ifdef __cplusplus
}
#endif
