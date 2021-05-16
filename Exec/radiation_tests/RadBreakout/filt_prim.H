#ifndef filt_prim_H
#define filt_prim_H

#include <AMReX.H>
#include <filter.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void filt_prim (int i, int j, int k,
                Array4<Real> const& Stmp,
                Array4<Real> const& Snew,
                int filt_T, int S,
                const GeometryData& geomdata,
                Real time)
{
#ifndef AMREX_USE_GPU
    // This filtering is only implemented in 1D today
    AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == 1);
#endif

    using namespace filter;

    bool filtered = false;

    if (time < problem::filter_timemax || Snew(i,j,k,URHO) < problem::filter_rhomax) {
        if (filt_T == 1) {
            filtered = true;
            Snew(i,j,k,URHO ) = ff1(0) *  Stmp(i  ,j,k,URHO ) +
                                ff1(1) * (Stmp(i-1,j,k,URHO ) + Stmp(i+1,j,k,URHO ));
            Snew(i,j,k,UMX  ) = ff1(0) *  Stmp(i  ,j,k,UMX  ) +
                                ff1(1) * (Stmp(i-1,j,k,UMX  ) + Stmp(i+1,j,k,UMX  ));
            Snew(i,j,k,UTEMP) = ff1(0) *  Stmp(i  ,j,k,UTEMP) +
                                ff1(1) * (Stmp(i-1,j,k,UTEMP) + Stmp(i+1,j,k,UTEMP));
        }
        else if (filt_T == 2) {
            filtered = true;
            Snew(i,j,k,URHO ) = ff2(0,S) *  Stmp(i  ,j,k,URHO ) +
                                ff2(1,S) * (Stmp(i-1,j,k,URHO ) + Stmp(i+1,j,k,URHO )) +
                                ff2(2,S) * (Stmp(i-2,j,k,URHO ) + Stmp(i+2,j,k,URHO ));
            Snew(i,j,k,UMX  ) = ff2(0,S) *  Stmp(i  ,j,k,UMX  ) +
                                ff2(1,S) * (Stmp(i-1,j,k,UMX  ) + Stmp(i+1,j,k,UMX  )) +
                                ff2(2,S) * (Stmp(i-2,j,k,UMX  ) + Stmp(i+2,j,k,UMX  ));
            Snew(i,j,k,UTEMP) = ff2(0,S) *  Stmp(i  ,j,k,UTEMP) +
                                ff2(1,S) * (Stmp(i-1,j,k,UTEMP) + Stmp(i+1,j,k,UTEMP)) +
                                ff2(2,S) * (Stmp(i-2,j,k,UTEMP) + Stmp(i+2,j,k,UTEMP));
        }
        else if (filt_T == 3) {
            filtered = true;
            Snew(i,j,k,URHO ) = ff3(0,S) *  Stmp(i  ,j,k,URHO ) +
                                ff3(1,S) * (Stmp(i-1,j,k,URHO ) + Stmp(i+1,j,k,URHO )) +
                                ff3(2,S) * (Stmp(i-2,j,k,URHO ) + Stmp(i+2,j,k,URHO )) +
                                ff3(3,S) * (Stmp(i-3,j,k,URHO ) + Stmp(i+3,j,k,URHO ));
            Snew(i,j,k,UMX  ) = ff3(0,S) *  Stmp(i  ,j,k,UMX  ) +
                                ff3(1,S) * (Stmp(i-1,j,k,UMX  ) + Stmp(i+1,j,k,UMX  )) +
                                ff3(2,S) * (Stmp(i-2,j,k,UMX  ) + Stmp(i+2,j,k,UMX  )) +
                                ff3(3,S) * (Stmp(i-3,j,k,UMX  ) + Stmp(i+3,j,k,UMX  ));
            Snew(i,j,k,UTEMP) = ff3(0,S) *  Stmp(i  ,j,k,UTEMP) +
                                ff3(1,S) * (Stmp(i-1,j,k,UTEMP) + Stmp(i+1,j,k,UTEMP)) +
                                ff3(2,S) * (Stmp(i-2,j,k,UTEMP) + Stmp(i+2,j,k,UTEMP)) +
                                ff3(3,S) * (Stmp(i-3,j,k,UTEMP) + Stmp(i+3,j,k,UTEMP));
        }
        else if (filt_T == 4) {
            filtered = true;
            Snew(i,j,k,URHO ) = ff4(0,S) *  Stmp(i  ,j,k,URHO ) +
                                ff4(1,S) * (Stmp(i-1,j,k,URHO ) + Stmp(i+1,j,k,URHO )) +
                                ff4(2,S) * (Stmp(i-2,j,k,URHO ) + Stmp(i+2,j,k,URHO )) +
                                ff4(3,S) * (Stmp(i-3,j,k,URHO ) + Stmp(i+3,j,k,URHO )) +
                                ff4(4,S) * (Stmp(i-4,j,k,URHO ) + Stmp(i+4,j,k,URHO ));
            Snew(i,j,k,UMX  ) = ff4(0,S) *  Stmp(i  ,j,k,UMX  ) +
                                ff4(1,S) * (Stmp(i-1,j,k,UMX  ) + Stmp(i+1,j,k,UMX  )) +
                                ff4(2,S) * (Stmp(i-2,j,k,UMX  ) + Stmp(i+2,j,k,UMX  )) +
                                ff4(3,S) * (Stmp(i-3,j,k,UMX  ) + Stmp(i+3,j,k,UMX  )) +
                                ff4(4,S) * (Stmp(i-4,j,k,UMX  ) + Stmp(i+4,j,k,UMX  ));
            Snew(i,j,k,UTEMP) = ff4(0,S) *  Stmp(i  ,j,k,UTEMP) +
                                ff4(1,S) * (Stmp(i-1,j,k,UTEMP) + Stmp(i+1,j,k,UTEMP)) +
                                ff4(2,S) * (Stmp(i-2,j,k,UTEMP) + Stmp(i+2,j,k,UTEMP)) +
                                ff4(3,S) * (Stmp(i-3,j,k,UTEMP) + Stmp(i+3,j,k,UTEMP)) +
                                ff4(4,S) * (Stmp(i-4,j,k,UTEMP) + Stmp(i+4,j,k,UTEMP));
        }
    }

    if (filtered) {

        Snew(i,j,k,URHO ) = amrex::max(Snew(i,j,k,URHO ), castro::small_dens);
        Snew(i,j,k,UTEMP) = amrex::max(Snew(i,j,k,UTEMP), castro::small_temp);

        eos_t eos_state;

        eos_state.rho = Snew(i,j,k,URHO);
        eos_state.T = Snew(i,j,k,UTEMP);

        Real rhotmpInv = 1.0_rt / Stmp(i,j,k,URHO);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = Stmp(i,j,k,UFS+n) * rhotmpInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = Stmp(i,j,k,UFX+n) * rhotmpInv;
        }
#endif

        eos(eos_input_rt, eos_state);

        Snew(i,j,k,UEINT) = Snew(i,j,k,URHO) * eos_state.e;
        Snew(i,j,k,UEDEN) = Snew(i,j,k,UEINT) + 0.5_rt *
                                                (Snew(i,j,k,UMX) * Snew(i,j,k,UMX)) /
                                                 Snew(i,j,k,URHO);
        for (int n = 0; n < NumAdv; ++n) {
            Snew(i,j,k,UFA+n) = Stmp(i,j,k,UFA+n) * rhotmpInv * Snew(i,j,k,URHO);
        }
        for (int n = 0; n < NumSpec; ++n) {
            Snew(i,j,k,UFS+n) = eos_state.xn[n] * Snew(i,j,k,URHO);
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; ++n) {
            Snew(i,j,k,UFX+n) = eos_state.aux[n] * Snew(i,j,k,URHO);
        }
#endif

    }
}

#endif
