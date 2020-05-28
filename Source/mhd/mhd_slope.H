#ifndef CASTRO_MHD_SLOPE_H
#define CASTRO_MHD_SLOPE_H

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void minmod(Real& dW, const Real WR,
                                                     const Real WL) {

    dW = 0.0;

    if (WL * WR > 0.0) {
        if (std::abs(WR) < std::abs(WL)) {
            dW = WR;
        } else {
            dW = WL;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void vanleer(Real& dW, const Real WR,
                                                      const Real WL) {

    dW = 0.0;

    if (WR * WL > 0.0) {
        dW = 2.0_rt * WR * WL / (WR + WL);
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void centerdif(Real& dW, const Real WR,
                                                        const Real WL) {
    dW = (WR + WL) / 2.0_rt;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void secondMC(Real& dW, const Real WR,
                                                       const Real WL) {

    Real dlim = 0.0;

    if (WR * WL >= 0.0_rt) {
        dlim = 2.0_rt * amrex::min(std::abs(WR), std::abs(WL));
    }

    dW = amrex::min(0.5_rt * std::abs(WR + WL), dlim) *
         std::copysign(1.0_rt, WR + WL);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void slope(Real& dW, const Real WR,
                                                    const Real WL,
                                                    const Real flat) {

    if (mhd_plm_slope == 0) {
        dW = 0.0;

    } else if (mhd_plm_slope == 1) {
        vanleer(dW, WR, WL);

    } else if (mhd_plm_slope == 2) {
        centerdif(dW, WR, WL);

    } else if (mhd_plm_slope == 3) {
        secondMC(dW, WR, WL);
    }

    if (use_flattening == 1) {
        dW = flat * dW;
    }
}

#endif
