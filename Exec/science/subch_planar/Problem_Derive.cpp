#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>

using namespace amrex;

using RealVector = amrex::Gpu::ManagedVector<amrex::Real>;

void ca_dergradpoverp(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // Compute grad p . U / (p |U|) -- this is what we do in the shock
    // detection algorithm

    const auto dx = geomdata.CellSizeArray();
    const int coord_type = geomdata.Coord();

    auto const dat = datfab.array();
    auto const der = derfab.array();

#if AMREX_SPACEDIM == 3
    amrex::Error("3D not supported");
#elif AMREX_SPACEDIM == 1
    return; // Skip for 1D
#endif

    Real dxinv = 1.0_rt / dx[0];
    Real dyinv = 1.0_rt / dx[1];

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real div_u = 0.0_rt;

        // get the pressure
        eos_t eos_state;

        Real rhoInv = 1.0 / dat(i,j,k,URHO);

        eos_state.rho = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }
#endif

        if (eos_state.e <= 0.0_rt) {
            eos(eos_input_rt, eos_state);
        } else {
            eos(eos_input_re, eos_state);
        }
        Real p_zone = eos_state.p;

        Real up = dat(i+1,j,k,UMX) / dat(i+1,j,k,URHO);
        Real um = dat(i-1,j,k,UMX) / dat(i-1,j,k,URHO);
        Real u0 = dat(i,j,k,UMX) / dat(i,j,k,URHO);

        Real vp = dat(i,j+1,k,UMY) / dat(i,j+1,k,URHO);
        Real vm = dat(i,j-1,k,UMY) / dat(i,j-1,k,URHO);
        Real v0 = dat(i,j,k,UMY) / dat(i,j,k,URHO);

        // construct div{U}
        if (coord_type == 0) {

            // Cartesian
            div_u += 0.5_rt * (up - um) * dxinv;
            div_u += 0.5_rt * (vp - vm) * dyinv;

        } else if (coord_type == 1) {

            // r-z
            Real rc = (i + 0.5_rt) * dx[0];
            Real rm = (i - 1 + 0.5_rt) * dx[0];
            Real rp = (i + 1 + 0.5_rt) * dx[0];

            div_u += 0.5_rt * (rp * up - rm * um) / (rc * dx[0]) +
                0.5_rt * (vp - vm) * dyinv;

#ifndef AMREX_USE_GPU
        } else {
            amrex::Error("ERROR: invalid coord_type in shock");
#endif
        }


        // we need to compute p in the full stencil

        Real p_ip1{};
        {
            Real lrhoInv = 1.0 / dat(i+1,j,k,URHO);

            eos_state.rho = dat(i+1,j,k,URHO);
            eos_state.T = dat(i+1,j,k,UTEMP);
            eos_state.e = dat(i+1,j,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i+1,j,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i+1,j,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_ip1 = eos_state.p;
        }

        Real p_im1{};
        {
            Real lrhoInv = 1.0 / dat(i-1,j,k,URHO);

            eos_state.rho = dat(i-1,j,k,URHO);
            eos_state.T = dat(i-1,j,k,UTEMP);
            eos_state.e = dat(i-1,j,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i-1,j,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i-1,j,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_im1 = eos_state.p;
      }

        Real dP_x = 0.5_rt * (p_ip1 - p_im1);

        Real p_jp1{};
        {
            Real lrhoInv = 1.0 / dat(i,j+1,k,URHO);

            eos_state.rho = dat(i,j+1,k,URHO);
            eos_state.T = dat(i,j+1,k,UTEMP);
            eos_state.e = dat(i,j+1,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i,j+1,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i,j+1,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_jp1 = eos_state.p;
        }

        Real p_jm1{};
        {
            Real lrhoInv = 1.0 / dat(i,j-1,k,URHO);

            eos_state.rho = dat(i,j-1,k,URHO);
            eos_state.T = dat(i,j-1,k,UTEMP);
            eos_state.e = dat(i,j-1,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i,j-1,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i,j-1,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_jm1 = eos_state.p;
        }

        Real dP_y = 0.5_rt * (p_jp1 - p_jm1);

        //Real gradPdx_over_P = std::sqrt(dP_x * dP_x + dP_y * dP_y + dP_z * dP_z) / dat(i,j,k,QPRES);
        Real vel = std::sqrt(u0 * u0 + v0 * v0);

        Real gradPdx_over_P{0.0_rt};
        if (vel != 0.0) {
            gradPdx_over_P = std::abs(dP_x * u0 + dP_y * v0) / vel;
        }
        gradPdx_over_P /= p_zone;

        der(i,j,k,0) = gradPdx_over_P;

    });

}


void ca_dergradpoverp1(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // compute grad p projected onto div{U} as an alternative to the shock detection

    const auto dx = geomdata.CellSizeArray();
    const int coord_type = geomdata.Coord();

    auto const dat = datfab.array();
    auto const der = derfab.array();

#if AMREX_SPACEDIM == 3
    amrex::Error("3D not supported");
#elif AMREX_SPACEDIM == 1
    return; // Skip for 1D
#endif

    Real dxinv = 1.0_rt / dx[0];
    Real dyinv = 1.0_rt / dx[1];


    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real div_u = 0.0_rt;

        // get the pressure
        eos_t eos_state;

        Real rhoInv = 1.0 / dat(i,j,k,URHO);

        eos_state.rho = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }
#endif

        if (eos_state.e <= 0.0_rt) {
            eos(eos_input_rt, eos_state);
        } else {
            eos(eos_input_re, eos_state);
        }
        Real p_zone = eos_state.p;

        Real up = dat(i+1,j,k,UMX) / dat(i+1,j,k,URHO);
        Real um = dat(i-1,j,k,UMX) / dat(i-1,j,k,URHO);
        Real u0 = dat(i,j,k,UMX) / dat(i,j,k,URHO);

        Real vp = dat(i,j+1,k,UMY) / dat(i,j+1,k,URHO);
        Real vm = dat(i,j-1,k,UMY) / dat(i,j-1,k,URHO);
        Real v0 = dat(i,j,k,UMY) / dat(i,j,k,URHO);

        Real du_x{};
        Real dv_y{};

        // construct div{U}
        if (coord_type == 0) {

            // Cartesian
            du_x = 0.5_rt * (up - um) * dxinv;
            dv_y = 0.5_rt * (vp - vm) * dyinv;

        } else if (coord_type == 1) {

            // r-z
            Real rc = (i + 0.5_rt) * dx[0];
            Real rm = (i - 1 + 0.5_rt) * dx[0];
            Real rp = (i + 1 + 0.5_rt) * dx[0];

            du_x = 0.5_rt * (rp * up - rm * um) / (rc * dx[0]);
            dv_y = 0.5_rt * (vp - vm) * dyinv;

#ifndef AMREX_USE_GPU
        } else {
            amrex::Error("ERROR: invalid coord_type in shock");
#endif
        }

        div_u = du_x + dv_y;

        // we need to compute p in the full stencil

        Real p_ip1{};
        {
            Real lrhoInv = 1.0 / dat(i+1,j,k,URHO);

            eos_state.rho = dat(i+1,j,k,URHO);
            eos_state.T = dat(i+1,j,k,UTEMP);
            eos_state.e = dat(i+1,j,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i+1,j,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i+1,j,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_ip1 = eos_state.p;
        }

        Real p_im1{};
        {
            Real lrhoInv = 1.0 / dat(i-1,j,k,URHO);

            eos_state.rho = dat(i-1,j,k,URHO);
            eos_state.T = dat(i-1,j,k,UTEMP);
            eos_state.e = dat(i-1,j,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i-1,j,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i-1,j,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_im1 = eos_state.p;
        }

        Real dP_x = 0.5_rt * (p_ip1 - p_im1);

        Real p_jp1{};
        {
            Real lrhoInv = 1.0 / dat(i,j+1,k,URHO);

            eos_state.rho = dat(i,j+1,k,URHO);
            eos_state.T = dat(i,j+1,k,UTEMP);
            eos_state.e = dat(i,j+1,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i,j+1,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i,j+1,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_jp1 = eos_state.p;
        }

        Real p_jm1{};
        {
            Real lrhoInv = 1.0 / dat(i,j-1,k,URHO);

            eos_state.rho = dat(i,j-1,k,URHO);
            eos_state.T = dat(i,j-1,k,UTEMP);
            eos_state.e = dat(i,j-1,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i,j-1,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i,j-1,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_jm1 = eos_state.p;
        }

        Real dP_y = 0.5_rt * (p_jp1 - p_jm1);

        //Real gradPdx_over_P = std::sqrt(dP_x * dP_x + dP_y * dP_y + dP_z * dP_z) / dat(i,j,k,QPRES);
        Real cdu_x = std::min(du_x, 0.0);
        Real cdv_y = std::min(dv_y, 0.0);

        Real divu_mag = std::sqrt(cdu_x * cdu_x + cdv_y * cdv_y + 1.e-30);

        Real gradPdx_over_P = std::abs(dP_x * cdu_x + dP_y * cdv_y) / divu_mag;
        gradPdx_over_P /= p_zone;

        der(i,j,k,0) = gradPdx_over_P;

  });

}


void ca_dergradpx(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    amrex::ignore_unused(geomdata);

    // Compute grad p / p . xhat

    auto const dat = datfab.array();
    auto const der = derfab.array();

#if AMREX_SPACEDIM == 3
    amrex::Error("3D not supported");
#elif AMREX_SPACEDIM == 1
    return; // Skip for 1D
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // get the pressure
        eos_t eos_state;

        Real rhoInv = 1.0 / dat(i,j,k,URHO);

        eos_state.rho = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }
#endif

        if (eos_state.e <= 0.0_rt) {
            eos(eos_input_rt, eos_state);
        } else {
            eos(eos_input_re, eos_state);
        }
        Real p_zone = eos_state.p;


        // we need to compute p in the full stencil

        Real p_ip1{};
        {
            Real lrhoInv = 1.0 / dat(i+1,j,k,URHO);

            eos_state.rho = dat(i+1,j,k,URHO);
            eos_state.T = dat(i+1,j,k,UTEMP);
            eos_state.e = dat(i+1,j,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i+1,j,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i+1,j,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_ip1 = eos_state.p;
        }

        Real p_im1{};
        {
            Real lrhoInv = 1.0 / dat(i-1,j,k,URHO);

            eos_state.rho = dat(i-1,j,k,URHO);
            eos_state.T = dat(i-1,j,k,UTEMP);
            eos_state.e = dat(i-1,j,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i-1,j,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i-1,j,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_im1 = eos_state.p;
      }

        Real dp_x = 0.5_rt * (p_ip1 - p_im1);

        der(i,j,k,0) = std::abs(dp_x) / p_zone;

    });

}


void ca_dergradpy(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    amrex::ignore_unused(geomdata);

    // compute grad p / p . yhat

    auto const dat = datfab.array();
    auto const der = derfab.array();

#if AMREX_SPACEDIM == 3
    amrex::Error("3D not supported");
#elif AMREX_SPACEDIM == 1
    return; // Skip for 1D
#endif

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        // get the pressure
        eos_t eos_state;

        Real rhoInv = 1.0 / dat(i,j,k,URHO);

        eos_state.rho = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }
#endif

        if (eos_state.e <= 0.0_rt) {
            eos(eos_input_rt, eos_state);
        } else {
            eos(eos_input_re, eos_state);
        }
        Real p_zone = eos_state.p;



        // we need to compute p in the full stencil

        Real p_jp1{};
        {
            Real lrhoInv = 1.0 / dat(i,j+1,k,URHO);

            eos_state.rho = dat(i,j+1,k,URHO);
            eos_state.T = dat(i,j+1,k,UTEMP);
            eos_state.e = dat(i,j+1,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i,j+1,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i,j+1,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_jp1 = eos_state.p;
        }

        Real p_jm1{};
        {
            Real lrhoInv = 1.0 / dat(i,j-1,k,URHO);

            eos_state.rho = dat(i,j-1,k,URHO);
            eos_state.T = dat(i,j-1,k,UTEMP);
            eos_state.e = dat(i,j-1,k,UEINT) * lrhoInv;
            for (int n = 0; n < NumSpec; n++) {
                eos_state.xn[n] = dat(i,j-1,k,UFS+n) * lrhoInv;
            }
#if NAUX_NET > 0
            for (int n = 0; n < NumAux; n++) {
                eos_state.aux[n] = dat(i,j-1,k,UFX+n) * lrhoInv;
            }
#endif

            if (eos_state.e <= 0.0_rt) {
                eos(eos_input_rt, eos_state);
            } else {
                eos(eos_input_re, eos_state);
            }
            p_jm1 = eos_state.p;
        }

        Real dp_y = 0.5_rt * (p_jp1 - p_jm1);

        der(i,j,k,0) = std::abs(dp_y) / p_zone;

  });

}


void ca_dergradrhooverrho(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                          const FArrayBox& datfab, const Geometry& geomdata,
                          Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    // compute grad rho / rho

    const auto dx = geomdata.CellSizeArray();
    const int coord_type = geomdata.Coord();

    auto const dat = datfab.array();
    auto const der = derfab.array();

#if AMREX_SPACEDIM == 3
    amrex::Error("3D not supported");
#elif AMREX_SPACEDIM == 1
    return; // Skip for 1D
#endif

    Real dxinv = 1.0_rt / dx[0];
    Real dyinv = 1.0_rt / dx[1];

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        Real div_u = 0.0_rt;

        Real up = dat(i+1,j,k,UMX) / dat(i+1,j,k,URHO);
        Real um = dat(i-1,j,k,UMX) / dat(i-1,j,k,URHO);
        Real u0 = dat(i,j,k,UMX) / dat(i,j,k,URHO);

        Real vp = dat(i,j+1,k,UMY) / dat(i,j+1,k,URHO);
        Real vm = dat(i,j-1,k,UMY) / dat(i,j-1,k,URHO);
        Real v0 = dat(i,j,k,UMY) / dat(i,j,k,URHO);

        // construct div{U}
        if (coord_type == 0) {

            // Cartesian
            div_u += 0.5_rt * (up - um) * dxinv;
            div_u += 0.5_rt * (vp - vm) * dyinv;

        } else if (coord_type == 1) {

            // r-z
            Real rc = (i + 0.5_rt) * dx[0];
            Real rm = (i - 1 + 0.5_rt) * dx[0];
            Real rp = (i + 1 + 0.5_rt) * dx[0];

            div_u += 0.5_rt * (rp * up - rm * um) / (rc * dx[0]) +
                0.5_rt * (vp - vm) * dyinv;

#ifndef AMREX_USE_GPU
        } else {
            amrex::Error("ERROR: invalid coord_type in shock");
#endif
        }


        Real drho_x = 0.5_rt * (dat(i+1,j,k,URHO) - dat(i-1,j,k,URHO));
        Real drho_y = 0.5_rt * (dat(i,j+1,k,URHO) - dat(i,j-1,k,URHO));

        Real vel = std::sqrt(u0 * u0 + v0 * v0);

        Real gradrhodx_over_rho{0.0_rt};
        if (vel != 0.0) {
            gradrhodx_over_rho = std::abs(drho_x * u0 + drho_y * v0) / vel;
        }
        gradrhodx_over_rho /= dat(i,j,k,URHO);

        der(i,j,k,0) = gradrhodx_over_rho;

    });

}
