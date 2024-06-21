#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>

using namespace amrex;

using RealVector = amrex::Gpu::ManagedVector<amrex::Real>;

void ca_dergradpoverp(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

  // derive the dynamic pressure

  const auto dx = geomdata.CellSizeArray();
  const auto problo = geomdata.ProbLoArray();
  const int coord_type = geomdata.Coord();

  auto const dat = datfab.array();
  auto const der = derfab.array();

#if AMREX_SPACEDIM == 3
  amrex::Error("3D not supported");
#endif
#if AMREX_SPACEDIM == 1
  amrex::Error("1D not supported");
#endif

  Real dxinv = 1.0_rt / dx[0];
  Real dyinv = 1.0_rt / dx[1];

  // Compute grad p . U / (p |U|) -- this is what we do in the shock
  // detection algorithm

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
          eos_state.xn[n] = dat(i,j,k,UFS+n) / dat(i,j,k,URHO);
      }
#if NAUX_NET > 0
      for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) / dat(i,j,k,URHO);
      }
#endif

      if (eos_state.e <= 0.0_rt) {
          eos(eos_input_rt, eos_state);
      } else {
          eos(eos_input_re, eos_state);
      }
      Real p_zone = eos_state.p;

      // construct div{U}
      if (coord_type == 0) {

          // Cartesian
          div_u += 0.5_rt * (dat(i+1,j,k,UMX) / dat(i+1,j,k,URHO) -
                             dat(i-1,j,k,UMX) / dat(i-1,j,k,URHO)) * dxinv;
#if (AMREX_SPACEDIM >= 2)
          div_u += 0.5_rt * (dat(i,j+1,k,UMY) / dat(i,j+1,k,URHO) -
                             dat(i,j-1,k,UMY) / dat(i,j-1,k,URHO)) * dyinv;
#endif

      } else if (coord_type == 1) {

          // r-z
          Real rc = (i + 0.5_rt) * dx[0];
          Real rm = (i - 1 + 0.5_rt) * dx[0];
          Real rp = (i + 1 + 0.5_rt) * dx[0];

          div_u += 0.5_rt * (rp * dat(i+1,j,k,UMX) / dat(i+1,j,k,URHO) -
                             rm * dat(i-1,j,k,UMX) / dat(i-1,j,k,URHO)) / (rc * dx[0]);
          div_u += 0.5_rt * (rp * dat(i+1,j,k,UMX) / dat(i+1,j,k,URHO) -
                             rm * dat(i-1,j,k,UMX) / dat(i-1,j,k,URHO)) / (rc * dx[0]) +
              0.5_rt * (dat(i,j+1,k,UMY) / dat(i,j+1,k,URHO) -
                        dat(i,j-1,k,UMY) / dat(i,j-1,k,URHO)) * dyinv;

#ifndef AMREX_USE_GPU

      } else {
          amrex::Error("ERROR: invalid coord_type in shock");
#endif
      }


      // now compute (grad P - rho g) . dx

      // we need to compute p in the full stencil

      Real p_ip1{};
      {
          Real rhoInv = 1.0 / dat(i+1,j,k,URHO);

          eos_state.rho = dat(i+1,j,k,URHO);
          eos_state.T = dat(i+1,j,k,UTEMP);
          eos_state.e = dat(i+1,j,k,UEINT) * rhoInv;
          for (int n = 0; n < NumSpec; n++) {
              eos_state.xn[n] = dat(i+1,j,k,UFS+n) * rhoInv;
          }
#if NAUX_NET > 0
          for (int n = 0; n < NumAux; n++) {
              eos_state.aux[n] = dat(i+1,j,k,UFX+n) / rhoInv
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
          Real rhoInv = 1.0 / dat(i-1,j,k,URHO);

          eos_state.rho = dat(i-1,j,k,URHO);
          eos_state.T = dat(i-1,j,k,UTEMP);
          eos_state.e = dat(i-1,j,k,UEINT) * rhoInv;
          for (int n = 0; n < NumSpec; n++) {
              eos_state.xn[n] = dat(i-1,j,k,UFS+n) * rhoInv;
          }
#if NAUX_NET > 0
          for (int n = 0; n < NumAux; n++) {
              eos_state.aux[n] = dat(i-1,j,k,UFX+n) / rhoInv
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
          Real rhoInv = 1.0 / dat(i,j+1,k,URHO);

          eos_state.rho = dat(i,j+1,k,URHO);
          eos_state.T = dat(i,j+1,k,UTEMP);
          eos_state.e = dat(i,j+1,k,UEINT) * rhoInv;
          for (int n = 0; n < NumSpec; n++) {
              eos_state.xn[n] = dat(i,j+1,k,UFS+n) * rhoInv;
          }
#if NAUX_NET > 0
          for (int n = 0; n < NumAux; n++) {
              eos_state.aux[n] = dat(i,j+1,k,UFX+n) / rhoInv
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
          Real rhoInv = 1.0 / dat(i,j-1,k,URHO);

          eos_state.rho = dat(i,j-1,k,URHO);
          eos_state.T = dat(i,j-1,k,UTEMP);
          eos_state.e = dat(i,j-1,k,UEINT) * rhoInv;
          for (int n = 0; n < NumSpec; n++) {
              eos_state.xn[n] = dat(i,j-1,k,UFS+n) * rhoInv;
          }
#if NAUX_NET > 0
          for (int n = 0; n < NumAux; n++) {
              eos_state.aux[n] = dat(i,j-1,k,UFX+n) / rhoInv
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
      Real vel = std::sqrt(dat(i,j,k,UMX) * dat(i,j,k,UMX) +
                           dat(i,j,k,UMY) * dat(i,j,k,UMY) +
                           dat(i,j,k,UMZ) * dat(i,j,k,UMZ)) / dat(i,j,k,URHO);

      Real gradPdx_over_P{0.0_rt};
      if (vel != 0.0) {
          gradPdx_over_P = std::abs(dP_x * dat(i,j,k,UMX) / dat(i,j,k,URHO) +
                                    dP_y * dat(i,j,k,UMY) / dat(i,j,k,URHO)) / vel;
      }
      gradPdx_over_P /= p_zone;

      der(i,j,k,0) = gradPdx_over_P;

  });

}
