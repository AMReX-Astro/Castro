#include "AMReX_REAL.H"

#include "Derive.H"
#include "Derive_F.H"
#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

#ifdef __cplusplus
extern "C"
{
#endif

    // Note that in the following routines, we are NOT passing
    // several variables to Fortran that would be unused.

    // These routines are called in an MFIter loop, so we do not
    // need to explicitly synchronize after GPU kernels.

    void ca_derpres(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        eos(eos_input_re, eos_state);

        der(i,j,k,0) = eos_state.p;
      });
    }

    void ca_dereint1(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& geomdata,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real rhoInv = 1.0_rt/dat(i,j,k,URHO);
        Real ux = dat(i,j,k,UMX)*rhoInv;
        Real uy = dat(i,j,k,UMY)*rhoInv;
        Real uz = dat(i,j,k,UMZ)*rhoInv;

        der(i,j,k,0) = dat(i,j,k,UEDEN)*rhoInv -
          0.5_rt * (ux*ux + uy*uy + uz*uz);
      });
    }

    void ca_dereint2(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& geomdata,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = dat(i,j,k,UEINT) / dat(i,j,k,URHO);
      });
    }

    void ca_derlogden(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
        der(i,j,k,0) = std::log10(dat(i,j,k,0));
      });
    }

    void ca_deruplusc(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {


        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        eos(eos_input_re, eos_state);

        der(i,j,k,0) = dat(i,j,k,UMX) / dat(i,j,k,URHO) + eos_state.cs;

      });
    }

    void ca_deruminusc(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {


        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        eos(eos_input_re, eos_state);

        der(i,j,k,0) = dat(i,j,k,UMX) / dat(i,j,k,URHO) - eos_state.cs;

      });
    }

    void ca_dersoundspeed(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                          const FArrayBox& datfab, const Geometry& geomdata,
                          Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {


        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        eos(eos_input_re, eos_state);

        der(i,j,k,0) = eos_state.cs;

      });
    }


    void ca_dergamma1(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {


        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        eos(eos_input_re, eos_state);

        der(i,j,k,0) = eos_state.gam1;

      });
    }

    void ca_dermachnumber(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                          const FArrayBox& datfab, const Geometry& geomdata,
                          Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {


        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        der(i,j,k,0) = std::sqrt(dat(i,j,k,UMX)*dat(i,j,k,UMX) +
                                 dat(i,j,k,UMY)*dat(i,j,k,UMY) +
                                 dat(i,j,k,UMZ)*dat(i,j,k,UMZ)) /
          (dat(i,j,k,URHO) * eos_state.cs);

      });
    }

    void ca_derentropy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                          const FArrayBox& datfab, const Geometry& geomdata,
                          Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {


        Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

        eos_t eos_state;
        eos_state.rho  = dat(i,j,k,URHO);
        eos_state.T = dat(i,j,k,UTEMP);
        eos_state.e = dat(i,j,k,UEINT) * rhoInv;
        for (int n = 0; n < NumSpec; n++) {
          eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
        }
        for (int n = 0; n < NumAux; n++) {
          eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
        }

        der(i,j,k,0) = eos_state.s;
      });
    }

#ifdef DIFFUSION
    void ca_dercond(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                    const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                    const int* lo, const int* hi,
                    const int* domain_lo, const int* domain_hi,
                    const Real* delta, const Real* xlo,
                    const Real* time, const Real* dt, const int* bcrec, 
                    const int* level, const int* grid_no)
    {

        IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
        IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

        Box bx(ilo, ihi);

#pragma gpu box(bx)
        dercond(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                AMREX_REAL_ANYD(delta));

    }

    void ca_derdiffcoeff(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                         const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                         const int* lo, const int* hi,
                         const int* domain_lo, const int* domain_hi,
                         const Real* delta, const Real* xlo,
                         const Real* time, const Real* dt, const int* bcrec, 
                         const int* level, const int* grid_no)
    {

        IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
        IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

        Box bx(ilo, ihi);

#pragma gpu box(bx)
        derdiffcoeff(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                     der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                     data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                     AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                     AMREX_REAL_ANYD(delta));

    }

    void ca_derdiffterm(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                        const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const Real* delta, const Real* xlo,
                        const Real* time, const Real* dt, const int* bcrec,
                        const int* level, const int* grid_no)
    {

        // Create an array for storing cell-centered conductivity data.
        // It needs to have a ghost zone for the next step.

        IntVect ilo(D_DECL(lo[0], lo[1], lo[2]));
        IntVect ihi(D_DECL(hi[0], hi[1], hi[2]));

        const Box bx(ilo, ihi);
        const Box& obx = amrex::grow(bx, 1);

        FArrayBox coeff_cc;
        coeff_cc.resize(obx, 1);
        Elixir elix_coeff_cc = coeff_cc.elixir();
        Array4<Real> const coeff_arr = coeff_cc.array();

        FArrayBox coeffs[3];
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            coeffs[dir].resize(amrex::surroundingNodes(bx, dir), 1);
        }

        Elixir elix_coeffs_x = coeffs[0].elixir();
        Elixir elix_coeffs_y = coeffs[1].elixir();
        Elixir elix_coeffs_z = coeffs[2].elixir();

#pragma gpu box(obx)
        ca_fill_temp_cond(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                          data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
                          BL_TO_FORTRAN_ANYD(coeff_cc));

        // Now average the data to zone edges.

        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            Array4<Real> const edge_coeff_arr = (coeffs[idir]).array();

            AMREX_PARALLEL_FOR_3D(nbx, i, j, k,
            {

              if (idir == 0) {
                edge_coeff_arr(i,j,k) = 0.5_rt * (coeff_arr(i,j,k) + coeff_arr(i-1,j,k));
              } else if (idir == 1) {
                       edge_coeff_arr(i,j,k) = 0.5_rt * (coeff_arr(i,j,k) + coeff_arr(i,j-1,k));
              } else {
                edge_coeff_arr(i,j,k) = 0.5_rt * (coeff_arr(i,j,k) + coeff_arr(i,j,k-1));
              }
            });

        }

#pragma gpu box(bx)
        derdiffterm(AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                    der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                    data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                    BL_TO_FORTRAN_ANYD(coeffs[0]),
#if AMREX_SPACEDIM >= 2
                    BL_TO_FORTRAN_ANYD(coeffs[1]),
#endif
#if AMREX_SPACEDIM == 3
                    BL_TO_FORTRAN_ANYD(coeffs[2]),
#endif
                    AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                    AMREX_REAL_ANYD(delta));

    }
#endif

#ifdef REACTIONS
    void ca_derenuctimescale(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geomdata,
                             Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto dx = geomdata.CellSizeArray();

      Real dd = 0.0_rt;
#if AMREX_SPACEDIM == 1
      dd = dx[0];
#elif AMREX_SPACEDIM == 2
      dd = amrex::min(dx[0], dx[1]);
#else
      dd = amrex::min(dx[0], dx[1], dx[2]);
#endif

      int enuc_comp = datfab.nComp()-1;

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        // the nuclear energy (rho H_nuc) is tacked onto the end of
        // the input state, after the NUM_STATE conserved state
        // quantities
        Real enuc = std::abs(dat(i,j,k,enuc_comp));

        if (enuc > 1.e-100_rt) {

          Real rhoInv = 1.0_rt / dat(i,j,k,URHO);

          // calculate the sound speed
          eos_t eos_state;
          eos_state.rho  = dat(i,j,k,URHO);
          eos_state.T = dat(i,j,k,UTEMP);
          eos_state.e = dat(i,j,k,UEINT) * rhoInv;
          for (int n = 0; n < NumSpec; n++) {
            eos_state.xn[n] = dat(i,j,k,UFS+n) * rhoInv;
          }
          for (int n = 0; n < NumAux; n++) {
            eos_state.aux[n] = dat(i,j,k,UFX+n) * rhoInv;
          }

          eos(eos_input_re, eos_state);

          Real t_e = eos_state.e / enuc;
          Real t_s = dd / eos_state.cs;

          der(i,j,k,0) = t_s/t_e;

        } else {
          der(i,j,k,0) = 0.0_rt;
        }

      });
    }
#endif

    void ca_dervel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
             der(i,j,k,0) = dat(i,j,k,1) / dat(i,j,k,0);
      });
    }


    void ca_dermagvel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real deninv = 1.0_rt/dat(i,j,k,0);

        der(i,j,k,0) = std::sqrt( (dat(i,j,k,1) * dat(i,j,k,1) +
                                   dat(i,j,k,2) * dat(i,j,k,2) +
                                   dat(i,j,k,3) * dat(i,j,k,3)) ) * deninv;
      });
    }


    void ca_dermaggrav(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        der(i,j,k,0) = std::sqrt(dat(i,j,k,0)*dat(i,j,k,0) +
                                 dat(i,j,k,1)*dat(i,j,k,1) +
                                 dat(i,j,k,2)*dat(i,j,k,2));

      });
    }

    void ca_derradialvel(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto dx = geomdata.CellSizeArray();

      // center calculated like advection_utils.cpp
      GpuArray<Real, 3> center;
      ca_get_center(center.begin());

      auto problo = geomdata.ProbLoArray();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real x = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];
#if AMREX_SPACEDIM >= 2
        Real y = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        Real y = 0.0_rt;
#endif
#if AMREX_SPACEDIM == 3
        Real z = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        Real z = 0.0_rt;
#endif

        Real r = std::sqrt(x*x + y*y + z*z);

        der(i,j,k,0) = (dat(i,j,k,1)*x +
                        dat(i,j,k,2)*y +
                        dat(i,j,k,3)*z) / ( dat(i,j,k,0)*r );
      });
    }


  void ca_dermagmom(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                    const FArrayBox& datfab, const Geometry& geomdata,
                    Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {

                         der(i,j,k,0) = std::sqrt(dat(i,j,k,0)*dat(i,j,k,0) +
                                                  dat(i,j,k,1)*dat(i,j,k,1) +
                                                  dat(i,j,k,2)*dat(i,j,k,2));

                       });
  }

  void ca_derangmomx (const Box& bx, FArrayBox& Lfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    int idir = 0;
    auto dx     = geomdata.CellSizeArray();
    auto problo = geomdata.ProbLoArray();

    // center calculated like advection_utils.cpp
    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    auto const dat = datfab.array();
    auto const L = Lfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {
                         Real loc[3];

                         //loc calculated like sum_utils.cpp
                         //This might be equivalent and more modular: position(i, j, k, geomdata, loc);
                         loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
                         loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
                         loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                         loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
                         loc[2] = 0.0_rt;
#endif

                         for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                           loc[dir] -= center[dir];
                         }

                         // Explicitly computing only the required cross-product as in inertial_to_rotational_velocity_c
                         if (idir == 0) { // cross_product(loc, mom): ang_mom(1)->x)
                           L(i,j,k,0) = loc[1] * dat(i,j,k,3) - loc[2] * dat(i,j,k,2);
                         }
                         else if (idir == 1) { // cross_product(loc, mom): ang_mom(2)->y)
                           L(i,j,k,0) = loc[2] * dat(i,j,k,1) - loc[0] * dat(i,j,k,3);
                         }
                         else { // cross_product(loc, mom): ang_mom(3)->z)
                           L(i,j,k,0) = loc[0] * dat(i,j,k,2) - loc[1] * dat(i,j,k,1);
                         }

                       });

  }

  void ca_derangmomy (const Box& bx, FArrayBox& Lfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    int idir = 1;
    auto dx     = geomdata.CellSizeArray();
    auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    auto const dat = datfab.array();
    auto const L = Lfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {
                         Real loc[3];

                         loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
                         loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
                         loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                         loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
                         loc[2] = 0.0_rt;
#endif
                         for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                           loc[dir] -= center[dir];
                         }

                         if (idir == 0) { // cross_product(loc, mom): ang_mom(1)->x)
                           L(i,j,k,0) = loc[1] * dat(i,j,k,3) - loc[2] * dat(i,j,k,2);
                         }
                         else if (idir == 1) { // cross_product(loc, mom): ang_mom(2)->y)
                           L(i,j,k,0) = loc[2] * dat(i,j,k,1) - loc[0] * dat(i,j,k,3);
                         }
                         else { // cross_product(loc, mom): ang_mom(3)->z)
                           L(i,j,k,0) = loc[0] * dat(i,j,k,2) - loc[1] * dat(i,j,k,1);
                         }

                       });

  }

  void ca_derangmomz (const Box& bx, FArrayBox& Lfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    int idir = 2;
    auto dx     = geomdata.CellSizeArray();
    auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    auto const dat = datfab.array();
    auto const L = Lfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {
                         Real loc[3];

                         loc[0] = problo[0] + (0.5_rt + i) * dx[0];

#if AMREX_SPACEDIM >= 2
                         loc[1] = problo[1] + (0.5_rt + j) * dx[1];
#else
                         loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
                         loc[2] = problo[2] + (0.5_rt + k) * dx[2];
#else
                         loc[2] = 0.0_rt;
#endif

                         for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                           loc[dir] -= center[dir];
                         }

                         if (idir == 0) { // cross_product(loc, mom): ang_mom(1)->x)
                           L(i,j,k,0) = loc[1] * dat(i,j,k,3) - loc[2] * dat(i,j,k,2);
                         }
                         else if (idir == 1) { // cross_product(loc, mom): ang_mom(2)->y)
                           L(i,j,k,0) = loc[2] * dat(i,j,k,1) - loc[0] * dat(i,j,k,3);
                         }
                         else { // cross_product(loc, mom): ang_mom(3)->z)
                           L(i,j,k,0) = loc[0] * dat(i,j,k,2) - loc[1] * dat(i,j,k,1);
                         }

                       });

  }

  void ca_derkineng (const Box& bx, FArrayBox& kinengfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& /*geomdata*/,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    auto const dat = datfab.array();
    auto const kineng = kinengfab.array();

    amrex::ParallelFor(bx,
                       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                       {
                         kineng(i,j,k,0) = 0.5_rt / dat(i,j,k,0) * ( dat(i,j,k,1)*dat(i,j,k,1) +
                                                                     dat(i,j,k,2)*dat(i,j,k,2) +
                                                                     dat(i,j,k,3)*dat(i,j,k,3) );
                       });

  }

  void ca_dernull(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                  const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                  const int* lo, const int* hi,
                  const int* domain_lo, const int* domain_hi,
                  const Real* delta, const Real* xlo,
                  const Real* time, const Real* dt, const int* bcrec, 
                  const int* level, const int* grid_no)
  {

    // This routine is used by particle_count.  Yes it does nothing.

  }

  void ca_derspec(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
             der(i,j,k,0) = dat(i,j,k,1) / dat(i,j,k,0);
      });
    }


  void ca_derabar(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real sum = 0.0_rt;
        Real xn;
        for (int n = 0; n < NumSpec; n++) {
          xn = dat(i,j,k,1+n)/dat(i,j,k,0);
          sum += xn/aion[n];
        }
        der(i,j,k,0) = 1.0_rt / sum;
      });
    }

  void ca_dermagvort(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& geomdata,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto dx = geomdata.CellSizeArray();

      const int coord_type = geomdata.Coord();

      auto problo = geomdata.ProbLoArray();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        // Calculate vorticity.
        if (coord_type == 0) {
          // Cartesian

          // dv/dx and dw/dx
          Real vx = 0.5_rt * (dat(i+1,j,k,2) / dat(i+1,j,k,0) -
                              dat(i-1,j,k,2) / dat(i-1,j,k,0)) / dx[0];
          Real wx = 0.5_rt * (dat(i+1,j,k,3) / dat(i+1,j,k,0) -
                              dat(i-1,j,k,3) / dat(i-1,j,k,0)) / dx[0];

          Real uy = 0.0;
          Real wy = 0.0;

#if AMREX_SPACEDIM >= 2
          // du/dy and dw/dy
          uy = 0.5_rt * (dat(i,j+1,k,1) / dat(i,j+1,k,0) -
                         dat(i,j-1,k,1) / dat(i,j-1,k,0)) / dx[1];
          wy = 0.5_rt * (dat(i,j+1,k,3) / dat(i,j+1,k,0) -
                         dat(i,j-1,k,3) / dat(i,j-1,k,0)) / dx[1];
#endif

          Real uz = 0.0;
          Real vz = 0.0;

#if AMREX_SPACEDIM == 3
          // du/dz and dv/dz
          uz = 0.5_rt * (dat(i,j,k+1,1) / dat(i,j,k+1,0) -
                         dat(i,j,k-1,1) / dat(i,j,k-1,0)) / dx[2];
          vz = 0.5_rt * (dat(i,j,k+1,2) / dat(i,j,k+1,0) -
                         dat(i,j,k-1,2) / dat(i,j,k-1,0)) / dx[2];
#endif

          // curl in Cartesian coords
         Real v1 = wy - vz;
         Real v2 = uz - wx;
         Real v3 = vx - uy;
         der(i,j,k,0) = std::sqrt(v1*v1 + v2*v2 + v3*v3);

        } else if (coord_type == 1) {
          // 2-d axisymmetric -- the coordinate ordering is r, z, phi

          Real r = (static_cast<Real>(i) + 0.5_rt)*dx[0] + problo[0];
          Real rm1 = (static_cast<Real>(i) - 0.5_rt)*dx[0] + problo[0];
          Real rp1 = (static_cast<Real>(i) + 1.5_rt)*dx[0] + problo[0];

          // dv_r/dz
          Real vr_z = 0.5_rt * (dat(i,j+1,k,1) / dat(i,j+1,k,0) -
                                dat(i,j-1,k,1) / dat(i,j-1,k,0)) / dx[1];

          // dv_phi/dz
          Real vphi_z = 0.5_rt * (dat(i,j+1,k,3) / dat(i,j+1,k,0) -
                                  dat(i,j-1,k,3) / dat(i,j-1,k,0)) / dx[1];

          // d (r v_phi)/dr
          Real rvphi_r = 0.5_rt * (rp1 * dat(i+1,j,k,3) / dat(i+1,j,k,0) -
                                   rm1 * dat(i-1,j,k,3) / dat(i-1,j,k,0)) / dx[0];

          // dv_z/dr
          Real vz_r = 0.5_rt * (dat(i+1,j,k,2) / dat(i+1,j,k,0) -
                                dat(i-1,j,k,2) / dat(i-1,j,k,0)) / dx[0];

          der(i,j,k,0) = std::sqrt(vphi_z*vphi_z +
                                   (vr_z - vz_r)*(vr_z - vz_r) +
                                   (rvphi_r/r)*(rvphi_r/r));

        } else if (coord_type == 2) {
          // 1-d spherical -- we don't really have a vorticity in this
          // case
          der(i,j,k,0) = 0.0;

        }
      });
    }

  void ca_derdivu(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

      auto const dat = datfab.array();
      auto const der = derfab.array();

      auto dx = geomdata.CellSizeArray();

      // center calculated like advection_utils.cpp
      GpuArray<Real, 3> center;
      ca_get_center(center.begin());

      auto problo = geomdata.ProbLoArray();

      const int coord_type = geomdata.Coord();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        Real uhi = dat(i+1,j,k,1) / dat(i+1,j,k,0);
        Real ulo = dat(i-1,j,k,1) / dat(i-1,j,k,0);

        Real vhi = dat(i,j+dg1,k,2) / dat(i,j+dg1,k,0);
        Real vlo = dat(i,j-dg1,k,2) / dat(i,j-dg1,k,0);

        Real whi = dat(i,j,k+dg2,3) / dat(i,j,k+dg2,0);
        Real wlo = dat(i,j,k-dg2,3) / dat(i,j,k-dg2,0);

        if (coord_type == 0) {
          // Cartesian divergence

          der(i,j,k,0) = 0.5_rt * (uhi - ulo) / dx[0];
#if AMREX_SPACEDIM >= 2
          der(i,j,k,0) += 0.5_rt * (vhi - vlo) / dx[1];
#endif
#if AMREX_SPACEDIM == 3
          der(i,j,k,0) += 0.5_rt * (whi - wlo) / dx[2];
#endif
        } else if (coord_type == 1) {
          // axisymmetric divergence -- defined only for 2-d axisymmetric

          Real r = (static_cast<Real>(i) + 0.5_rt)*dx[0] + problo[0];
          Real rm1 = (static_cast<Real>(i) - 0.5_rt)*dx[0] + problo[0];
          Real rp1 = (static_cast<Real>(i) + 1.5_rt)*dx[0] + problo[0];

          der(i,j,k,0) = 0.5_rt * (rp1*uhi - rm1*ulo) / (r*dx[0]) +
            0.5_rt*(vhi - vlo) / dx[1];

        } else if (coord_type == 2) {

          Real r = (static_cast<Real>(i) + 0.5_rt)*dx[0] + problo[0];
          Real rm1 = (static_cast<Real>(i) - 0.5_rt)*dx[0] + problo[0];
          Real rp1 = (static_cast<Real>(i) + 1.5_rt)*dx[0] + problo[0];

          der(i,j,k,0) = 0.5_rt * (rp1*rp1*uhi - rm1*rm1*ulo) /
            (r*r * dx[0]);
        }

      });
  }

  void ca_derstate(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const FArrayBox& datfab, const Geometry& geomdata,
                   Real /*time*/, const int* /*bcrec*/, int /*level*/)
  {

    auto const dat = datfab.array();
      auto const der = derfab.array();

      amrex::ParallelFor(bx,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {

        // density
        der(i,j,k,0) = dat(i,j,k,0);

        // temperature
        der(i,j,k,1) = dat(i,j,k,1);

        // (rho X)_1 = X_1
        der(i,j,k,2) = dat(i,j,k,2) / dat(i,j,k,0);

      });

    }

#ifdef __cplusplus
}
#endif
