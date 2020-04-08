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

    void ca_derpres(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derpres(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                AMREX_REAL_ANYD(delta));

    }

    void ca_dereint1(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dereint1(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                 data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                 AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                 AMREX_REAL_ANYD(delta));

    }

    void ca_dereint2(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dereint2(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                 data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                 AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                 AMREX_REAL_ANYD(delta));

    }

    void ca_derlogden(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derlogden(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(delta));

    }

    void ca_deruplusc(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        deruplusc(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(delta));

    }

    void ca_deruminusc(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        deruminusc(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_dersoundspeed(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dersoundspeed(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                      data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                      AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                      AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                      AMREX_REAL_ANYD(delta));

    }

    void ca_dergamma1(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dergamma1(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(delta));

    }

    void ca_dermachnumber(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dermachnumber(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                      data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                      AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                      AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                      AMREX_REAL_ANYD(delta));

    }

    void ca_derentropy(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derentropy(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

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
    void ca_derenuctimescale(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derenuctimescale(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                         data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                         AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                         AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                         AMREX_REAL_ANYD(delta));

    }
#endif

    void ca_dervel(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dervel(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
               data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
               AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
               AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
               AMREX_REAL_ANYD(delta));

    }

    void ca_dermagvel(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dermagvel(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(delta));

    }

    void ca_dermaggrav(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dermaggrav(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_derradialvel(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derradialvel(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                     data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                     AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                     AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                     AMREX_REAL_ANYD(delta));

    }

    void ca_dermagmom(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dermagmom(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(delta));

    }

    void ca_derangmomx(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derangmomx(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_derangmomy(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derangmomy(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_derangmomz(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derangmomz(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_derkineng (const Box& bx, FArrayBox& kinengfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& /*geomdata*/,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
    {

        Box bx(ilo, ihi);
		auto const dat = datfab.array();
		auto const kineng = kinengfab.array();

		amrex::ParallelFor(bx,
		[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
		{
		    kineng(i,j,k,0) = 0.5_rt / dat(i,j,k,0) * ( dat(i,j,k,1)*dat(i,j,k,1) + 
														dat(i,j,k,2)*dat(i,j,k,2) + 
														dat(i,j,k,3)*dat(i,j,k,3) )
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

    void ca_derspec(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derspec(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                AMREX_REAL_ANYD(delta));

    }


    void ca_derabar(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derabar(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                AMREX_REAL_ANYD(delta));

    }

    void ca_dermagvort(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        dermagvort(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_derdivu(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derdivu(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                AMREX_REAL_ANYD(delta));

    }

    void ca_derstate(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
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
        derstate(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                 data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                 AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                 AMREX_REAL_ANYD(delta));

    }

#ifdef __cplusplus
}
#endif
