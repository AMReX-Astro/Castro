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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

        FArrayBox coeffs[3];
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            coeffs[dir].resize(amrex::surroundingNodes(bx, dir), 1);
        }

        Elixir elix_coeffs_x = coeffs[0].elixir();
        Elixir elix_coeffs_y = coeffs[1].elixir();
        Elixir elix_coeffs_z = coeffs[2].elixir();

#pragma gpu
        ca_fill_temp_cond(AMREX_INT_ANYD(obx.loVect()), AMREX_INT_ANYD(obx.hiVect()),
                          data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi),
                          BL_TO_FORTRAN_ANYD(coeff_cc));

        // Now average the data to zone edges.

        for (int idir = 0; idir < AMREX_SPACEDIM; ++idir) {

            const Box& nbx = amrex::surroundingNodes(bx, idir);

            const int idir_f = idir + 1;

#pragma gpu
            ca_average_coef_cc_to_ec(AMREX_INT_ANYD(nbx.loVect()), AMREX_INT_ANYD(nbx.hiVect()),
                                     BL_TO_FORTRAN_ANYD(coeff_cc),
                                     BL_TO_FORTRAN_ANYD(coeffs[idir]),
                                     idir_f);

        }

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
        derangmomz(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(delta));

    }

    void ca_derkineng(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                      const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const Real* delta, const Real* xlo,
                      const Real* time, const Real* dt, const int* bcrec, 
                      const int* level, const int* grid_no)
    {

#pragma gpu
        derkineng(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(delta));

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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
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

#pragma gpu
        derstate(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                 data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                 AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                 AMREX_REAL_ANYD(delta));

    }

#ifdef __cplusplus
}
#endif
