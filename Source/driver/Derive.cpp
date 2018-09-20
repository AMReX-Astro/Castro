#include "AMReX_REAL.H"

#include "Derive.H"
#include "Derive_F.H"

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
                AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                 AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                 AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                  AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                  AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                      AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                  AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                      AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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

        dercond(der, ARLIM_3D(der_lo), ARLIM_3D(der_hi), *nvar,
                data, ARLIM_3D(data_lo), ARLIM_3D(data_hi), *ncomp,
                ARLIM_3D(lo), ARLIM_3D(hi),
                ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
                ZFILL(delta), ZFILL(xlo));

    }

    void ca_derdiffcoeff(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                         const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                         const int* lo, const int* hi,
                         const int* domain_lo, const int* domain_hi,
                         const Real* delta, const Real* xlo,
                         const Real* time, const Real* dt, const int* bcrec, 
                         const int* level, const int* grid_no)
    {

        derdiffcoeff(der, ARLIM_3D(der_lo), ARLIM_3D(der_hi), *nvar,
                     data, ARLIM_3D(data_lo), ARLIM_3D(data_hi), *ncomp,
                     ARLIM_3D(lo), ARLIM_3D(hi),
                     ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
                     ZFILL(delta), ZFILL(xlo));

    }

    void ca_derdiffterm(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                        const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const Real* delta, const Real* xlo,
                        const Real* time, const Real* dt, const int* bcrec, 
                        const int* level, const int* grid_no)
    {

        derdiffterm(der, ARLIM_3D(der_lo), ARLIM_3D(der_hi), *nvar,
                    data, ARLIM_3D(data_lo), ARLIM_3D(data_hi), *ncomp,
                    ARLIM_3D(lo), ARLIM_3D(hi),
                    ARLIM_3D(domain_lo), ARLIM_3D(domain_hi),
                    ZFILL(delta), ZFILL(xlo));

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
                         AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
               AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                  AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                     AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                  AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                  AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                   AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

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
                 AMREX_REAL_ANYD(delta), AMREX_REAL_ANYD(xlo));

    }

#ifdef __cplusplus
}
#endif
