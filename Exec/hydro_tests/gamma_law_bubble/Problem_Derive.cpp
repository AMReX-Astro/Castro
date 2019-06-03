#include "AMReX_REAL.H"

#include "Derive.H"
#include "Problem_Derive_F.H"
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

    void ca_derpi(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                  const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                  const int* lo, const int* hi,
                  const int* domain_lo, const int* domain_hi,
                  const Real* dx, const Real* xlo,
                  const Real* time, const Real* dt, const int* bcrec, 
                  const int* level, const int* grid_no)
    {

#pragma gpu
        derpi(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
              data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
              AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
              AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
              AMREX_REAL_ANYD(dx), *time);

    }


    void ca_derpioverp0(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                        const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const Real* dx, const Real* xlo,
                        const Real* time, const Real* dt, const int* bcrec, 
                        const int* level, const int* grid_no)
    {

#pragma gpu
        derpioverp0(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                    data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                    AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                    AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                    AMREX_REAL_ANYD(dx), *time);

    }


    void ca_derrhopert(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                       const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                       const int* lo, const int* hi,
                       const int* domain_lo, const int* domain_hi,
                       const Real* dx, const Real* xlo,
                       const Real* time, const Real* dt, const int* bcrec, 
                       const int* level, const int* grid_no)
    {

#pragma gpu
        derrhopert(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                   data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                   AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                   AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                   AMREX_REAL_ANYD(dx), *time);

    }


    void ca_dertpert(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                     const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const Real* dx, const Real* xlo,
                     const Real* time, const Real* dt, const int* bcrec, 
                     const int* level, const int* grid_no)
    {

#pragma gpu
        dertpert(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                 data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                 AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                 AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                 AMREX_REAL_ANYD(dx), *time);

    }


#ifdef __cplusplus
}
#endif
