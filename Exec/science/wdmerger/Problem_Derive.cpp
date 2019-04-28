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

    void ca_derinertialmomentumx(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                                 const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                                 const int* lo, const int* hi,
                                 const int* domain_lo, const int* domain_hi,
                                 const Real* dx, const Real* xlo,
                                 const Real* time, const Real* dt, const int* bcrec, 
                                 const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialmomentumx(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                             data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                             AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                             AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                             AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialmomentumy(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                                 const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                                 const int* lo, const int* hi,
                                 const int* domain_lo, const int* domain_hi,
                                 const Real* dx, const Real* xlo,
                                 const Real* time, const Real* dt, const int* bcrec, 
                                 const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialmomentumy(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                             data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                             AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                             AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                             AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialmomentumz(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                                 const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                                 const int* lo, const int* hi,
                                 const int* domain_lo, const int* domain_hi,
                                 const Real* dx, const Real* xlo,
                                 const Real* time, const Real* dt, const int* bcrec, 
                                 const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialmomentumz(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                             data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                             AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                             AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                             AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialangmomx(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                               const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const Real* xlo,
                               const Real* time, const Real* dt, const int* bcrec, 
                               const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialangmomx(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                           data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                           AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                           AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialangmomy(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                               const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const Real* xlo,
                               const Real* time, const Real* dt, const int* bcrec, 
                               const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialangmomy(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                           data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                           AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                           AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialangmomz(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                               const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const Real* xlo,
                               const Real* time, const Real* dt, const int* bcrec, 
                               const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialangmomz(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                           data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                           AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                           AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialradmomx(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                               const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const Real* xlo,
                               const Real* time, const Real* dt, const int* bcrec, 
                               const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialradmomx(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                           data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                           AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                           AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialradmomy(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                               const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const Real* xlo,
                               const Real* time, const Real* dt, const int* bcrec, 
                               const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialradmomy(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                           data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                           AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                           AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derinertialradmomz(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                               const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                               const int* lo, const int* hi,
                               const int* domain_lo, const int* domain_hi,
                               const Real* dx, const Real* xlo,
                               const Real* time, const Real* dt, const int* bcrec, 
                               const int* level, const int* grid_no)
    {

#pragma gpu
        derinertialradmomz(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                           data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                           AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                           AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                           AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derphieff(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                      const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const Real* dx, const Real* xlo,
                      const Real* time, const Real* dt, const int* bcrec, 
                      const int* level, const int* grid_no)
    {

#pragma gpu
        derphieff(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                  data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                  AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                  AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                  AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derphieffpm_p(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                          const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* time, const Real* dt, const int* bcrec, 
                          const int* level, const int* grid_no)
    {

#pragma gpu
        derphieffpm_p(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                      data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                      AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                      AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                      AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derphieffpm_s(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                          const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* time, const Real* dt, const int* bcrec, 
                          const int* level, const int* grid_no)
    {

#pragma gpu
        derphieffpm_s(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                      data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                      AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                      AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                      AMREX_REAL_ANYD(dx), *time);

    }


    void ca_derrhophiGrav(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                          const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const Real* dx, const Real* xlo,
                          const Real* time, const Real* dt, const int* bcrec, 
                          const int* level, const int* grid_no)
    {

#pragma gpu
        derrhophiGrav(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                      data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                      AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                      AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                      AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derrhophiRot(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                         const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                         const int* lo, const int* hi,
                         const int* domain_lo, const int* domain_hi,
                         const Real* dx, const Real* xlo,
                         const Real* time, const Real* dt, const int* bcrec, 
                         const int* level, const int* grid_no)
    {

#pragma gpu
        derrhophiRot(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                     data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                     AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                     AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                     AMREX_REAL_ANYD(dx), *time);

    }

    void ca_derprimarymask(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                           const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                           const int* lo, const int* hi,
                           const int* domain_lo, const int* domain_hi,
                           const Real* dx, const Real* xlo,
                           const Real* time, const Real* dt, const int* bcrec, 
                           const int* level, const int* grid_no)
    {

#pragma gpu
        derprimarymask(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                       data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                       AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                       AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                       AMREX_REAL_ANYD(dx), *time);

    }

    void ca_dersecondarymask(Real* der, const int* der_lo, const int* der_hi, const int* nvar,
                             const Real* data, const int* data_lo, const int* data_hi, const int* ncomp,
                             const int* lo, const int* hi,
                             const int* domain_lo, const int* domain_hi,
                             const Real* dx, const Real* xlo,
                             const Real* time, const Real* dt, const int* bcrec, 
                             const int* level, const int* grid_no)
    {

#pragma gpu
        dersecondarymask(der, AMREX_INT_ANYD(der_lo), AMREX_INT_ANYD(der_hi), *nvar,
                         data, AMREX_INT_ANYD(data_lo), AMREX_INT_ANYD(data_hi), *ncomp,
                         AMREX_INT_ANYD(lo), AMREX_INT_ANYD(hi),
                         AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi),
                         AMREX_REAL_ANYD(dx), *time);

    }
    
#ifdef __cplusplus
}
#endif
