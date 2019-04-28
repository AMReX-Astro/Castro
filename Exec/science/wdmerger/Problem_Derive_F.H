#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void derinertialmomentumx(BL_FORT_FAB_ARG_3D(der), const int nvar,
                            const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                            const int* lo, const int* hi,
                            const int* domain_lo, const int* domain_hi,
                            const amrex::Real* dx, const amrex::Real time);

  void derinertialmomentumy(BL_FORT_FAB_ARG_3D(der), const int nvar,
                            const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                            const int* lo, const int* hi,
                            const int* domain_lo, const int* domain_hi,
                            const amrex::Real* dx, const amrex::Real time);

  void derinertialmomentumz(BL_FORT_FAB_ARG_3D(der), const int nvar,
                            const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                            const int* lo, const int* hi,
                            const int* domain_lo, const int* domain_hi,
                            const amrex::Real* dx, const amrex::Real time);

  void derinertialangmomx(BL_FORT_FAB_ARG_3D(der), const int nvar,
                          const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real time);

  void derinertialangmomy(BL_FORT_FAB_ARG_3D(der), const int nvar,
                          const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real time);

  void derinertialangmomz(BL_FORT_FAB_ARG_3D(der), const int nvar,
                          const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real time);

  void derinertialradmomx(BL_FORT_FAB_ARG_3D(der), const int nvar,
                          const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real time);

  void derinertialradmomy(BL_FORT_FAB_ARG_3D(der), const int nvar,
                          const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real time);

  void derinertialradmomz(BL_FORT_FAB_ARG_3D(der), const int nvar,
                          const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                          const int* lo, const int* hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* dx, const amrex::Real time);

  void derphieff(BL_FORT_FAB_ARG_3D(der), const int nvar,
                 const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                 const int* lo, const int* hi,
                 const int* domain_lo, const int* domain_hi,
                 const amrex::Real* dx, const amrex::Real time);

  void derphieffpm_p(BL_FORT_FAB_ARG_3D(der), const int nvar,
                     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const amrex::Real* dx, const amrex::Real time);

  void derphieffpm_s(BL_FORT_FAB_ARG_3D(der), const int nvar,
                     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const amrex::Real* dx, const amrex::Real time);

  void derrhophiGrav(BL_FORT_FAB_ARG_3D(der), const int nvar,
                     const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                     const int* lo, const int* hi,
                     const int* domain_lo, const int* domain_hi,
                     const amrex::Real* dx, const amrex::Real time);

  void derrhophiRot(BL_FORT_FAB_ARG_3D(der), const int nvar,
                    const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                    const int* lo, const int* hi,
                    const int* domain_lo, const int* domain_hi,
                    const amrex::Real* dx, const amrex::Real time);

  void derprimarymask(BL_FORT_FAB_ARG_3D(der), const int nvar,
                      const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                      const int* lo, const int* hi,
                      const int* domain_lo, const int* domain_hi,
                      const amrex::Real* dx, const amrex::Real time);

  void dersecondarymask(BL_FORT_FAB_ARG_3D(der), const int nvar,
                        const BL_FORT_FAB_ARG_3D(data), const int ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const amrex::Real* dx, const amrex::Real time);

#ifdef __cplusplus
}
#endif
