
BL_FORT_PROC_DECL(CA_DERPI,ca_derpi)
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERPIOVERP0,ca_derpioverp0)
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERRHOPERT,ca_derrhopert)
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

BL_FORT_PROC_DECL(CA_DERTPERT,ca_dertpert)
    (BL_FORT_FAB_ARG_3D(der),const int* nvar,
     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
     const int* lo, const int* hi,
     const int* domain_lo, const int* domain_hi,
     const Real* delta, const Real* xlo,
     const Real* time, const Real* dt, const int* bcrec, 
     const int* level, const int* grid_no);

