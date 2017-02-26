extern "C"
{

  void ca_derinertialmomentumx(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			       const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			       const int* lo, const int* hi,
			       const int* domain_lo, const int* domain_hi,
			       const Real* delta, const Real* xlo,
			       const Real* time, const Real* dt, const int* bcrec, 
			       const int* level, const int* grid_no);

  void ca_derinertialmomentumy(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			       const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			       const int* lo, const int* hi,
			       const int* domain_lo, const int* domain_hi,
			       const Real* delta, const Real* xlo,
			       const Real* time, const Real* dt, const int* bcrec, 
			       const int* level, const int* grid_no);

  void ca_derinertialmomentumz(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			       const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			       const int* lo, const int* hi,
			       const int* domain_lo, const int* domain_hi,
			       const Real* delta, const Real* xlo,
			       const Real* time, const Real* dt, const int* bcrec, 
			       const int* level, const int* grid_no);

  void ca_derinertialangmomx(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			     const int* lo, const int* hi,
			     const int* domain_lo, const int* domain_hi,
			     const Real* delta, const Real* xlo,
			     const Real* time, const Real* dt, const int* bcrec, 
			     const int* level, const int* grid_no);

  void ca_derinertialangmomy(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			     const int* lo, const int* hi,
			     const int* domain_lo, const int* domain_hi,
			     const Real* delta, const Real* xlo,
			     const Real* time, const Real* dt, const int* bcrec, 
			     const int* level, const int* grid_no);

  void ca_derinertialangmomz(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			     const int* lo, const int* hi,
			     const int* domain_lo, const int* domain_hi,
			     const Real* delta, const Real* xlo,
			     const Real* time, const Real* dt, const int* bcrec, 
			     const int* level, const int* grid_no);

  void ca_derinertialradmomx(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			     const int* lo, const int* hi,
			     const int* domain_lo, const int* domain_hi,
			     const Real* delta, const Real* xlo,
			     const Real* time, const Real* dt, const int* bcrec, 
			     const int* level, const int* grid_no);

  void ca_derinertialradmomy(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			     const int* lo, const int* hi,
			     const int* domain_lo, const int* domain_hi,
			     const Real* delta, const Real* xlo,
			     const Real* time, const Real* dt, const int* bcrec, 
			     const int* level, const int* grid_no);

  void ca_derinertialradmomz(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			     const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			     const int* lo, const int* hi,
			     const int* domain_lo, const int* domain_hi,
			     const Real* delta, const Real* xlo,
			     const Real* time, const Real* dt, const int* bcrec, 
			     const int* level, const int* grid_no);

  void ca_derphieff(BL_FORT_FAB_ARG_3D(der),const int* nvar,
		    const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
		    const int* lo, const int* hi,
		    const int* domain_lo, const int* domain_hi,
		    const Real* delta, const Real* xlo,
		    const Real* time, const Real* dt, const int* bcrec,
		    const int* level, const int* grid_no);

  void ca_derphieffpm_p(BL_FORT_FAB_ARG_3D(der),const int* nvar,
                        const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			const int* lo, const int* hi,
			const int* domain_lo, const int* domain_hi,
			const Real* delta, const Real* xlo,
			const Real* time, const Real* dt, const int* bcrec,
			const int* level, const int* grid_no);

  void ca_derphieffpm_s(BL_FORT_FAB_ARG_3D(der),const int* nvar,
                        const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			const int* lo, const int* hi,
			const int* domain_lo, const int* domain_hi,
			const Real* delta, const Real* xlo,
			const Real* time, const Real* dt, const int* bcrec,
			const int* level, const int* grid_no);

  void ca_derrhophiGrav(BL_FORT_FAB_ARG_3D(der),const int* nvar,
                        const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                        const int* lo, const int* hi,
                        const int* domain_lo, const int* domain_hi,
                        const Real* delta, const Real* xlo,
                        const Real* time, const Real* dt, const int* bcrec,
                        const int* level, const int* grid_no);

  void ca_derrhophiRot(BL_FORT_FAB_ARG_3D(der),const int* nvar,
                       const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
                       const int* lo, const int* hi,
                       const int* domain_lo, const int* domain_hi,
                       const Real* delta, const Real* xlo,
                       const Real* time, const Real* dt, const int* bcrec,
                       const int* level, const int* grid_no);

  void ca_derprimarymask(BL_FORT_FAB_ARG_3D(der),const int* nvar,
                         const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			 const int* lo, const int* hi,
			 const int* domain_lo, const int* domain_hi,
			 const Real* delta, const Real* xlo,
			 const Real* time, const Real* dt, const int* bcrec,
			 const int* level, const int* grid_no);

  void ca_dersecondarymask(BL_FORT_FAB_ARG_3D(der),const int* nvar,
			   const BL_FORT_FAB_ARG_3D(data),const int* ncomp,
			   const int* lo, const int* hi,
			   const int* domain_lo, const int* domain_hi,
			   const Real* delta, const Real* xlo,
			   const Real* time, const Real* dt, const int* bcrec,
			   const int* level, const int* grid_no);

}
