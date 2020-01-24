
  void pm_compute_delta_mass
    (const int* lo, const int* hi, amrex::Real* pm_delta_mass,
     const BL_FORT_FAB_ARG_3D(state),
     const BL_FORT_FAB_ARG_3D(stateout),
     const BL_FORT_FAB_ARG_3D(volume),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);


  void pm_fix_solution
    (const int* lo, const int* hi,
     const BL_FORT_FAB_ARG_3D(stateold),
     BL_FORT_FAB_ARG_3D(statenew),
     const amrex::Real* problo, const amrex::Real* dx,
     const amrex::Real time, const amrex::Real dt);


  void set_pointmass(const amrex::Real* pointmass);


  void ca_gsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(s_old),
     const BL_FORT_FAB_ARG_3D(phi),
     const BL_FORT_FAB_ARG_3D(grav),
     BL_FORT_FAB_ARG_3D(source),
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time);

  void ca_corrgsrc
    (const int* lo, const int* hi,
     const int* domlo, const int* domhi,
     const BL_FORT_FAB_ARG_3D(S_old),
     const BL_FORT_FAB_ARG_3D(S_new),
     const BL_FORT_FAB_ARG_3D(phi_old),
     const BL_FORT_FAB_ARG_3D(phi_new),
     const BL_FORT_FAB_ARG_3D(grav_old),
     const BL_FORT_FAB_ARG_3D(grav_new),
     const BL_FORT_FAB_ARG_3D(volume),
     const BL_FORT_FAB_ARG_3D(xflux),
     const BL_FORT_FAB_ARG_3D(yflux),
     const BL_FORT_FAB_ARG_3D(zflux),
     const BL_FORT_FAB_ARG_3D(source),
     const amrex::Real* dx, const amrex::Real dt, const amrex::Real time);
