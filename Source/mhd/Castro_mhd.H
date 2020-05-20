    void just_the_mhd(amrex::Real time, amrex::Real dt);

    void
    check_for_mhd_cfl_violation(const amrex::Box& bx, const amrex::Real dt,
                                amrex::Array4<amrex::Real const> const& q_arr,
                                amrex::Array4<amrex::Real const> const& qaux_arr);

    void
    consup_mhd(const amrex::Box& bx,
               amrex::Array4<amrex::Real> const& update,
               amrex::Array4<amrex::Real const> const& flux0,
               amrex::Array4<amrex::Real const> const& flux1,
               amrex::Array4<amrex::Real const> const& flux2);

    void
    PrimToCons(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const& q_arr,
               amrex::Array4<amrex::Real> const& u_arr);

    void
    prim_half(const amrex::Box& bx,
              amrex::Array4<amrex::Real> const& q2D,
              amrex::Array4<amrex::Real const> const& q_arr,
              amrex::Array4<amrex::Real const> const& flxx,
              amrex::Array4<amrex::Real const> const& flxy,
              amrex::Array4<amrex::Real const> const& flxz,
              const amrex::Real dt);


    void
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    qflux(amrex::Real* qflx, amrex::Real* flx, amrex::Real* q_zone);
