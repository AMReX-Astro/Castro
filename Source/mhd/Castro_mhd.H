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

