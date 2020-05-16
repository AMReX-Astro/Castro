    void just_the_mhd(amrex::Real time, amrex::Real dt);

    void
    mhd_speeds(const amrex::Box& bx,
               amrex::Array4<amrex::Real const> const& Bx,
               amrex::Array4<amrex::Real const> const& By,
               amrex::Array4<amrex::Real const> const& Bz,
               amrex::Array4<amrex::Real const> const& q_arr,
               amrex::Array4<amrex::Real const> const& qaux_arr,
               amrex::Array4<amrex::Real> const& cx,
               amrex::Array4<amrex::Real> const& cy,
               amrex::Array4<amrex::Real> const& cz);
    void
    consup_mhd(const amrex::Box& bx,
               amrex::Array4<amrex::Real> const& update,
               amrex::Array4<amrex::Real const> const& flux0,
               amrex::Array4<amrex::Real const> const& flux1,
               amrex::Array4<amrex::Real const> const& flux2);

