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

    void
    electric_edge_x(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& q_arr,
                    amrex::Array4<amrex::Real> const& E,
                    amrex::Array4<amrex::Real const> const& flxy,
                    amrex::Array4<amrex::Real const> const& flxz);

    void
    electric_edge_y(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& q_arr,
                    amrex::Array4<amrex::Real> const& E,
                    amrex::Array4<amrex::Real const> const& flxx,
                    amrex::Array4<amrex::Real const> const& flxz);

    void
    electric_edge_z(const amrex::Box& bx,
                    amrex::Array4<amrex::Real const> const& q_arr,
                    amrex::Array4<amrex::Real> const& E,
                    amrex::Array4<amrex::Real const> const& flxx,
                    amrex::Array4<amrex::Real const> const& flxy);

    void
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    electric(amrex::Real* q_zone, amrex::Real& E_zone, const int comp);


    void
    corner_couple(const amrex::Box& bx,
                  amrex::Array4<amrex::Real> const& qr_out,
                  amrex::Array4<amrex::Real> const& ql_out,
                  amrex::Array4<amrex::Real const> const& ur,
                  amrex::Array4<amrex::Real const> const& ul,
                  amrex::Array4<amrex::Real const> const& flxd2,
                  amrex::Array4<amrex::Real const> const& Ed1,
                  amrex::Array4<amrex::Real const> const& Ed3,
                  const int d1, const int d2, const int d3,
                  const amrex::Real dt);

    void
    half_step(const amrex::Box& bx,
              amrex::Array4<amrex::Real> const& qr_out,
              amrex::Array4<amrex::Real> const& ql_out,
              amrex::Array4<amrex::Real const> const& ur,
              amrex::Array4<amrex::Real const> const& ul,
              amrex::Array4<amrex::Real const> const& flxd1,
              amrex::Array4<amrex::Real const> const& flxd2,
              amrex::Array4<amrex::Real const> const& Ed,
              amrex::Array4<amrex::Real const> const& Ed1,
              amrex::Array4<amrex::Real const> const& Ed2,
              const int d, const int d1, const int d2, const amrex::Real dt);

    void
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    ConsToPrim(amrex::Real* q_zone, amrex::Real* U_zone);


    void
    hlld(const amrex::Box& bx,
         amrex::Array4<amrex::Real const> const& qleft,
         amrex::Array4<amrex::Real const> const& qright,
         amrex::Array4<amrex::Real> const& flx,
         const int dir);


    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    PToC(amrex::Array1D<amrex::Real, 0, NQ>& q,
         amrex::Array1D<amrex::Real, 0, NUM_STATE+2> U,
         amrex::Real& gam1);


