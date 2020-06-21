
///
/// Construct rotation source terms at old time
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_rotation_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);


///
/// Construct rotation source terms at new time.
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old State
/// @param state_new    New State
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_rotation_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Fill rotational potential and acceleration
///
/// @param phi      Rotational potential
/// @param rot      Rotational acceleration
/// @param state    Current state
/// @param time     current time
///
    void fill_rotation_field(amrex::MultiFab& phi, amrex::MultiFab& rot, amrex::MultiFab& state, amrex::Real time);


    AMREX_GPU_HOST_DEVICE
    void inertial_to_rotational_velocity_c(const int i, const int j, const int k,
                                           const amrex::GeometryData& geomdata,
                                           const Real* center,
                                           const Real* omega,
                                           const amrex::Real time, amrex::Real* v);

    void
    rsrc(const Box& box,
         Array4<Real const> const& phi,
         Array4<Real const> const& rot,
         Array4<Real const> const& uold,
         Array4<Real> const& source, 
         Array4<Real const> const& vol,
         const Real dt, const Real time);

    void
    corrrsrc(const Box& bx,
             Array4<Real const> const& phi_old,
             Array4<Real const> const& phi_new,
             Array4<Real const> const& rold,
             Array4<Real const> const& rnew,
             Array4<Real const> const& uold,
             Array4<Real const> const& unew,
             Array4<Real> const& source,
             Array4<Real const> const& flux1,
             Array4<Real const> const& flux2,
             Array4<Real const> const& flux3,
             const Real dt, const Real time,
             Array4<Real const> const& vol);

    void
    rotational_acceleration(GpuArray<Real, 3>& r, GpuArray<Real, 3>& v,
                            GpuArray<Real, 3> const& omega,
                            const bool coriolis, Real* Sr);

    Real
    rotational_potential(GpuArray<Real, 3>& r, GpuArray<Real, 3> const& omega);

    void
    fill_rotational_potential(const Box& bx,
                              Array4<Real> const& phi,
                              const Real time);

    void
    fill_rotational_acceleration(const Box& bx,
                                 Array4<Real> const& rot,
                                 Array4<Real const> const& state,
                                 const Real time);

    void
    fill_rotational_psi(const Box& bx,
                        Array4<Real> const& psi,
                        const Real time);










