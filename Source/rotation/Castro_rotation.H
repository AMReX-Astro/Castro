
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

///
/// Add the rotation source for input state uold to the full hydrodynamic source
/// term source
///
/// @param bx        the box to operate over
/// @param rot       the rotational acceleration
/// @param source    the full hydrodynamical source term
/// @param dt        current timestep
///
    void
    rsrc(const Box& bx,
         Array4<Real const> const& rot,
         Array4<Real const> const& uold,
         Array4<Real> const& source, 
         const Real dt);

///
/// Compute the correction source term for rotation
///
/// @param bx       the box to operate over
/// @param phi_old  rotational potential at old time
/// @param phi_new  rotational potential at new time
/// @param rold     rotational acceleration at old time
/// @param rnew     rotational acceleration at new time
/// @param uold     old time state
/// @param unew     new time state
/// @param source   full hydrodynamic state
/// @param flux0    mass flux in x coord dir
/// @param flux1    mass flux in y coord dir
/// @param flux2    mass flux in z coord dir
/// @param dt       current timestep
/// @param vol      cell volume
///
    void
    corrrsrc(const Box& bx,
             Array4<Real const> const& phi_old,
             Array4<Real const> const& phi_new,
             Array4<Real const> const& rold,
             Array4<Real const> const& rnew,
             Array4<Real const> const& uold,
             Array4<Real const> const& unew,
             Array4<Real> const& source,
             Array4<Real const> const& flux0,
             Array4<Real const> const& flux1,
             Array4<Real const> const& flux2,
             const Real dt,
             Array4<Real const> const& vol);

///
/// Compute the rotational acceleration over a box
///
/// @param bx         the box to operate over
/// @param rot        the rotational acceleration
/// @param state      the input state
////
    void
    fill_rotational_acceleration(const Box& bx,
                                 Array4<Real> const& rot,
                                 Array4<Real const> const& state);

///
/// Compute the rotational acceleration for a single zone
/// (Coriolis and centrifugal)
///
/// @param r         distance from origin of rotation vector
/// @param v         velocity
/// @param omega     rotation vector
/// @param coriolis  do we include the Coriolis force
/// @param Sr        rotational acceleration
///
    AMREX_GPU_HOST_DEVICE
    void
    rotational_acceleration(GpuArray<Real, 3>& r, GpuArray<Real, 3>& v,
                            GpuArray<Real, 3> const& omega,
                            const bool coriolis, Real* Sr);

    void
    fill_rotational_potential(const Box& bx,
                              Array4<Real> const& phi,
                              const Real time);


    void
    fill_rotational_psi(const Box& bx,
                        Array4<Real> const& psi,
                        const Real time);
