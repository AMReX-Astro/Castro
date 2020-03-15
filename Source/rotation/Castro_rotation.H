
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


    void inertial_to_rotational_velocity(const int i, const int j, const int k,
                                         const amrex::GeometryData& geomdata,
                                         Real* center,
                                         Real* omega,
                                         const amrex::Real time, amrex::Real* v);
