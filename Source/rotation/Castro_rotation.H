#ifndef CASTRO_ROTATION_H
#define CASTRO_ROTATION_H

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
    void fill_rotation_field(amrex::MultiFab& phi, amrex::MultiFab& state, amrex::Real time);

///
/// Add the rotation source for input state uold to the full hydrodynamic source
/// term source
///
/// @param bx        the box to operate over
/// @param source    the full hydrodynamical source term
/// @param dt        current timestep
///
    void
    rsrc(const Box& bx,
         Array4<Real const> const& uold,
         Array4<Real> const& source,
         const Real dt);

///
/// Compute the correction source term for rotation
///
/// @param bx       the box to operate over
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
             Array4<Real const> const& uold,
             Array4<Real const> const& unew,
             Array4<Real> const& source,
             Array4<Real const> const& flux0,
             Array4<Real const> const& flux1,
             Array4<Real const> const& flux2,
             const Real dt,
             Array4<Real const> const& vol);

    void
    fill_rotational_potential(const Box& bx,
                              Array4<Real> const& phi,
                              const Real time);


    void
    fill_rotational_psi(const Box& bx,
                        Array4<Real> const& psi,
                        const Real time);
#endif
