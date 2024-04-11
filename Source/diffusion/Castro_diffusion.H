#ifndef CASTRO_DIFFUSION_H
#define CASTRO_DIFFUSION_H

///
/// Construct diffusion source at old time
///
/// @param source   MultiFab to save source terms to
/// @param state    Old state
/// @param time     current time
/// @param dt       timestep
///
    void construct_old_diff_source(amrex::MultiFab& source, amrex::MultiFab& state, amrex::Real time, amrex::Real dt);

///
/// Construct diffusion source at new time
///
/// @param source       MultiFab to save source terms to
/// @param state_old    Old state
/// @param state_new    New state
/// @param time         current time
/// @param dt           timestep
///
    void construct_new_diff_source(amrex::MultiFab& source, amrex::MultiFab& state_old, amrex::MultiFab& state_new, amrex::Real time, amrex::Real dt);


///
/// Get thermal conductivity diffusion term at given time
///
/// @param time     current time
/// @param state    Current state
/// @param DiffTerm MultiFab to save term to
///
void getTempDiffusionTerm (amrex::Real time, amrex::MultiFab& state, amrex::MultiFab& DiffTerm);


///
/// Calculate temperature or enthalpty diffusion terms and add to ``ext_src`` (multiplied by ``mult_factor``).
///
/// @param ext_src      Source terms to add diffusion sources to
/// @param state        Current state
/// @param DiffTerm     MultiFab to save diffusion sources to
/// @param t            Current time
/// @param mult_factor  Real, factor to multiply diffusion
///                     sources by before adding to other sources.
///
void add_temp_diffusion_to_source (amrex::MultiFab& ext_src, amrex::MultiFab& state,
                                   amrex::MultiFab& DiffTerm, amrex::Real t,
                                   amrex::Real mult_factor = 1.0);


#endif
