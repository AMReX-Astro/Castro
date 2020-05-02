///
/// React state through a single reaction timestep, saving reaction source terms to ``reactions`` MultiFab. This is the CTU version.
///
/// @param state        Current state
/// @param reactions    MultiFab to save reaction sources to
/// @param time         current time
/// @param dt           reaction timestep
///
    bool react_state(amrex::MultiFab& state,
                     amrex::MultiFab& reactions,
                     amrex::Real time,
                     amrex::Real dt);

///
/// Simplified SDC version of react_state. Reacts the current state through a single timestep.
///
/// @param time     current time
/// @param dt       timestep
///
    bool react_state(amrex::Real time, amrex::Real dt, const int num_sdc_iters);

///
/// Are there any zones in ``State`` that can burn?
///
/// @param State    State MultiFab
///
    bool valid_zones_to_burn(amrex::MultiFab& State);
