///
/// React state through a single reaction timestep, saving reaction source terms to ``reactions`` MultiFab. This is the CTU version.
///
/// @param state        Current state
/// @param reactions    MultiFab to save reaction sources to
/// @param mask
/// @param weights
/// @param time         current time
/// @param dt_react     reaction timestep
/// @param strang_half  Is this a Strang-split half timestep?
/// @param ngrow        Number of ghost cells
///
    bool react_state(amrex::MultiFab& state,
		     amrex::MultiFab& reactions,
		     const amrex::iMultiFab& mask,
		     amrex::MultiFab& weights,
		     amrex::Real time,
		     amrex::Real dt_react,
		     int strang_half, int  ngrow = 0);

///
/// Simplified SDC version of react_state. Reacts the current state through a single timestep.
///
/// @param time     current time
/// @param dt       timestep
///
    bool react_state(amrex::Real time, amrex::Real dt);

///
/// React through first half timestep
///
/// @param time     current time
/// @param dt       timestep
///
    bool strang_react_first_half(amrex::Real time, amrex::Real dt);


///
/// React through second half stimestep
///
/// @param time     current time
/// @param dt       timestep
///
    bool strang_react_second_half(amrex::Real time, amrex::Real dt);

///
/// Are there any zones in ``State`` that can burn?
///
/// @param State    State MultiFab
///
    bool valid_zones_to_burn(amrex::MultiFab& State);
