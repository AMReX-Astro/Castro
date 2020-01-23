
///
/// Calculate primitive variables from conserved variables (uses StateData)
///
/// @param time     current time
///
    void cons_to_prim(const amrex::Real time);

///
/// Calculate primitive variables from given conserved variables
///
/// @param u        MultiFab of conserved variables
/// @param q        MultiFab to save primitive variables to
/// @param qaux     MultiFab of auxiliary variables
/// @param time     Current time
///
    void cons_to_prim(amrex::MultiFab& u, amrex::MultiFab& q, amrex::MultiFab& qaux, const amrex::Real time);

///
/// convert the conservative state cell averages to primitive cell
/// averages with 4th order accuracy
///
/// @param time     current time
///
    void cons_to_prim_fourth(const amrex::Real time);

///
/// Check to see if the CFL condition has been violated
///
/// @param dt timestep
///
    void check_for_cfl_violation(const amrex::Real dt);

///
/// this constructs the hydrodynamic source (essentially the flux
/// divergence) using the CTU framework for unsplit hydrodynamics
///
/// @param time     current time
/// @param dt       timestep
///
    void construct_ctu_hydro_source(amrex::Real time, amrex::Real dt);

///
/// this constructs the hydrodynamic source (essentially the flux
/// divergence) using method of lines integration.  The output, is the
/// divergence of the fluxes, A = -div{F(U)}
///
/// @param time     current time
/// @param dt   timestep
/// @param A_update   divergence of the fluxes
///
    void construct_mol_hydro_source(amrex::Real time, amrex::Real dt, amrex::MultiFab& A_update);

    void src_to_prim(const amrex::Box& bx,
                     amrex::Array4<amrex::Real> const q, amrex::Array4<amrex::Real> const qaux,
                     amrex::Array4<amrex::Real> const src, amrex::Array4<amrex::Real> const srcQ);


    void add_hybrid_advection_source_c(const int i, const int j, const int k,
                                       const amrex::Real dt, const amrex::Real* dx,
                                       amrex::Array4<amrex::Real> const update,
                                       amrex::Array4<amrex::Real> const qx,
                                       amrex::Array4<amrex::Real> const qy,
                                       amrex::Array4<amrex::Real> const qz);
