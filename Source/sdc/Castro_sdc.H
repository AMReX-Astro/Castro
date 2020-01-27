#ifndef AMREX_USE_CUDA
    void do_sdc_update(int m1, int m2, amrex::Real dt);
#endif

#ifdef REACTIONS
/// Take an input conserved state U_state, evaluate the instantaneous
/// reaction rates, and store the result in R_source/
///
/// Special consideration is needed for fourth_order accuracy:
///
///  * If we come in with U_state representing an average
///    (input_is_average = true), then, for fourth_order, we will
///    convert U_state to cell-centers (consuming one ghost cell),
///    burn on the cell-centers (including one ghost cell) and then
///    average back to cell-centers (consuming another ghosst cell).
///    This requires >= 2 ghost cells for U_state, which we assume
///    are valid.
///
///  * If we come in with U_state representing cell centers
///    (input_is_average = false), then we just return R_source
///    at the cell-centers.
///
    void construct_old_react_source(amrex::MultiFab& U_state,
                                    amrex::MultiFab& R_source,
                                    const bool input_is_average);
#endif

