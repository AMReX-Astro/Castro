#ifndef CASTRO_SDC_H
#define CASTRO_SDC_H

#ifndef AMREX_USE_GPU
void do_sdc_update(int m_start, int m_end, amrex::Real dt);
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

void ca_sdc_update_advection_o2_lobatto(const amrex::Box& bx,
                                        amrex::Real dt_m, amrex::Real dt,
                                        amrex::Array4<const amrex::Real> const& k_m,
                                        amrex::Array4<amrex::Real> const& k_n,
                                        amrex::Array4<const amrex::Real> const& A_m,
                                        amrex::Array4<const amrex::Real> const& A_0_old,
                                        amrex::Array4<const amrex::Real> const& A_1_old,
                                        int m_start);

void ca_sdc_update_advection_o2_radau(const amrex::Box& bx,
                                      amrex::Real dt_m, amrex::Real dt,
                                      amrex::Array4<const amrex::Real> const& k_m,
                                      amrex::Array4<amrex::Real> const& k_n,
                                      amrex::Array4<const amrex::Real> const& A_m,
                                      amrex::Array4<const amrex::Real> const& A_0_old,
                                      amrex::Array4<const amrex::Real> const& A_1_old,
                                      amrex::Array4<const amrex::Real> const& A_2_old,
                                      int m_start);

void ca_sdc_update_advection_o4_lobatto(const amrex::Box& bx,
                                        amrex::Real dt_m, amrex::Real dt,
                                        amrex::Array4<const amrex::Real> const& k_m,
                                        amrex::Array4<amrex::Real> const& k_n,
                                        amrex::Array4<const amrex::Real> const& A_m,
                                        amrex::Array4<const amrex::Real> const& A_0_old,
                                        amrex::Array4<const amrex::Real> const& A_1_old,
                                        amrex::Array4<const amrex::Real> const& A_2_old,
                                        int m_start);

void ca_sdc_update_advection_o4_radau(const amrex::Box& bx,
                                      amrex::Real dt_m, amrex::Real dt,
                                      amrex::Array4<const amrex::Real> const& k_m,
                                      amrex::Array4<amrex::Real> const& k_n,
                                      amrex::Array4<const amrex::Real> const& A_m,
                                      amrex::Array4<const amrex::Real> const& A_0_old,
                                      amrex::Array4<const amrex::Real> const& A_1_old,
                                      amrex::Array4<const amrex::Real> const& A_2_old,
                                      amrex::Array4<const amrex::Real> const& A_3_old,
                                      int m_start);
#ifdef REACTIONS
void ca_sdc_compute_C2_lobatto(const amrex::Box& bx,
                               amrex::Real dt_m, amrex::Real dt,
                               amrex::Array4<const amrex::Real> const& A_m,
                               amrex::Array4<const amrex::Real> const& A_0_old,
                               amrex::Array4<const amrex::Real> const& A_1_old,
                               amrex::Array4<const amrex::Real> const& R_0_old,
                               amrex::Array4<const amrex::Real> const& R_1_old,
                               amrex::Array4<amrex::Real> const& C,
                               int m_start);

void ca_sdc_compute_C2_radau(const amrex::Box& bx,
                             amrex::Real dt_m, amrex::Real dt,
                             amrex::Array4<const amrex::Real> const& A_m,
                             amrex::Array4<const amrex::Real> const& A_0_old,
                             amrex::Array4<const amrex::Real> const& A_1_old,
                             amrex::Array4<const amrex::Real> const& A_2_old,
                             amrex::Array4<const amrex::Real> const& R_0_old,
                             amrex::Array4<const amrex::Real> const& R_1_old,
                             amrex::Array4<const amrex::Real> const& R_2_old,
                             amrex::Array4<amrex::Real> const& C,
                             int m_start);

void ca_sdc_compute_C4_lobatto(const amrex::Box& bx,
                               amrex::Real dt_m, amrex::Real dt,
                               amrex::Array4<const amrex::Real> const& A_m,
                               amrex::Array4<const amrex::Real> const& A_0_old,
                               amrex::Array4<const amrex::Real> const& A_1_old,
                               amrex::Array4<const amrex::Real> const& A_2_old,
                               amrex::Array4<const amrex::Real> const& R_0_old,
                               amrex::Array4<const amrex::Real> const& R_1_old,
                               amrex::Array4<const amrex::Real> const& R_2_old,
                               amrex::Array4<amrex::Real> const& C,
                               int m_start);

void ca_sdc_compute_C4_radau(const amrex::Box& bx,
                             amrex::Real dt_m, amrex::Real dt,
                             amrex::Array4<const amrex::Real> const& A_m,
                             amrex::Array4<const amrex::Real> const& A_0_old,
                             amrex::Array4<const amrex::Real> const& A_1_old,
                             amrex::Array4<const amrex::Real> const& A_2_old,
                             amrex::Array4<const amrex::Real> const& A_3_old,
                             amrex::Array4<const amrex::Real> const& R_0_old,
                             amrex::Array4<const amrex::Real> const& R_1_old,
                             amrex::Array4<const amrex::Real> const& R_2_old,
                             amrex::Array4<const amrex::Real> const& R_3_old,
                             amrex::Array4<amrex::Real> const& C,
                             int m_start);

void ca_sdc_conservative_update(const amrex::Box& bx, amrex::Real const dt_m,
                                amrex::Array4<const amrex::Real> const& U_old,
                                amrex::Array4<amrex::Real> const& U_new,
                                amrex::Array4<const amrex::Real> const& C,
                                amrex::Array4<const amrex::Real> const& R_new);
#endif

void ca_sdc_compute_initial_guess(const amrex::Box& bx,
                                  amrex::Array4<const amrex::Real> const& U_old,
                                  amrex::Array4<const amrex::Real> const& U_new,
                                  amrex::Array4<const amrex::Real> const& A_old,
                                  amrex::Array4<const amrex::Real> const& R_old,
                                  amrex::Array4<amrex::Real> const& U_guess,
                                  amrex::Real const dt_m, int const sdc_iteration);

#ifdef REACTIONS
void ca_store_reaction_state(const amrex::Box& bx,
                             amrex::Array4<const amrex::Real> const& R_old,
                             amrex::Array4<amrex::Real> const& R_store);
#endif

#endif
