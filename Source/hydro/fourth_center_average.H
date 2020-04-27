#ifndef _FOURTH_CENTER_AVERAGE_H_
#define _FOURTH_CENTER_AVERAGE_H_

Real compute_laplacian(int i, int j, int k, int n,
                       Array4<Real const> const a,
                       GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

Real trans_laplacian(int i, int j, int k, int n,
                     int idir,
                     Array4<Real const> const a,
                     GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

void make_cell_center(const Box& bx,
                      Array4<Real const> const& U,
                      Array4<Real> const& U_cc,
                      GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

void make_cell_center_in_place(const Box& bx,
                               Array4<Real> const& U,
                               Array4<Real> const& tmp,
                               GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

void compute_lap_term(const Box& bx,
                      Array4<Real const> const& U,
                      Array4<Real> const& lap, const int ncomp,
                      GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

void make_fourth_average(const Box& bx,
                         Array4<Real> const& q,
                         Array4<Real const> const& q_bar,
                         GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

void make_fourth_in_place(const Box& bx,
                          Array4<Real> const& q,
                          Array4<Real> const& tmp,
                          GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

void make_fourth_in_place_n(const Box& bx,
                            Array4<Real> const& q, const int ncomp,
                            Array4<Real> const& tmp,
                            GpuArray<Real, 3>& domlo, GpuArray<Real, 3>& domhi);

#endif
