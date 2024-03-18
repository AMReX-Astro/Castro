#ifndef problem_initialize_mhd_data_H
#define problem_initialize_mhd_data_H

#include <prob_parameters.H>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real A_x (int i, int j, int k,
          const GeometryData& geomdata)
{
    GpuArray<Real, 3> loc;

    bool ccx = false;
    bool ccy = true;
    bool ccz = true;
    position(i, j, k, geomdata, loc, ccx, ccy, ccz);

    Real x_P = loc[0] - problem::center_P_initial[0];
    Real y_P = loc[0] - problem::center_P_initial[1];
    Real z_P = loc[0] - problem::center_P_initial[2];

    Real x_S = loc[0] - problem::center_S_initial[0];
    Real y_S = loc[0] - problem::center_S_initial[1];
    Real z_S = loc[0] - problem::center_S_initial[2];

    Real dist_P = std::sqrt(std::pow(x_P, 2) + std::pow(y_P, 2) + std::pow(z_P, 2));
    Real dist_S = std::sqrt(std::pow(x_S, 2) + std::pow(y_S, 2) + std::pow(z_S, 2));

    Real theta_P = std::atan2(std::sqrt(x_P * x_P + y_P * y_P), z_P);
    Real phi_P = std::atan2(y_P, x_P);

    Real theta_S = std::atan2(std::sqrt(x_S * x_S + y_S * y_S), z_S);
    Real phi_S = std::atan2(y_S, x_S);

    Real A;

    if (problem::mass_P > 0.0_rt && (dist_P < problem::radius_P)) {

        // Vector potential inside the primary star

        A = -problem::m_0 * dist_P * std::sin(theta_P) * std::sin(phi_P) / std::pow(problem::radius_P, 3);

    }
    else if (problem::mass_S > 0.0_rt && (dist_S < problem::radius_S)) {

        // Vector potential inside the secondary star

        A = -problem::m_0 * dist_S * std::sin(theta_S) * std::sin(phi_S) / std::pow(problem::radius_S, 3);

    }
    else  {

        // Dipole outside the primary star

        A = -problem::m_0 * std::sin(theta_P) * std::sin(phi_P) / (dist_P * dist_P);

        // Add contribution of dipole outside the secondary star

        A += -problem::m_0 * std::sin(theta_S) * std::sin(phi_S) / (dist_S * dist_S);

    }

    return A;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real A_y (int i, int j, int k,
          const GeometryData& geomdata)
{
    GpuArray<Real, 3> loc;

    bool ccx = true;
    bool ccy = false;
    bool ccz = true;
    position(i, j, k, geomdata, loc, ccx, ccy, ccz);

    Real x_P = loc[0] - problem::center_P_initial[0];
    Real y_P = loc[0] - problem::center_P_initial[1];
    Real z_P = loc[0] - problem::center_P_initial[2];

    Real x_S = loc[0] - problem::center_S_initial[0];
    Real y_S = loc[0] - problem::center_S_initial[1];
    Real z_S = loc[0] - problem::center_S_initial[2];

    Real dist_P = std::sqrt(std::pow(x_P, 2) + std::pow(y_P, 2) + std::pow(z_P, 2));
    Real dist_S = std::sqrt(std::pow(x_S, 2) + std::pow(y_S, 2) + std::pow(z_S, 2));

    Real theta_P = std::atan2(std::sqrt(x_P * x_P + y_P * y_P), z_P);
    Real phi_P = std::atan2(y_P, x_P);

    Real theta_S = std::atan2(std::sqrt(x_S * x_S + y_S * y_S), z_S);
    Real phi_S = std::atan2(y_S, x_S);

    Real A;

    if (problem::mass_P > 0.0_rt && (dist_P < problem::radius_P)) {

        // Vector potential inside the primary star

        A = -problem::m_0 * dist_P * std::sin(theta_P) * std::cos(phi_P) / std::pow(problem::radius_P, 3);

    }
    else if (problem::mass_S > 0.0_rt && (dist_S < problem::radius_S)) {

        // Vector potential inside the secondary star

        A = -problem::m_0 * dist_S * std::sin(theta_S) * std::cos(phi_S) / std::pow(problem::radius_S, 3);

    }
    else  {

        // Dipole outside the primary star

        A = -problem::m_0 * std::sin(theta_P) * std::cos(phi_P) / (dist_P * dist_P);

        // Add contribution of dipole outside the secondary star

        A += -problem::m_0 * std::sin(theta_P) * std::cos(phi_P) / (dist_S * dist_S);

    }

    return A;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_initialize_mhd_data (int i, int j, int k,
                                  Array4<Real> const& B, const int idir,
                                  const GeometryData& geomdata)
{
    // A is of the form A = (A_x, A_y, 0)
    // B_x = -dA_x/dz
    // B_y = dA_y/dz
    // B_z = dA_y/dx - dA_x/dy

    const Real* dx = geomdata.CellSize();

    if (idir == 0) {
        B(i,j,k) = (A_y(i,j,k,geomdata) - A_y(i,j,k+1,geomdata)) / dx[2];
    }
    else if (idir == 1) {
        B(i,j,k) = (A_x(i,j,k+1,geomdata) - A_x(i,j,k,geomdata)) / dx[2];
    }
    else {
        B(i,j,k) = (A_y(i+1,j,k,geomdata) - A_y(i,j,k,geomdata)) / dx[0] -
                   (A_x(i,j+1,k,geomdata) - A_x(i,j,k,geomdata)) / dx[1];
    }
}

#endif
