#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <Castro_F.H>
#include <fundamental_constants.H>
#include <prob_parameters.H>
#include <wdmerger_util.H>
#include <wdmerger_data.H>
#ifdef ROTATION
#include <Castro_rotation_F.H>
#endif

using namespace amrex;

void ca_derinertialmomentumx(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geomdata,
                             Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel = inertial_velocity(loc, vel, omega);
        der(i,j,k,0) = rho * inertial_vel[0];
    });
}

void ca_derinertialmomentumy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geomdata,
                             Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel = inertial_velocity(loc, vel, omega);
        der(i,j,k,0) = rho * inertial_vel[1];
    });
}

void ca_derinertialmomentumz(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geomdata,
                             Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel = inertial_velocity(loc, vel, omega);
        der(i,j,k,0) = rho * inertial_vel[2];
    });
}

void ca_derinertialangmomx(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geomdata,
                           Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive angular momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel = inertial_velocity(loc, vel, omega);

        GpuArray<Real, 3> angular_vel;
        cross_product(loc, inertial_vel, angular_vel);

        der(i,j,k,0) = rho * angular_vel[0];
    });
}

void ca_derinertialangmomy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geomdata,
                           Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive angular momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel = inertial_velocity(loc, vel, omega);

        GpuArray<Real, 3> angular_vel;
        cross_product(loc, inertial_vel, angular_vel);

        der(i,j,k,0) = rho * angular_vel[1];
    });
}

void ca_derinertialangmomz(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geomdata,
                           Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive angular momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel = inertial_velocity(loc, vel, omega);

        GpuArray<Real, 3> angular_vel;
        cross_product(loc, inertial_vel, angular_vel);

        der(i,j,k,0) = rho * angular_vel[2];
    });
}

void ca_derinertialradmomx(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geomdata,
                           Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive radial momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        GpuArray<Real, 3> mom{dat(i,j,k,1), dat(i,j,k,2), dat(i,j,k,3)};
        GpuArray<Real, 3> inertial_mom = inertial_velocity(loc, mom, omega);

        Real radInv = 1.0_rt / std::sqrt(loc[1] * loc[1] + loc[2] * loc[2]);

        der(i,j,k,0) = loc[1] * radInv * inertial_mom[1] + loc[2] * radInv * inertial_mom[2];
    });
}

void ca_derinertialradmomy(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geomdata,
                           Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive radial momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        GpuArray<Real, 3> mom{dat(i,j,k,1), dat(i,j,k,2), dat(i,j,k,3)};
        GpuArray<Real, 3> inertial_mom = inertial_velocity(loc, mom, omega);

        Real radInv = 1.0_rt / std::sqrt(loc[0] * loc[0] + loc[2] * loc[2]);

        der(i,j,k,0) = loc[0] * radInv * inertial_mom[0] + loc[2] * radInv * inertial_mom[2];
    });
}

void ca_derinertialradmomz(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geomdata,
                           Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive radial momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    GpuArray<Real, 3> center;
    ca_get_center(center.begin());

    GpuArray<Real, 3> omega;
    get_omega(omega.begin());

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - center[2];
#else
        loc[2] = 0.0_rt;
#endif

        GpuArray<Real, 3> mom{dat(i,j,k,1), dat(i,j,k,2), dat(i,j,k,3)};
        GpuArray<Real, 3> inertial_mom = inertial_velocity(loc, mom, omega);

        Real radInv = 1.0_rt / std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]);

        der(i,j,k,0) = loc[0] * radInv * inertial_mom[0] + loc[1] * radInv * inertial_mom[1];
    });
}

void ca_derphieff(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geomdata,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive the effective potential phiEff = phiGrav + phiRot

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,0) = dat(i,j,k,0) + dat(i,j,k,1);
    });
}

void ca_derphieffpm_p(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive an approximation to the effective potential of the primary only,
    // by treating it as a point-mass at its center of mass.
    // The u array contains the rotational potential, so we only need to calculate
    // the gravitational potential from the point-mass.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,0) = 0.0_rt;

        // Don't do anything here if the star no longer exists,
        // or if it never existed.

        if (mass_P <= 0.0_rt) return;

        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real r = std::sqrt((loc[0] - com_P[0]) * (loc[0] - com_P[0]) +
                           (loc[1] - com_P[1]) * (loc[1] - com_P[1]) +
                           (loc[2] - com_P[2]) * (loc[2] - com_P[2]));

        der(i,j,k,0) = -C::Gconst * mass_P / r + dat(i,j,k,0);
    });
}

void ca_derphieffpm_s(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Same as above, but for the secondary.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,0) = 0.0_rt;

        // Don't do anything here if the star no longer exists,
        // or if it never existed.

        if (mass_S <= 0.0_rt) return;

        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real r = std::sqrt((loc[0] - com_S[0]) * (loc[0] - com_S[0]) +
                           (loc[1] - com_S[1]) * (loc[1] - com_S[1]) +
                           (loc[2] - com_S[2]) * (loc[2] - com_S[2]));

        der(i,j,k,0) = -C::Gconst * mass_S / r + dat(i,j,k,0);
    });
}

void ca_derrhophiGrav(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geomdata,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,0) = dat(i,j,k,0) * dat(i,j,k,1);
    });
}

void ca_derrhophiRot(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& geomdata,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k,0) = dat(i,j,k,0) * dat(i,j,k,1);
    });
}

void ca_derprimarymask(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geomdata,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Create a mask for all zones considered to be within the primary star.
    // It uses the same prescription as above for the effective potential of the
    // star, and uses the stellar density threshold input parameter to determine
    // what parts of the domain should be considered stellar material.
    // The convention will be that the mask is positive (1) for zones inside the
    // star and negative (-1) for zones outside the star.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        // By default, assume we're not inside the star.

        der(i,j,k,0) = -1.0_rt;

        // Don't do anything here if the star no longer exists,
        // or if it never existed.

        if (mass_P <= 0.0_rt) return;

        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
        loc[2] = 0.0_rt;
#endif

        // Ignore zones whose density is too low.

        if (dat(i,j,k,0) < stellar_density_threshold) return;

        Real r_P = std::sqrt((loc[0] - com_P[0]) * (loc[0] - com_P[0]) +
                             (loc[1] - com_P[1]) * (loc[1] - com_P[1]) +
                             (loc[2] - com_P[2]) * (loc[2] - com_P[2]));

        Real r_S = std::sqrt((loc[0] - com_S[0]) * (loc[0] - com_S[0]) +
                             (loc[1] - com_S[1]) * (loc[1] - com_S[1]) +
                             (loc[2] - com_S[2]) * (loc[2] - com_S[2]));

        Real phi_p = -C::Gconst * mass_P / r_P + dat(i,j,k,1);
        Real phi_s = -C::Gconst * mass_S / r_S + dat(i,j,k,1);

        if (phi_p < 0.0_rt && phi_p < phi_s) {
            der(i,j,k,0) = 1.0_rt;
        }
    });
}

void ca_dersecondarymask(const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& datfab, const Geometry& geomdata,
                         Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Same as above, but for the secondary.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geomdata.CellSizeArray();
    const auto problo = geomdata.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k) noexcept
    {
        // By default, assume we're not inside the star.

        der(i,j,k,0) = -1.0_rt;

        // Don't do anything here if the star no longer exists,
        // or if it never existed.

        if (mass_S <= 0.0_rt) return;

        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2];
#else
        loc[2] = 0.0_rt;
#endif

        // Ignore zones whose density is too low.

        if (dat(i,j,k,0) < stellar_density_threshold) return;

        Real r_P = std::sqrt((loc[0] - com_P[0]) * (loc[0] - com_P[0]) +
                             (loc[1] - com_P[1]) * (loc[1] - com_P[1]) +
                             (loc[2] - com_P[2]) * (loc[2] - com_P[2]));

        Real r_S = std::sqrt((loc[0] - com_S[0]) * (loc[0] - com_S[0]) +
                             (loc[1] - com_S[1]) * (loc[1] - com_S[1]) +
                             (loc[2] - com_S[2]) * (loc[2] - com_S[2]));

        Real phi_p = -C::Gconst * mass_P / r_P + dat(i,j,k,1);
        Real phi_s = -C::Gconst * mass_S / r_S + dat(i,j,k,1);

        if (phi_s < 0.0_rt && phi_s < phi_p) {
            der(i,j,k,0) = 1.0_rt;
        }
    });
}
