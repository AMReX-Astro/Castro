#include <AMReX_REAL.H>

#include <Derive.H>
#include <Castro.H>
#include <fundamental_constants.H>
#include <prob_parameters.H>
#include <wdmerger_util.H>
#include <wdmerger_data.H>
#ifdef ROTATION
#include <Rotation.H>
#endif

using namespace amrex;

void ca_derinertialmomentumx(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geom,
                             Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif
        der(i,j,k,0) = rho * inertial_vel[0];
    });
}

void ca_derinertialmomentumy(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geom,
                             Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif
        der(i,j,k,0) = rho * inertial_vel[1];
    });
}

void ca_derinertialmomentumz(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                             const FArrayBox& datfab, const Geometry& geom,
                             Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif
        der(i,j,k,0) = rho * inertial_vel[2];
    });
}

void ca_derinertialangmomx(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geom,
                           Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive angular momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif

        GpuArray<Real, 3> angular_vel;
        cross_product(loc, inertial_vel, angular_vel);

        der(i,j,k,0) = rho * angular_vel[0];
    });
}

void ca_derinertialangmomy(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geom,
                           Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive angular momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif

        GpuArray<Real, 3> angular_vel;
        cross_product(loc, inertial_vel, angular_vel);

        der(i,j,k,0) = rho * angular_vel[1];
    });
}

void ca_derinertialangmomz(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geom,
                           Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive angular momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        Real rho = dat(i,j,k,0);
        GpuArray<Real, 3> vel{dat(i,j,k,1) / rho, dat(i,j,k,2) / rho, dat(i,j,k,3) / rho};
        GpuArray<Real, 3> inertial_vel{vel};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_vel);
#endif

        GpuArray<Real, 3> angular_vel;
        cross_product(loc, inertial_vel, angular_vel);

        der(i,j,k,0) = rho * angular_vel[2];
    });
}

void ca_derinertialradmomx(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geom,
                           Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive radial momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        GpuArray<Real, 3> mom{dat(i,j,k,1), dat(i,j,k,2), dat(i,j,k,3)};
        GpuArray<Real, 3> inertial_mom{mom};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_mom);
#endif

        Real radInv = 1.0_rt / std::sqrt(loc[1] * loc[1] + loc[2] * loc[2]);

        der(i,j,k,0) = loc[1] * radInv * inertial_mom[1] + loc[2] * radInv * inertial_mom[2];
    });
}

void ca_derinertialradmomy(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geom,
                           Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive radial momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        GpuArray<Real, 3> mom{dat(i,j,k,1), dat(i,j,k,2), dat(i,j,k,3)};
        GpuArray<Real, 3> inertial_mom{mom};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_mom);
#endif

        Real radInv = 1.0_rt / std::sqrt(loc[0] * loc[0] + loc[2] * loc[2]);

        der(i,j,k,0) = loc[0] * radInv * inertial_mom[0] + loc[2] * radInv * inertial_mom[2];
    });
}

void ca_derinertialradmomz(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                           const FArrayBox& datfab, const Geometry& geom,
                           Real time, const int* /*bcrec*/, int /*level*/)
{
    // Derive radial momentum, given the grid momenta.

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        GpuArray<Real, 3> mom{dat(i,j,k,1), dat(i,j,k,2), dat(i,j,k,3)};
        GpuArray<Real, 3> inertial_mom{mom};
#ifdef ROTATION
        rotational_to_inertial_velocity(i, j, k, geomdata, time, inertial_mom);
#endif

        Real radInv = 1.0_rt / std::sqrt(loc[0] * loc[0] + loc[1] * loc[1]);

        der(i,j,k,0) = loc[0] * radInv * inertial_mom[0] + loc[1] * radInv * inertial_mom[1];
    });
}

void ca_derphieff(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                  const FArrayBox& datfab, const Geometry& geom,
                  Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive the effective potential phiEff = phiGrav + phiRot

    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        der(i,j,k,0) = dat(i,j,k,0) + rotational_potential(loc);
    });
}

void ca_derphieffpm_p(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geom,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Derive an approximation to the effective potential of the primary only,
    // by treating it as a point-mass at its center of mass.
    // The u array contains the rotational potential, so we only need to calculate
    // the gravitational potential from the point-mass.

    amrex::ignore_unused(datfab);

    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        der(i,j,k,0) = 0.0_rt;

        // Don't do anything here if the star no longer exists,
        // or if it never existed.

        if (problem::mass_P <= 0.0_rt) return;

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

        Real r = std::sqrt((loc[0] - problem::com_P[0]) * (loc[0] - problem::com_P[0]) +
                           (loc[1] - problem::com_P[1]) * (loc[1] - problem::com_P[1]) +
                           (loc[2] - problem::com_P[2]) * (loc[2] - problem::com_P[2]));

        for (int iloc = 0; iloc < AMREX_SPACEDIM; ++iloc) {
            loc[iloc] -= problem::center[iloc];
        }

        der(i,j,k,0) = -C::Gconst * problem::mass_P / r + rotational_potential(loc);
    });
}

void ca_derphieffpm_s(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geom,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    // Same as above, but for the secondary.

    amrex::ignore_unused(datfab);

    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        der(i,j,k,0) = 0.0_rt;

        // Don't do anything here if the star no longer exists,
        // or if it never existed.

        if (problem::mass_S <= 0.0_rt) return;

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

        Real r = std::sqrt((loc[0] - problem::com_S[0]) * (loc[0] - problem::com_S[0]) +
                           (loc[1] - problem::com_S[1]) * (loc[1] - problem::com_S[1]) +
                           (loc[2] - problem::com_S[2]) * (loc[2] - problem::com_S[2]));

        for (int iloc = 0; iloc < AMREX_SPACEDIM; ++iloc) {
            loc[iloc] -= problem::center[iloc];
        }

        der(i,j,k,0) = -C::Gconst * problem::mass_S / r + rotational_potential(loc);
    });
}

void ca_derrhophiGrav(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                      const FArrayBox& datfab, const Geometry& geom,
                      Real /*time*/, const int* /*bcrec*/, int /*level*/)
{

    amrex::ignore_unused(geom);

    auto const dat = datfab.array();
    auto const der = derfab.array();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        der(i,j,k,0) = dat(i,j,k,0) * dat(i,j,k,1);
    });
}

void ca_derrhophiRot(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                     const FArrayBox& datfab, const Geometry& geom,
                     Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    auto const dat = datfab.array();
    auto const der = derfab.array();

    const auto dx = geom.CellSizeArray();
    const auto problo = geom.ProbLoArray();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {

        GpuArray<Real, 3> loc;
        loc[0] = problo[0] + (static_cast<Real>(i) + 0.5_rt) * dx[0] - problem::center[0];

#if AMREX_SPACEDIM >= 2
        loc[1] = problo[1] + (static_cast<Real>(j) + 0.5_rt) * dx[1] - problem::center[1];
#else
        loc[1] = 0.0_rt;
#endif

#if AMREX_SPACEDIM == 3
        loc[2] = problo[2] + (static_cast<Real>(k) + 0.5_rt) * dx[2] - problem::center[2];
#else
        loc[2] = 0.0_rt;
#endif

        der(i,j,k,0) = dat(i,j,k,0) * rotational_potential(loc);
    });
}

void ca_derprimarymask(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                       const FArrayBox& datfab, const Geometry& geom,
                       Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    Array4<Real const> const rho = datfab.array();
    Array4<Real> const mask = derfab.array();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        const bool is_primary = true;
        mask(i,j,k) = stellar_mask(i, j, k, geomdata, rho(i,j,k), is_primary);
    });
}

void ca_dersecondarymask(const Box& bx, FArrayBox& derfab, int /*dcomp*/, int /*ncomp*/,
                         const FArrayBox& datfab, const Geometry& geom,
                         Real /*time*/, const int* /*bcrec*/, int /*level*/)
{
    Array4<Real const> const rho = datfab.array();
    Array4<Real> const mask = derfab.array();

    const auto geomdata = geom.data();

    amrex::ParallelFor(bx,
    [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
    {
        bool is_primary = false;
        mask(i,j,k) = stellar_mask(i, j, k, geomdata, rho(i,j,k), is_primary);
    });
}
