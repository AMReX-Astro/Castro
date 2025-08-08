#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_MultiFab.H>

#include <global.H>

#include <eos.H>

#ifdef DIFFUSION
#include <conductivity.H>
#endif

#ifdef MHD
#include <mhd_util.H>
#endif

#ifdef ROTATION
#include <Rotation.H>
#endif

#ifdef REACTIONS
#include <actual_network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#else
#include <actual_rhs.H>
#endif
#endif

#ifdef RADIATION
#include <Radiation.H>
#endif

namespace timestep {

    enum estdt_type { cfl = 0,
                      mhd,
                      burning,
                      radiation,
                      diffusion
    };

    // Estimate timestep given various constraints

    template <estdt_type T, typename... MFs>
    AMREX_INLINE
    amrex::ValLocPair<amrex::Real, amrex::IntVect>
    estdt (const amrex::GeometryData& geomdata,
           const amrex::MultiFab& maskMF,
           const amrex::MultiFab& stateMF,
           const MFs&... mfs)
    {
        amrex::ignore_unused(maskMF);

        constexpr std::size_t N_mfs = sizeof...(MFs);

        const auto* dx = geomdata.CellSize();
        const auto* problo = geomdata.ProbLo();
        const auto coord = geomdata.Coord();
        amrex::ignore_unused(problo, coord);

        auto const& ua = stateMF.const_arrays();

        bool mask_covered_zones = maskMF.isDefined();

        MultiArray4<Real const> empty_arr{};
        const auto& ma = mask_covered_zones ? maskMF.const_arrays() : empty_arr;

        if constexpr (T == cfl) {

            // Validate we don't have any extra MultiFabs.
            static_assert(N_mfs == 0);

            // Courant-condition limited timestep

            auto r = amrex::ParReduce(TypeList<ReduceOpMin>{}, TypeList<ValLocPair<Real, IntVect>>{}, stateMF,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> GpuTuple<ValLocPair<Real, IntVect>>
            {
                Array4<Real const> const& u = ua[box_no];

                IntVect idx(AMREX_D_DECL(i,j,k));

                Real rhoInv = 1.0_rt / u(i,j,k,URHO);

                eos_rep_t eos_state;
                eos_state.rho = u(i,j,k,URHO);
                eos_state.T = u(i,j,k,UTEMP);
                eos_state.e = u(i,j,k,UEINT) * rhoInv;
                for (int n = 0; n < NumSpec; n++) {
                    eos_state.xn[n] = u(i,j,k,UFS+n) * rhoInv;
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                    eos_state.aux[n] = u(i,j,k,UFX+n) * rhoInv;
                }
#endif

                eos(eos_input_re, eos_state);

                // Compute velocity and then calculate CFL timestep.

                Real ux = u(i,j,k,UMX) * rhoInv;
#if AMREX_SPACEDIM >= 2
                Real uy = u(i,j,k,UMY) * rhoInv;
#endif
#if AMREX_SPACEDIM == 3
                Real uz = u(i,j,k,UMZ) * rhoInv;
#endif

                Real c = eos_state.cs;

                Real dt1 = dx[0]/(c + std::abs(ux));

                Real dt2;
#if AMREX_SPACEDIM >= 2
                dt2 = dx[1]/(c + std::abs(uy));
                if (coord == 2) {
                    // dx[1] in Spherical2D is just dtheta, need rdtheta for physical length
                    // so just multiply by the smallest r
                    dt2 *= problo[0] + 0.5_rt * dx[0];
                }
#else
                dt2 = dt1;
#endif

                Real dt3;
#if AMREX_SPACEDIM == 3
                dt3 = dx[2]/(c + std::abs(uz));
#else
                dt3 = dt1;
#endif

                // The CTU method has a less restrictive timestep than MOL-based
                // schemes (including the true SDC).  Since the simplified SDC
                // solver is based on CTU, we can use its timestep.
                if (castro::time_integration_method == 0 || castro::time_integration_method == 3) {
                    Real dt = amrex::min(dt1, dt2, dt3);
                    return {ValLocPair<Real, IntVect>{dt, idx}};

                } else {
                    // method of lines-style constraint is tougher
                    Real dt_tmp = 1.0_rt/dt1;
#if AMREX_SPACEDIM >= 2
                    dt_tmp += 1.0_rt/dt2;
#endif
#if AMREX_SPACEDIM == 3
                    dt_tmp += 1.0_rt/dt3;
#endif

                    return {ValLocPair<Real, IntVect>{1.0_rt/dt_tmp, idx}};
                }
            });

            return r;

        } // cfl
#ifdef MHD
        else if constexpr (T == mhd) {

            // MHD timestep limiter

            // Validate we have three extra MultiFabs (Bx, By, and Bz).

            static_assert(N_mfs == 3);

            auto const& [bx, by, bz] = std::tie(mfs...);

            auto const& bxa = bx.const_arrays();
            auto const& bya = by.const_arrays();
            auto const& bza = bz.const_arrays();

            auto r = amrex::ParReduce(TypeList<ReduceOpMin>{}, TypeList<ValLocPair<Real, IntVect>>{}, stateMF,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> GpuTuple<ValLocPair<Real, IntVect>>
            {

                Array4<Real const> const& u_arr = ua[box_no];

                Array4<Real const> const& bx_arr = bxa[box_no];
                Array4<Real const> const& by_arr = bya[box_no];
                Array4<Real const> const& bz_arr = bza[box_no];

                IntVect idx(AMREX_D_DECL(i,j,k));

                Real rhoInv = 1.0_rt / u_arr(i,j,k,URHO);
                Real bcx = 0.5_rt * (bx_arr(i+1,j,k) + bx_arr(i,j,k));
                Real bcy = 0.5_rt * (by_arr(i,j+1,k) + by_arr(i,j,k));
                Real bcz = 0.5_rt * (bz_arr(i,j,k+1) + bz_arr(i,j,k));

                Real ux = u_arr(i,j,k,UMX) * rhoInv;
                Real uy = u_arr(i,j,k,UMY) * rhoInv;
                Real uz = u_arr(i,j,k,UMZ) * rhoInv;

                eos_rep_t eos_state;
                eos_state.rho = u_arr(i,j,k,URHO);
                eos_state.e = u_arr(i,j,k,UEINT) * rhoInv;
                eos_state.T = u_arr(i,j,k,UTEMP);
                for (int n = 0; n < NumSpec; n++) {
                  eos_state.xn[n] = u_arr(i,j,k,UFS+n) * rhoInv;
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; n++) {
                  eos_state.aux[n] = u_arr(i,j,k,UFX+n) * rhoInv;
                }
#endif

                eos(eos_input_re, eos_state);

                Real e  = u_arr(i,j,k,UEINT) * rhoInv;

                Real as = eos_state.gam1 * eos_state.p * rhoInv;
                Real ca = (bcx*bcx + bcy*bcy + bcz*bcz) * rhoInv;

                Real cx;
                Real cy;
                Real cz;

                if (e > 0_rt) {
                  Real cad = bcx*bcx * rhoInv;
                  eos_soundspeed_mhd(cx, as, ca, cad);

                  cad = bcy*bcy * rhoInv;
                  eos_soundspeed_mhd(cy, as, ca, cad);

                  cad = bcz*bcz * rhoInv;
                  eos_soundspeed_mhd(cz, as, ca, cad);

                } else {
                  cx = 0.0_rt;
                  cy = 0.0_rt;
                  cz = 0.0_rt;
                }

                Real dt1 = dx[0]/(cx + std::abs(ux));

                Real dt2;
#if AMREX_SPACEDIM >= 2
                dt2 = dx[1]/(cy + std::abs(uy));
                if (coord == 2) {
                    dt2 *= problo[0] + 0.5_rt * dx[0];
                }
#else
                dt2 = dt1;
#endif

                Real dt3;
#if AMREX_SPACEDIM == 3
                dt3 = dx[2]/(cz + std::abs(uz));
#else
                dt3 = dt1;
#endif

                return {ValLocPair<Real, IntVect>{amrex::min(dt1, dt2, dt3), idx}};

            });

            return r;

        } // mhd
#endif // MHD
#ifdef DIFFUSION
        else if constexpr (T == diffusion) {

            // Diffusion-limited timestep
            //
            // dt < 0.5 dx**2 / D
            // where D = k/(rho c_v), and k is the conductivity

            const Real ldiffuse_cutoff_density = castro::diffuse_cutoff_density;
            const Real lmax_dt = castro::max_dt;
            const Real lcfl = castro::cfl;

            auto r = amrex::ParReduce(TypeList<ReduceOpMin>{}, TypeList<ValLocPair<Real, IntVect>>{}, stateMF,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> GpuTuple<ValLocPair<Real, IntVect>>
            {
                Array4<Real const> const& ustate = ua[box_no];

                IntVect idx(AMREX_D_DECL(i,j,k));

                if (ustate(i,j,k,URHO) > ldiffuse_cutoff_density) {

                    Real rho_inv = 1.0_rt/ustate(i,j,k,URHO);

                    // we need cv
                    eos_t eos_state;

                    eos_state.rho = ustate(i,j,k,URHO);
                    eos_state.T = ustate(i,j,k,UTEMP);
                    eos_state.e = ustate(i,j,k,UEINT) * rho_inv;
                    for (int n = 0; n < NumSpec; n++) {
                        eos_state.xn[n]  = ustate(i,j,k,UFS+n) * rho_inv;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; n++) {
                        eos_state.aux[n] = ustate(i,j,k,UFX+n) * rho_inv;
                    }
#endif

                    eos(eos_input_re, eos_state);

                    // we also need the conductivity
                    conductivity(eos_state);

                    // maybe we should check (and take action) on negative cv here?
                    Real D = eos_state.conductivity * rho_inv / eos_state.cv;

                    Real dt1 = 0.5_rt * dx[0]*dx[0] / D;

                    Real dt2;
#if AMREX_SPACEDIM >= 2
                    dt2 = 0.5_rt * dx[1]*dx[1] / D;
                    if (coord == 2) {
                        Real rc = problo[0] + 0.5_rt * dx[0];
                        dt2 *= rc * rc;
                    }
#else
                    dt2 = dt1;
#endif

                    Real dt3;
#if AMREX_SPACEDIM >= 3
                    dt3 = 0.5_rt * dx[2]*dx[2] / D;
#else
                    dt3 = dt1;
#endif

                    return {ValLocPair<Real, IntVect>{amrex::min(dt1, dt2, dt3), idx}};

                } else {
                    return {ValLocPair<Real, IntVect>{lmax_dt/lcfl, idx}};
                }
            });

            return r;

        } // diffusion
#endif // DIFFUSION
#ifdef REACTIONS
        else if constexpr (T == burning) {

            if (castro::dtnuc_e > 1.e199_rt && castro::dtnuc_X > 1.e199_rt) {
                IntVect idx(AMREX_D_DECL(0,0,0));
                return {ValLocPair<Real, IntVect>{1.e200_rt, idx}};
            }

            auto r = amrex::ParReduce(TypeList<ReduceOpMin>{}, TypeList<ValLocPair<Real, IntVect>>{}, stateMF,
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) -> GpuTuple<ValLocPair<Real, IntVect>>
            {
                Real dt_tmp = 1.e200_rt;

                Array4<Real const> const& S = ua[box_no];
                Array4<Real const> const& mask = ma[box_no];

                IntVect idx(AMREX_D_DECL(i,j,k));

                // Don't consider zones that are masked out.

                if (mask_covered_zones && mask.contains(i,j,k)) {
                    if (mask(i,j,k) == 0.0_rt) {
                        return {ValLocPair<Real, IntVect>{dt_tmp, idx}};
                    }
                }

                // Set a floor on the minimum size of a derivative. This floor
                // is small enough such that it will result in no timestep limiting.

                const Real derivative_floor = 1.e-50_rt;

                // We want to limit the timestep so that it is not larger than
                // dtnuc_e * (e / (de/dt)).  If the timestep factor dtnuc is
                // equal to 1, this says that we don't want the
                // internal energy to change by any more than its current
                // magnitude in the next timestep.
                //
                // If dtnuc is less than one, it controls the fraction we will
                // allow the internal energy to change in this timestep due to
                //  nuclear burning, provided that our instantaneous estimate
                // of the energy release is representative of the full timestep.
                //
                // We also do the same thing for the species, using a timestep
                // limiter dtnuc_X * (X_k / (dX_k/dt)). To prevent changes
                // due to trace isotopes that we probably are not interested in,
                // only apply the limiter to species with an abundance greater
                // than a user-specified threshold.
                //
                // To estimate de/dt and dX/dt, we are going to call the RHS of the
                // burner given the current state data. We need to do an EOS
                // call before we do the RHS call so that we have accurate
                // values for the thermodynamic data like abar, zbar, etc.
                // But we will call in (rho, T) mode, which is inexpensive.

                Real rhoInv = 1.0_rt / S(i,j,k,URHO);

                burn_t burn_state;

#if AMREX_SPACEDIM == 1
                burn_state.dx = dx[0];
#else
                Real dx1 = dx[1];
#if AMREX_SPACEDIM >= 2
                if (coord == 2) {
                    dx1 *= problo[0] + 0.5_rt * dx[0];
                }
#endif
                burn_state.dx = amrex::min(AMREX_D_DECL(dx[0], dx1, dx[2]));
#endif

                burn_state.rho = S(i,j,k,URHO);
                burn_state.T   = S(i,j,k,UTEMP);
                burn_state.e   = S(i,j,k,UEINT) * rhoInv;
                for (int n = 0; n < NumSpec; ++n) {
                    burn_state.xn[n] = S(i,j,k,UFS+n) * rhoInv;
                }
#if NAUX_NET > 0
                for (int n = 0; n < NumAux; ++n) {
                    burn_state.aux[n] = S(i,j,k,UFX+n) * rhoInv;
                }
#endif

                // Don't consider zones that are outside our burning criteria.

                if (burn_state.T < castro::react_T_min || burn_state.T > castro::react_T_max ||
                    burn_state.rho < castro::react_rho_min || burn_state.rho > castro::react_rho_max) {
                    return {ValLocPair<Real, IntVect>{dt_tmp, idx}};
                }

                Real e = burn_state.e;
                Real X[NumSpec];
                for (int n = 0; n < NumSpec; ++n) {
                    X[n] = amrex::max(burn_state.xn[n], small_x);
                }

                eos(eos_input_rt, burn_state);

                Array1D<Real, 1, neqs> ydot;
                actual_rhs(burn_state, ydot);

                Real dedt = ydot(net_ienuc);
                Real dXdt[NumSpec];
                for (int n = 0; n < NumSpec; ++n) {
                    dXdt[n] = ydot(n+1) * aion[n];
                }

                // Apply a floor to the derivatives. This ensures that we don't
                // divide by zero; it also gives us a quick method to disable
                // the timestep limiting, because the floor is small enough
                // that the implied timestep will be very large, and thus
                // ignored compared to other limiters.

                dedt = amrex::max(std::abs(dedt), derivative_floor);

                for (int n = 0; n < NumSpec; ++n) {
                    if (X[n] >= castro::dtnuc_X_threshold) {
                        dXdt[n] = amrex::max(std::abs(dXdt[n]), derivative_floor);
                    } else {
                        dXdt[n] = derivative_floor;
                    }
                }

#ifdef NSE

#ifdef SIMPLIFIED_SDC
                // if we are doing simplified-SDC + NSE, then the `in_nse()`
                // check will use burn_state.y[], so we need to ensure that
                // those are initialized
                for (int n = 0; n < NumSpec; ++n) {
                    burn_state.y[SFS+n] = burn_state.rho * burn_state.xn[n];
                }

                burn_state.y[SEINT] = burn_state.rho * burn_state.e;

#endif // SIMPLIFIED_SDC

#ifdef NSE_NET
                burn_state.mu_p = S(i,j,k,UMUP);
                burn_state.mu_n = S(i,j,k,UMUN);
#endif

                if (!in_nse(burn_state)) {
#endif // NSE
                    dt_tmp = castro::dtnuc_e * e / dedt;
#ifdef NSE
                }
#endif // NSE
                for (int n = 0; n < NumSpec; ++n) {
                    dt_tmp = amrex::min(dt_tmp, castro::dtnuc_X * (X[n] / dXdt[n]));
                }

                return {ValLocPair<Real, IntVect>{dt_tmp, idx}};
            });

            return r;

        } // burning
#endif // REACTIONS

    }

#ifdef RADIATION
    ///
    /// Radiation hydro timestep
    ///
    amrex::Real estdt_rad (const amrex::MultiFab& stateMF, const amrex::MultiFab& radMF,
                           const amrex::GeometryData& geomdata);
#endif

} // namespace timestep

#endif // TIMESTEP_H
