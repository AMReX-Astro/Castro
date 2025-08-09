#include <iomanip>

#include <Castro.H>

#ifdef GRAVITY
#include <Gravity.H>
#endif

#include <problem_diagnostics.H>

using namespace amrex;

void
Castro::sum_integrated_quantities ()
{
    if (verbose <= 0) {
        return;
    }

    BL_PROFILE("Castro::sum_integrated_quantities()");

    bool local_flag = true;

    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real dt          = parent->dtLevel(0);
    if (time == 0.0) {
        dt = 0.0; // dtLevel returns the next timestep for t = 0, so overwrite
    }
    int timestep     = parent->levelSteps(0);
    Real mass        = 0.0;
    Real mom[3]      = { 0.0 };
    Real ang_mom[3]  = { 0.0 };
#ifdef HYBRID_MOMENTUM
    Real hyb_mom[3]  = { 0.0 };
#endif
    Real com[3]      = { 0.0 };
    Real com_vel[3]  = { 0.0 };
    Real rho_e       = 0.0;
    Real rho_K       = 0.0;
    Real rho_E       = 0.0;
#ifdef GRAVITY
    Real rho_phi     = 0.0;
    Real total_energy = 0.0;
#endif

    Real T_max       = 0.0;
    Real rho_max     = 0.0;
    Real ts_te_max   = 0.0;

    int datprecision = 16;

    int datwidth     = 25; // Floating point data in scientific notation
    int fixwidth     = 25; // Floating point data not in scientific notation
    int intwidth     = 12; // Integer data

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);
        const MultiFab& S_new = ca_lev.get_new_data(State_Type);
#ifdef GRAVITY
        MultiFab& phi_new = ca_lev.get_new_data(PhiGrav_Type);
#endif
#ifdef REACTIONS
        MultiFab& R_new = ca_lev.get_new_data(Reactions_Type);
#endif

        mass   += ca_lev.volWgtSum(S_new, URHO, local_flag);
        mom[0] += ca_lev.volWgtSum(S_new, UMX, local_flag);
        mom[1] += ca_lev.volWgtSum(S_new, UMY, local_flag);
        mom[2] += ca_lev.volWgtSum(S_new, UMZ, local_flag);

        ang_mom[0] += ca_lev.volWgtSum("angular_momentum_x", time, local_flag);
        ang_mom[1] += ca_lev.volWgtSum("angular_momentum_y", time, local_flag);
        ang_mom[2] += ca_lev.volWgtSum("angular_momentum_z", time, local_flag);

#ifdef HYBRID_MOMENTUM
        hyb_mom[0] += ca_lev.volWgtSum(S_new, UMR, local_flag);
        hyb_mom[1] += ca_lev.volWgtSum(S_new, UML, local_flag);
        hyb_mom[2] += ca_lev.volWgtSum(S_new, UMP, local_flag);
#endif

        com[0] += ca_lev.locWgtSum(S_new, URHO, 0, local_flag);
        com[1] += ca_lev.locWgtSum(S_new, URHO, 1, local_flag);
        com[2] += ca_lev.locWgtSum(S_new, URHO, 2, local_flag);

        rho_e += ca_lev.volWgtSum(S_new, UEINT, local_flag);
        rho_K += ca_lev.volWgtSum("kineng", time, local_flag);
        rho_E += ca_lev.volWgtSum(S_new, UEDEN, local_flag);
#ifdef GRAVITY
        if (gravity->get_gravity_type() == "PoissonGrav") {
            rho_phi += ca_lev.volProductSum(S_new, phi_new, URHO, 0, local_flag);
        }
#endif

        // Compute extrema

#ifdef REACTIONS
        auto dx = ca_lev.geom.CellSizeArray();

        Real dd = 0.0_rt;
#if AMREX_SPACEDIM == 1
        dd = dx[0];
#elif AMREX_SPACEDIM == 2
        dd = amrex::min(dx[0], dx[1]);
#else
        dd = amrex::min(dx[0], dx[1], dx[2]);
#endif
#endif

        bool mask_available = true;
        if (lev == parent->finestLevel()) {
            mask_available = false;
        }

        MultiFab tmp_mf;
        const MultiFab& mask_mf = mask_available ? getLevel(lev+1).build_fine_mask() : tmp_mf;

        ReduceOps<ReduceOpMax, ReduceOpMax, ReduceOpMax> reduce_op;
        ReduceData<Real, Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();

            auto U = S_new[mfi].array();
#ifdef REACTIONS
            auto R = R_new[mfi].array();
#endif

            const auto & level_mask = mask_available ? mask_mf[mfi].array() : Array4<Real>{};

            reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                Real maskFactor = 1.0;
                if (mask_available) {
                    maskFactor = level_mask(i,j,k);
                }

                Real T = U(i,j,k,UTEMP) * maskFactor;
                Real rho = U(i,j,k,URHO) * maskFactor;
                Real ts_te = 0.0_rt;

#ifdef REACTIONS
                Real enuc = std::abs(R(i,j,k,0)) / U(i,j,k,URHO);

                if (enuc > 1.e-100_rt && maskFactor == 1.0) {

                    Real rhoInv = 1.0_rt / rho;

                    // Calculate sound speed
                    eos_rep_t eos_state;
                    eos_state.rho = rho;
                    eos_state.T   = T;
                    eos_state.e   = U(i,j,k,UEINT) * rhoInv;
                    for (int n = 0; n < NumSpec; ++n) {
                        eos_state.xn[n] = U(i,j,k,UFS+n) * rhoInv;
                    }
#if NAUX_NET > 0
                    for (int n = 0; n < NumAux; ++n) {
                        eos_state.aux[n] = U(i,j,k,UFX+n) * rhoInv;
                    }
#endif

                    eos(eos_input_re, eos_state);

                    Real t_e = eos_state.e / enuc;
                    Real t_s = dd / eos_state.cs;

                    ts_te = t_s / t_e;
                }
#endif

                return {T, rho, ts_te};
            });

        }

        ReduceTuple hv = reduce_data.value();

        T_max = amrex::max(T_max, amrex::get<0>(hv));
        rho_max = amrex::max(rho_max, amrex::get<1>(hv));
#ifdef REACTIONS
        ts_te_max = amrex::max(ts_te_max, amrex::get<2>(hv));
#endif

    }

    if (verbose > 0)
    {

#ifdef HYBRID_MOMENTUM
#ifdef GRAVITY
       const int nfoo = 17;
#else
       const int nfoo = 16;
#endif
#else
#ifdef GRAVITY
       const int nfoo = 14;
#else
       const int nfoo = 13;
#endif
#endif

        Real foo[nfoo] = {mass, mom[0], mom[1], mom[2],
                          com[0], com[1], com[2],
                          ang_mom[0], ang_mom[1], ang_mom[2],
#ifdef HYBRID_MOMENTUM
                          hyb_mom[0], hyb_mom[1], hyb_mom[2],
#endif
#ifdef GRAVITY
                          rho_e, rho_K, rho_E, rho_phi};
#else
                          rho_e, rho_K, rho_E};
#endif

        const int nfoo_max = 3;

        Real foo_max[nfoo_max] = {T_max, rho_max, ts_te_max};

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif

        ParallelDescriptor::ReduceRealSum(foo, nfoo, ParallelDescriptor::IOProcessorNumber());

        ParallelDescriptor::ReduceRealMax(foo_max, nfoo_max, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor()) {

            int i = 0;
            mass       = foo[i++];
            mom[0]     = foo[i++];
            mom[1]     = foo[i++];
            mom[2]     = foo[i++];
            com[0]     = foo[i++];
            com[1]     = foo[i++];
            com[2]     = foo[i++];
            ang_mom[0] = foo[i++];
            ang_mom[1] = foo[i++];
            ang_mom[2] = foo[i++];
#ifdef HYBRID_MOMENTUM
            hyb_mom[0] = foo[i++];
            hyb_mom[1] = foo[i++];
            hyb_mom[2] = foo[i++];
#endif
            rho_e      = foo[i++];
            rho_K      = foo[i++];
            rho_E      = foo[i++];
#ifdef GRAVITY
            rho_phi    = foo[i++];

            // Total energy is 1/2 * rho * phi + rho * E for self-gravity,
            // and rho * phi + rho * E for externally-supplied gravity.
            std::string gravity_type = gravity->get_gravity_type();
            if (gravity_type == "PoissonGrav" || gravity_type == "MonopoleGrav") {
                total_energy = 0.5 * rho_phi + rho_E;
            }
            else {
                total_energy = rho_phi + rho_E;
            }
#endif

            for (int idir = 0; idir < 3; idir++) {
                com[idir]     = com[idir] / mass;
                com_vel[idir] = mom[idir] / mass;
            }

            i = 0;
            T_max     = foo_max[i++];
            rho_max   = foo_max[i++];
            ts_te_max = foo_max[i++];    // NOLINT(clang-analyzer-deadcode.DeadStores)

            std::cout << '\n';
            std::cout << "TIME= " << time << " MASS        = "   << mass      << '\n';
            std::cout << "TIME= " << time << " XMOM        = "   << mom[0]    << '\n';
            std::cout << "TIME= " << time << " YMOM        = "   << mom[1]    << '\n';
            std::cout << "TIME= " << time << " ZMOM        = "   << mom[2]    << '\n';
            std::cout << "TIME= " << time << " ANG MOM X   = "   << ang_mom[0] << '\n';
            std::cout << "TIME= " << time << " ANG MOM Y   = "   << ang_mom[1] << '\n';
            std::cout << "TIME= " << time << " ANG MOM Z   = "   << ang_mom[2] << '\n';
#ifdef HYBRID_MOMENTUM
            std::cout << "TIME= " << time << " HYB MOM R   = "   << hyb_mom[0] << '\n';
            std::cout << "TIME= " << time << " HYB MOM L   = "   << hyb_mom[1] << '\n';
            std::cout << "TIME= " << time << " HYB MOM P   = "   << hyb_mom[2] << '\n';
#endif
            std::cout << "TIME= " << time << " RHO*e       = "   << rho_e     << '\n';
            std::cout << "TIME= " << time << " RHO*K       = "   << rho_K     << '\n';
            std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
#ifdef GRAVITY
            std::cout << "TIME= " << time << " RHO*PHI     = "   << rho_phi   << '\n';
            std::cout << "TIME= " << time << " TOTAL ENERGY= "   << total_energy << '\n';
#endif
            std::cout << "TIME= " << time << " CENTER OF MASS X-LOC = " << com[0]     << '\n';
            std::cout << "TIME= " << time << " CENTER OF MASS X-VEL = " << com_vel[0] << '\n';

            std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC = " << com[1]     << '\n';
            std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL = " << com_vel[1] << '\n';

            std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC = " << com[2]     << '\n';
            std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL = " << com_vel[2] << '\n';

            std::cout << "TIME= " << time << " MAXIMUM TEMPERATURE  = " << T_max << '\n';
            std::cout << "TIME= " << time << " MAXIMUM DENSITY      = " << rho_max << '\n';
#ifdef REACTIONS
            std::cout << "TIME= " << time << " MAXIMUM T_S / T_E    = " << ts_te_max << '\n';
#endif

            std::ostream& data_log1 = *Castro::data_logs[0];

            if (data_log1.good()) {

               if (time == 0.0) {

                   std::ostringstream header;

                   int n = 0;

                   header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                   header << std::setw(fixwidth) << "                     TIME"; ++n;
                   header << std::setw(datwidth) << "                     MASS"; ++n;
                   header << std::setw(datwidth) << "                     XMOM"; ++n;
                   header << std::setw(datwidth) << "                     YMOM"; ++n;
                   header << std::setw(datwidth) << "                     ZMOM"; ++n;
                   header << std::setw(datwidth) << "              ANG. MOM. X"; ++n;
                   header << std::setw(datwidth) << "              ANG. MOM. Y"; ++n;
                   header << std::setw(datwidth) << "              ANG. MOM. Z"; ++n;
#if (AMREX_SPACEDIM == 3)
#ifdef HYBRID_MOMENTUM
                   header << std::setw(datwidth) << "              HYB. MOM. R"; ++n;
                   header << std::setw(datwidth) << "              HYB. MOM. L"; ++n;
                   header << std::setw(datwidth) << "              HYB. MOM. P"; ++n;
#endif
#endif
                   header << std::setw(datwidth) << "              KIN. ENERGY"; ++n;
                   header << std::setw(datwidth) << "              INT. ENERGY"; ++n;
                   header << std::setw(datwidth) << "               GAS ENERGY"; ++n;
#ifdef GRAVITY
                   header << std::setw(datwidth) << "             GRAV. ENERGY"; ++n;
                   header << std::setw(datwidth) << "             TOTAL ENERGY"; ++n;
#endif
                   header << std::setw(datwidth) << "     CENTER OF MASS X-LOC"; ++n;
                   header << std::setw(datwidth) << "     CENTER OF MASS Y-LOC"; ++n;
                   header << std::setw(datwidth) << "     CENTER OF MASS Z-LOC"; ++n;
                   header << std::setw(datwidth) << "     CENTER OF MASS X-VEL"; ++n;
                   header << std::setw(datwidth) << "     CENTER OF MASS Y-VEL"; ++n;
                   header << std::setw(datwidth) << "     CENTER OF MASS Z-VEL"; ++n;
                   header << std::setw(datwidth) << "      MAXIMUM TEMPERATURE"; ++n;
                   header << std::setw(datwidth) << "          MAXIMUM DENSITY"; ++n;
#ifdef REACTIONS
                   header << std::setw(datwidth) << "        MAXIMUM T_S / T_E"; ++n;
#endif

                   header << std::endl;

                   data_log1 << std::setw(intwidth) << "#   COLUMN 1";
                   data_log1 << std::setw(fixwidth) << 2;

                   for (int icol = 3; icol <= n; ++icol) {
                       data_log1 << std::setw(datwidth) << icol;
                   }

                   data_log1 << std::endl;

                   data_log1 << header.str();
               }

               // Write the quantities at this time
               data_log1 << std::fixed;

               data_log1 << std::setw(intwidth) <<  timestep;

               if (time < 1.e-4_rt || time > 1.e4_rt) {
                   data_log1 << std::scientific;
               } else {
                   data_log1 << std::fixed;
               }

               data_log1 << std::setw(fixwidth) <<  std::setprecision(datprecision) << time;

               data_log1 << std::scientific;

               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mass;
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[0];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[1];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << mom[2];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << ang_mom[0];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << ang_mom[1];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << ang_mom[2];
#ifdef HYBRID_MOMENTUM
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << hyb_mom[0];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << hyb_mom[1];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << hyb_mom[2];
#endif
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_K;
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_e;
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_E;
#ifdef GRAVITY
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_phi;
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << total_energy;
#endif
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << com[0];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << com[1];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << com[2];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << com_vel[0];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << com_vel[1];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << com_vel[2];
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << T_max;
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << rho_max;
#ifdef REACTIONS
               data_log1 << std::setw(datwidth) <<  std::setprecision(datprecision) << ts_te_max;
#endif

               data_log1 << std::endl;

            }
        }
#ifdef BL_LAZY
        });
#endif
    }

#ifdef GRAVITY
    // Gravity diagnostics
    {
        // Gravitational wave amplitudes

        Real h_plus_1  = 0.0;
        Real h_cross_1 = 0.0;

        Real h_plus_2  = 0.0;
        Real h_cross_2 = 0.0;

        Real h_plus_3  = 0.0;
        Real h_cross_3 = 0.0;

        for (int lev = 0; lev <= finest_level; lev++)
        {
            Castro& ca_lev = getLevel(lev);

#if (AMREX_SPACEDIM > 1)
            // Gravitational wave signal. This is designed to add to these quantities so we can send them directly.
            ca_lev.gwstrain(time, h_plus_1, h_cross_1, h_plus_2, h_cross_2, h_plus_3, h_cross_3, local_flag);
#endif

        }

        const int nfoo_sum = 6;

        amrex::Vector<Real> foo_sum(nfoo_sum);

        foo_sum[0] = h_plus_1;
        foo_sum[1] = h_cross_1;
        foo_sum[2] = h_plus_2;
        foo_sum[3] = h_cross_2;
        foo_sum[4] = h_plus_3;
        foo_sum[5] = h_cross_3;

        amrex::ParallelDescriptor::ReduceRealSum(foo_sum.dataPtr(), nfoo_sum);

        h_plus_1   = foo_sum[0];
        h_cross_1  = foo_sum[1];
        h_plus_2   = foo_sum[2];
        h_cross_2  = foo_sum[3];
        h_plus_3   = foo_sum[4];
        h_cross_3  = foo_sum[5];

        if (ParallelDescriptor::IOProcessor()) {

            std::ostream& log = *Castro::data_logs[1];

            // Write header row

            if (time == 0.0) {

                int n = 0;

                std::ostringstream header;

                header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                header << std::setw(fixwidth) << "                     TIME"; ++n;

                header << std::setw(datwidth) << "                  h_+ (x)"; ++n;
                header << std::setw(datwidth) << "                  h_x (x)"; ++n;
                header << std::setw(datwidth) << "                  h_+ (y)"; ++n;
                header << std::setw(datwidth) << "                  h_x (y)"; ++n;
                header << std::setw(datwidth) << "                  h_+ (z)"; ++n;
                header << std::setw(datwidth) << "                  h_x (z)"; ++n;

                header << std::endl;

                log << std::setw(intwidth) << "#   COLUMN 1";
                log << std::setw(fixwidth) << 2;

                for (int i = 3; i <= n; ++i) {
                    log << std::setw(datwidth) << i;
                }

                log << std::endl;

                log << header.str();

            }

            log << std::fixed;

            log << std::setw(intwidth)                                    << timestep;

            if (time < 1.e-4_rt || time > 1.e4_rt) {
                log << std::scientific;
            }
            else {
                log << std::fixed;
            }

            log << std::setw(fixwidth) << std::setprecision(datprecision) << time;

            log << std::scientific;

            log << std::setw(datwidth) << std::setprecision(datprecision) << h_plus_1;
            log << std::setw(datwidth) << std::setprecision(datprecision) << h_cross_1;
            log << std::setw(datwidth) << std::setprecision(datprecision) << h_plus_2;
            log << std::setw(datwidth) << std::setprecision(datprecision) << h_cross_2;
            log << std::setw(datwidth) << std::setprecision(datprecision) << h_plus_3;
            log << std::setw(datwidth) << std::setprecision(datprecision) << h_cross_3;

            log << std::endl;

        }

    }
#endif

    // Species

    {
        std::vector<Real> species_mass(NumSpec);
        std::vector<std::string> species_names(NumSpec);

        // Species names

        for (int i = 0; i < NumSpec; i++) {
            species_names[i] = desc_lst[State_Type].name(UFS+i);
            species_names[i] = species_names[i].substr(4,std::string::npos);
            species_mass[i] = 0.0;
        }

        // Integrated mass of all species on the domain

        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab& S_new = getLevel(lev).get_new_data(State_Type);
            for (int i = 0; i < NumSpec; ++i) {
                species_mass[i] += getLevel(lev).volWgtSum(S_new, UFS + i, local_flag) / C::M_solar;
            }
        }

        int nfoo_sum = NumSpec;

        amrex::Vector<Real> foo_sum(nfoo_sum);

        for (int i = 0; i < NumSpec; ++i) {
            foo_sum[i] = species_mass[i];
        }

        amrex::ParallelDescriptor::ReduceRealSum(foo_sum.dataPtr(), nfoo_sum);

        for (int i = 0; i < NumSpec; ++i) {
            species_mass[i] = foo_sum[i];
        }

        if (ParallelDescriptor::IOProcessor()) {

            std::ostream& log = *Castro::data_logs[2];

            if (time == 0.0) {

                int n = 0;

                std::ostringstream header;

                header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                header << std::setw(fixwidth) << "                     TIME"; ++n;

                for (int i = 0; i < NumSpec; i++) {
                    header << std::setw(datwidth) << ("Mass " + species_names[i]); ++n;
                }

                header << std::endl;

                log << std::setw(intwidth) << "#   COLUMN 1";
                log << std::setw(fixwidth) << 2;

                for (int i = 3; i <= n; ++i) {
                    log << std::setw(datwidth) << i;
                }

                log << std::endl;

                log << header.str();

            }

            log << std::fixed;

            log << std::setw(intwidth)                                    << timestep;

            if (time < 1.e-4_rt || time > 1.e4_rt) {
                log << std::scientific;
            } else {
                log << std::fixed;
            }

            log << std::setw(fixwidth) << std::setprecision(datprecision) << time;

            log << std::scientific;

            for (int i = 0; i < NumSpec; i++) {
                log << std::setw(datwidth) << std::setprecision(datprecision) << species_mass[i];
            }

            log << std::endl;

        }

    }

    // Information about the AMR driver.

    {
        // Calculate wall time for the step.

        Real wall_time = 0.0;

        if (time > 0.0) {
            wall_time = amrex::ParallelDescriptor::second() - wall_time_start;
        }

        // Calculate GPU memory consumption.

#ifdef AMREX_USE_GPU
        Long gpu_size_free_MB = Gpu::Device::freeMemAvailable() / (1024 * 1024);
        ParallelDescriptor::ReduceLongMin(gpu_size_free_MB, ParallelDescriptor::IOProcessorNumber());

        Long gpu_size_used_MB = (Gpu::Device::totalGlobalMem() - Gpu::Device::freeMemAvailable()) / (1024 * 1024);
        ParallelDescriptor::ReduceLongMax(gpu_size_used_MB, ParallelDescriptor::IOProcessorNumber());
#endif

        // Calculate maximum number of advance subcycles across all levels.

        int max_num_subcycles = 0;
        if (time > 0.0_rt) {
            for (int lev = 0; lev <= parent->finestLevel(); ++lev) {
                max_num_subcycles = std::max(max_num_subcycles, getLevel(lev).num_subcycles_taken);
            }
        }

        if (ParallelDescriptor::IOProcessor()) {

            std::ostream& log = *Castro::data_logs[3];

            if (time == 0.0) {

                int n = 0;

                std::ostringstream header;

                header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                header << std::setw(fixwidth) << "                     TIME"; ++n;
                header << std::setw(fixwidth) << "                       DT"; ++n;
                header << std::setw(intwidth) << "  FINEST LEV";              ++n;
                header << std::setw(fixwidth) << "  MAX NUMBER OF SUBCYCLES"; ++n;
                header << std::setw(datwidth) << " COARSE TIMESTEP WALLTIME"; ++n;
#ifdef AMREX_USE_GPU
                header << std::setw(datwidth) << "  MAXIMUM GPU MEMORY USED"; ++n;
                header << std::setw(datwidth) << "  MINIMUM GPU MEMORY FREE"; ++n;
#endif

                header << std::endl;

                log << std::setw(intwidth) << "#   COLUMN 1";
                log << std::setw(fixwidth) << 2;

                for (int i = 3; i < 4; ++i) {
                    log << std::setw(fixwidth) << i;
                }

                log << std::setw(intwidth) << 4; // Handle the finest lev column
                log << std::setw(fixwidth) << 5; // Handle the subcycle count column

                for (int i = 6; i <= n; ++i) {
                    log << std::setw(datwidth) << i;
                }

                log << std::endl;

                log << header.str();

            }

            log << std::fixed;

            log << std::setw(intwidth)                                    << timestep;

            if (time < 1.e-4_rt || time > 1.e4_rt) {
                log << std::scientific;
            } else {
                log << std::fixed;
            }

            log << std::setw(fixwidth) << std::setprecision(datprecision) << time;

            log << std::fixed;

            log << std::setw(fixwidth) << std::setprecision(datprecision) << dt;
            log << std::setw(intwidth)                                    << parent->finestLevel();
            log << std::setw(fixwidth)                                    << max_num_subcycles;
            log << std::setw(datwidth) << std::setprecision(datprecision) << wall_time;
#ifdef AMREX_USE_GPU
            log << std::setw(datwidth)                                    << gpu_size_used_MB;
            log << std::setw(datwidth)                                    << gpu_size_free_MB;
#endif

            log << std::endl;

        }

    }

    problem_diagnostics();
}
