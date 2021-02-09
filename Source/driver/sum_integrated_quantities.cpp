#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

#ifdef GRAVITY
#include <Gravity.H>
#endif

#include <problem_diagnostics.H>

using namespace amrex;

void
Castro::sum_integrated_quantities ()
{

    if (verbose <= 0) return;

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

    int datprecision = 16;

    int datwidth     = 25; // Floating point data in scientific notation
    int fixwidth     = 25; // Floating point data not in scientific notation
    int intwidth     = 12; // Integer data

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass   += ca_lev.volWgtSum("density", time, local_flag);
        mom[0] += ca_lev.volWgtSum("xmom", time, local_flag);
        mom[1] += ca_lev.volWgtSum("ymom", time, local_flag);
        mom[2] += ca_lev.volWgtSum("zmom", time, local_flag);

        ang_mom[0] += ca_lev.volWgtSum("angular_momentum_x", time, local_flag);
        ang_mom[1] += ca_lev.volWgtSum("angular_momentum_y", time, local_flag);
        ang_mom[2] += ca_lev.volWgtSum("angular_momentum_z", time, local_flag);

#ifdef HYBRID_MOMENTUM
        hyb_mom[0] += ca_lev.volWgtSum("rmom", time, local_flag);
        hyb_mom[1] += ca_lev.volWgtSum("lmom", time, local_flag);
        hyb_mom[2] += ca_lev.volWgtSum("zmom", time, local_flag);
#endif

        if (show_center_of_mass) {
           com[0] += ca_lev.locWgtSum("density", time, 0, local_flag);
           com[1] += ca_lev.locWgtSum("density", time, 1, local_flag);
           com[2] += ca_lev.locWgtSum("density", time, 2, local_flag);
        }

       rho_e += ca_lev.volWgtSum("rho_e", time, local_flag);
       rho_K += ca_lev.volWgtSum("kineng", time, local_flag);
       rho_E += ca_lev.volWgtSum("rho_E", time, local_flag);
#ifdef GRAVITY
        if (gravity->get_gravity_type() == "PoissonGrav")
               rho_phi += ca_lev.volProductSum("density", "phiGrav", time, local_flag);
#endif

    }

    if (verbose > 0)
    {

#ifdef HYBRID_MOMENTUM
#ifdef GRAVITY
       const int nfoo = 14;
#else
       const int nfoo = 13;
#endif
#else
#ifdef GRAVITY
       const int nfoo = 11;
#else
       const int nfoo = 10;
#endif
#endif

        Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], ang_mom[0], ang_mom[1], ang_mom[2],
#ifdef HYBRID_MOMENTUM
                          hyb_mom[0], hyb_mom[1], hyb_mom[2],
#endif
#ifdef GRAVITY
                          rho_e, rho_K, rho_E, rho_phi};
#else
                          rho_e, rho_K, rho_E};
#endif

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif

        ParallelDescriptor::ReduceRealSum(foo, nfoo, ParallelDescriptor::IOProcessorNumber());

        if (show_center_of_mass)
            ParallelDescriptor::ReduceRealSum(com, 3, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor()) {

            int i = 0;
            mass       = foo[i++];
            mom[0]     = foo[i++];
            mom[1]     = foo[i++];
            mom[2]     = foo[i++];
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

            // Total energy is -1/2 * rho * phi + rho * E for self-gravity,
            // and -rho * phi + rho * E for externally-supplied gravity.
            std::string gravity_type = gravity->get_gravity_type();
            if (gravity_type == "PoissonGrav" || gravity_type == "MonopoleGrav")
              total_energy = -0.5 * rho_phi + rho_E;
            else
              total_energy = -rho_phi + rho_E;
#endif

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

            std::ostream& data_log1 = *Castro::data_logs[0];

            if (data_log1.good()) {

               if (time == 0.0) {
                   data_log1 << std::setw(datwidth) <<  "          time";
                   data_log1 << std::setw(datwidth) <<  "          mass";
                   data_log1 << std::setw(datwidth) <<  "          xmom";
                   data_log1 << std::setw(datwidth) <<  "          ymom";
                   data_log1 << std::setw(datwidth) <<  "          zmom";
                   data_log1 << std::setw(datwidth) <<  "     ang mom x";
                   data_log1 << std::setw(datwidth) <<  "     ang mom y";
                   data_log1 << std::setw(datwidth) <<  "     ang mom z";
#ifdef HYBRID_MOMENTUM
                   data_log1 << std::setw(datwidth) <<  "     hyb mom r";
                   data_log1 << std::setw(datwidth) <<  "     hyb mom l";
                   data_log1 << std::setw(datwidth) <<  "     hyb mom p";
#endif
                   data_log1 << std::setw(datwidth) <<  "         rho_K";
                   data_log1 << std::setw(datwidth) <<  "         rho_e";
                   data_log1 << std::setw(datwidth) <<  "         rho_E";
#ifdef GRAVITY
                   data_log1 << std::setw(datwidth) <<  "       rho_phi";
                   data_log1 << std::setw(datwidth) <<  "  total energy";
#endif
                   data_log1 << std::endl;
               }

               // Write the quantities at this time
               data_log1 << std::setw(datwidth) <<  time;
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
               data_log1 << std::endl;

            }

            if (show_center_of_mass) {
                for (int idir = 0; idir <= 2; idir++) {
                  com[idir]     = com[idir] / mass;
                  com_vel[idir] = mom[idir] / mass;
                }

                std::cout << "TIME= " << time << " CENTER OF MASS X-LOC = " << com[0]     << '\n';
                std::cout << "TIME= " << time << " CENTER OF MASS X-VEL = " << com_vel[0] << '\n';

                std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC = " << com[1]     << '\n';
                std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL = " << com_vel[1] << '\n';

                std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC = " << com[2]     << '\n';
                std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL = " << com_vel[2] << '\n';
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

                log << std::setw(intwidth) << "#   COLUMN 1";
                log << std::setw(fixwidth) << "                         2";
                log << std::setw(fixwidth) << "                         3";
                log << std::setw(fixwidth) << "                         4";
                log << std::setw(fixwidth) << "                         5";
                log << std::setw(fixwidth) << "                         6";
                log << std::setw(fixwidth) << "                         7";

                std::ostringstream header;

                header << std::setw(intwidth) << "#   TIMESTEP";
                header << std::setw(fixwidth) << "                     TIME";

                header << std::setw(datwidth) << "             h_+ (x)";
                header << std::setw(datwidth) << "             h_x (x)";
                header << std::setw(datwidth) << "             h_+ (y)";
                header << std::setw(datwidth) << "             h_x (y)";
                header << std::setw(datwidth) << "             h_+ (z)";
                header << std::setw(datwidth) << "             h_x (z)";

                header << std::endl;

                log << std::endl;

                log << header.str();

            }

            log << std::fixed;

            log << std::setw(intwidth)                                    << timestep;
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
            for (int i = 0; i < NumSpec; ++i) {
                species_mass[i] += getLevel(lev).volWgtSum("rho_" + species_names[i], time, local_flag) / C::M_solar;
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

                header << std::setw(intwidth) << "#   TIMESTEP";           ++n;
                header << std::setw(fixwidth) << "                  TIME"; ++n;

                // We need to be careful here since the species names have differing numbers of characters

                for (int i = 0; i < NumSpec; i++) {
                    std::string outString  = "";
                    std::string massString = "Mass ";
                    std::string specString = species_names[i];
                    while (outString.length() + specString.length() + massString.length() < datwidth) outString += " ";
                    outString += massString;
                    outString += specString;
                    header << std::setw(datwidth) << outString; ++n;
                }

                header << std::endl;

                log << std::setw(intwidth) << "#   COLUMN 1";
                log << std::setw(fixwidth) << "                        2";

                for (int i = 3; i <= n; ++i)
                    log << std::setw(datwidth) << i;

                log << std::endl;

                log << header.str();

            }

            log << std::fixed;

            log << std::setw(intwidth)                                    << timestep;
            log << std::setw(fixwidth) << std::setprecision(datprecision) << time;

            log << std::scientific;

            for (int i = 0; i < NumSpec; i++)
                log << std::setw(datwidth) << std::setprecision(datprecision) << species_mass[i];

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

        if (ParallelDescriptor::IOProcessor()) {

            std::ostream& log = *Castro::data_logs[3];

            if (time == 0.0) {

                int n = 0;

                std::ostringstream header;

                header << std::setw(intwidth) << "#   TIMESTEP";              ++n;
                header << std::setw(fixwidth) << "                     TIME"; ++n;
                header << std::setw(fixwidth) << "                       DT"; ++n;
                header << std::setw(intwidth) << "  FINEST LEV";              ++n;
                header << std::setw(fixwidth) << " COARSE TIMESTEP WALLTIME"; ++n;
#ifdef AMREX_USE_GPU
                header << std::setw(fixwidth) << "  MAXIMUM GPU MEMORY USED"; ++n;
                header << std::setw(fixwidth) << "  MINIMUM GPU MEMORY FREE"; ++n;
#endif

                header << std::endl;

                log << std::setw(intwidth) << "#   COLUMN 1";
                log << std::setw(fixwidth) << "                        2";

                for (int i = 3; i < 4; ++i) {
                    log << std::setw(datwidth) << i;
                }

                log << std::setw(intwidth) << 4; // Handle the finest lev column

                for (int i = 5; i <= n; ++i) {
                    log << std::setw(datwidth) << i;
                }

                log << std::endl;

                log << header.str();

            }

            log << std::fixed;

            log << std::setw(intwidth)                                    << timestep;
            log << std::setw(fixwidth) << std::setprecision(datprecision) << time;
            log << std::setw(fixwidth) << std::setprecision(datprecision) << dt;
            log << std::setw(intwidth)                                    << parent->finestLevel();
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
