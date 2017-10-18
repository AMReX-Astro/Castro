#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

using namespace amrex;

void
Castro::sum_integrated_quantities ()
{
#ifndef SGS
    if (verbose <= 0) return;
#endif

    bool local_flag = true;

    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
    Real prev_time   = state[State_Type].prevTime();
    Real mass        = 0.0;
    Real mom[3]      = { 0.0 };
    Real com[3]      = { 0.0 };
    Real com_vel[3]  = { 0.0 };
    Real rho_e       = 0.0;
    Real rho_K       = 0.0;
    Real rho_E       = 0.0;
    Real dt_crse     = parent->dtLevel(0);
#ifdef SGS
    Real Etot        = 0.0;
    Real delta_E     = 0.0;
    Real delta_K     = 0.0;
    Real prod_sgs    = 0.0;
    Real diss_sgs    = 0.0;
    Real turb_src    = 0.0;
    Real rms_mach    = 0.0;
#endif

    Real T_max = -1.0e200;
    Real T_min =  1.0e200;
    Real grad_T_max = -1.0e200;
    Real rho_fuel_dot = 0.0;
    Real rho_fuel_initial = 1.0;

    Real flame_width = 0.0;
    Real flame_speed = 0.0;
    
    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass   += ca_lev.volWgtSum("density", time, local_flag);
        mom[0] += ca_lev.volWgtSum("xmom", time, local_flag);
	mom[1] += ca_lev.volWgtSum("ymom", time, local_flag);
	mom[2] += ca_lev.volWgtSum("zmom", time, local_flag);

	if (show_center_of_mass) {
	   com[0] += ca_lev.locWgtSum("density", time, 0, local_flag);
	   com[1] += ca_lev.locWgtSum("density", time, 1, local_flag);
	   com[2] += ca_lev.locWgtSum("density", time, 2, local_flag);
	}

       rho_e += ca_lev.volWgtSum("rho_e", time, local_flag);
       rho_K += ca_lev.volWgtSum("kineng", time, local_flag);
       rho_E += ca_lev.volWgtSum("rho_E", time, local_flag);

#ifdef SGS
        Real  cur_time = state[SGS_Type].curTime();
        Real prev_time = state[SGS_Type].prevTime();

        delta_E  += ca_lev.volWgtSum("rho_E", cur_time, local_flag);
        delta_E  -= ca_lev.volWgtSum("rho_E", prev_time, local_flag);

        delta_K  += ca_lev.volWgtSum("kineng", cur_time, local_flag);
        delta_K  -= ca_lev.volWgtSum("kineng", prev_time, local_flag);

        rms_mach  += ca_lev.volWgtSquaredSum("MachNumber", time, local_flag);

        prod_sgs += 0.5 * ca_lev.volWgtSum("prod_sgs", prev_time, local_flag) * dt_crse;
        prod_sgs += 0.5 * ca_lev.volWgtSum("prod_sgs",  cur_time, local_flag) * dt_crse;
        diss_sgs += 0.5 * ca_lev.volWgtSum("diss_sgs", prev_time, local_flag) * dt_crse;
        diss_sgs += 0.5 * ca_lev.volWgtSum("diss_sgs",  cur_time, local_flag) * dt_crse;
        turb_src += 0.5 * ca_lev.volWgtSum("turb_src", prev_time, local_flag) * dt_crse;
        turb_src += 0.5 * ca_lev.volWgtSum("turb_src",  cur_time, local_flag) * dt_crse;

        sum_turb_src = sum_turb_src + turb_src;
#endif

	ca_lev.flame_width_properties(time, T_max, T_min, grad_T_max);

	// For now we only do the flame speed calculation on the coarse grid since
	// it doesn't share the same timestep as the fine grids.

	if (lev == 0)
	  ca_lev.flame_speed_properties(time, rho_fuel_dot);
    }
 
    if (verbose > 0)
    {
#ifdef SGS
	const int nfoo = 14;
	Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], rho_e, rho_K, rho_E, Etot, delta_E, delta_K, 
			  prod_sgs, diss_sgs, turb_src, rms_mach};
#else
	const int nfoo = 7;
	Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], rho_e, rho_K, rho_E};
#endif

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif

	ParallelDescriptor::ReduceRealSum(foo, nfoo, ParallelDescriptor::IOProcessorNumber());

	if (show_center_of_mass)
	    ParallelDescriptor::ReduceRealSum(com, 3, ParallelDescriptor::IOProcessorNumber());

	ParallelDescriptor::ReduceRealMax(T_max);
	ParallelDescriptor::ReduceRealMin(T_min);
	ParallelDescriptor::ReduceRealMax(grad_T_max);

	flame_width = (T_max - T_min) / grad_T_max;

	ParallelDescriptor::ReduceRealSum(rho_fuel_dot);

	// Note that rho_fuel_dot has already been multiplied by dx, so
	// the dimensionality here checks out.

	flame_speed = -rho_fuel_dot / rho_fuel_initial;

	if (ParallelDescriptor::IOProcessor()) {

	    int i = 0;
	    mass     = foo[i++];
	    mom[0]   = foo[i++];
            mom[1]   = foo[i++];
            mom[2]   = foo[i++];
	    rho_e    = foo[i++];
	    rho_K    = foo[i++];
            rho_E    = foo[i++];
#ifdef SGS
	    Etot     = foo[i++];
	    delta_E  = foo[i++];
	    delta_K  = foo[i++];
	    prod_sgs = foo[i++];
	    diss_sgs = foo[i++];
	    turb_sgs = foo[i++];
	    rms_mach = foo[i++];
#endif

	    std::cout << '\n';
	    std::cout << "TIME= " << time << " MASS        = "   << mass      << '\n';
	    std::cout << "TIME= " << time << " XMOM        = "   << mom[0]    << '\n';
	    std::cout << "TIME= " << time << " YMOM        = "   << mom[1]    << '\n';
	    std::cout << "TIME= " << time << " ZMOM        = "   << mom[2]    << '\n';
	    std::cout << "TIME= " << time << " RHO*e       = "   << rho_e     << '\n';
	    std::cout << "TIME= " << time << " RHO*K       = "   << rho_K     << '\n';
	    std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
#ifdef SGS
	    Etot     = rho_E + rho_K;
	    std::cout << "TIME= " << time << " TOTAL E     = "   << Etot      << '\n';
	    std::cout << "TIME= " << time << " DELTA E     = "   << delta_E   << '\n';
	    std::cout << "TIME= " << time << " DELTA K     = "   << delta_K   << '\n';
	    std::cout << "TIME= " << time << " DELTA TOT   = "   << delta_K+delta_E   << '\n';
	    std::cout << "TIME= " << time << " PROD_SGS    = "   << prod_sgs  << '\n';
	    std::cout << "TIME= " << time << " DISS_SGS    = "   << diss_sgs  << '\n';
	    std::cout << "TIME= " << time << " TURB_SRC    = "   << turb_src  << '\n';
	    std::cout << "TIME= " << time << " DE+DK-TURB_SRC = "   << delta_E+delta_K-turb_src  << '\n';
#endif
	    std::cout << "TIME= " << time << " FLAME WIDTH = "   << flame_width << '\n';
	    std::cout << "TIME= " << time << " FLAME SPEED = "   << flame_speed << '\n';
	    
	    if (parent->NumDataLogs() > 0 ) {

	       std::ostream& data_log1 = parent->DataLog(0);

	       if (data_log1.good()) {

		  if (time == 0.0) {
		      data_log1 << std::setw(14) <<  "      time    ";
		      data_log1 << std::setw(14) <<  "         mass ";
		      data_log1 << std::setw(14) <<  "         xmom ";
		      data_log1 << std::setw(14) <<  "         ymom ";
		      data_log1 << std::setw(14) <<  "         zmom ";
		      data_log1 << std::setw(14) <<  "        rho_K ";
		      data_log1 << std::setw(14) <<  "        rho_e ";
		      data_log1 << std::setw(14) <<  "        rho_E ";
		      data_log1 << std::setw(14) <<  "  flame width ";
		      data_log1 << std::setw(14) <<  "  flame speed " << std::endl;
		  }

		      // Write the quantities at this time
		  data_log1 << std::setw(14) <<  time;
		  data_log1 << std::setw(14) <<  std::setprecision(6) << mass;
		  data_log1 << std::setw(14) <<  std::setprecision(6) << mom[0];
		  data_log1 << std::setw(14) <<  std::setprecision(6) << mom[1];
		  data_log1 << std::setw(14) <<  std::setprecision(6) << mom[2];
		  data_log1 << std::setw(14) <<  std::setprecision(6) << rho_K;
		  data_log1 << std::setw(14) <<  std::setprecision(6) << rho_e;
		  data_log1 << std::setw(14) <<  std::setprecision(6) << rho_E;
		  data_log1 << std::setw(14) <<  std::setprecision(6) << flame_width;
		  data_log1 << std::setw(14) <<  std::setprecision(6) << flame_speed << std::endl;

	       }

	    }

#ifdef SGS
	    if (parent->NumDataLogs() > 1) { 

	       std::ostream& data_log2 = parent->DataLog(1);

	       if (data_log2.good()) {

		  // Write the quantities that represent changes from prev_time to cur_time
		  if (time == 0.0) {
		      data_log2 << std::setw(14) <<  "      time    ";
		      data_log2 << std::setw(14) <<  "         Etot ";
		      data_log2 << std::setw(16) <<  "  Etot-sum_turb  ";
		      data_log2 << std::setw(14) <<  "      rms_mach";
		      data_log2 << std::setw(14) <<  "      delta_E ";
		      data_log2 << std::setw(14) <<  "      delta_K ";
		      data_log2 << std::setw(14) <<  "      prod_sgs";
		      data_log2 << std::setw(14) <<  "      diss_sgs";
		      data_log2 << std::setw(14) <<  "      turb_src" << std::endl;
		  }

		  data_log2 << std::setw(14) <<  std::setprecision(6) << time;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << Etot;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << delta_E;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << delta_K;
		  data_log2 << std::setw(16) <<  std::setprecision(10) << Etot-sum_turb_src;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << rms_mach;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << prod_sgs;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << diss_sgs;
		  data_log2 << std::setw(14) <<  std::setprecision(6) << turb_src << std::endl;
	       }
	    }
#endif
	    
	    if (show_center_of_mass) {
	        for (int i = 0; i <= 2; i++) {
		  com[i]     = com[i] / mass;
		  com_vel[i] = mom[i] / mass;
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
}
