#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

#ifdef GRAVITY
#include <Gravity.H>
#endif

void
Castro::sum_integrated_quantities ()
{
#ifndef SGS
    if (verbose <= 0) return;
#endif

    bool local_flag = true;

    int finest_level = parent->finestLevel();
    Real time        = state[State_Type].curTime();
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
    Real rho_phi     = 0.0;
#ifdef GRAVITY
    Real total_energy = 0.0;
#endif
#ifdef SGS
    Real dt_crse     = parent->dtLevel(0);
    Real Etot        = 0.0;
    Real delta_E     = 0.0;
    Real delta_K     = 0.0;
    Real prod_sgs    = 0.0;
    Real diss_sgs    = 0.0;
    Real turb_src    = 0.0;
    Real rms_mach    = 0.0;
#endif

    int datwidth     = 14;
    int datprecision = 6;

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
       rho_phi += ca_lev.volProductSum("density", "phiGrav", time, local_flag);
#endif

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
    }
 
    if (verbose > 0)
    {
#ifdef SGS
#ifdef HYBRID_MOMENTUM
	const int nfoo = 21;
#else
	const int nfoo = 18;
#endif
	Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], ang_mom[0], ang_mom[1], ang_mom[2],
#ifdef HYBRID_MOMENTUM
			  hyb_mom[0], hyb_mom[1], hyb_mom[2],
#endif
			  rho_e, rho_K, rho_E, rho_phi, Etot, delta_E, delta_K, 
			  prod_sgs, diss_sgs, turb_src, rms_mach};
#else
#ifdef HYBRID_MOMENTUM
	const int nfoo = 14;
#else
	const int nfoo = 11;
#endif
	Real foo[nfoo] = {mass, mom[0], mom[1], mom[2], ang_mom[0], ang_mom[1], ang_mom[2],
#ifdef HYBRID_MOMENTUM
			  hyb_mom[0], hyb_mom[1], hyb_mom[2],
#endif
			  rho_e, rho_K, rho_E, rho_phi};
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
	    rho_phi    = foo[i++];
#ifdef SGS
	    Etot       = foo[i++];
	    delta_E    = foo[i++];
	    delta_K    = foo[i++];
	    prod_sgs   = foo[i++];
	    diss_sgs   = foo[i++];
	    turb_sgs   = foo[i++];
	    rms_mach   = foo[i++];
#endif

#ifdef GRAVITY
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
	    if (parent->NumDataLogs() > 0 ) {

	       std::ostream& data_log1 = parent->DataLog(0);

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

	    }

#ifdef SGS
	    if (parent->NumDataLogs() > 1) { 

	       std::ostream& data_log2 = parent->DataLog(1);

	       if (data_log2.good()) {

		  // Write the quantities that represent changes from prev_time to cur_time
		  if (time == 0.0) {
		      data_log2 << std::setw(datwidth) <<  "      time    ";
		      data_log2 << std::setw(datwidth) <<  "         Etot ";
		      data_log2 << std::setw(datwidth) <<  " Etot-sum_turb";
		      data_log2 << std::setw(datwidth) <<  "      rms_mach";
		      data_log2 << std::setw(datwidth) <<  "      delta_E ";
		      data_log2 << std::setw(datwidth) <<  "      delta_K ";
		      data_log2 << std::setw(datwidth) <<  "      prod_sgs";
		      data_log2 << std::setw(datwidth) <<  "      diss_sgs";
		      data_log2 << std::setw(datwidth) <<  "      turb_src" << std::endl;
		  }

		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << time;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << Etot;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << delta_E;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << delta_K;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << Etot-sum_turb_src;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << rms_mach;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << prod_sgs;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << diss_sgs;
		  data_log2 << std::setw(datwidth) <<  std::setprecision(datprecision) << turb_src << std::endl;
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
