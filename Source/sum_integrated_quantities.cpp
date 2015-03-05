#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

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
    Real xmom        = 0.0;
    Real rho_E       = 0.0;
#ifdef SGS
    Real dt_crse     = parent->dtLevel(0);
    Real rho_e       = 0.0;
    Real rho_K       = 0.0;
    Real Etot        = 0.0;
    Real delta_E     = 0.0;
    Real delta_K     = 0.0;
    Real prod_sgs    = 0.0;
    Real diss_sgs    = 0.0;
    Real turb_src    = 0.0;
    Real rms_mach    = 0.0;
#endif
    Real com_xloc    = 0.0;
#if (BL_SPACEDIM>=2)
    Real ymom        = 0.0;
    Real com_yloc    = 0.0;
#endif
#if (BL_SPACEDIM==3)
    Real zmom        = 0.0;
    Real com_zloc    = 0.0;
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass     += ca_lev.volWgtSum("density", time, local_flag);
        xmom     += ca_lev.volWgtSum("xmom", time, local_flag);
#if (BL_SPACEDIM == 2)
       if (Geometry::IsRZ()) 
          xmom = 0.;
#endif

       if (show_center_of_mass) 
	   com_xloc += ca_lev.locWgtSum("density", time, 0, local_flag);
#if (BL_SPACEDIM>=2)
       ymom     += ca_lev.volWgtSum("ymom", time, local_flag);
       if (show_center_of_mass) 
	   com_yloc += ca_lev.locWgtSum("density", time, 1, local_flag);
#endif
#if (BL_SPACEDIM==3)
       zmom     += ca_lev.volWgtSum("zmom", time, local_flag);
       if (show_center_of_mass) 
	   com_zloc += ca_lev.locWgtSum("density", time, 2, local_flag);
#endif
       rho_E    += ca_lev.volWgtSum("rho_E", time, local_flag);

#ifdef SGS
        Real  cur_time = state[SGS_Type].curTime();
        Real prev_time = state[SGS_Type].prevTime();

        rho_e    += ca_lev.volWgtSum("rho_e", time, local_flag);
        rho_K    += ca_lev.volWgtSum("rho_K", time, local_flag);

        delta_E  += ca_lev.volWgtSum("rho_E", cur_time, local_flag);
        delta_E  -= ca_lev.volWgtSum("rho_E", prev_time, local_flag);

        delta_K  += ca_lev.volWgtSum("rho_K", cur_time, local_flag);
        delta_K  -= ca_lev.volWgtSum("rho_K", prev_time, local_flag);

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
	const int nfoo = 10+BL_SPACEDIM;
	Real foo[nfoo] = {mass, D_DECL(xmom,ymom,zmom), rho_E, rho_K, Etot, delta_E, delta_K, 
			  prod_sgs, diss_sgs, turb_src, rms_mach};
#else
	const int nfoo = 2+BL_SPACEDIM;
	Real foo[nfoo] = {mass, D_DECL(xmom,ymom,zmom), rho_E};
#endif

	Real coms[BL_SPACEDIM] = {D_DECL(com_xloc, com_yloc, com_zloc)};

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif

	ParallelDescriptor::ReduceRealSum(foo, nfoo, ParallelDescriptor::IOProcessorNumber());

	if (show_center_of_mass)
	    ParallelDescriptor::ReduceRealSum(coms, BL_SPACEDIM, ParallelDescriptor::IOProcessorNumber());

	if (ParallelDescriptor::IOProcessor()) {

	    int i = 0;
	    mass     = foo[i++];
	    xmom     = foo[i++];
#if (BL_SPACEDIM>=2)
            ymom     = foo[i++];
#endif
#if (BL_SPACEDIM==3)
            zmom     = foo[i++];
#endif
            rho_E    = foo[i++];
#ifdef SGS
	    rho_K    = foo[i++];
	    Etot     = foo[i++];
	    delta_E  = foo[i++];
	    delta_K  = foo[i++];
	    prod_sgs = foo[i++];
	    diss_sgs = foo[i++];
	    turb_sgs = foo[i++];
	    rms_mach = foo[i++];
#endif

	    std::cout << '\n';
	    std::cout << "TIME= " << time << " MASS        = "   << mass  << '\n';
	    std::cout << "TIME= " << time << " XMOM        = "   << xmom     << '\n';
#if (BL_SPACEDIM>=2)
	    std::cout << "TIME= " << time << " YMOM        = "   << ymom     << '\n';
#endif
#if (BL_SPACEDIM==3)
	    std::cout << "TIME= " << time << " ZMOM        = "   << zmom     << '\n';
#endif
	    std::cout << "TIME= " << time << " RHO*E       = "   << rho_E     << '\n';
#ifdef SGS
	    std::cout << "TIME= " << time << " RHO*K       = "   << rho_K     << '\n';
	    Etot     = rho_E + rho_K;
	    std::cout << "TIME= " << time << " TOTAL E     = "   << Etot      << '\n';
	    std::cout << "TIME= " << time << " DELTA E     = "   << delta_E   << '\n';
	    std::cout << "TIME= " << time << " DELTA K     = "   << delta_K   << '\n';
	    std::cout << "TIME= " << time << " DELTA TOT   = "   << delta_K+delta_E   << '\n';
	    std::cout << "TIME= " << time << " PROD_SGS    = "   << prod_sgs  << '\n';
	    std::cout << "TIME= " << time << " DISS_SGS    = "   << diss_sgs  << '\n';
	    std::cout << "TIME= " << time << " TURB_SRC    = "   << turb_src  << '\n';
	    std::cout << "TIME= " << time << " DE+DK-TURB_SRC = "   << delta_E+delta_K-turb_src  << '\n';

   	    std::ostream& data_log1 = parent->DataLog(0);
	    
	    if (time == 0.0) {
		data_log1 << std::setw(14) <<  "      time    ";
		data_log1 << std::setw(14) <<  "        rho_E ";
		data_log1 << std::setw(14) <<  "        rho_K ";
		data_log1 << std::setw(14) <<  "        rho_e ";
		data_log1 << std::setw(16) <<  "  Etot-sum_turb  ";
		data_log1 << std::setw(14) <<  "      rms_mach" << std::endl;
	    }

		// Write the quantities at this time
	    data_log1 << std::setw(14) <<  time;
	    data_log1 << std::setw(14) <<  std::setprecision(6) << rho_E;
	    data_log1 << std::setw(14) <<  std::setprecision(6) << rho_K;
	    data_log1 << std::setw(14) <<  std::setprecision(6) << rho_e;
	    data_log1 << std::setw(16) <<  std::setprecision(10) << Etot-sum_turb_src;
	    data_log1 << std::setw(14) <<  std::setprecision(6) << rms_mach << std::endl;

	    std::ostream& data_log2 = parent->DataLog(1);

		// Write the quantities that represent changes from prev_time to cur_time
	    if (time == 0.0) {
		data_log2 << std::setw(14) <<  "      time    ";
		data_log2 << std::setw(14) <<  "      delta_E ";
		data_log2 << std::setw(14) <<  "      delta_K ";
		data_log2 << std::setw(14) <<  "      prod_sgs";
		data_log2 << std::setw(14) <<  "      diss_sgs";
		data_log2 << std::setw(14) <<  "      turb_src" << std::endl;
	    }
	    
	    data_log2 << std::setw(14) <<  std::setprecision(6) << time;
	    data_log2 << std::setw(14) <<  std::setprecision(6) << delta_E;
	    data_log2 << std::setw(14) <<  std::setprecision(6) << delta_K;
	    data_log2 << std::setw(14) <<  std::setprecision(6) << prod_sgs;
	    data_log2 << std::setw(14) <<  std::setprecision(6) << diss_sgs;
	    data_log2 << std::setw(14) <<  std::setprecision(6) << turb_src << std::endl;
#endif
	    
	    if (show_center_of_mass) {
		com_xloc = com_xloc / mass;
		Real com_xvel = xmom / mass;
		std::cout << "TIME= " << time << " CENTER OF MASS X-LOC = " << com_xloc  << '\n';
		std::cout << "TIME= " << time << " CENTER OF MASS X-VEL = " << com_xvel  << '\n';
#if (BL_SPACEDIM>=2)
		com_yloc = com_yloc / mass;
		Real com_yvel = ymom / mass;
		std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC = " << com_yloc  << '\n';
		std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL = " << com_yvel  << '\n';
#endif
#if (BL_SPACEDIM==3)
		com_zloc = com_zloc / mass;
		Real com_zvel = zmom / mass;
		std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC = " << com_zloc  << '\n';
		std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL = " << com_zvel  << '\n';
#endif
	    }
	}
#ifdef BL_LAZY
	});
#endif
    }
}
