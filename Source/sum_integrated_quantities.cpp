#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

void
Castro::sum_integrated_quantities ()
{
#ifndef SGS
    if (verbose <= 0) return;
#endif

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
    Real com_xvel    = 0.0;
#if (BL_SPACEDIM>=2)
    Real ymom        = 0.0;
    Real com_yloc    = 0.0;
    Real com_yvel    = 0.0;
#endif
#if (BL_SPACEDIM==3)
    Real zmom        = 0.0;
    Real com_zloc    = 0.0;
    Real com_zvel    = 0.0;
#endif

    for (int lev = 0; lev <= finest_level; lev++)
    {
        Castro& ca_lev = getLevel(lev);

        mass     += ca_lev.volWgtSum("density", time);
        xmom     += ca_lev.volWgtSum("xmom", time);
#if (BL_SPACEDIM == 2)
       if (Geometry::IsRZ()) 
          xmom = 0.;
#endif

       if (show_center_of_mass) 
          com_xloc += ca_lev.locWgtSum("density", time, 0);
#if (BL_SPACEDIM>=2)
       ymom     += ca_lev.volWgtSum("ymom", time);
       if (show_center_of_mass) 
          com_yloc += ca_lev.locWgtSum("density", time, 1);
#endif
#if (BL_SPACEDIM==3)
       zmom     += ca_lev.volWgtSum("zmom", time);
       if (show_center_of_mass) 
          com_zloc += ca_lev.locWgtSum("density", time, 2);
#endif
        rho_E    += ca_lev.volWgtSum("rho_E", time);

#ifdef SGS
        Real  cur_time = state[SGS_Type].curTime();
        Real prev_time = state[SGS_Type].prevTime();

        rho_e    += ca_lev.volWgtSum("rho_e", time);
        rho_K    += ca_lev.volWgtSum("rho_K", time);

        delta_E  += ca_lev.volWgtSum("rho_E", cur_time);
        delta_E  -= ca_lev.volWgtSum("rho_E", prev_time);

        delta_K  += ca_lev.volWgtSum("rho_K", cur_time);
        delta_K  -= ca_lev.volWgtSum("rho_K", prev_time);

        rms_mach  += ca_lev.volWgtSquaredSum("MachNumber", time);

        prod_sgs += 0.5 * ca_lev.volWgtSum("prod_sgs", prev_time) * dt_crse;
        prod_sgs += 0.5 * ca_lev.volWgtSum("prod_sgs",  cur_time) * dt_crse;
        diss_sgs += 0.5 * ca_lev.volWgtSum("diss_sgs", prev_time) * dt_crse;
        diss_sgs += 0.5 * ca_lev.volWgtSum("diss_sgs",  cur_time) * dt_crse;
        turb_src += 0.5 * ca_lev.volWgtSum("turb_src", prev_time) * dt_crse;
        turb_src += 0.5 * ca_lev.volWgtSum("turb_src",  cur_time) * dt_crse;

        sum_turb_src = sum_turb_src + turb_src;
#endif
    }
 
    if (verbose > 0 && ParallelDescriptor::IOProcessor())
    {
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
	std::ostream& data_log2 = parent->DataLog(1);

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
           com_xvel = xmom / mass;
           std::cout << "TIME= " << time << " CENTER OF MASS X-LOC = " << com_xloc  << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS X-VEL = " << com_xvel  << '\n';
#if (BL_SPACEDIM>=2)
           com_yloc = com_yloc / mass;
           com_yvel = ymom / mass;
           std::cout << "TIME= " << time << " CENTER OF MASS Y-LOC = " << com_yloc  << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS Y-VEL = " << com_yvel  << '\n';
#endif
#if (BL_SPACEDIM==3)
           com_zloc = com_zloc / mass;
           com_zvel = zmom / mass;
           std::cout << "TIME= " << time << " CENTER OF MASS Z-LOC = " << com_zloc  << '\n';
           std::cout << "TIME= " << time << " CENTER OF MASS Z-VEL = " << com_zvel  << '\n';
#endif
        }
	std::cout<<'\n';
    }
}

Real
Castro::volWgtSum (const std::string& name,
                   Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMMASS,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSquaredSum (const std::string& name,
                          Real               time)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if(BL_SPACEDIM < 3) 
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if(BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 2)
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s,rad,irlo,irhi);
#elif(BL_SPACEDIM == 3)
	BL_FORT_PROC_CALL(CA_SUMSQUARED,ca_summass)
            (BL_TO_FORTRAN(fab),lo,hi,dx,&s);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSum (const std::string& name,
                   Real               time,
                   int                idir)
{
    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    MultiFab*   mf      = derive(name,time,0);

    BL_ASSERT(mf != 0);

    BoxArray baf;

    if (level < parent->finestLevel())
    {
        baf = parent->boxArray(level+1);
        baf.coarsen(fine_ratio);
    }

    for (MFIter mfi(*mf); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        if (level < parent->finestLevel())
        {
            std::vector< std::pair<int,Box> > isects = baf.intersections(grids[mfi.index()]);

            for (int ii = 0; ii < isects.size(); ii++)
            {
                fab.setVal(0,isects[ii].second,0,fab.nComp());
            }
        }
        Real s;
        const Box& box  = mfi.validbox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();
#if (BL_SPACEDIM < 3)
        const Real* rad = radius[mfi.index()].dataPtr();
        int irlo        = lo[0]-radius_grow;
        int irhi        = hi[0]+radius_grow;
#endif

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
#if (BL_SPACEDIM == 1) 
	BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,rad,irlo,irhi,idir);
#elif (BL_SPACEDIM == 2)
        int geom_flag = Geometry::IsRZ() ? 1 : 0;
        if (idir == 0 && geom_flag == 1) {
            s = 0.0;
        } else {
	   BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
               (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,rad,irlo,irhi,idir);
        }
#else
	BL_FORT_PROC_CALL(CA_SUMLOCMASS,ca_sumlocmass)
            (BL_TO_FORTRAN(fab),lo,hi,geom.ProbLo(),dx,&s,idir);
#endif
        sum += s;
    }

    delete mf;

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}
