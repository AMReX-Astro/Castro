
#include <iostream>
#include "AMReX_DataServices.H"
#include <Sedov_F.H>

using namespace amrex;

std::string inputs_name = "";

//
// Prototypes
//
void GetInputArgs (const int argc, char** argv,
                   string& pltfile, string& slcfile,
                   double& xctr, double& yctr, double & zctr,
                   bool& sphr);

void PrintHelp ();


int main(int argc, char* argv[])
{

	amrex::Initialize(argc, argv, false);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

	// wallclock time
	const Real strt_total = amrex::second();

	{
		// Input arguments
		string pltfile;
		string slcfile;
		double xctr = 0.0;
		double yctr = 0.0;
		double zctr = 0.0;
		bool sphr = false;

		GetInputArgs (argc, argv, pltfile, slcfile, xctr, yctr, zctr, sphr);

		// Start dataservices (no clue why we need to do this)
		DataServices::SetBatchMode();

		// Define the type of file
		Amrvis::FileType fileType(Amrvis::NEWPLT);
		DataServices dataServices (pltfile, fileType);

		if (!dataServices.AmrDataOk())
			DataServices::Dispatch(DataServices::ExitRequest, NULL);

		// get data from plot file
		AmrData& data = dataServices.AmrDataRef();

		int finestLevel = data.FinestLevel();

		const int factor = 2;

		// get variable names
		const Vector<string>& varNames = data.PlotVarNames();

		// get the index bounds and dx.
		Box domain = data.ProbDomain()[finestLevel];

		Vector<Real> dx(AMREX_SPACEDIM);
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			dx[i] = data.ProbSize()[i] / domain.length(i);

		// const Vector<Vector<Real> >& dx = data.DxLevel();
		const Vector<Real>& problo = data.ProbLo();
		const Vector<Real>& probhi = data.ProbHi();

		// compute the size of the radially-binned array -- we'll do it to
		// the furtherest corner of the domain
#if (AMREX_SPACEDIM == 1)
		double maxdist = fabs(probhi[0] - problo[0]);
#elif (AMREX_SPACEDIM == 2)
		double x_maxdist = max(fabs(probhi[0] - xctr), fabs(problo[0] - xctr));
		double y_maxdist = max(fabs(probhi[1] - yctr), fabs(problo[1] - yctr));
		double maxdist = sqrt(x_maxdist*x_maxdist + y_maxdist*y_maxdist);
#else
		double maxdist = max(abs(probhi[0] - problo[0]), abs(probhi[1] - problo[1]), abs(probhi[2] - problo[2]));
#endif

		double dx_fine = *(std::min_element(dx.begin(), dx.end()));

		int nbins = int(maxdist / dx_fine);

		Vector<Real> r(nbins);

		for (auto i = 0; i < nbins; i++)
			r[i] = (i + 0.5) * dx_fine;

		// find variable indices
		auto dens_comp = data.StateNumber("density");
		auto xmom_comp = data.StateNumber("xmom");
#if (AMREX_SPACEDIM >= 2)
		auto ymom_comp = data.StateNumber("ymom");
#endif
#if (AMREX_SPACEDIM == 3)
		auto zmom_comp = data.StateNumber("zmom");
#endif
		auto pres_comp = data.StateNumber("pressure");
		auto rhoe_comp = data.StateNumber("rho_e");

		if (dens_comp < 0 || xmom_comp < 0 || pres_comp < 0 || rhoe_comp < 0)
			Abort("ERROR: variable(s) not found");

#if (AMREX_SPACEDIM == 3)
		if (ymom_comp < 0 || zmom_comp < 0)
			Abort("ERROR: variable(s) not found");
#endif

		// allocate storage for data
		Vector<Real> dens_bin(nbins, 0.);
		Vector<Real> vel_bin(nbins, 0.);
		Vector<Real> pres_bin(nbins, 0.);
		Vector<Real> e_bin(nbins, 0.);
		Vector<int> ncount(nbins, 0);

		auto r1 = 1.0;

		Vector<int> rr = data.RefRatio();

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++) {
			fill_comps[i] = i;
		}

		const BoxArray& ba_fine = data.boxArray(finestLevel);
		const DistributionMapping& dm_fine = data.DistributionMap(finestLevel);

		MultiFab data_mf(ba_fine, dm_fine, data.NComp(), data.NGrow());
		data.FillVar(data_mf, finestLevel, varNames, fill_comps);

		// ! imask will be set to false if we've already output the data.
		// ! Note, imask is defined in terms of the finest level.  As we loop
		// ! over levels, we will compare to the finest level index space to
		// ! determine if we've already output here
		int mask_size = nbins;
		for (auto i = 0; i < finestLevel - 1; i++)
			mask_size *= rr[i];


#if (AMREX_SPACEDIM == 1)
		Vector<int> imask(mask_size);
#elif (AMREX_SPACEDIM >=2)
		Vector<int> imask(pow(mask_size, AMREX_SPACEDIM));
#endif

		for (auto it=imask.begin(); it!=imask.end(); ++it)
			*it = 1;

		// extract the 1d data
		for (int l = finestLevel; l >= 0; l--) {

			int refratio = 1;

			for (int ll = 0; ll < l; ll++) refratio *= rr[ll];

			Vector<Real> level_dx = data.DxLevel()[l];

			const BoxArray& ba = data.boxArray(l);
			const DistributionMapping& dm = data.DistributionMap(l);

			MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
			data.FillVar(lev_data_mf, l, varNames, fill_comps);

			for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
				const Box& bx = mfi.tilebox();

#if (AMREX_SPACEDIM == 1)
				fextract1d(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				           BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				           nbins, (*dens_bin).dataPtr(),
				           vel_bin.dataPtr(), pres_bin.dataPtr(),
				           e_bin.dataPtr(), imask.dataPtr(), mask_size, r1,
				           dens_comp, xmom_comp, pres_comp, rhoe_comp);
#elif (AMREX_SPACEDIM == 2)
				if (sphr) {

                    fextract2d_sph(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
					               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
					               nbins, dens_bin.dataPtr(),
					               vel_bin.dataPtr(), pres_bin.dataPtr(),
					               e_bin.dataPtr(), ncount.dataPtr(),
					               imask.dataPtr(), mask_size, r1,
					               dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp,
					               dx_fine, level_dx.dataPtr(),
					               xctr, yctr);

				} else {
					fextract2d_cyl(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
					               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
					               nbins, dens_bin.dataPtr(),
					               vel_bin.dataPtr(), pres_bin.dataPtr(),
					               e_bin.dataPtr(), ncount.dataPtr(),
					               imask.dataPtr(), mask_size, r1,
					               dens_comp, xmom_comp, ymom_comp, pres_comp, rhoe_comp,
					               dx_fine, level_dx.dataPtr(),
					               xctr, yctr);
				}
#else
				if (sphr) {
                    fextract3d_sph(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
					               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
					               nbins, dens_bin.dataPtr(),
					               vel_bin.dataPtr(), pres_bin.dataPtr(),
                                   e_bin.dataPtr(), ncount.dataPtr(),
					               imask.dataPtr(), mask_size, r1,
					               dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp,
					               dx_fine, level_dx.dataPtr(), xctr, yctr, zctr);

				} else {
					fextract3d_cyl(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
					               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
					               nbins, dens_bin.dataPtr(),
					               vel_bin.dataPtr(), pres_bin.dataPtr(),
                                   ncount.dataPtr(),
					               imask.dataPtr(), mask_size, r1,
					               dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp,
					               dx_fine, level_dx.dataPtr(), refratio, xctr, yctr);
				}
#endif
			}

			// adjust r1 for the next lowest level
			if (l != 0) r1 *= rr[l-1];
		}

#if (AMREX_SPACEDIM >= 2)
		//normalize
		for (int i = 0; i < nbins; i++) {
			if (ncount[i] != 0) {
				dens_bin[i] /= ncount[i];
				vel_bin[i] /= ncount[i];
				pres_bin[i] /= ncount[i];
				e_bin[i] /= ncount[i];
			}
		}
#endif

		// now open the slicefile and write out the data
		std::ofstream slicefile;
		slicefile.open(slcfile);
		slicefile.setf(std::ios::scientific);
		slicefile.precision(12);
		const auto w = 24;

		// write the header
		slicefile << std::setw(w) << "x" << std::setw(w) << "density" << std::setw(w) << "velocity" << std::setw(w) << "pressure" << std::setw(w) << "int. energy" << std::endl;

		// write the data in columns
		const auto SMALL = 1.e-20;
		for (auto i = 0; i < nbins; i++) {
			if (fabs(dens_bin[i]) < SMALL) dens_bin[i] = 0.0;
			if (fabs( vel_bin[i]) < SMALL) vel_bin[i] = 0.0;
			if (fabs(pres_bin[i]) < SMALL) pres_bin[i] = 0.0;
			if (fabs(   e_bin[i]) < SMALL) e_bin[i] = 0.0;

			slicefile << std::setw(w) << r[i] << std::setw(w) << dens_bin[i] << std::setw(w) << vel_bin[i] << std::setw(w) << pres_bin[i] << std::setw(w) << e_bin[i] << std::endl;

		}

		slicefile.close();
	}


	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}

//
// Parse command line arguments
//
void GetInputArgs ( const int argc, char** argv,
                    string& pltfile, string& slcfile,
                    double& xctr, double& yctr, double& zctr,
                    bool &sphr)
{

	int i = 1; // skip program name

	while ( i < argc)
	{

		if ( !strcmp(argv[i], "-p") || !strcmp(argv[i],"--pltfile") )
		{
			pltfile = argv[++i];
		}
		else if ( !strcmp(argv[i], "-s") || !strcmp(argv[i],"--slicefile") )
		{
			slcfile = argv[++i];
		}
		else if ( !strcmp(argv[i],"--xctr") )
		{
			xctr = std::atof(argv[++i]);
		}
		else if ( !strcmp(argv[i],"--yctr") )
		{
			yctr = std::atof(argv[++i]);
		}
		else if ( !strcmp(argv[i],"--zctr") )
		{
			zctr = std::atof(argv[++i]);
		}
		else if ( !strcmp(argv[i],"--sphr") )
		{
			sphr = true;
		}
		else
		{
			std::cout << "\n\nOption " << argv[i] << " not recognized" << std::endl;
			PrintHelp ();
			exit ( EXIT_FAILURE );
		}

		// Go to the next parameter name
		++i;
	}

	if (pltfile.empty() && slcfile.empty())
	{
		PrintHelp();
		Abort("Missing input file");
	}

	Print() << "\npltfile = \"" << pltfile << "\"" << std::endl;
	Print() << "slcfile = \"" << slcfile << "\"" << std::endl;
	Print() << std::endl;
}



//
// Print usage info
//
void PrintHelp ()
{
	Print() << "\nusage: executable_name args"
	        << "\nargs [-p|--pltfile]     plotfile : plot file directory (required)"
	        << "\n     [-s|--slicefile] slice file : slice file          (required)"
	        << "\n     [--xctr]               xctr : central x coord     (non-cartesian only)"
	        << "\n     [--yctr]               yctr : central y coord     (non-cartesian only)"
	        << "\n     [--zctr]               zctr : central z coord     (non-cartesian only)"
	        << "\n     [--sphr]          spherical : spherical problem"
	        << "\n\n" << std::endl;

}
