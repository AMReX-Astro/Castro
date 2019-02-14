
#include <iostream>
#include "AMReX_DataServices.H"
#include <Sedov_F.H>

using namespace amrex;

std::string inputs_name = "";

//
// Prototypes
//
void GetInputArgs (const int argc, char** argv,
                   string& pltfile, string& slcfile);

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

		GetInputArgs (argc, argv, pltfile, slcfile);

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

		// get variable names
		const Vector<string>& varNames = data.PlotVarNames ();

		// get the index bounds and dx.
		const Vector<Vector<Real> >& dxLevel = data.DxLevel();
		const Vector<Real>& problo = data.ProbLo();
		const Vector<Real>& probhi = data.ProbHi();

		// compute the size of the radially-binned array -- we'll do it to
		// the furtherest corner of the domain
		double maxdist = abs(probhi[0] - problo[0]);
		double dx_fine = *(std::min_element(dxLevel[finestLevel].begin(), dxLevel[finestLevel].end()));

		auto n_bins = int(maxdist / dx_fine);

		Vector<Real> r(n_bins);
		for (auto i = 0; i < n_bins; i++)
			r[i] = (i + 0.5) * dx_fine;

		// find variable indices
		auto dens_comp = -1;
		auto xmom_comp = -1;
		auto pres_comp = -1;
		auto rhoe_comp = -1;

		for (auto i = 0; i < data.NComp(); ++i) {
			if (varNames[i] == "density") dens_comp = i;
			else if (varNames[i] == "xmom") xmom_comp = i;
			else if (varNames[i] == "pressure") pres_comp = i;
			else if (varNames[i] == "rho_e") rhoe_comp = i;
		}

		if (dens_comp < 0 || xmom_comp < 0 || pres_comp < 0 || rhoe_comp < 0) {
			Abort("ERROR: variable(s) not found");
		}

		// allocate storage for data
		Vector<Real> dens_bin(n_bins);
		Vector<Real> vel_bin(n_bins);
		Vector<Real> pres_bin(n_bins);
		Vector<Real> e_bin(n_bins);

		for (auto i = 0; i < n_bins; ++i) {
			dens_bin[i] = 0.0;
			vel_bin[i] = 0.0;
			pres_bin[i] = 0.0;
			e_bin[i] = 0.0;
		}

		auto r1 = 1.0;

		auto rr = data.RefRatio();

		// fill a multifab with the data
		Vector<int> fill_comps(data.NComp());
		for (auto i = 0; i < data.NComp(); i++) {
			fill_comps[i] = i;
        }

		MultiFab data_mf;
		data.FillVar(data_mf, finestLevel, varNames, fill_comps);

		// ! imask will be set to false if we've already output the data.
		// ! Note, imask is defined in terms of the finest level.  As we loop
		// ! over levels, we will compare to the finest level index space to
		// ! determine if we've already output here
        int mask_size = n_bins;
        for (auto i = 0; i < finestLevel - 1; i++) {
            mask_size *= rr[i];
        }

		Vector<int> imask(mask_size);

        for (auto it=imask.begin(); it!=imask.end(); ++it)
            *it = 1;

        // extract the 1d data
		for (auto l = finestLevel-1; l >= 0; l--) {

            MultiFab lev_data_mf;
            data.FillVar(lev_data_mf, l, varNames, fill_comps);

			for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
				const Box& bx = mfi.tilebox();

				fextract1d(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				           BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
				           n_bins, dens_bin.dataPtr(),
				           vel_bin.dataPtr(), pres_bin.dataPtr(),
				           e_bin.dataPtr(), imask.dataPtr(), mask_size, r1,
				           dens_comp, xmom_comp, pres_comp, rhoe_comp);
			}

            // adjust r1 for the next lowest level
            if (l != 1) r1 *= rr[l-1];
		}

        // now open the slicefile and write out the data
        std::ofstream slicefile;
        slicefile.open(slcfile);
        slicefile.precision(9);

        // write the header
        slicefile << std::setw(12) << "x" << std::setw(12) << "density" << std::setw(12) << "velocity" << std::setw(12) << "pressure" << std::setw(12) << "int. energy" << std::endl;

        // write the data in columns
        const auto SMALL = 1.e-99;
        for (auto i = 0; i < n_bins; i++) {
            if (abs(dens_bin[i]) < SMALL) dens_bin[i] = 0.0;
            if (abs( vel_bin[i]) < SMALL)  vel_bin[i] = 0.0;
            if (abs(pres_bin[i]) < SMALL) pres_bin[i] = 0.0;
            if (abs(   e_bin[i]) < SMALL)    e_bin[i] = 0.0;

            slicefile << std::setw(12) << r[i] << std::setw(12) << dens_bin[i] << std::setw(12) << vel_bin[i] << std::setw(12) << pres_bin[i] << std::setw(12) << e_bin[i] << std::endl;

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
                    string& pltfile, string& slcfile)
{

	int i = 1; // skip program name

	while ( i < argc) // Skip last two: those are the files names
	{

		if ( !strcmp(argv[i],"-p") || !strcmp(argv[i],"--pltfile") )
		{
			pltfile = argv[++i];
		}
		else if ( strcmp(argv[i],"-s") || strcmp(argv[i],"--slicefile") )
		{
			slcfile = argv[++i];
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
	        << "\nargs [-p|--pltfile]   plotfile   : plot file directory (required)"
	        << "\n     [-s|--slicefile] slice file : slice file          (required)"
	        << "\n\n" << std::endl;

}
