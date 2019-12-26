//
// Process a plotfile to find the location where the timestep is being limited
//
#include <iostream>
#include <fstream>
// #include <stringstream>
#include <regex>
#include "AMReX_DataServices.H"
#include <AMReX_ParmParse.H>
#include <AMReX_BaseFab.H>
#include <Limiter_F.H>
#include "Castro_F.H"

using namespace amrex;

std::string inputs_name = "";

//
// Prototypes
//
void GetInputArgs (const int argc, char** argv,
                   string& pltfile);

void ProcessJobInfo(string job_info_file, string inputs_file_name);


int main(int argc, char* argv[])
{

    int dummy = 0;

	amrex::Initialize(dummy, argv);

	// timer for profiling
	BL_PROFILE_VAR("main()", pmain);

    // Tell AMReX to be quiet
    amrex::system::verbose = 0;

	// Input arguments
	string pltfile;

	GetInputArgs (argc, argv, pltfile);

    string job_info = pltfile + "/job_info";
    string inputs_file = "inputs.txt";

    ProcessJobInfo(job_info, inputs_file);

    // this will process the input parameters from the job_info file 
    amrex::ParmParse::Initialize(0,0,inputs_file.c_str());

    // Read in the input values to Fortran.
    int NUM_GROW;

    ca_get_method_params(&NUM_GROW);
    ca_set_castro_method_params();

	// Start dataservices
	DataServices::SetBatchMode();

	// Define the type of file
	Amrvis::FileType fileType(Amrvis::NEWPLT);
	DataServices dataServices (pltfile, fileType);

	if (!dataServices.AmrDataOk())
		DataServices::Dispatch(DataServices::ExitRequest, NULL);

    // initialize microphysics stuff 
    auto probin_name = "probin";
    auto probin_file_length = 6;
    Vector<int> probin_file(probin_file_length);
    for (auto i = 0; i < probin_file_length; i++) {
        probin_file[i] = probin_name[i];
    }

    microphysics_initialize(probin_file.dataPtr(), &probin_file_length);

	// get data from plot file
	AmrData& data = dataServices.AmrDataRef();

	// get variable names
	const Vector<string>& varNames = data.PlotVarNames();

    int dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp, spec_comp, temp_comp;
    int time_integration_method = 3;

	// find variable indices
	dens_comp = data.StateNumber("density");
	xmom_comp = data.StateNumber("xmom");
#if (AMREX_SPACEDIM >= 2)
	ymom_comp = data.StateNumber("ymom");
#endif
#if (AMREX_SPACEDIM == 3)
	zmom_comp = data.StateNumber("zmom");
#endif
	pres_comp = data.StateNumber("pressure");
	rhoe_comp = data.StateNumber("rho_e");
	temp_comp = data.StateNumber("Temp");

	if (dens_comp < 0 || xmom_comp < 0 || pres_comp < 0 || rhoe_comp < 0)
		Abort("ERROR: variable(s) not found");

#if (AMREX_SPACEDIM == 3)
	if (ymom_comp < 0 || zmom_comp < 0)
		Abort("ERROR: variable(s) not found");
#endif

    // we're going to find the spec comp by looking for the first variable name
    // that begins with 'X'
    std::string first_spec_name = "";
    for (auto &it : varNames) {
        if (it[0] == 'X') {
            first_spec_name = it;
            break;
        }
    }

    if (first_spec_name == "")
        Abort("ERROR: no species were found");

    spec_comp = data.StateNumber(first_spec_name);

	// fill a multifab with the data
	Vector<int> fill_comps(data.NComp());
	for (auto i = 0; i < data.NComp(); i++)
		fill_comps[i] = i;

    Vector<Real> dt_loc = {0.,0.,0.};
    Vector<Real> burning_dt_loc = {0.,0.,0.};
    Vector<Real> diffusion_dt_loc = {0.,0.,0.};
    Real dt = 1.e99;
    Real burning_dt = 1.e99;
    Real diffusion_dt = 1.e99;

	// loop over the data
	for (int lev=0; lev <= data.FinestLevel(); lev++) {

        Vector<Real> level_dx = data.DxLevel()[lev];

		const BoxArray& ba = data.boxArray(lev);
		const DistributionMapping& dm = data.DistributionMap(lev);

		MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());

        for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
			const Box& bx = mfi.tilebox();

            // we only want to load the data a Fab at the time as otherwise 
            // we can experience memory issues. 
            // Therefore, we'll iterate over the variables using FillVar
            // to fill a temporary fab, then copy this over the to MultiFab
            Array4<Real> lev_data_arr = lev_data_mf.array(mfi);
            FArrayBox tmp_fab(bx);

            for (auto n = 0; n < data.NComp(); n++) {
                data.FillVar(&tmp_fab, bx, lev, varNames[n], ParallelDescriptor::IOProcessorNumber());

                Array4<Real> tmp_arr = tmp_fab.array();

                AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                    lev_data_arr(i,j,k,n) = tmp_arr(i,j,k);
                });
            }

            find_timestep_limiter(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                            BL_TO_FORTRAN_FAB(lev_data_mf[mfi]), 
                            dens_comp, xmom_comp, ymom_comp, zmom_comp, pres_comp, rhoe_comp, spec_comp,
                            time_integration_method,
                            level_dx.dataPtr(), &dt, dt_loc.dataPtr());

            find_timestep_limiter_burning(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
                            BL_TO_FORTRAN_FAB(lev_data_mf[mfi]), 
                            dens_comp, temp_comp, rhoe_comp, spec_comp,
                            level_dx.dataPtr(), &burning_dt, burning_dt_loc.dataPtr());
#ifdef DIFFUSION
				find_timestep_limiter_diffusion(ARLIM_3D(bx.loVect()), ARLIM_3D(bx.hiVect()),
				               BL_TO_FORTRAN_FAB(lev_data_mf[mfi]), 
				               dens_comp, temp_comp, rhoe_comp, spec_comp,
				               level_dx.dataPtr(), &diffusion_dt, diffusion_dt_loc.dataPtr());
#endif
		}
	}
    Print() << std::endl;

    Print() << "dt = " << dt << " at location";
    for (auto i = 0; i < AMREX_SPACEDIM; i++) {
        Print() << ' ' << dt_loc[i];
    }
    Print() << std::endl;

    Print() << "burning_dt = " << burning_dt << " at location";
    for (auto i = 0; i < AMREX_SPACEDIM; i++) {
        Print() << ' ' << burning_dt_loc[i];
    }
    Print() << std::endl;

#ifdef DIFFUSION
    Print() << "diffusion_dt = " << diffusion_dt << " at location";
    for (auto i = 0; i < AMREX_SPACEDIM; i++) {
        Print() << ' ' << diffusion_dt_loc[i];
    }
    Print() << std::endl;
#endif

    Print() << std::endl;

    // finalize microphysics stuff 
    microphysics_finalize();

	// destroy timer for profiling
	BL_PROFILE_VAR_STOP(pmain);

	amrex::Finalize();
}

//
// Parse command line arguments
//
void GetInputArgs ( const int argc, char** argv,
                    string& pltfile)
{
    pltfile = argv[1];

	if (pltfile.empty())
	{
		Abort("Missing input file");
	}

	Print() << "Finding limiting timestep in plotfile  = \"" << pltfile << "\"" << std::endl;
}

//
// Reads in a job_info file, extracts the inputs parameters and saves 
// them to a new file 
//
void ProcessJobInfo(string job_info_file, string inputs_file_name)
{
    std::ifstream job_info (job_info_file);
    std::ofstream inputs_file (inputs_file_name);

    std::regex inputs_rgx("Inputs File Parameters");

    if (job_info.is_open() && inputs_file.is_open()) {
        string line;
        bool found_inputs = false;
        while (getline(job_info, line)) {
            if (found_inputs) {
                // if the last character is a space, then there is no variable
                if (line.back() != ' ') {
                    if (line.front() == '[') {
                        inputs_file << line.substr(3) << std::endl;
                    } else {
                        inputs_file << line << std::endl;
                    }
                } 
            } else {
                if (std::regex_search(line, inputs_rgx)) {
                    // found a match! 
                    found_inputs = true;

                    // skip a line 
                    getline(job_info, line);
                }
            }
        }

        job_info.close();
        inputs_file.close();
    }
}