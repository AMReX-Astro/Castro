//
// Analysis routine for RHD_shocktube
//
#include <iostream>
#include <regex>
#include <AMReX_DataServices.H>
#include <Radiation_F.H>
#include <Radiation_utils.H>

using namespace amrex;

std::string inputs_name = "";

int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv, false);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    if (AMREX_SPACEDIM != 1)
        Abort("ERROR: rhd_shocktube diagnostic only works for DIM=1");

    // Input arguments
    string pltfile, groupfile, slcfile;
    int j = 1;         // skip program name

    while ( j < argc)
    {

        if ( !strcmp(argv[j], "-p") || !strcmp(argv[j],"--pltfile") )
        {
            pltfile = argv[++j];
        }
        else if ( !strcmp(argv[j], "-g") || !strcmp(argv[j],"--groupfile") )
        {
            groupfile = argv[++j];
        }
        else if ( !strcmp(argv[j], "-s") || !strcmp(argv[j],"--slicefile") )
        {
            slcfile = argv[++j];
        }
        else
        {
            std::cout << "\n\nOption " << argv[j] << " not recognized" << std::endl;
            PrintHelp();
            exit ( EXIT_FAILURE );
        }

        // Go to the next parameter name
        ++j;
    }

    if (pltfile.empty() || groupfile.empty() || slcfile.empty())
    {
        PrintHelp();
        Abort("Missing input file");
    }

    Print() << "\nplotfile  = \"" << pltfile << "\"" << std::endl;
    Print() <<   "groupfile = \"" << groupfile << "\"" << std::endl;
    Print() <<   "slicefile = \"" << slcfile << "\"" << std::endl;
    Print() << std::endl;

    // open the group file and read in the number of groups
    std::ifstream group_file;
    group_file.open(groupfile);

    string header_line;
    string line;

    std::getline(group_file, line);
    const std::regex re("=\\s*([0-9]*)");
    std::smatch m;
    std::regex_search(line, m, re);

    int ngroups = stoi(m[1].str());

    group_file.close();

    // Start dataservices
    DataServices::SetBatchMode();

    // Define the type of file
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices (pltfile, fileType);

    if (!dataServices.AmrDataOk())
        DataServices::Dispatch(DataServices::ExitRequest, NULL);

    // get data from plot file
    AmrData& data = dataServices.AmrDataRef();

    if (data.FinestLevel() != 0)
        Abort("ERROR: rhd_shocktube only works for single level");

    // get variable names
    const Vector<string>& varNames = data.PlotVarNames();

    // get the index bounds and dx.
    Box domain = data.ProbDomain()[0];
    Vector<Real> dx = data.CellSize(0);
    const Vector<Real>& problo = data.ProbLo();
    const Vector<Real>& probhi = data.ProbHi();

    double dx_coarse = *(std::min_element(dx.begin(), dx.end()));

    int nbins = domain.length()[0];

    // x-coord
    Vector<Real> x(nbins);

    for (auto i = 0; i < nbins; i++)
        x[i] = (i + 0.5) * dx_coarse;

    // find variable indices
    Vector<std::string> compVarNames = {"density", "x_velocity", "pressure", "rad0"};
    auto varComps = GetComponents(data, compVarNames);

    auto dens_comp = varComps[0];
    auto velx_comp = varComps[1];
    auto pres_comp = varComps[2];
    auto rad0_comp = varComps[3];

    // allocate storage for data
    Vector<Real> dens_bin(nbins, 0.);
    Vector<Real> vel_bin(nbins, 0.);
    Vector<Real> pres_bin(nbins, 0.);
    Vector<Real> rad_bin(nbins, 0.);

    Vector<int> fill_comps(data.NComp());
    for (auto i = 0; i < data.NComp(); i++)
        fill_comps[i] = i;

    // extract the 1d data
    auto l = 0;

    Vector<Real> level_dx = data.DxLevel()[l];

    const BoxArray& ba = data.boxArray(l);
    const DistributionMapping& dm = data.DistributionMap(l);

    MultiFab lev_data_mf(ba, dm, data.NComp(), data.NGrow());
    data.FillVar(lev_data_mf, l, varNames, fill_comps);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(lev_data_mf, true); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        frhdshocktube(AMREX_ARLIM_3D(bx.loVect()), AMREX_ARLIM_3D(bx.hiVect()),
                      BL_TO_FORTRAN_FAB(lev_data_mf[mfi]),
                      nbins, dens_bin.dataPtr(), vel_bin.dataPtr(),
                      pres_bin.dataPtr(), rad_bin.dataPtr(),
                      dens_comp, velx_comp, pres_comp, rad0_comp,
                      ngroups);
    }

    // write slicefile
    Vector<Vector<Real> > vars = {dens_bin, vel_bin, pres_bin, rad_bin};
    Vector<std::string> slcvarNames = {"density", "velocity", "pressure", "rad"};
    WriteSlicefile(nbins, x, slcvarNames, vars, slcfile);

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
