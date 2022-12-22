//
// Process a sedov problem to produce rho, u, and p as a
// function of r, for comparison to the analytic solution.
//
#include <iostream>
// #include <stringstream>
#include <regex>
#include <string>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

std::string inputs_name = "";

//
// Prototypes
//
void GetInputArgs (const int argc, char** argv,
                   std::string& pltfile, std::string& slcfile,
                   bool& sphr);

std::string GetVarFromJobInfo (const std::string pltfile, const std::string varname);

Vector<Real> GetCenter (const std::string pltfile);

void PrintHelp ();


std::pair<Real, Real> get_coord_info(const Array<Real, AMREX_SPACEDIM>& p,
                                     const Vector<Real>& center,
                                     const Array<Real, AMREX_SPACEDIM>& dx_level,
                                     const int coord, const bool sphr) {


    // compute radial coord and index

    Real r_zone{0.0_rt};
    Real vol = std::numeric_limits<Real>::quiet_NaN();

#if AMREX_SPACEDIM == 1
    // 1-d spherical geometry / spherical Sedov explosion

    AMREX_ASSERT(coord == 2);

    r_zone = p[0] - center[0];
    Real r_r = problo[0]+static_cast<Real>(i+1)*dx_level[0];
    Real r_l = problo[0]+static_cast<Real>(i)*dx_level[0];
    vol = (4.0_rt/3.0_rt) * M_PI * dx_level[0] *
        (r_r*r_r + r_l*r_r + r_l*r_l);

#elif AMREX_SPACEDIM == 2
    if (sphr) {
        // 2-d axisymmetric geometry / spherical Sedov explosion

        AMREX_ASSERT(coord == 1);

        // axisymmetric V = pi (r_r**2 - r_l**2) * dz
        //                = pi dr * dz * (r_r + r_l)
        //                = 2 pi r dr dz

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]));
        vol = 2 * M_PI * p[0] * dx_level[0] * dx_level[1];

    } else {
        // 2-d Cartesian geometry / cylindrical Sedov explosion

        AMREX_ASSERT(coord == 0);

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]));
        vol = dx_level[0] * dx_level[1];

    }

#else

    AMREX_ASSERT(coord == 0);

    vol = dx_level[0] * dx_level[1] * dx_level[2];

    if (sphr) {
        // 3-d Cartesian geometry / spherical Sedov explosion

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]) +
                           (p[2] - center[2]) * (p[2] - center[2]));

    } else {
        // 3-d Cartesian geometry / cylindrical Sedov explosion

        r_zone = std::sqrt((p[0] - center[0]) * (p[0] - center[0]) +
                           (p[1] - center[1]) * (p[1] - center[1]));

    }
#endif

    return {r_zone, vol};

}

int main(int argc, char* argv[])
{

    amrex::Initialize(argc, argv, false);

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // Input arguments
    std::string pltfile, slcfile;
    bool sphr = false;

    GetInputArgs (argc, argv, pltfile, slcfile, sphr);

    auto center = GetCenter(pltfile);
    double xctr = center[0];
    double yctr = center[1];
    double zctr = center[2];

    PlotFileData pf(pltfile);

    int fine_level = pf.finestLevel();
    const int dim = pf.spaceDim();

    // get the index bounds and dx.
    Box domain = pf.probDomain(fine_level);
    auto dx = pf.cellSize(fine_level);
    auto problo = pf.probLo();
    auto probhi = pf.probHi();

    // compute the size of the radially-binned array -- we'll do it to
    // the furtherest corner of the domain

#if (AMREX_SPACEDIM == 1)
    double maxdist = std::abs(probhi[0] - problo[0]);

#elif (AMREX_SPACEDIM == 2)
    double x_maxdist = amrex::max(std::abs(probhi[0] - xctr),
                                  std::abs(problo[0] - xctr));
    double y_maxdist = amrex::max(std::abs(probhi[1] - yctr),
                                  std::abs(problo[1] - yctr));
    double maxdist = std::sqrt(x_maxdist*x_maxdist +
                               y_maxdist*y_maxdist);

#else
    double x_maxdist = amrex::max(std::abs(probhi[0] - xctr),
                                  std::fabs(problo[0] - xctr));
    double y_maxdist = amrex::max(std::abs(probhi[1] - yctr),
                                  std::abs(problo[1] - yctr));
    double z_maxdist = amrex::max(std::abs(probhi[2] - zctr),
                                  std::abs(problo[2] - zctr));
    double maxdist = std::sqrt(x_maxdist*x_maxdist +
                               y_maxdist*y_maxdist +
                               z_maxdist*z_maxdist);
#endif

    double dx_fine = *(std::min_element(dx.begin(), dx.end()));

    int nbins = int(maxdist / dx_fine);

    // radial coordinate
    Vector<Real> r(nbins);

    for (auto i = 0; i < nbins; i++)
        r[i] = (i + 0.5) * dx_fine;

    // find variable indices
    const Vector<std::string>& var_names_pf = pf.varNames();

    int dens_comp = std::distance(var_names_pf.cbegin(),
                              std::find(var_names_pf.cbegin(), var_names_pf.cend(), "density"));

    int xmom_comp = std::distance(var_names_pf.cbegin(),
                              std::find(var_names_pf.cbegin(), var_names_pf.cend(), "xmom"));

    int ymom_comp = std::distance(var_names_pf.cbegin(),
                              std::find(var_names_pf.cbegin(), var_names_pf.cend(), "ymom"));

    int zmom_comp = std::distance(var_names_pf.cbegin(),
                              std::find(var_names_pf.cbegin(), var_names_pf.cend(), "zmom"));

    int pres_comp = std::distance(var_names_pf.cbegin(),
                              std::find(var_names_pf.cbegin(), var_names_pf.cend(), "pressure"));

    int rhoe_comp = std::distance(var_names_pf.cbegin(),
                              std::find(var_names_pf.cbegin(), var_names_pf.cend(), "rho_e"));

    int coord = pf.coordSys();

    // allocate storage for data
    Vector<Real> dens_bin(nbins, 0.);
    Vector<Real> vel_bin(nbins, 0.);
    Vector<Real> pres_bin(nbins, 0.);
    Vector<Real> e_bin(nbins, 0.);
    Vector<Real> volcount(nbins, 0);

    // we will use a mask that tells us if a zone on the current level
    // is covered by data on a finer level.

    for (int ilev = 0; ilev <= fine_level; ++ilev) {

        Array<Real, AMREX_SPACEDIM> dx_level = pf.cellSize(ilev);

        if (ilev < fine_level) {
            IntVect ratio{pf.refRatio(ilev)};
            for (int idim = dim; idim < AMREX_SPACEDIM; ++idim) {
                ratio[idim] = 1;
            }
            const iMultiFab mask = makeFineMask(pf.boxArray(ilev), pf.DistributionMap(ilev),
                                                pf.boxArray(ilev+1), ratio);

            const MultiFab& lev_data_mf = pf.get(ilev);

            for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.ok()) {
                    const auto& m = mask.array(mfi);
                    const auto& fab = lev_data_mf.array(mfi);
                    const auto lo = amrex::lbound(bx);
                    const auto hi = amrex::ubound(bx);

                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {
                                if (m(i,j,k) == 0) { // not covered by fine

                                    Array<Real,AMREX_SPACEDIM> p
                                        = {AMREX_D_DECL(problo[0]+static_cast<Real>(i+0.5)*dx_level[0],
                                                        problo[1]+static_cast<Real>(j+0.5)*dx_level[1],
                                                        problo[2]+static_cast<Real>(k+0.5)*dx_level[2])};


                                    const auto& [r_zone, vol] = get_coord_info(p, center, dx_level, coord, sphr);

                                    int index = static_cast<int>(r_zone / dx_fine);

                                    // add to the bin, weighting by the size

                                    dens_bin[index] += fab(i,j,k,dens_comp) * vol;

                                    vel_bin[index] +=
                                        std::sqrt(std::pow(fab(i,j,k,xmom_comp), 2) +
                                                  std::pow(fab(i,j,k,ymom_comp), 2) +
                                                  std::pow(fab(i,j,k,zmom_comp), 2)) /
                                        fab(i,j,k,dens_comp) * vol;

                                    pres_bin[index] += fab(i,j,k,pres_comp) * vol;

                                    e_bin[index] += (fab(i,j,k,rhoe_comp) / fab(i,j,k,dens_comp)) * vol;

                                    volcount[index] += vol;

                                } // mask

                            }
                        }
                    }

                } // bx.ok()

            } // MFIter

        } else {
            // this is the finest level

            const MultiFab& lev_data_mf = pf.get(ilev);

            for (MFIter mfi(lev_data_mf); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                if (bx.ok()) {
                    const auto& fab = lev_data_mf.array(mfi);
                    const auto lo = amrex::lbound(bx);
                    const auto hi = amrex::ubound(bx);

                    for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                            for (int i = lo.x; i <= hi.x; ++i) {

                                Array<Real,AMREX_SPACEDIM> p
                                    = {AMREX_D_DECL(problo[0]+static_cast<Real>(i+0.5)*dx_level[0],
                                                    problo[1]+static_cast<Real>(j+0.5)*dx_level[1],
                                                    problo[2]+static_cast<Real>(k+0.5)*dx_level[2])};

                                const auto& [r_zone, vol] = get_coord_info(p, center, dx_level, coord, sphr);

                                int index = static_cast<int>(r_zone / dx_fine);

                                // add to the bin, weighting by the size

                                dens_bin[index] += fab(i,j,k,dens_comp) * vol;

                                vel_bin[index] +=
                                    std::sqrt(std::pow(fab(i,j,k,xmom_comp), 2) +
                                              std::pow(fab(i,j,k,ymom_comp), 2) +
                                              std::pow(fab(i,j,k,zmom_comp), 2)) /
                                    fab(i,j,k,dens_comp) * vol;

                                pres_bin[index] += fab(i,j,k,pres_comp) * vol;

                                e_bin[index] += (fab(i,j,k,rhoe_comp) / fab(i,j,k,dens_comp)) * vol;

                                volcount[index] += vol;
                                
                            }
                        }
                    }

                } // bx.ok()

            } // MFIter


        }

    } // level loop


    //normalize
    for (int i = 0; i < nbins; i++) {
        if (volcount[i] != 0.0) {
            dens_bin[i] /= volcount[i];
            vel_bin[i] /= volcount[i];
            pres_bin[i] /= volcount[i];
            e_bin[i] /= volcount[i];
        }
    }

    // now open the slicefile and write out the data
    std::ofstream slicefile;
    slicefile.open(slcfile);
    slicefile.setf(std::ios::scientific);
    slicefile.precision(12);
    const auto w = 24;

    // write the header
    slicefile << "# " << std::setw(w) << "x" << std::setw(w) << "density" << std::setw(w) << "velocity" << std::setw(w) << "pressure" << std::setw(w) << "int. energy" << std::endl;

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

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}

//
// Parse command line arguments
//
void GetInputArgs ( const int argc, char** argv,
                    std::string& pltfile, std::string& slcfile,
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

#if (AMREX_SPACEDIM == 1)
    Print() << "Extracting slice from 1d problem" << std::endl;
#elif (AMREX_SPACEDIM >= 2)
    if (sphr) {
        Print() << "Extracting slice from " << AMREX_SPACEDIM << "d spherical problem" << std::endl;
    } else {
        Print() << "Extracting slice from " << AMREX_SPACEDIM << "d cylindrical problem" << std::endl;
    }
#endif

    Print() << "\nplotfile  = \"" << pltfile << "\"" << std::endl;
    Print() << "slicefile = \"" << slcfile << "\"" << std::endl;
    Print() << std::endl;
}

///
/// Gets the variable ``varname`` from the ``job_info`` file and returns as a
/// string
///
std::string GetVarFromJobInfo (const std::string pltfile, const std::string varname) {
    std::string filename = pltfile + "/job_info";
    std::regex re("(?:[ \\t]*)" + varname + "\\s*:\\s*(.*)\\s*\\n");

    std::smatch m;

    std::ifstream jobfile(filename);
    if (jobfile.is_open()) {
        std::stringstream buf;
        buf << jobfile.rdbuf();
        std::string file_contents = buf.str();

        if (std::regex_search(file_contents, m, re)) {
            return m[1];
        } else {
            Print() << "Unable to find " << varname << " in job_info file!" << std::endl;
        }
    } else {
        Print() << "Could not open job_info file!" << std::endl;
    }

    return "";
}

// Get the center from the job info file and return as a Real Vector
Vector<Real> GetCenter (const std::string pltfile) {
    auto center_str = GetVarFromJobInfo(pltfile, "center");

    // split string
    std::istringstream iss {center_str};
    Vector<Real> center;

    std::string s;
    while (std::getline(iss, s, ','))
        center.push_back(stod(s));

    return center;
}

//
// Print usage info
//
void PrintHelp ()
{
    Print() << "\nusage: executable_name args"
            << "\nargs [-p|--pltfile]     plotfile : plot file directory (required)"
            << "\n     [-s|--slicefile] slice file : slice file          (required)"
#if AMREX_SPACEDIM >= 2
            << "\n     [--sphr]          spherical : spherical problem"
#endif
            << "\n\n" << std::endl;

}
