#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <limits>
#include <iterator>
#include <fstream>
#include <cmath>
#include <iterator>
#include <algorithm>

// find the thermodynamic state corresponding to the larged abs(enuc)
// and output it

using namespace amrex;

void main_main()
{
    const int narg = amrex::command_argument_count();

    std::string pltfile;
    std::string varname_arg;

    // the executable name is the first arg

    int farg{1};

    if (narg != 1) {
        amrex::Print()
            << "\n"
            << " Output the thermodynamic state corresponding to the largest abs(enuc)\n"
            << " we only consider the finest level\n"
            << " Usage:\n"
            << "    fburn_weight_avg plotfile\n"
            << "\n"
            << std::endl;
        return;
    }

    int ntime = narg;

    const std::string& filename = amrex::get_command_argument(farg);

    PlotFileData pf(filename);

    const Vector<std::string>& var_names_pf = pf.varNames();

    // we need rho, T, X, and enuc

    auto ienuc = std::distance(var_names_pf.cbegin(), std::find(var_names_pf.cbegin(), var_names_pf.cend(), "enuc"));

    int fine_level = pf.finestLevel();

    Vector<Real> lstate(var_names_pf.size(), 0.0);

    int level = fine_level;

    Real enuc_max = std::numeric_limits<Real>::lowest();

    const MultiFab& mf = pf.get(level);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        if (bx.ok()) {
            const auto& fab = mf.array(mfi);
            const auto lo = amrex::lbound(bx);
            const auto hi = amrex::ubound(bx);
            for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                    for (int i = lo.x; i <= hi.x; ++i) {
                        if (std::abs(fab(i,j,k,ienuc)) > enuc_max) {
                            enuc_max = std::abs(fab(i,j,k,ienuc));
                            for (int ivar = 0; ivar < var_names_pf.size(); ++ivar) {
                                lstate[ivar] = fab(i,j,k,ivar);
                            }
                        }
                    }
                }
            }
        }
    }

    std::cout << "enuc_max = " << enuc_max << std::endl;

    //ParallelDescriptor::ReduceRealSum(lstate.data(), lstate.size());

    // output the header
    for (int ivar = 0; ivar < var_names_pf.size(); ++ivar) {
        std::cout << std::setw(25) << std::left << var_names_pf[ivar] << std::setw(25) << lstate[ivar] << std::endl;
    }
    std::cout << std::endl;

}

int main (int argc, char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
}
