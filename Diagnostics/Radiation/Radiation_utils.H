#ifndef _Radiation_utils_H_
#define _Radiation_utils_H_
#include <AMReX_BLFort.H>
#include <AMReX_DataServices.H>
#include <regex>

using namespace amrex;

template <typename T>
Vector<size_t> sort_indexes(const Vector<T> &v);

void GetInputArgs ( const int argc, char** argv,
                    string& pltfile, string& slcfile, int& dir);

Vector<int> GetComponents(AmrData& data, const Vector<std::string> varNames);

void PrintHelp ();

string GetVarFromJobInfo (const string pltfile, const string varname);

Vector<Real> GetCenter (const string pltfile);


template <typename T>
Vector<size_t> sort_indexes(const Vector<T> &v) {

    // initialize original index locations
    Vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {
        return v[i1] < v[i2];
    });

    return idx;
}

//
// Parse command line arguments
//
void GetInputArgs ( const int argc, char** argv,
                    string& pltfile, string& slcfile, int& dir)
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
        else if ( !strcmp(argv[i], "-d") || !strcmp(argv[i],"--direction") )
        {
            dir = std::atoi(argv[++i]);
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

    Print() << "\nplotfile  = \"" << pltfile << "\"" << std::endl;
    Print() << "slicefile = \"" << slcfile << "\"" << std::endl;
    Print() << std::endl;
}

//
// Get the indices of the components
//
Vector<int> GetComponents(AmrData& data, const Vector<std::string> varNames) {

    Vector<int> varComps;

    for (auto it=varNames.begin(); it!=varNames.end(); ++it)
        varComps.push_back(data.StateNumber(*it));

    return varComps;

}

//
// Write the slice file
//
void WriteSlicefile(const int nbins, const Vector<Real> r,
                    const Vector<std::string> varNames,
                    Vector<Vector<Real> > vars,
                    const std::string slcfile) {

    // now open the slicefile and write out the data
    std::ofstream slicefile;
    slicefile.open(slcfile);
    slicefile.setf(std::ios::scientific);
    slicefile.precision(12);
    const auto w = 24;

    // write the header
    slicefile << std::setw(w) << "x";
    for (auto it=varNames.begin(); it!=varNames.end(); ++it)
        slicefile << std::setw(w) << *it;

    slicefile << std::endl;

    // write the data in columns
    const auto SMALL = 1.e-20;
    for (auto i = 0; i < nbins; i++) {

        slicefile << std::setw(w) << r[i];

        for (auto it=vars.begin(); it!=vars.end(); ++it) {
            if(std::abs((*it)[i]) < SMALL) (*it)[i] = 0.0;
            slicefile << std::setw(w) << (*it)[i];
        }

        slicefile << std::endl;

    }

    slicefile.close();
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
            << "\n\n" << std::endl;

}

///
/// Gets the variable ``varname`` from the ``job_info`` file and returns as a
/// string
///
string GetVarFromJobInfo (const string pltfile, const string varname) {
    string filename = pltfile + "/job_info";
    std::regex re("(?:[ \\t]*)" + varname + "\\s*:\\s*(.*)\\s*\\n");

    std::smatch m;

    std::ifstream jobfile(filename);
    if (jobfile.is_open()) {
        std::stringstream buf;
        buf << jobfile.rdbuf();
        string file_contents = buf.str();

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
Vector<Real> GetCenter (const string pltfile) {
    auto center_str = GetVarFromJobInfo(pltfile, "center");

    // split string
    std::istringstream iss {center_str};
    Vector<Real> center;

    std::string s;
    while (std::getline(iss, s, ','))
        center.push_back(stod(s));

    return center;
}

#endif
