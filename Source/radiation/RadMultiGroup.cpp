#ifndef RADIATION
#define RADIATION
#endif

#include <Radiation.H>

#include <RAD_F.H>

#include <AMReX_ParmParse.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

void Radiation::get_groups(int verbose)
{
  group_print_factor   = 1.0;
  group_units          = " (units are Hz)";

  Vector<Real> dlognugroup;

  xnu.resize(nGroups+1, 0.0); // Bounds of the frequency group
  nugroup.resize(nGroups, 1.0); // Geometric center of the frequency group
  dnugroup.resize(nGroups, 0.0); // Width of the frequency group
  dlognugroup.resize(nGroups, 0.0); // Log of the width of the frequency group

  ParmParse pp("radiation");

  if (nGroups > 1) {
    Real lowest;
    pp.get("lowestGroupHz", lowest);

    Real groupGrowFactor;
    if (pp.query("groupGrowFactor", groupGrowFactor)) {
      xnu[0] = lowest;

      Real firstGroupWidthHz;
      pp.get("firstGroupWidthHz", firstGroupWidthHz);
      dnugroup[0] = firstGroupWidthHz;

      xnu[1] = xnu[0] + dnugroup[0];

      if (lowest == 0.0) {
        nugroup[0] = 0.5*dnugroup[0];
        dlognugroup[0] = 2.0 * (log(xnu[1]) - log(nugroup[0]));
      }
      else {
        nugroup[0] = sqrt(xnu[0]*xnu[1]);
        dlognugroup[0] = log(xnu[1]) - log(xnu[0]);
      }

      for (int i=1; i<nGroups; i++) {
        dnugroup[i] = dnugroup[i-1] * groupGrowFactor;
        xnu[i+1] = xnu[i] + dnugroup[i];
        nugroup[i] = sqrt(xnu[i]*xnu[i+1]);
        dlognugroup[i] = log(xnu[i+1]) - log(xnu[i]);
      }
    }
    else {
      Real highest;
      pp.get("highestGroupHz", highest);

      Real loglowest = log10(lowest);
      Real loghighest = log10(highest);
      Real dlognu = (loghighest - loglowest) / Real(nGroups);

      for (int i=0; i<nGroups; i++) {
        xnu[i] = pow(10.0, loglowest+i*dlognu);
        nugroup[i] = pow(10.0, loglowest+(i+0.5)*dlognu);
      }
      xnu[nGroups] = highest;

      for (int i=0; i<nGroups; i++) {
        dnugroup[i] = xnu[i+1] - xnu[i];
        dlognugroup[i] = log(xnu[i+1]) - log(xnu[i]);
      }
    }
  }

  int nG0 = 0, nG1 = 0;

  if (SolverType == MGFLDSolver) { 
    BL_FORT_PROC_CALL(CA_INITGROUPS3,ca_initgroups3)
      (nugroup.dataPtr(), dnugroup.dataPtr(), dlognugroup.dataPtr(), xnu.dataPtr(), 
       nGroups, nG0, nG1);
  }
  else if (xnu.size() > 0) {
    BL_FORT_PROC_CALL(CA_INITGROUPS2,ca_initgroups2)
      (nugroup.dataPtr(), dnugroup.dataPtr(), xnu.dataPtr(), nGroups);
  }
  else {
    BL_FORT_PROC_CALL(CA_INITGROUPS,ca_initgroups)
      (nugroup.dataPtr(), dnugroup.dataPtr(), nGroups, nG0, nG1);
  }

  if (ParallelDescriptor::IOProcessor()) {

    std::ofstream groupfile;
    groupfile.open("group_structure.dat");

    groupfile << "# total number of groups = " << nGroups << std::endl;
    groupfile << "# group center, group weight" << group_units << std::endl;

    groupfile.precision(10);
    for (int i = 0; i < nGroups; i++) {
      groupfile.width(15);
      groupfile <<  nugroup[i] * group_print_factor << "\t";
      groupfile.width(15);
      groupfile << dnugroup[i] * group_print_factor << std::endl;
    }
    groupfile << std::endl;

    if (xnu.size() > 0) {
      groupfile << "# group lower boundaries" << std::endl;
      for (int i = 0; i < xnu.size(); i++) {
        groupfile << "group(" << i << ") = "
                  << xnu[i] * group_print_factor << std::endl;
      }
    }

    if (dlognugroup.size() > 0) {
      groupfile << std::endl;
      groupfile << "# group width in log space" << std::endl;
      for (int i = 0; i < dlognugroup.size(); i++) {
        groupfile << "group(" << i << ") = "
                  << dlognugroup[i] << std::endl;
      }      
    }

    groupfile.close();

    if (verbose >= 1) {
      write_groups(std::cout);
    }
  }
}

void Radiation::write_groups(ostream& os)
{
  os << "# total number of groups = " << nGroups << std::endl;

  if (nGroups > 1) {
    os << "# group center, group weight" << group_units << std::endl;
    int oldprec = os.precision(10);
    for (int i = 0; i < nGroups; i++) {
      os.width(3);
      os << i << ": ";
      os.width(15);
      os <<  nugroup[i] * group_print_factor << ", ";
      os.width(15);
      os << dnugroup[i] * group_print_factor << std::endl;
    }
    os.precision(oldprec);
  }

  if (xnu.size() > 0) {
    os << "# group lower boundaries" << std::endl;
    for (int i = 0; i < xnu.size(); i++) {
      os << "group(" << i << ") = "
         << xnu[i] * group_print_factor << std::endl;
    }
  }
}
