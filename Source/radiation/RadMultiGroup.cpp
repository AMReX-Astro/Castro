#ifndef RADIATION
#define RADIATION
#endif

#include "Radiation.H"

#include "RAD_F.H"

#include "RadTests.H"

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

  if (RadTests::do_rad_sphere) {

    xnu.resize(nGroups+1, 0.0);  

    // fundamental constants in CGS
    Real ev2erg = 1.e-6*convert_MeV_erg;
    
    Real E_min = 0.5e0*ev2erg;
    Real E_max = 306e3*ev2erg;
    
    // alpha is the geometric spacing between groups
    Real alpha = pow((E_max/E_min), 1.e0/static_cast<Real>(nGroups));
    
    xnu[0] = E_min/hPlanck;
    
    for(int i = 1; i <= nGroups; i++) {
        xnu[i] = alpha*xnu[i-1];
    }
  }
  else if (SolverType == MGFLDSolver) {

    // xnu is being used by MGFLDSolver

    xnu.resize(nGroups+1, 0.0);
    nugroup.resize(nGroups, 1.0);
    dnugroup.resize(nGroups, 0.0);
    dlognugroup.resize(nGroups, 0.0);

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
  }
  else {
    std::cout << "What groups to use?" << std::endl;
    exit(1);
  }

  if (nugroup.size() == 0) {  // if not initialized above
    nugroup.resize(nGroups, 0.0);

    for (int i = 0; i < nGroups; i++) {
        nugroup[i] = sqrt(xnu[i] * xnu[i+1]);
    }
  }

  if (dnugroup.size() == 0) {  // if not initialized above
    dnugroup.resize(nGroups, 0.0);
    dlognugroup.resize(nGroups, 0.0);

    if (SolverType == MGFLDSolver && xnu.size() > 0) {
      for (int i=0; i<nGroups; i++) {
	dnugroup[i] = xnu[i+1] - xnu[i];
	dlognugroup[i] = log(xnu[i+1]) - log(xnu[i]);
      }
    }
    else {
      int i = 0;
      dnugroup[i] = 0.5 * (nugroup[i+1] - nugroup[i]);
      for (i = 1; i < nGroups - 1; i++) {
        dnugroup[i] = 0.5 * (nugroup[i+1] - nugroup[i-1]);
      }
      i = nGroups - 1;
      dnugroup[i] = 0.5 * (nugroup[i] - nugroup[i-1]);
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
