#ifndef RADIATION
#define RADIATION
#endif

#include "Radiation.H"

#include "RAD_F.H"
#include "LHH.H"

#include <AMReX_ParmParse.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

void Radiation::get_groups(int verbose)
{
  group_print_factor   = 1.0;
  group_units          = " (units are Hz)";

  Array<Real> dlognugroup;

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
  else if (SolverType == MGFLDSolver && radiation_type == Neutrino) {

    xnu.resize(nGroups+3, 0.0); 
    nugroup.resize(nGroups, 0.0);
    dnugroup.resize(nGroups, 0.0);
    dlognugroup.resize(nGroups, 0.0);

    Array<Real> lowest, highest;
    ParmParse pp("radiation");
    pp.getarr( "lowestGroupMeV",  lowest, 0, nNeutrinoSpecies);
    pp.getarr("highestGroupMeV", highest, 0, nNeutrinoSpecies);

    Real hPlanckMeV = hPlanck / convert_MeV_erg; // Planck in (MeV sec)

    group_print_factor = hPlanckMeV;
    group_units        = " (units are MeV)";

    int ibase = 0;
    for (int ispec = 0; ispec < nNeutrinoSpecies; ispec++) {
      if (nNeutrinoGroups[ispec] > 0) {
        BL_ASSERT(nNeutrinoGroups[ispec] > 1); // no "gray" species
        Real lowlog =  log( lowest[ispec]);
        Real faclog = ((log(highest[ispec]) - lowlog)
                       / (nNeutrinoGroups[ispec] - 1));
        for (int igroup = 0; igroup < nNeutrinoGroups[ispec]; igroup++) {
          int i = ibase + igroup;
          nugroup[i] = exp(lowlog + igroup * faclog) / hPlanckMeV;
        }
      }
      ibase += nNeutrinoGroups[ispec];
    }

    // set up xnu & dnu & dlognu
    int ibase_nu = 0;
    int ibase_xnu = 0;
    for (int ispec = 0; ispec < nNeutrinoSpecies; ispec++) {
      if (nNeutrinoGroups[ispec] > 0) {
	int igroup, inu, ixnu;
        for (igroup = 1; igroup < nNeutrinoGroups[ispec]; igroup++) {
          inu = ibase_nu + igroup;
	  ixnu = ibase_xnu + igroup;
	  xnu[ixnu] = sqrt(nugroup[inu-1] * nugroup[inu]);
        }
	igroup = 0;
	inu = ibase_nu + igroup;
	ixnu = ibase_xnu + igroup;
	xnu[ixnu] = nugroup[inu]*nugroup[inu]/xnu[ixnu+1];
	igroup = nNeutrinoGroups[ispec];
	inu = ibase_nu + igroup;
	ixnu = ibase_xnu + igroup;
	xnu[ixnu] = nugroup[inu-1]*nugroup[inu-1]/xnu[ixnu-1];

	for (igroup=0; igroup < nNeutrinoGroups[ispec]; igroup++) {
          inu = ibase_nu + igroup;
	  ixnu = ibase_xnu + igroup;
	  dnugroup[inu] = xnu[ixnu+1] - xnu[ixnu]; 
	  dlognugroup[inu] = log(xnu[ixnu+1]) - log(xnu[ixnu]); 
	}
      }
      ibase_nu  += nNeutrinoGroups[ispec];
      ibase_xnu += nNeutrinoGroups[ispec]+1;
    }
    
  }
  else if (SolverType == MGFLDSolver) {

    // xnu is being used by MGFLDSolver

    xnu.resize(nGroups+3, 0.0);   // large enough for three neutrino species
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

    if (nNeutrinoSpecies == 0) {
      for (int i = 0; i < nGroups; i++) {
        nugroup[i] = sqrt(xnu[i] * xnu[i+1]);
      }
    }
    else {
      int ibase = 0;
      for (int ispec = 0; ispec < nNeutrinoSpecies; ispec++) {
        for (int igroup = 0; igroup < nNeutrinoGroups[ispec]; igroup++) {
          int i = ibase + igroup;
          nugroup[i] = 0.5 * (xnu[igroup] + xnu[igroup+1]);
        }
        ibase += nNeutrinoGroups[ispec];
      }
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
    else if (nNeutrinoSpecies == 0) {
      int i = 0;
      dnugroup[i] = 0.5 * (nugroup[i+1] - nugroup[i]);
      for (i = 1; i < nGroups - 1; i++) {
        dnugroup[i] = 0.5 * (nugroup[i+1] - nugroup[i-1]);
      }
      i = nGroups - 1;
      dnugroup[i] = 0.5 * (nugroup[i] - nugroup[i-1]);
    }
    else {
      int ibase = 0;
      for (int ispec = 0; ispec < nNeutrinoSpecies; ispec++) {
        if (nNeutrinoGroups[ispec] == 0) {
          // no groups for this species, do nothing
        }
        else if (nNeutrinoGroups[ispec] == 1) {
          // should we support "gray" species?
          std::cout << "Species with single energy group not supported" << std::endl;
          exit(1);
        }
        else {
          int i = ibase;
          dnugroup[i] = 0.5 * (nugroup[i+1] - nugroup[i]);
          for (i = ibase + 1; i < ibase + nNeutrinoGroups[ispec] - 1; i++) {
            dnugroup[i] = 0.5 * (nugroup[i+1] - nugroup[i-1]);
          }
          i = ibase + nNeutrinoGroups[ispec] - 1;
          dnugroup[i] = 0.5 * (nugroup[i] - nugroup[i-1]);
        }

        ibase += nNeutrinoGroups[ispec];
      }
    }
  }

  int nG0 = 0, nG1 = 0;
  if (nNeutrinoSpecies >= 2) {
    nG0 = nNeutrinoGroups[0];
    nG1 = nNeutrinoGroups[1];
  }
  else if (nNeutrinoSpecies == 1) {
    nG0 = nNeutrinoGroups[0];
    nG1 = 0;
  }

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
    if (nNeutrinoSpecies > 0) {
      groupfile << "# " << nNeutrinoSpecies
                << " neutrino species, numbers of groups are: ";
      for (int n = 0; n < nNeutrinoSpecies; n++) {
        if (n > 0) {
          groupfile << ", ";
        }
        groupfile << nNeutrinoGroups[n];
      }
      groupfile << std::endl;
    }
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
  if (nNeutrinoSpecies > 0) {
    os << "# " << nNeutrinoSpecies
       << " neutrino species, numbers of groups are: ";
    for (int n = 0; n < nNeutrinoSpecies; n++) {
      if (n > 0) {
        os << ", ";
      }
      os << nNeutrinoGroups[n];
    }
    os << std::endl;
  }

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
