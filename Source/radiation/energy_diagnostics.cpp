#include <iomanip>

#include <cstdio>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "Castro.H"
#include "Radiation.H"
#include "Castro_F.H"

using namespace amrex;

// RHOYLTEST just prints out extra diagnostics to check the Yl implementation.

#define RHOYLTEST

void
Castro::do_energy_diagnostics()
{
  int finest_level = parent->finestLevel();

  int v = verbose;
  v = (do_radiation && radiation->verbose > v) ? radiation->verbose : v;

  if (do_radiation && (v > 2 || (level == 0 && v > 0))) {
    Real prev_time = state[State_Type].prevTime();
    Real dt = parent->dtLevel(level);

    Real m = 0.0, s = 0.0, r = 0.0, rr = 0.0, y = 0.0, ry = 0.0, rry = 0.0;
    Real r_nu[3] = {0.0};
#ifdef RHOYLTEST
    Real yl = 0.0;
#endif
    Real rhomax = 0.0;

    for (int lev = 0; lev <= finest_level; lev++) {
      m += getLevel(lev).volWgtSum("density", prev_time + dt);
      s += getLevel(lev).volWgtSum("rho_E",   prev_time + dt);
      if (Radiation::nNeutrinoSpecies == 0 ||
          Radiation::nNeutrinoGroups[0] == 0) {
        // This is either a regular photon problem, or a photon
        // problem being run through the neutrino solver
        if (!Radiation::do_multigroup) {
          r += getLevel(lev).volWgtSum("rad", prev_time + dt);
        }
        else {
          char rad_name[10];
          for (int igroup = 0; igroup < Radiation::nGroups; igroup++) {
            sprintf(rad_name, "rad%d", igroup);
            r += getLevel(lev).volWgtSum(rad_name, prev_time + dt);
          }
        }
        if (lev < finest_level) {
          // If using deferred sync, also include flux register energy
          FluxRegister* sync_flux = radiation->consRegister(lev + 1);
          if (sync_flux) {
            for (int k = 0; k < Radiation::nGroups; k++) {
              // Minus sign in following relates to FluxRegister defn:
              Real tmp = -sync_flux->SumReg(k) * parent->dtLevel(lev);
              rr  += tmp;
            }
          }
        }
      }
      else {
        // This is a real neutrino problem
	MultiFab & S_new = getLevel(lev).get_new_data(State_Type);
	Real levelrhomax = S_new.max(URHO);
	rhomax = (rhomax > levelrhomax) ? rhomax : levelrhomax;
#ifdef RHOYLTEST
        yl += getLevel(lev).volWgtSum("rho_Yl", prev_time + dt);
#endif
        y  += getLevel(lev).volWgtSum("rho_Ye", prev_time + dt);
        for (int j = 0, k = 0; j < Radiation::nNeutrinoSpecies; j++) {
          int species_sign = 0;
          if (j == 0) {
            species_sign =  1;
          }
          else if (j == 1) {
            species_sign = -1;
          }
          char rad_name[10];
          for (int i = 0; i < Radiation::nNeutrinoGroups[j]; i++, k++) {
            sprintf(rad_name, "rads%dg%d", j, i);
            Real tmp = getLevel(lev).volWgtSum(rad_name, prev_time + dt);
            r  += tmp;
            ry += tmp * species_sign / radiation->group_center(k);
	    r_nu[j] += tmp;
          }
        }
        if (lev < finest_level) {
          FluxRegister* sync_flux = radiation->consRegister(lev + 1);
          if (sync_flux) {
            for (int j=0, k=0; j < Radiation::nNeutrinoSpecies; j++) {
              int species_sign = 0;
              if (j == 0) {
                species_sign =  1;
              }
              else if (j == 1) {
                species_sign = -1;
              }
              for (int i=0; i < Radiation::nNeutrinoGroups[j]; i++,k++) {
                // minus sign in following relates to FluxRegister defn:
                Real tmp = -sync_flux->SumReg(k) * parent->dtLevel(lev);
                rr  += tmp;
                rry += tmp * species_sign / radiation->group_center(k);
              }
            }
          }
        }
      }
    }

    // Units conversions:

    if (Radiation::nNeutrinoSpecies > 0) {
      // convert neutrino zeroth moment into energy units
      // (using c instead of clight so that we can run certain contrived
      // single-group tests and get the right energy diagnostic)
      r   *= Radiation::radtoE; 
      r_nu[0] *= Radiation::radtoE; 
      r_nu[1] *= Radiation::radtoE; 
      r_nu[2] *= Radiation::radtoE; 
      // flux register mult by additional (c * dt) to get to energy units
      rr  *= Radiation::radfluxtoF;
      if (Radiation::nNeutrinoGroups[0] > 0) {
        // rho Y_e to lepton number
        y   *= Radiation::Avogadro;
        // same conversion as for r, also group center to group energy
        ry  *=  Radiation::radtoE /  Radiation::hPlanck;
        // flux register mult by additional (c * dt)
        rry *= Radiation::radfluxtoF / Radiation::hPlanck;
      }
    }

    // Geometric corrections:

#if (BL_SPACEDIM == 1)
    if (Geom().IsSPHERICAL()) {
      // Internal to the radiation class, fluxes in flux registers have
      // been weighted by r^2 in spherical coords.  We multiply by an
      // additional 4*pi here because we want to weight them as covering
      // the full surface of a spherical shell.  This is compatible
      // with what volWgtSum does for volumes.
      rr  *= (4.0 * M_PI);
      rry *= (4.0 * M_PI);
    }
#elif (BL_SPACEDIM == 2)
    if (Geom().IsRZ()) {
      // Internal to the radiation class, fluxes in flux registers have
      // been weighted by r in RZ coordinates.  We multiply by an
      // additional 2*pi here to be compatible with what volWgtSum
      // does for volumes.
      rr  *= (2.0 * M_PI);
      rry *= (2.0 * M_PI);
    }
#endif

    if (ParallelDescriptor::IOProcessor()) {
      int oldprec = cout.precision(20);
      cout << "Integrated  Fluid   Mass  is " << m << '\n';
      cout << "Integrated  Fluid  Energy is " << s << '\n';
      cout << "Integrated Radiant Energy is " << r << '\n';
      cout << "     Flux Register Energy is " << rr << '\n';
      const char* nunames[3] = {"    Electron Neutrinos", 
				"Electron Antineutrinos" , 
				"Muon and Tau Neutrinos"};
      if (Radiation::nNeutrinoSpecies > 0 &&
          Radiation::nNeutrinoGroups[0] > 0) {
	for (int j=0; j<Radiation::nNeutrinoSpecies; j++) {
	  cout << "     " << nunames[j] << ": " << r_nu[j]<<endl;
	}
      }
      cout << "Integrated  Total  Energy is " << s + r + rr << endl;
      if (Radiation::nNeutrinoSpecies > 0 &&
          Radiation::nNeutrinoGroups[0] > 0) {
        cout << "Integrated  Fluid  Lepton Number is " << y << '\n';
        cout << "Integrated Radiant Lepton Number is " << ry << '\n';
        cout << "     Flux Register Lepton Number is " << rry << '\n';
        cout << "Integrated  Total  Lepton Number is " << y+ry+rry << endl;
#ifdef RHOYLTEST
        cout << "Integrated  Density * Yl         is " << yl << endl;
        cout << "Integrated  Density * Yl (Check) is "
             << (y + ry) / Radiation::Avogadro << endl;
#endif
	cout << "The maximum density is " << rhomax << endl;
      }
      cout.precision(oldprec);
    }
  }
  else if (v > 2 || (level == 0 && v > 0)) {
    Real prev_time = state[State_Type].prevTime();
    Real dt = parent->dtLevel(level);

    Real m = 0.0, s = 0.0;
    for (int lev = 0; lev <= finest_level; lev++) {
      m += getLevel(lev).volWgtSum("density", prev_time + dt);
      s += getLevel(lev).volWgtSum("rho_E",   prev_time + dt);
    }
    if (ParallelDescriptor::IOProcessor()) {
      int oldprec = cout.precision(20);
      cout << "Integrated  Fluid   Mass  is " << m << '\n';
      cout << "Integrated  Fluid  Energy is " << s << '\n';
      cout.precision(oldprec);
    }
  }
}
