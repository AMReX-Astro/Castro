#include <iomanip>

#include <cstdio>
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <Castro.H>
#include <Radiation.H>

using namespace amrex;

void
Castro::do_energy_diagnostics()
{
  int finest_level = parent->finestLevel();

  int v = verbose;
  v = (do_radiation && radiation->verbose > v) ? radiation->verbose : v;

  if (do_radiation && (v > 2 || (level == 0 && v > 0))) {
    Real prev_time = state[State_Type].prevTime();
    Real dt = parent->dtLevel(level);

    Real m = 0.0, s = 0.0, r = 0.0, rr = 0.0, rry = 0.0;

    for (int lev = 0; lev <= finest_level; lev++) {
      m += getLevel(lev).volWgtSum("density", prev_time + dt);
      s += getLevel(lev).volWgtSum("rho_E",   prev_time + dt);
      if (!Radiation::do_multigroup) {
          r += getLevel(lev).volWgtSum("rad", prev_time + dt);
      }
      else {
          std::string rad_name;
          for (int igroup = 0; igroup < Radiation::nGroups; igroup++) {
              rad_name = "rad" + std::to_string(igroup);
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

    // Geometric corrections:

#if (AMREX_SPACEDIM == 1)
    if (Geom().IsSPHERICAL()) {
      // Internal to the radiation class, fluxes in flux registers have
      // been weighted by r^2 in spherical coords.  We multiply by an
      // additional 4*pi here because we want to weight them as covering
      // the full surface of a spherical shell.  This is compatible
      // with what volWgtSum does for volumes.
      rr  *= (4.0 * M_PI);
      rry *= (4.0 * M_PI);
    }
#elif (AMREX_SPACEDIM == 2)
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
      cout << "Integrated  Total  Energy is " << s + r + rr << endl;
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
