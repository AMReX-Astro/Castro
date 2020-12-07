#ifndef wdmerger_util_H
#define wdmerger_util_H

#include <AMReX.H>

#include <math.H>
#include <fundamental_constants.H>

#include <Rotation.H>

#include <initial_model.H>

void freefall_velocity (Real mass, Real distance, Real& vel);

void kepler_third_law (Real radius_1, Real mass_1, Real radius_2, Real mass_2,
                       Real& period, Real eccentricity, Real phi, Real& a,
                       Real& r_1, Real& r_2, Real& v_1r, Real& v_2r, Real& v_1p, Real& v_2p);

void set_wd_composition (Real mass, Real& envelope_mass, Real core_comp[NumSpec], Real envelope_comp[NumSpec]);

void ensure_primary_mass_larger ();

void initialize_problem ();

void finalize_probdata ();

void binary_setup ();

void set_small ();

void set_star_data ();

void update_center (Real time);

#endif
