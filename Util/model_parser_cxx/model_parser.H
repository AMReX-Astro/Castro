#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <network.H>
#include <model_parser_data.H>
#include <AMReX_Print.H>
#include <castro_params.H>
#include <eos.H>
#include <ambient.H>

using namespace amrex;

///
/// read in an initial model and return arrays with the model data.
/// take care to match the species available in the model file to
/// those defined by the network
///
/// the model file is assumed to be of the follow form:
/// # npts = 896
/// # num of variables = 6
/// # density
/// # temperature
/// # pressure
/// # carbon-12
/// # oxygen-16
/// # magnesium-24
/// 195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
/// 585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0
///
/// we read in the number of variables and their order and use this to
/// map them into the model_state array.  We ignore anything other than
/// density, temperature, pressure and composition.
///
/// composition is assumed to be in terms of mass fractions

// remove whitespace -- from stackoverflow


namespace model_string
{
    inline std::string& ltrim(std::string& s)
    {
	auto it = std::find_if(s.begin(), s.end(),
                               [](int c) {
                                   return !std::isspace(c);
                               });
	s.erase(s.begin(), it);
	return s;
    }

    inline std::string& rtrim(std::string& s)
    {
	auto it = std::find_if(s.rbegin(), s.rend(),
                               [](int c) {
                                   return !std::isspace(c);
                               });
	s.erase(it.base(), s.end());
	return s;
    }
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
int
locate(const Real r, const int model_index) {

    int loc;

    if (r <= model::profile(model_index).r(0)) {
       loc = 0;

    } else if (r > model::profile(model_index).r(model::npts-2)) {
       loc = model::npts-1;

    } else {

        int ilo = 0;
        int ihi = model::npts-2;

        while (ilo+1 != ihi) {
            int imid = (ilo + ihi) / 2;

            if (r <= model::profile(model_index).r(imid)) {
                ihi = imid;
            } else {
                ilo = imid;
            }
        }

        loc = ihi;
    }

    return loc;
}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
Real
interpolate(const Real r, const int var_index, const int model_index=0) {

    // find the value of model_state component var_index at point r
    // using linear interpolation.  Eventually, we can do something
    // fancier here.

    int id = locate(r, model_index);

    Real slope;
    Real interp;

    if (id == 0) {

       slope = (model::profile(model_index).state(id+1, var_index) -
                model::profile(model_index).state(id, var_index)) /
           (model::profile(model_index).r(id+1) - model::profile(model_index).r(id));
       interp = slope * (r - model::profile(model_index).r(id)) + model::profile(model_index).state(id, var_index);

       // safety check to make sure interp lies within the bounding points
       Real minvar = amrex::min(model::profile(model_index).state(id+1, var_index),
                                model::profile(model_index).state(id, var_index));
       Real maxvar = amrex::max(model::profile(model_index).state(id+1, var_index),
                                model::profile(model_index).state(id, var_index));
       interp = amrex::max(interp, minvar);
       interp = amrex::min(interp, maxvar);

    } else if (id == model::npts-1) {

       slope = (model::profile(model_index).state(id, var_index) -
                model::profile(model_index).state(id-1, var_index)) /
           (model::profile(model_index).r(id) - model::profile(model_index).r(id-1));
       interp = slope * (r - model::profile(model_index).r(id)) + model::profile(model_index).state(id, var_index);


       // safety check to make sure interp lies within the bounding points
       Real minvar = amrex::min(model::profile(model_index).state(id-1, var_index),
                                model::profile(model_index).state(id, var_index));
       Real maxvar = amrex::max(model::profile(model_index).state(id-1, var_index),
                                model::profile(model_index).state(id, var_index));
       interp = amrex::max(interp, minvar);
       interp = amrex::min(interp, maxvar);

    } else {

        if (r >= model::profile(model_index).r(id)) {

            slope = (model::profile(model_index).state(id+1, var_index) -
                     model::profile(model_index).state(id, var_index)) /
                (model::profile(model_index).r(id+1) - model::profile(model_index).r(id));
            interp = slope * (r - model::profile(model_index).r(id)) + model::profile(model_index).state(id, var_index);

        } else {

            slope = (model::profile(model_index).state(id, var_index) -
                     model::profile(model_index).state(id-1, var_index)) /
                (model::profile(model_index).r(id) - model::profile(model_index).r(id-1));
            interp = slope * (r - model::profile(model_index).r(id)) + model::profile(model_index).state(id, var_index);

        }

    }

    return interp;

}

// Subsample the interpolation to get an averaged profile. For this we need to know the
// 3D coordinate (relative to the model center) and cell size.

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real interpolate_3d (const Real* loc, const Real* dx, int var_index, int nsub = 1, int model_index = 0)
{
    // We perform a sub-grid-scale interpolation, where
    // nsub determines the number of intervals we split the zone into.
    // If nsub = 1, we simply interpolate using the cell-center location.
    // As an example, if nsub = 3, then the sampled locations will be
    // k = 0 --> z = loc(3) - dx(3) / 3   (1/6 of way from left edge of zone)
    // k = 1 --> z = loc(3)               (halfway between left and right edge)
    // k = 2 --> z = loc(3) + dx(3) / 3   (1/6 of way from right edge of zone)

    Real interp = 0.0_rt;

    for (int k = 0; k < nsub; ++k) {
        Real z = loc[2] + (static_cast<Real>(k) + 0.5_rt * (1 - nsub)) * dx[2] / nsub;

        for (int j = 0; j < nsub; ++j) {
            Real y = loc[1] + (static_cast<Real>(j) + 0.5_rt * (1 - nsub)) * dx[1] / nsub;

            for (int i = 0; i < nsub; ++i) {
                Real x = loc[0] + (static_cast<Real>(i) + 0.5_rt * (1 - nsub)) * dx[0] / nsub;

                Real dist = std::sqrt(x * x + y * y + z * z);

                interp += interpolate(dist, var_index, model_index);
            }
        }
    }

    // Now normalize by the number of intervals.

    interp /= (nsub * nsub * nsub);

    return interp;
}

// Establish an isothermal initial model. The constraints are:
// dx: the spacing of the points
// temperature: uniform stellar temperature
// core_comp: uniform composition of the core
// envelope_mass: the mass of the envelope with composition core_comp (can be set to zero to ignore)
//
// The user can pass in either a target mass_want or central_density_want, and should leave the
// other one negative. Both mass_want and central_density_want will be updated with the final mass
// and central density of the model. The radius will also be returned, which is the distance at which
// we switch from the initial model to ambient material.

AMREX_INLINE
void establish_hse (Real& mass_want, Real& central_density_want, Real& radius,
                    const Real core_comp[NumSpec], Real temperature, Real dx,
                    Real envelope_mass, const Real envelope_comp[NumSpec], int model_index = 0)
{
    // Note that if central_density > 0, then this initial model generator will use it in calculating
    // the model. If mass is also provided in this case, we assume it is an estimate used for the purpose of
    // determining the envelope mass boundary.

    // Check to make sure we've specified at least one of them.

    if (mass_want < 0.0_rt && central_density_want < 0.0_rt) {
        amrex::Error("Error: Must specify either mass or central density in the initial model generator.");
    }

    // If we are specifying the mass, then we don't know what WD central density
    // will give the desired total mass, so we need to do a secant iteration
    // over central density. rho_c_old is the 'old' guess for the central
    // density and rho_c is the current guess.  After two loops, we can
    // start estimating the density required to yield our desired mass.

    // If instead we are specifying the central density, then we only need to do a
    // single HSE integration.

    const int max_hse_iter = 250;
    int max_mass_iter;

    Real rho_c, rho_c_old, drho_c;
    Real mass, mass_old;
    Real p_want, drho;

    if (central_density_want > 0.0_rt) {

        max_mass_iter = 1;

        rho_c_old = central_density_want;
        rho_c     = central_density_want;

    }
    else {

        max_mass_iter = max_hse_iter;

        rho_c_old = -1.0_rt;
        rho_c     = 1.e7_rt;     // A reasonable starting guess for moderate-mass WDs

    }

    // Check to make sure the initial temperature makes sense.

    if (temperature < castro::small_temp) {
        amrex::Error("Error: WD central temperature is less than small_temp. Aborting.");
    }

    bool mass_converged = false;

    model::initial_model_t& model = model::profile(model_index);

    for (int mass_iter = 1; mass_iter <= max_mass_iter; ++mass_iter) {

        bool fluff = false;

        // We start at the center of the WD and integrate outward.  Initialize
        // the central conditions.

        model.state(0, model::itemp) = temperature;
        model.state(0, model::idens) = rho_c;
        for (int n = 0; n < NumSpec; ++n) {
            model.state(0, model::ispec + n) = core_comp[n];
        }

        eos_t eos_state;
        eos_state.rho  = model.state(0, model::idens);
        eos_state.T    = model.state(0, model::itemp);
        for (int n = 0; n < NumSpec; ++n) {
            eos_state.xn[n] = model.state(0, model::ispec + n);
        }

        eos(eos_input_rt, eos_state);

        model.state(0, model::ipres) = eos_state.p;

        model.r(0) = 0.5 * dx;

        int icutoff = NPTS_MODEL;

        // Make the initial guess be completely uniform.

        for (int i = 1; i < NPTS_MODEL; ++i) {
            model.state(i, model::idens) = model.state(0, model::idens);
            model.state(i, model::itemp) = model.state(0, model::itemp);
            model.state(i, model::ipres) = model.state(0, model::ipres);
            for (int n = 0; n < NumSpec; ++n) {
                model.state(i, model::ispec + n) = model.state(0, model::ispec + n);
            }
            model.r(i) = model.r(i-1) + dx;
        }

        // Keep track of the mass enclosed below the current zone.

        Real rl = 0.0_rt;
        Real rr = rl + dx;

        Real M_enclosed = (4.0_rt / 3.0_rt) * M_PI * (std::pow(rr, 3) - std::pow(rl, 3)) * model.state(0, model::idens);
        mass = M_enclosed;

        //-------------------------------------------------------------------------
        // HSE solve
        //-------------------------------------------------------------------------
        for (int i = 1; i < NPTS_MODEL; ++i) {

            rl += dx;
            rr += dx;

            // As the initial guess for the density, use the underlying zone.

            model.state(i, model::idens) = model.state(i-1, model::idens);

            if (mass_want > 0.0_rt && M_enclosed >= mass_want - envelope_mass) {
                for (int n = 0; n < NumSpec; ++n) {
                    model.state(i, model::ispec + n) = envelope_comp[n];
                    eos_state.xn[n] = model.state(i, model::ispec + n);
                }
            }
            else {
                for (int n = 0; n < NumSpec; ++n) {
                    model.state(i, model::ispec + n) = core_comp[n];
                    eos_state.xn[n] = model.state(i, model::ispec + n);
                }
            }

            Real g = -C::Gconst * M_enclosed / (std::pow(rl, 2));


            //----------------------------------------------------------------------
            // Iteration loop
            //----------------------------------------------------------------------

            // Start off the Newton loop by assuming that the zone has not converged.

            bool converged_hse = false;

            for (int hse_iter = 1; hse_iter <= max_hse_iter; ++hse_iter) {

                if (fluff) {
                    model.state(i, model::idens) = ambient::ambient_state[URHO];
                    eos_state.rho = ambient::ambient_state[URHO];
                    break;
                }

                // The core is isothermal, so we just need to constrain
                // the density and pressure to agree with the EOS and HSE.

                // We difference HSE about the interface between the current
                // zone and the one just inside.

                Real rho_avg = 0.5_rt * (model.state(i, model::idens) + model.state(i-1, model::idens));
                p_want = model.state(i-1, model::ipres) + dx * rho_avg * g;

                eos(eos_input_rt, eos_state);

                drho = (p_want - eos_state.p) / (eos_state.dpdr - 0.5_rt * dx * g);

                model.state(i, model::idens) = amrex::max(0.9_rt * model.state(i, model::idens),
                                                     amrex::min(model.state(i, model::idens) + drho,
                                                                1.1_rt * model.state(i, model::idens)));
                eos_state.rho = model.state(i, model::idens);

                if (model.state(i, model::idens) < ambient::ambient_state[URHO]) {
                    icutoff = i;
                    fluff = true;
                }

                if (std::abs(drho) < model::hse_tol * model.state(i, model::idens)) {
                    converged_hse = true;
                    break;
                }

            }

            if (!converged_hse && (!fluff)) {

                std::cout << "Error: zone " <<  i << " did not converge in init_hse()" << std::endl;
                std::cout << model.state(i, model::idens) << " " << model.state(i, model::itemp) << std::endl;
                std::cout << p_want << " " << eos_state.p;
                std::cout << drho << " " << model::hse_tol * model.state(i, model::idens);
                amrex::Error("Error: HSE non-convergence.");

            }

            // Call the EOS to establish the final properties of this zone.

            eos(eos_input_rt, eos_state);

            model.state(i, model::ipres) = eos_state.p;

            // Discretize the mass enclosed as (4 pi / 3) * rho * dr * (rl**2 + rl * rr + rr**2).

            Real dM = (4.0_rt / 3.0_rt) * M_PI * model.state(i, model::idens) * dx *
                      (rr * rr + rl * rr + rl * rl);
            M_enclosed += dM;

            if (i <= icutoff) {
                // Also update the final WD mass if we're not in the ambient material.
                mass += dM;
            }

        } // End loop over zones

        radius = model.r(icutoff);

        if (rho_c_old < 0.0_rt) {

            // Not enough iterations yet -- use an arbitrary guess for the next iteration.

            rho_c_old = rho_c;
            rho_c = 0.5_rt * rho_c_old;

        }
        else {

            // Check if we have converged.

            if (std::abs(mass - mass_want) / mass_want < model::mass_tol) {
                mass_converged = true;
                break;
            }

            // Do a secant iteration:
            // M_tot = M(rho_c) + dM/drho |_rho_c x drho + ...

            drho_c = (mass_want - mass) / ((mass  - mass_old) / (rho_c - rho_c_old));

            rho_c_old = rho_c;
            rho_c = amrex::min(1.1e0_rt * rho_c_old, amrex::max((rho_c + drho_c), 0.9e0_rt * rho_c_old));

        }

        mass_old = mass;

    } // End mass constraint loop

    if (!mass_converged && max_mass_iter > 1) {
        amrex::Error("ERROR: WD mass did not converge.");
    }

    central_density_want = model.state(0, model::idens);
    mass_want = mass;

    model::initialized = true;
    model::npts = NPTS_MODEL;
}

AMREX_INLINE
void
read_model_file(std::string& model_file, const int model_index=0) {

    bool found_model, found_dens, found_temp, found_pres, found_velr;
    bool found_spec[NumSpec];
#if NAUX_NET > 0
    bool found_aux[NumAux];
#endif

    // read in the initial model

    std::ifstream initial_model_file;

    initial_model_file.open(model_file, std::ios::in);
    if (!initial_model_file.is_open()) {
        amrex::Error("Error opening the initial model");
    }

    std::string line;

    // first the header line -- this tells us the number of points

    getline(initial_model_file, line);
    std::string npts_string = line.substr(line.find("=")+1, line.length());
    model::npts = std::stoi(npts_string);

    if (model::npts > NPTS_MODEL) {
        amrex::Error("Error: model has more than NPTS_MODEL points,  Increase MAX_NPTS_MODEL");
    }

    // next line tells use the number of variables

    getline(initial_model_file, line);
    std::string num_vars_string = line.substr(line.find("=")+1, line.length());
    int nvars_model_file = std::stoi(num_vars_string);

    // now read in the names of the variables

    std::vector<std::string> varnames_stored;
    for (int n = 0; n < nvars_model_file; n++) {
        getline(initial_model_file, line);
        std::string var_string = line.substr(line.find("#")+1, line.length());
        varnames_stored.push_back(model_string::ltrim(model_string::rtrim(var_string)));
    }

    // allocate storage for the model data

    amrex::Print() << "reading initial model" << std::endl;
    amrex::Print() << model::npts << " points found in the initial model" << std::endl;
    amrex::Print() << nvars_model_file << " variables found in the initial model file" << std::endl;

    // start reading in the data

    amrex::Vector<Real> vars_stored;
    vars_stored.resize(nvars_model_file);

    for (int i = 0; i < model::npts; i++) {
        initial_model_file >> model::profile(model_index).r(i);

        for (int j = 0; j < nvars_model_file; j++) {
            initial_model_file >> vars_stored[j];
        }

        for (int j = 0; j < model::nvars; j++) {
            model::profile(model_index).state(i,j) = 0.0_rt;
        }

        // make sure that each of the variables we care about is found
        found_dens = false;
        found_temp = false;
        found_pres = false;
        found_velr = false;
        for (int n = 0; n < NumSpec; n++) {
            found_spec[n] = false;
        }
#if NAUX_NET > 0
        for (int n = 0; n < NumAux; n++) {
            found_aux[n] = false;
        }
#endif

        for (int j = 0; j < nvars_model_file; j++) {

            // keep track of whether the current variable from the model
            // file is one that we care about
            found_model = false;

            if (varnames_stored[j] == "density") {
                model::profile(model_index).state(i,model::idens) = vars_stored[j];
                found_model = true;
                found_dens = true;

            } else if (varnames_stored[j] == "temperature") {
                model::profile(model_index).state(i,model::itemp) = vars_stored[j];
                found_model = true;
                found_temp = true;

            } else if (varnames_stored[j] == "pressure") {
                model::profile(model_index).state(i,model::ipres) = vars_stored[j];
                found_model = true;
                found_pres = true;

            } else if (varnames_stored[j] == "velocity") {
                model::profile(model_index).state(i,model::ivelr) = vars_stored[j];
                found_model = true;
                found_velr = true;

            } else {
                for (int comp = 0; comp < NumSpec; comp++) {
                    if (varnames_stored[j] == spec_names_cxx[comp]) {
                        model::profile(model_index).state(i,model::ispec+comp) = vars_stored[j];
                        found_model = true;
                        found_spec[comp] = true;
                        break;
                    }
                }
#if NAUX_NET > 0
                if (!found_model) {
                    for (int comp = 0; comp < NumAux; comp++) {
                        if (varnames_stored[j] == aux_names_cxx[comp]) {
                            model::profile(model_index).state(i,model::iaux+comp) = vars_stored[j];
                            found_model = true;
                            found_aux[comp] = true;
                            break;
                        }
                    }
                }
#endif
            }

            // yell if we didn't find the current variable

            if (!found_model && i == 0) {
                amrex::Print() << "Warning variable not found: " << varnames_stored[j] << std::endl;
            }

        } // end loop over nvars_model_file

        //  were all the variables we care about provided?

        if (i == 0) {
            if (!found_dens) {
                amrex::Print() << "WARNING: density not provided in inputs file" << std::endl;
            }

            if (!found_temp) {
                amrex::Print() << "WARNING: temperature not provided in inputs file" << std::endl;
            }

            if (!found_pres) {
                amrex::Print() << "WARNING: pressure not provided in inputs file" << std::endl;
            }

            if (!found_velr) {
                amrex::Print() << "WARNING: velocity not provided in inputs file" << std::endl;
            }

            for (int comp = 0; comp < NumSpec; comp++) {
                if (!found_spec[comp]) {
                    amrex::Print() << "WARNING: " << spec_names_cxx[comp] << " not provided in inputs file" << std::endl;
                }
            }

#if NAUX_NET > 0
            for (int comp = 0; comp < NumAux; comp++) {
                if (!found_aux[comp]) {
                    amrex::Print() << "WARNING: " << aux_names_cxx[comp] << " not provided in inputs file" << std::endl;
                }
            }
#endif
        }

    }  // end of loop over lines in the model file

    initial_model_file.close();

    model::initialized = true;
}


#endif
