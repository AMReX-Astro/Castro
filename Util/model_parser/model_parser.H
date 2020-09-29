#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>


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


// integer keys for indexing the model_state array
constexpr int nvars_model = 3 + NumSpec;
constexpr int idens_model = 0;
constexpr int itemp_model = 1;
constexpr int ipres_model = 2;
constexpr int ispec_model = 3;

extern AMREX_GPU_MANAGED int npts_model;
extern AMREX_GPU_MANAGED bool model_initialized;

extern amrex::Gpu::ManagedVector<Real> model_state;
extern amrex::Gpu::ManagedVector<Real> model_r;

#define MODELINDEX(I, J) npts_model*I+J;


AMREX_INLINE
void
read_model_file(std::string& model_file) {

    bool found_model, found_dens, found_temp, found_pres;
    bool found_spec[NumSpec];

    // read in the initial model

    std::ifstream initial_model_file;

    initial_model_file.open("model.hse.cool.coulomb", std::ios::in);
    if (!initial_model_file.is_open()) {
        amrex::Error("Error opening the initial model");
    }

    std::string line;

    // first the header line -- this tells us the number of points

    getline(initial_model_file, line);
    std::string npts_string = line.substr(line.find("=")+1, line.length());
    npts_model = std::stoi(npts_string);

    // next line tells use the number of variables

    getline(initial_model_file, line);
    std::string num_vars_string = line.substr(line.find("=")+1, line.length());
    int nvars_model_file = std::stoi(num_vars_string);

    // now read in the names of the variables

    std::vector<std::string> varnames_stored;
    for (int n = 0; n < nvars_model_file; n++) {
        getline(initial_model_file, line);
        std::string var_string = line.substr(line.find("#")+1, line.length());
        varnames_stored.push_back(ltrim(var_string));
    }


    // allocate storage for the model data

    model_state.resize(npts_model * nvars_model);
    model_r.resize(npts_model);

    amrex::Print() << "reading initial model" << std::endl;
    amrex::Print() << npts_model << " points found in the initial model" << std::endl;
    amrex::Print() << nvars_model_file << " variables found in the initial model file" << std::endl;

    // start reading in the data
    Real vars_stored[nvars_model_file];

    for (int i = 0; i < npts_model; i++) {
        initial_model_file >> model_r[i];

        for (int j = 0; j < nvars_model_file; j++) {
            initial_model_file >> vars_stored[j];
        }

        for (int j = 0; j < nvars_model; j++) {
            model_state(MODELINDEX(i,j)) = 0.0_rt;
        }

        // make sure that each of the variables we care about is found
        found_dens = false;
        found_temp = false;
        found_pres = false;
        for (int n = 0; n < NumSpec; n++) {
            found_spec[n] = false;
        }

        for (int j = 0; j < nvars_model_file; j++) {

            // keep track of whether the current variable from the model
            // file is one that we care about
            found_model = false;

            if (varnames_stored[j] == "density") {
                model_state(MODELINDEX(i,idens_model)) = vars_stored[j];
                found_model = true;
                found_dens = true;

            } else if (varnames_stored[j] == "temperature") {
                model_state(MODELINDEX(i,itemp_model)) = vars_stored[j];
                found_model = true;
                found_temp = true;

            } else if (varnames_stored[j] == "pressure") {
                model_state(MODELINDEX(i,ipres_model)) = vars_stored[j];
                found_model = true;
                found_pres = true;

            } else {
                for (int comp = 0; comp < NumSpec; comp++) {
                    if (varnames_stored(j) == spec_names_cxx[comp]) {
                        model_state(MODELINDEX(i,ispec_model+comp)) = vars_stored[j];
                        found_model = true;
                        found_spec[comp] = true;
                        break;
                    }
                }
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

            for (int comp = 0; comp < NumSpec; comp++) {
                if (!found_spec[comp]) {
                    amrex::Print() << "WARNING: " << spec_names_cxx[comp] << " not provided in inputs file" << std::end;
                }
            }
        }

    }

    initial_model_file.close();

    model_initialized = true;
}

AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
interpolate(const Read r, const int var_index) {

    // find the value of model_state component var_index at point r
    // using linear interpolation.  Eventually, we can do something
    // fancier here.

    int id = locate(r);

    Real slope;
    Real interp;

    if (id == 0) {

       slope = (model_state(MODELINDEX(id+1, var_index)) -
                model_state(MODELINDEX(id, var_index))) /
           (model_r(id+1) - model_r(id));
       interp = slope * (r - model_r(id)) + model_state(MODELINDEX(id, var_index));

       // safety check to make sure interp lies within the bounding points
       Real minvar = amrex::min(model_state(MODELINDEX(id+1, var_index)),
                                model_state(MODELINDEX(id, var_index)));
       Real maxvar = amrex::max(model_state(MODELINDEX(id+1, var_index)),
                                model_state(MODELINDEX(id, var_index)));
       interp = amrex::max(interp, minvar);
       interp = amrex::min(interp, maxvar);

    } else if (id == npts_model-1) {

       slope = (model_state(MODELINDEX(id, var_index)) -
                model_state(MODELINDEX(id-1, var_index))) /
           (model_r(id) - model_r(id-1));
       interp = slope * (r - model_r(id)) + model_state(MODELINDEX(id, var_index));


       // safety check to make sure interp lies within the bounding points
       Real minvar = amrex::min(model_state(MODELINDEX(id-1, var_index)),
                                model_state(MODELINDEX(id, var_index)));
       Real maxvar = amrex::max(model_state(MODELINDEX(id-1, var_index)),
                                model_state(MODELINDEX(id, var_index)));
       interp = amrex::max(interp, minvar);
       interp = amrex::min(interp, maxvar);

    } else {

       if (r >= model_r[id]) {

          slope = (model_state(MODELINDEX(id+1, var_index)) -
                   model_state(MODELINDEX(id, var_index))) /
              (model_r(id+1) - model_r(id));
          interp = slope * (r - model_r[id]) + model_state(MODELINDEX(id, var_index));

       } else {

           slope = (model_state(MODELINDEX(id, var_index)) -
                    model_state(MODELINDEX(id-1, var_index))) /
               (model_r(id) - model_r(id-1));
           interp = slope * (r - model_r(id)) + model_state(MODELINDEX(id, var_index));

       }

    }

    return interp;

}


AMREX_INLINE AMREX_GPU_HOST_DEVICE
void
locate_sub(const Real r) {

    int loc;

    if (r <= model_r[0]) {
       loc = 0;

    } else if (r > model_r[npts_model-2]) {
       loc = npts_model-1;

    } else {

        int ilo = 0;
        int ihi = npts_model-2;

        while (ilo+1 != ihi) {
            int imid = (ilo + ihi) / 2;

            if (r <= model_r[imid]) {
                ihi = imid;
            } else {
                ilo = imid;
            }
        }

        loc = ihi;
    }

    return loc;
}

#endif
