#ifndef interpolate_H
#define interpolate_H

#include <AMReX.H>

using namespace amrex;

namespace interpolate {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
int locate (Real x, int n, const Real* xs)
{
    int loc;

    int ilo, ihi, imid;

    if (x <= xs[0]) {
        loc = 0;
    }
    else if (x > xs[n-2]) {
        loc = n-1;
    }
    else {
        ilo = 0;
        ihi = n-2;

        while (ilo + 1 != ihi) {
            imid = (ilo + ihi) / 2;

            if (x <= xs[imid]) {
                ihi = imid;
            }
            else {
                ilo = imid;
            }
        }

        loc = ihi;
    }

    return loc;
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
Real interpolate(Real r, int npts_model, const Real* model_r, const Real* model_var, int iloc = -1)
{
    // given the array of model coordinates (model_r), and variable (model_var),
    // find the value of model_var at point r (var_r) using linear interpolation.
    // Eventually, we can do something fancier here.

    // find the location in the coordinate array where we want to interpolate
    int id;
    if (iloc >= 0) {
        id = iloc;
    }
    else {
        id = locate(r, npts_model, model_r);
    }

    Real interp;

    if (id == 0) {

        Real slope = (model_var[id+1] - model_var[id])/  (model_r[id+1] - model_r[id]);
        interp = slope * (r - model_r[id]) + model_var[id];

        // safety check to make sure interp lies within the bounding points
        Real minvar = amrex::min(model_var[id+1], model_var[id]);
        Real maxvar = amrex::max(model_var[id+1], model_var[id]);
        interp = amrex::max(interp, minvar);
        interp = amrex::min(interp, maxvar);

    }
    else if (id == npts_model-1) {

        Real slope = (model_var[id] - model_var[id-1]) / (model_r[id] - model_r[id-1]);
        interp = slope * (r - model_r[id]) + model_var[id];

        // safety check to make sure interp lies within the bounding points
        Real minvar = amrex::min(model_var[id], model_var[id-1]);
        Real maxvar = amrex::max(model_var[id], model_var[id-1]);
        interp = amrex::max(interp, minvar);
        interp = amrex::min(interp, maxvar);

    }
    else {

        if (r >= model_r[id]) {

            Real slope = (model_var[id+1] - model_var[id]) / (model_r[id+1] - model_r[id]);
            interp = slope * (r - model_r[id]) + model_var[id];

            // // safety check to make sure interp lies within the bounding points
            // Real minvar = amrex::min(model_var[id+1], model_var[id]);
            // Real maxvar = amrex::max(model_var[id+1], model_var[id]);
            // interp = amrex::max(interp, minvar);
            // interp = amrex::min(interp, maxvar);

        }
        else {

            Real slope = (model_var[id] - model_var[id-1]) / (model_r[id] - model_r[id-1]);
            interp = slope * (r - model_r[id]) + model_var[id];

            // // safety check to make sure interp lies within the bounding points
            // Real minvar = amrex:: min(model_var[id], model_var[id-1]);
            // Real maxvar = amrex::max(model_var[id], model_var[id-1]);
            // interp = amrex::max(interp, minvar);
            // interp = amrex::min(interp, maxvar);

        }

    }

    return interp;
}

}

#endif
