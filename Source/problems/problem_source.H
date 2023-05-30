#ifndef problem_source_H
#define problem_source_H

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void problem_source (int i, int j, int k,
                     GeometryData const& geomdata,
                     Array4<const Real> const& state,
                     Array4<Real> const& src,
                     const Real dt, const Real time) {

    amrex::ignore_unused(i);
    amrex::ignore_unused(j);
    amrex::ignore_unused(k);
    amrex::ignore_unused(geomdata);
    amrex::ignore_unused(state);
    amrex::ignore_unused(src);
    amrex::ignore_unused(dt);
    amrex::ignore_unused(time);

}

#endif
