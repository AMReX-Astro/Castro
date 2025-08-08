#ifndef TIMESTEP_H
#define TIMESTEP_H

#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_MultiFab.H>

namespace timestep {

    ///
    /// Compute the CFL timestep
    ///
    amrex::ValLocPair<amrex::Real, amrex::IntVect> estdt_cfl (const amrex::MultiFab& stateMF,
                                                              const amrex::GeometryData& geomdata);


#ifdef MHD
    ///
    /// Compute the CFL timestep
    ///
    amrex::ValLocPair<amrex::Real, amrex::IntVect>  estdt_mhd (const amrex::MultiFab& stateMF, const amrex::MultiFab& bx,
                                                               const amrex::MultiFab& by, const amrex::MultiFab& bz,
                                                               const amrex::GeometryData& geomdata);
#endif

    ///
    /// Diffusion-limited timestep
    ///
    amrex::ValLocPair<amrex::Real, amrex::IntVect> estdt_temp_diffusion (const amrex::MultiFab& stateMF,
                                                                         const amrex::GeometryData& geomdata);

#ifdef REACTIONS
    ///
    /// Reactions-limited timestep
    ///
    amrex::ValLocPair<amrex::Real, amrex::IntVect> estdt_burning (const amrex::MultiFab& stateMF,
                                                                  const amrex::MultiFab& maskMF,
                                                                  const amrex::GeometryData& geomdata);
#endif

#ifdef RADIATION
    ///
    /// Radiation hydro timestep
    ///
    amrex::Real estdt_rad (const amrex::MultiFab& stateMF, const amrex::MultiFab& radMF,
                           const amrex::GeometryData& geomdata);
#endif

}

#endif // TIMESTEP_H
