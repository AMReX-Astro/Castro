---
title: 'CASTRO: A Massively Parallel Compressible Astrophysics Simulation Code'
tags:
  - C++
  - Fortran90
  - convection
  - hydrodynamics
  - nuclear reactions
  - nucleosynthesis
  - abundances
  - supernovae
authors:
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    affiliation: 1
  - name: Maria Barrios Sazo
    orcid: 0000-0002-3185-9809
    affiliation: 2
  - name: John Bell
    orcid: 0000-0002-5749-334X
    affiliation: 1
  - name: Alice Harpole
    orcid: 0000-0002-1530-781X
    affiliation: 2
  - name: Max Katz
    orcid: 0000-0003-0439-4556
    affiliation: 3
  - name: Donald Willcox
    orcid: 0000-0003-2300-5165
    affiliation: 1
  - name: Weiqun Zhang
    orcid: 0000-0001-8092-1974
    affiliation: 1
  - name: Michael Zingale
    orcid: 0000-0001-8401-030X
    affiliation: 2
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Department of Physics and Astronomy, Stony Brook University
    index: 2
  - name: NVIDIA Corp
    index: 3
date: 17 February 2020
bibliography: paper.bib
---

# Summary Castro is a highly parallel, adaptive mesh, multiphysics
simulation code for compressible astrophysical flows.  It has been
used to simulate different progenitor models of Type Ia supernovae,
X-ray bursts, core-collapse and electron capture supernovae, and
dynamics in exoplanets.  Together, Castro, the low Mach number code
MAESTROeX [@maestroex], and the cosmology code Nyx [@nyx] make up the
AMReX-Astrophysics Suite of open-source, adaptive mesh, performance
portable astrophysical simulation codes.

The core hydrodynamics solver in Castro [@castro] is based on a the
directionally unsplit corner transport upwind method of [@ctu] with
piecewise parabolic reconstruction [@ppm].  Modeling reactive flows in
stellar environments is a core capability of Castro.  Astrophysical
reaction networks are stiff and require implicit integration
techniques for accurate and stable solutions.  In Castro, we have
several modes of coupling the reactions to hydro.  The simplest method
is the traditional operator splitting approach, using Strang splitting
to achieve second-order in time.  However, when the reactions are
energetic this coupling can break down, and we have two different
implementations based on spectral deferred corrections (SDC), a method
that aims prevent the hydro and reactions from becoming decoupled.  The
simplified SDC method uses the CTU PPM hydro together with an
iterative scheme to fully couple the reactions and hydro, still to
second order [@simple_sdc].  Both of these methods have a retry
scheme, where a timestep will be rejected if the burning solve fails
to meet its tolerance, negative densities are generated, or we violate
one of the timestepping criterion.  Alternately, we have implemented a
traditional SDC method that couples hydro and reactions to both second
and fourth-order in space and time [@castro_sdc] (at present, this
method is single-level only).

In addition to reactive hydrodynamics, Castro includes full
self-gravity with isolated boundary conditions and rotation, both
implemented in an energy-conserving fashion, explicit thermal
diffusion, and gray [@castroII] and multigroup [@castroIII] flux
limited diffusion radiation hydrodynamics.  An MHD solver has been
implemented and will be merged into the master branch in the near
future.  Castro can use an arbitrary equation of state and reaction
network and these microphysics routines are provided by the StarKiller
project [@starkiller].

Castro is built on the AMReX [@AMReX] adaptive mesh refinement (AMR)
library and is largely written in C++ with Fortran compute kernels.
AMR levels are advanced at their own timestep (subcycling) and jumps
by factors of 2 and 4 are supported between levels.  We use MPI to
distribute AMR grids across nodes and using logical tiling with OpenMP
to divide a grid across threads for manycore machines or CUDA to
spread the work across GPU threads on GPU-based machines.  All of the
core physics runs on GPUs (CTU PPM hydro, self-gravity, diffusion,
reactions) and has been shown to scale well to 1000s of GPUs
[@castro_2019].  For performance portability, the same compute kernel
is used for CPUs or GPUs, using either a custom preprocessor pragrma
for Fortran or lambda-capturing for C++ kernels.





# Acknowledgements

The work at Stony Brook was supported by DOE/Office of Nuclear Physics
grant DE-FG02-87ER40317 and NSF award AST-1211563.  MZ acknowledges
support from the Simons Foundation.  This research was supported by
the Exascale Computing Project (17-SC-20-SC), a collaborative effort
of the U.S. Department of Energy Office of Science and the National
Nuclear Security Administration.  The work at LBNL was supported by
U.S. Department of Energy under contract No. DE-AC02-05CH11231.  We
also thank NVIDIA Corporation for the donation of a Titan X Pascal and
Titan V used in this research.  The GPU development of \castro\
benefited greatly from numerous GPU hackathons arranged by OLCF.

# References

