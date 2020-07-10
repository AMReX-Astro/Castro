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
  - name: Jean Sexton
    orcid: 0000-0003-2551-1678
    affiliation: 1
  - name: Donald Willcox
    orcid: 0000-0003-2300-5165
    affiliation: 1
  - name: Weiqun Zhang
    orcid: 0000-0001-8092-1974
    affiliation: 1
  - name: Michael Zingale
    orcid: 0000-0001-8401-030X
    affiliation: "2, 4"
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Department of Physics and Astronomy, Stony Brook University
    index: 2
  - name: NVIDIA Corporation
    index: 3
  - name: Center for Computational Astrophysics, Flatiron Institute
    index: 4
date: 02 July 2020
bibliography: paper.bib
---

# Summary 
Castro is a highly parallel, adaptive mesh, multiphysics
simulation code for compressible astrophysical flows.  It has been
used to simulate different progenitor models of Type Ia supernovae,
X-ray bursts, core-collapse and electron capture supernovae, and
dynamics in exoplanets.  Together, Castro, the low Mach number code
MAESTROeX [@maestroex], and the cosmology code Nyx [@nyx] make up the
AMReX-Astrophysics Suite of open-source, adaptive mesh, performance
portable astrophysical simulation codes.

The core hydrodynamics solver in Castro [@castro] is based on the
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
that aims to prevent the hydro and reactions from becoming decoupled.  The
simplified SDC method uses the CTU PPM hydro together with an
iterative scheme to fully couple the reactions and hydro, still to
second order [@simple_sdc].  Alternatively, we have implemented a
traditional SDC method that couples hydro and reactions to both second
and fourth-order in space and time [@castro_sdc] (at present, this
method is single-level only).  The Strang splitting and simplified SDC
methods have a retry scheme, where a timestep will be rejected and retried
at a smaller, subcycled timestep if the burning solve fails to meet its
tolerance, negative densities are generated, or we violate one of the
timestepping criteria.

In addition to reactive hydrodynamics, Castro includes full
self-gravity with isolated boundary conditions and rotation, both
implemented in an energy-conserving fashion, explicit thermal
diffusion, and gray [@castroII] and multigroup [@castroIII] flux
limited diffusion radiation hydrodynamics.  A constrained transport
MHD solver based on the CTU algorithm is also available, and can use
the same physics source terms.  Castro can use an arbitrary equation of
state and reaction network, and these microphysics routines are
provided by the StarKiller project [@starkiller].

Castro is built on the AMReX [@AMReX] adaptive mesh refinement (AMR)
library and is largely written in C++ with a few Fortran compute
kernels.  AMR levels are advanced at their own timestep (subcycling)
and jumps by factors of 2 and 4 are supported between levels.  We use
MPI to distribute AMR grids across nodes and use logical tiling with
OpenMP to divide a grid across threads for multi-core CPU machines
(exposing coarse-grained parallelism) or CUDA to spread the work across
GPU threads on GPU-based machines (fine-grained parallelism).  All of
the core physics can run on GPUs and has been shown to scale well to
thousands of GPUs [@castro_2019] and hundreds of thousands of CPU cores
[@castro_2017].  For performance portability, we use the same source code
for both CPUs and GPUs, and implement our parallel loops in an abstraction
layer provided by AMReX. An abstract parallel for loop accepts as arguments
a range of indices and the body of the loop to execute for a given index,
and the AMReX backend dispatches the work appropriately (e.g. one zone per
GPU thread). This strategy is similar to the way the Kokkos [@Kokkos] and
RAJA [@RAJA] abstraction models provide performance portability in C++.

While there are a number of astrophysical hydrodynamics simulation codes, Castro
offers a few unique features.  The original motivation for developing
Castro was to build a simulation code based on a modern,
well-supported AMR library (BoxLib which evolved to AMReX), using
unsplit integration techniques and targeting problems in nuclear
astrophysics.  The radiation solver was a key design consideration in
the early development.  The large developer community contributing to AMReX
(representing a large number of application codes across various domains)
results in Castro continually gaining optimizations for new
architectures.  As Castro evolved, we adopted a fully open development
model (as does the Enzo [@enzo] code, for example).  We pride ourselves in
making all of the science problems available in the Castro git repository as
we are developing them, and the infrastructure we use for running our problems
and writing our science papers is publicly available in the AMReX-Astro organization.
Other simulation codes, like Flash [@flash], also work with a general equation of
state and reaction network, but Castro is unique in focusing on
spectral deferred correction techniques for coupling the hydro and
reactions.  Finally, while some astrophysics codes have performance portable forks
(like K-Athena [@kathena], which uses Kokkos), Castro's current design -- which targets both
CPUs and GPUs for all solvers -- achieves performance portability as a core
design principle, avoiding the need for a fractured development model.




# Acknowledgements

The work at Stony Brook was supported by DOE/Office of Nuclear Physics
grant DE-FG02-87ER40317 and NSF award AST-1211563.  MZ acknowledges
support from the Simons Foundation.  This research was supported by
the Exascale Computing Project (17-SC-20-SC), a collaborative effort
of the U.S. Department of Energy Office of Science and the National
Nuclear Security Administration.  The work at LBNL was supported by
U.S. Department of Energy under contract No. DE-AC02-05CH11231.  We
also thank NVIDIA Corporation for the donation of a Titan X Pascal and
Titan V used in this research.  The GPU development of Castro
benefited greatly from numerous GPU hackathons arranged by OLCF.

# References

