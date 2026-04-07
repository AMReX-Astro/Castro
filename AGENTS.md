# Castro Agents Guide

Use this guide whenever you orchestrate explorers/workers inside the Castro repository. It covers both Castro developers (PR reviews, bug hunts, new features, and documentation) and Castro users who ask agents for help learning or building with Castro. Castro itself is a C++ repository for computational astrophysics simulations. Castro relies on AMReX to provide block-structured adaptive mesh refinement (AMR) targeting large-scale PDE simulations on CPU and GPU architectures (CUDA, HIP, SYCL).

## Purpose & Personas

- **Castro developers** – structure every agent task (reviews, fixes, features, documentation) so it is scoped, reproducible, and merged with confidence.
- **Castro users** – route questions about capabilities, docs, tutorials, builds, or troubleshooting through the authoritative resources already shipped with the repo.

## Repository Layout at a Glance

- `Source/` – Primary C++ implementation, organized by physics module:
  - `driver/` – Main `Castro` class (`Castro.H`, `Castro.cpp`), time advancement (`Castro_advance*.cpp`), I/O, setup, and diagnostics.
  - `hydro/` – Unsplit finite-volume hydrodynamics (2nd/4th order, 1D/2D/3D).
  - `radiation/` – Multigroup flux-limited diffusion radiation hydrodynamics.
  - `gravity/` – Multiple gravity solvers, including a full Poisson gravity solver with isolated boundary conditions.
  - `mhd/` – Constrained transport ideal MHD.
  - `reactions/` – Nuclear reaction networks and burning.
  - `sdc/` – Spectral Deferred Corrections time integration.
  - `diffusion/` – Explicit thermal diffusion.
  - `rotation/` – Rotating frame physics.
  - `scf/` – Self-consistent field gravity module.
  - `sources/` – Source term handling.
  - `particles/` – Tracer particle transport.
- `Exec/` – Problem directories, each containing its own `GNUmakefile` and source files. Organized into groups: `hydro_tests/`, `radiation_tests/`, `mhd_tests/`, `reacting_tests/`, `gravity_tests/`, `scf_tests/`, `unit_tests/`, `science/`.
- `Docs/` – Sphinx documentation source (`.rst` files). Published at https://amrex-astro.github.io/Castro/docs/. Edit these when updating guides.
- `Util/` – Helper scripts and tools.
- `external/` – Git submodule dependencies:
  - `amrex/` – AMReX adaptive mesh refinement library (mesh, parallelism, data structures).
  - `Microphysics/` – Equations of state, nuclear reaction networks, conductivity models.

## Build System

Castro uses **GNUmake exclusively** (not CMake). You always build from inside a problem directory under `Exec/`, never from the repo root or `Source/`.

```bash
cd Exec/hydro_tests/Sod
make clean   # run this if the existing executable is significantly older than recent source changes
make -j4
./Castro1d.gnu.ex inputs-sod-x
```

The executable name encodes dimension and compiler: `Castro<Nd>.<comp>[.MPI][.OMP].ex`.

### Key GNUmakefile variables

| Variable | Values | Purpose |
|---|---|---|
| `DIM` | `1`, `2`, `3` | Spatial dimensionality |
| `COMP` | `gnu`, `llvm`, `intel`, `cray` | Compiler suite |
| `DEBUG` | `TRUE`/`FALSE` | Debug vs. optimized build |
| `USE_MPI` | `TRUE`/`FALSE` | MPI parallelism |
| `USE_OMP` | `TRUE`/`FALSE` | OpenMP threading |
| `USE_CUDA` | `TRUE`/`FALSE` | NVIDIA GPU support |
| `USE_HIP` | `TRUE`/`FALSE` | AMD GPU support |
| `EOS_DIR` | e.g., `gamma_law`, `helmholtz` | EOS from `Microphysics/EOS/` |
| `NETWORK_DIR` | e.g., `general_null`, `aprox13` | Reaction network from `Microphysics/networks/` |
| `CASTRO_HOME` | path to repo root | Set automatically by most problem makefiles |

Variables can be set in the `GNUmakefile` or passed on the command line:
```bash
make -j4 DEBUG=TRUE USE_MPI=FALSE
```

### Physics feature flags

These are set in the problem's `GNUmakefile`:

| Flag | Physics enabled |
|---|---|
| `USE_RADIATION = TRUE` | Multigroup flux-limited diffusion (also requires `HYPRE_DIR`) |
| `USE_GRAVITY = TRUE` | Full Poisson gravity |
| `USE_MHD = TRUE` | Magnetohydrodynamics |
| `USE_REACTIONS = TRUE` | Nuclear burning |
| `USE_DIFFUSION = TRUE` | Thermal diffusion |
| `USE_SDC = TRUE` | Spectral Deferred Corrections integration |
| `USE_ROTATION = TRUE` | Rotating frame physics |

## Problem Anatomy

Each problem under `Exec/` is a self-contained directory. Required files:

- **`GNUmakefile`** – sets `DIM`, `COMP`, physics flags, `EOS_DIR`, `NETWORK_DIR`; ends with `include $(CASTRO_HOME)/Exec/Make.Castro`.
- **`problem_initialize.H`** – C++ header with per-grid problem initialization routine.
- **`problem_initialize_state_data.H`** – fills the initial state data on the mesh.
- **`Make.package`** – lists any problem-specific source files for the build system.
- **`inputs`** – runtime control file; sets AMReX and Castro parameters (see below).

Optional:
- **`_prob_params`** – defines problem-specific runtime parameters. Python scripts auto-generate C++ code from this file. Format: `name  datatype  default  namelist?  size`.

**To create a new problem**: copy an existing similar problem directory and modify it. Do not start from scratch.

## Runtime Parameters

All runtime configuration goes in the `inputs` file, not in source code. Parameter namespaces:

- `amr.*` – AMReX mesh and refinement control (e.g., `amr.max_level`, `amr.n_cell`)
- `castro.*` – Castro physics parameters (e.g., `castro.cfl`, `castro.use_flattening`)
- `geometry.*` – Domain geometry
- `problem.*` – Problem-specific parameters defined in `_prob_params`

Full parameter reference: https://amrex-astro.github.io/Castro/docs/inputs.html

## Code Conventions

- **C++20** standard; no Fortran in the codebase.
- Use **`amrex::Real`** for floating-point values, not raw `double` or `float`.
- GPU portability: compute kernels use `amrex::ParallelFor` with lambda captures. Do not write raw loops over mesh indices in physics code.
- State array component indices are named constants defined in header files — use them, not raw integers.
- Each `Source/` subdirectory has a `Make.package` listing its files; update it when adding new source files.

## Operating Principles

- **Development Model**: Work from short-lived branches based on the latest `development`, and never commit directly on the tracking `development` branch (see the “Development Model” section of `README.md`). The `main` branch is used for monthly releases.
- **CHANGES.md**: Add a line summarizing any bug fix or new feature, referencing the PR number.
- **Issue logging & hand-off**: Keep a personal, untracked scratchpad on each machine (we recommend `agent-notes/<NN>-<component>-<short-description>.md`). Use it to capture open questions, repro notes, or follow-ups. Include suggested patches whenever possible so the next agent can act quickly.
- **Learn from past bugs**: If you already keep a local `agent-notes/` notebook, skim it before diving into similar code to refresh common pitfalls.

## Developer Playbooks

### PR and Bug Reviews
1. **Sync & inspect** – Update the local branch, note the PR/issue scope, and identify which physics module(s) are affected.
2. **Reproduce** – Build the relevant problem directory with the same flags as the PR author. Check the `GNUmakefile` for `DIM`, `USE_*` flags, `EOS_DIR`, and `NETWORK_DIR`.
3. **Read the diff** – Confirm changes follow code conventions (`amrex::Real`, `ParallelFor` kernels, named state indices, `Make.package` updated for new files).
4. **Run** – Execute the problem with the appropriate `inputs` file and compare output to expected results if a baseline exists.
5. **Report** – Summarize findings (blocking issues first), cite files and line numbers, note which tests need to pass.
6. **Log follow-ups** – Record remaining work in `agent-notes/` or the PR description.

### Feature or Fix Implementation
1. **Understand scope** – Identify which `Source/` subdirectory owns the physics, and which `USE_*` flag gates it.
2. **Find a reference problem** – Locate a problem under `Exec/` that exercises the affected physics. Use it to build and test.

   ```bash
   cd Exec/hydro_tests/Sod       # example; pick the appropriate problem
   make -j4 DIM=1 COMP=gnu
   ./Castro1d.gnu.ex inputs-sod-x
   ```

3. **Implement** – Touch only files relevant to the task. If adding a new source file to `Source/<module>/`, add it to that module's `Make.package`.
4. **Test** – Rebuild and run. For physics changes, compare against a known-good baseline or analytic solution.
5. **Document** – Update `Docs/source/` RST files if behavior changes. Update `CHANGES.md` with a one-line summary referencing the PR.
6. **Hand off** – Record open questions, test commands, and outputs in `agent-notes/` or the PR description.

### Creating a New Problem
1. Copy an existing problem of similar physics into a new directory under the appropriate `Exec/` group.
2. Update `GNUmakefile`: set `DIM`, `EOS_DIR`, `NETWORK_DIR`, and any required `USE_*` flags.
3. Edit `problem_initialize.H` and `problem_initialize_state_data.H` for the new initial conditions.
4. Add or remove entries in `_prob_params` for any problem-specific runtime parameters.
5. Create or update the `inputs` file with appropriate `amr.*`, `castro.*`, and `problem.*` settings.
6. Build and run to confirm the setup is correct before writing physics.

### Documentation Updates
- RST source lives in `Docs/source/`. The published guide mirrors it exactly.
- When adding a new runtime parameter, document it in the appropriate `Docs/source/inputs.rst` section.
- When adding a new physics module or `USE_*` flag, add it to `Docs/source/build_system.rst`.

## Quick Checklist

1. Confirm you are on a task-specific branch that tracks `development` cleanly.
2. Identify the physics module (`Source/<module>/`) and a representative problem (`Exec/`) before writing any code.
3. Build with: `cd Exec/<group>/<problem> && make -j4` using the correct `DIM`, `COMP`, and `USE_*` flags. Run `make clean` first if the existing executable is significantly older than recent source changes.
4. Use `amrex::Real`, `amrex::ParallelFor`, and named state indices — not raw types or loops.
5. Update `Make.package` when adding new source files.
6. Update `CHANGES.md` and relevant `Docs/source/` RST files when the change is user-visible.
7. Capture unresolved work and test commands in `agent-notes/` so future agents can pick up where you left off.
