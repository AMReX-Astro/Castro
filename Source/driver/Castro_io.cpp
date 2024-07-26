
#ifndef WIN32
#include <unistd.h>
#endif

#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#include <AMReX_Utility.H>
#include <Castro.H>
#include <Castro_io.H>
#include <AMReX_ParmParse.H>

#ifdef RADIATION
#include <Radiation.H>
#endif

#ifdef GRAVITY
#include <Gravity.H>
#endif

#ifdef DIFFUSION
#include <Diffusion.H>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <problem_initialize_state_data.H>
#include <problem_checkpoint.H>
#include <problem_restart.H>
#include <AMReX_buildInfo.H>

using std::string;
using namespace amrex;

// Castro maintains an internal checkpoint version numbering system.
// This allows us to ensure that we don't attempt to restart from a
// checkpoint that is incompatible with the current code. The version
// number is stored in the CastroHeader file inside a checkpoint.
// The following were the changes made in updating version numbers:
// 0: all checkpoints before we began the numbering system (July 27, 2015)
// 1: PhiGrav_Type was added to the checkpoint
// 2: Source_Type was added to the checkpoint
// 3: A ReactHeader file was generated and the maximum de/dt was stored there
// 4: Reactions_Type added to checkpoint; ReactHeader functionality deprecated
// 5: Simplified_SDC_Source_Type and Simplified_SDC_React_Type added to checkpoint
// 6: Simplified_SDC_Source_Type removed from Castro
// 7: A weights field was added to Reactions_Type; number of ghost zones increased to NUM_GROW
// 8: Reactions_Type modified to use rho * omegadot instead of omegadot; rho * auxdot added
// 9: Rotation_Type was removed from Castro
// 10: Reactions_Type was removed from checkpoints
// 11: PhiRot_Type was removed from Castro
// 12: State_Type's additional ghost zone, used when radiation is enabled, has been removed

namespace
{
    int input_version = -1;
    int current_version = 12;
}

// I/O routines for Castro

void
Castro::restart (Amr&     papa,
                 istream& is,
                 bool     bReadSpecial)
{
    // Let's check Castro checkpoint version first;
    // trying to read from checkpoint; if nonexisting, set it to 0.
    if (input_version == -1) {
        if (ParallelDescriptor::IOProcessor()) {
            std::ifstream CastroHeaderFile;
            std::string FullPathCastroHeaderFile = papa.theRestartFile();
            FullPathCastroHeaderFile += "/CastroHeader";
            CastroHeaderFile.open(FullPathCastroHeaderFile.c_str(), std::ios::in);
            if (CastroHeaderFile.good()) {
                char foo[256];
                // first line: Checkpoint version: ?
                CastroHeaderFile.getline(foo, 256, ':');
                CastroHeaderFile >> input_version;
                CastroHeaderFile.close();
            } else {
                input_version = 0;
            }
        }
        ParallelDescriptor::Bcast(&input_version, 1, ParallelDescriptor::IOProcessorNumber());
    }

    // Check that the version number matches.

    if (input_version != current_version) {
        amrex::Error("Checkpoint format incompatible with current code");
    }

    // Check if there's a file in the header indicating that the
    // previous timestep was limited to hit a plot interval. If so,
    // read in the value, so that the timestep after this restart is
    // limited appropriately.

    if (ParallelDescriptor::IOProcessor()) {

        std::ifstream dtHeaderFile;
        std::string FullPathdtHeaderFile = papa.theRestartFile();
        FullPathdtHeaderFile += "/dtHeader";
        dtHeaderFile.open(FullPathdtHeaderFile.c_str(), std::ios::in);

        if (dtHeaderFile.good()) {

            lastDtPlotLimited = 1;
            dtHeaderFile >> lastDtBeforePlotLimiting;
            dtHeaderFile.close();

        }

    }

    // Check if we have the same state variables

    if (ParallelDescriptor::IOProcessor()) {

      std::ifstream StateListFile;
      std::string FullPathStateList = papa.theRestartFile();
      FullPathStateList += "/state_names.txt";
      StateListFile.open(FullPathStateList.c_str(), std::ios::in);

      if (StateListFile.good()) {
        std::string var;
        for (int n = 0; n < NUM_STATE; n++) {
          StateListFile >> var;
          if (var != desc_lst[State_Type].name(n)) {
            amrex::Error("state variables do not agree");
          }
        }
        StateListFile.close();
      }
    }

    ParallelDescriptor::Bcast(&lastDtPlotLimited, 1, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(&lastDtBeforePlotLimiting, 1, ParallelDescriptor::IOProcessorNumber());

    // also need to mod checkPoint function to store the new version in a text file

    AmrLevel::restart(papa,is,bReadSpecial);

    buildMetrics();

    initMFs();

    // get the elapsed CPU time to now;
    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
      // get elapsed CPU time
      std::ifstream CPUFile;
      std::string FullPathCPUFile = parent->theRestartFile();
      FullPathCPUFile += "/CPUtime";
      CPUFile.open(FullPathCPUFile.c_str(), std::ios::in);

      if (CPUFile.good()) {

          CPUFile >> previousCPUTimeUsed;
          CPUFile.close();

      }

      std::cout << "read CPU time: " << previousCPUTimeUsed << "\n";

    }

#ifdef GRAVITY
    if (use_point_mass && level == 0)
    {

        // get the current value of the point mass
        std::ifstream PMFile;
        std::string FullPathPMFile = parent->theRestartFile();
        FullPathPMFile += "/point_mass";
        PMFile.open(FullPathPMFile.c_str(), std::ios::in);

        if (PMFile.good()) {
            PMFile >> point_mass;
            PMFile.close();
        }

    }
#endif

#ifdef ROTATION
    if (do_rotation && level == 0)
    {
        // get current value of the rotation period
        std::ifstream RotationFile;
        std::string FullPathRotationFile = parent->theRestartFile();
        FullPathRotationFile += "/Rotation";
        RotationFile.open(FullPathRotationFile.c_str(), std::ios::in);

        if (RotationFile.is_open()) {
            RotationFile >> castro::rotational_period;
            amrex::Print() << "  Based on the checkpoint, setting the rotational period to "
                           << std::setprecision(7) << std::fixed << castro::rotational_period << " s.\n";
            RotationFile.close();
        }
    }
#endif

    if (level == 0)
    {
        // get problem-specific stuff -- note all processors do this,
        // eliminating the need for a broadcast
        std::string dir = parent->theRestartFile();

        problem_restart(dir);
    }

    if ( (grown_factor > 1) && (parent->maxLevel() < 1) )
    {
       std::cout << "grown_factor is " << grown_factor << std::endl;
       std::cout << "max_level is " << parent->maxLevel() << std::endl;
       amrex::Error("Must have max_level > 0 if doing special restart with grown_factor");
    }

    if (grown_factor > 1 && level == 0)
    {
       if (verbose && ParallelDescriptor::IOProcessor()) {
         std::cout << "Doing special restart with grown_factor " << grown_factor << std::endl;
       }
       MultiFab& S_new = get_new_data(State_Type);

       Box orig_domain;
       if (star_at_center == 0) {
          orig_domain = amrex::coarsen(geom.Domain(),grown_factor);
       } else if (star_at_center == 1) {

          Box domain(geom.Domain());
          int lo=0, hi=0;
          if (geom.IsRZ()) {
             if (grown_factor != 2) {
                amrex::Abort("Must have grown_factor = 2");
             }

             int d = 0;
             int dlen =  domain.size()[d];
             lo = 0;
             hi = dlen/2;
             orig_domain.setSmall(d,lo);
             orig_domain.setBig(d,hi);

             d = 1;
             dlen =  domain.size()[d];
             lo =   dlen/4    ;
             hi = 3*dlen/4 - 1;
             orig_domain.setSmall(d,lo);
             orig_domain.setBig(d,hi);

          } else {
             for (int d = 0; d < AMREX_SPACEDIM; d++)
             {
                int dlen =  domain.size()[d];
                if (grown_factor == 2) {
                   lo =   dlen/4    ;
                   hi = 3*dlen/4 - 1;
                } else if (grown_factor == 3) {
                   lo =   (dlen)/3    ;
                   hi = 2*(dlen)/3 - 1;
                } else {
                   amrex::Abort("Must have grown_factor = 2 or 3");
                }
                orig_domain.setSmall(d,lo);
                orig_domain.setBig(d,hi);
             }
          }
       } else {
          if (ParallelDescriptor::IOProcessor()) {
             std::cout << "... invalid value of star_at_center: " << star_at_center << std::endl;
          }
          amrex::Abort();
       }

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
       for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
           const Box& bx = mfi.tilebox();

           if (!orig_domain.contains(bx)) {
               auto s = S_new[mfi].array();
               auto geomdata = geom.data();

#ifdef RNG_STATE_INIT
               amrex::Error("Error: random initialization not yet supported with grown factor");
#else
               amrex::ParallelFor(bx,
               [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                   // C++ problem initialization; has no effect if not implemented
                   // by a problem setup (defaults to an empty routine).
                   problem_initialize_state_data(i, j, k, s, geomdata);
               });
#endif
           }
       }
    }

    if (grown_factor > 1 && level == 1) {
      getLevel(0).avgDown();
    }

#ifdef GRAVITY
    if (do_grav && level == 0) {
       BL_ASSERT(gravity == nullptr);
       gravity = new Gravity(parent,parent->finestLevel(),&phys_bc, URHO);
    }
#endif

#ifdef DIFFUSION
    if (level == 0) {
       BL_ASSERT(diffusion == 0);
       diffusion = new Diffusion(parent,&phys_bc);
    }
#endif

#ifdef RADIATION
    if (do_radiation) {
      if (radiation == 0) {
        // radiation is a static object, only alloc if not already there
        int rad_restart = 1; // disables quasi-steady initialization
        radiation = new Radiation(parent, this, rad_restart);
      }
      radiation->regrid(level, grids, dmap);
      radiation->restart(level, grids, dmap, parent->theRestartFile(), is);

      rad_solver.reset(new RadSolve(parent, level, grids, dmap));
    }
#endif

    // If we want, we can restart the checkpoint at a new time.

    if (reset_checkpoint_time > -1.e199) {

        if (!parent->RegridOnRestart()) {
          amrex::Error("reset_checkpoint_time only makes sense when amr.regrid_on_restart=1");
        }

        const Real cur_time = get_state_data(State_Type).curTime();
        const Real prev_time = get_state_data(State_Type).prevTime();
        const Real dt = cur_time - prev_time;

        parent->setStartTime(reset_checkpoint_time);
        parent->setCumTime(reset_checkpoint_time);

        for (int n = 0; n < num_state_type; ++n) {
            StateData& cur_state = get_state_data(n);
            cur_state.setOldTimeLevel(reset_checkpoint_time-dt);
            cur_state.setNewTimeLevel(reset_checkpoint_time   );
        }

    }

    if (reset_checkpoint_step > -1) {

        if (!parent->RegridOnRestart()) {
          amrex::Error("reset_checkpoint_step only makes sense when amr.regrid_on_restart=1");
        }

        parent->setLevelSteps(level, reset_checkpoint_step);
        parent->setLevelCount(level, reset_checkpoint_step);

    }
}

void
Castro::set_state_in_checkpoint (Vector<int>& state_in_checkpoint)
{
    for (int i=0; i<num_state_type; ++i) {
        state_in_checkpoint[i] = 1;
    }
}

void
Castro::checkPoint(const std::string& dir,
                   std::ostream&  os,
                   VisMF::How     how,
                   bool /*dump_old_default*/)
{

  const Real io_start_time = ParallelDescriptor::second();

  AmrLevel::checkPoint(dir, os, how, dump_old);

  const Real io_time = ParallelDescriptor::second() - io_start_time;

#ifdef RADIATION
  if (do_radiation) {
    radiation->checkPoint(level, dir, os, how);
  }
#endif

#ifdef AMREX_PARTICLES
  ParticleCheckPoint(dir);
#endif

  if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        {
            std::ofstream CastroHeaderFile;
            std::string FullPathCastroHeaderFile = dir;
            FullPathCastroHeaderFile += "/CastroHeader";
            CastroHeaderFile.open(FullPathCastroHeaderFile.c_str(), std::ios::out);

            CastroHeaderFile << "Checkpoint version: " << current_version << std::endl;
            CastroHeaderFile.close();

            writeJobInfo(dir, io_time);

            // output the list of state variables, so we can do a sanity check on restart
            std::ofstream StateListFile;
            std::string FullPathStateList = dir;
            FullPathStateList += "/state_names.txt";
            StateListFile.open(FullPathStateList.c_str(), std::ios::out);

            for (int n = 0; n < NUM_STATE; n++) {
              StateListFile << desc_lst[State_Type].name(n) << "\n";
            }
            StateListFile.close();
        }

        // If we have limited this last timestep to hit a plot interval,
        // store the timestep we took before limiting. After the restart,
        // this dt will be read in and used to limit the next timestep
        // appropriately, rather than the shortened timestep.

        if (lastDtPlotLimited == 1) {

            std::ofstream dtHeaderFile;
            std::string FullPathdtHeaderFile = dir;
            FullPathdtHeaderFile += "/dtHeader";
            dtHeaderFile.open(FullPathdtHeaderFile.c_str(), std::ios::out);

            dtHeaderFile << lastDtBeforePlotLimiting << std::endl;
            dtHeaderFile.close();

        }

        {
            // store elapsed CPU time
            std::ofstream CPUFile;
            std::string FullPathCPUFile = dir;
            FullPathCPUFile += "/CPUtime";
            CPUFile.open(FullPathCPUFile.c_str(), std::ios::out);

            CPUFile << std::setprecision(17) << getCPUTime();
            CPUFile.close();
        }

#ifdef GRAVITY
        if (use_point_mass) {

            // store current value of the point mass
            std::ofstream PMFile;
            std::string FullPathPMFile = dir;
            FullPathPMFile += "/point_mass";
            PMFile.open(FullPathPMFile.c_str(), std::ios::out);

            PMFile << std::setprecision(17) << point_mass << std::endl;

            PMFile.close();

        }
#endif

#ifdef ROTATION
        if (do_rotation) {
            // store current value of the rotation period
            std::ofstream RotationFile;
            std::string FullPathRotationFile = dir;
            FullPathRotationFile += "/Rotation";
            RotationFile.open(FullPathRotationFile.c_str(), std::ios::out);

            RotationFile << std::scientific;
            RotationFile.precision(19);

            RotationFile << std::setw(30) << castro::rotational_period << std::endl;

            RotationFile.close();
        }
#endif

        {
            // store any problem-specific stuff
            problem_checkpoint(dir);
        }
    }

}

std::string
Castro::thePlotFileType () const
{
    //
    // Increment this whenever the writePlotFile() format changes.
    //
    static const std::string the_plot_file_type("HyperCLaw-V1.1");

    return the_plot_file_type;
}

void
Castro::setPlotVariables ()
{
  AmrLevel::setPlotVariables();

  // Don't add the Source_Type data to the plotfile, we only
  // want to store it in the checkpoints.

  for (int i = 0; i < desc_lst[Source_Type].nComp(); i++) {
    parent->deleteStatePlotVar(desc_lst[Source_Type].name(i));
  }

#ifdef SIMPLIFIED_SDC
#ifdef REACTIONS
  if (time_integration_method == SimplifiedSpectralDeferredCorrections) {
      for (int i = 0; i < desc_lst[Simplified_SDC_React_Type].nComp(); i++) {
          parent->deleteStatePlotVar(desc_lst[Simplified_SDC_React_Type].name(i));
      }
  }
#endif
#endif

}



void
Castro::writeJobInfo (const std::string& dir, const Real io_time)
{

  // job_info file with details about the run
  std::ofstream jobInfoFile;
  std::string FullPathJobInfoFile = dir;
  FullPathJobInfoFile += "/job_info";
  jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

  std::string PrettyLine = std::string(78, '=') + "\n";
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // job information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Castro Job Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "job name: " << job_name << "\n\n";
  jobInfoFile << "inputs file: " << inputs_name << "\n\n";

  jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef AMREX_USE_OMP
  jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif
  jobInfoFile << "\n";
  jobInfoFile << "hydro tile size:         " << hydro_tile_size << "\n";

  jobInfoFile << "\n";
  jobInfoFile << "CPU time used since start of simulation (CPU-hours): " <<
    getCPUTime()/3600.0;

  jobInfoFile << "\n\n";

  // plotfile information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Plotfile Information\n";
  jobInfoFile << PrettyLine;

  const std::time_t now = time(nullptr);
  jobInfoFile << "output date / time: "
              << std::put_time(std::localtime(&now), "%c\n") << "\n";

  jobInfoFile << "output dir:         " << amrex::FileSystem::CurrentPath() << "\n";

  jobInfoFile << "I/O time (s):       " << io_time << "\n";

  jobInfoFile << "\n\n";

#ifdef AMREX_USE_GPU
  // This output assumes for simplicity that every rank uses the
  // same type of GPU.

  jobInfoFile << PrettyLine;
  jobInfoFile << "GPU Information:       " << "\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "GPU model name: " << Gpu::Device::deviceName() << "\n";
  jobInfoFile << "Number of GPUs used: " << Gpu::Device::numDevicesUsed() << "\n";

  jobInfoFile << "\n\n";
#endif

  // build information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Build Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
  jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
  jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
  jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
  jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  jobInfoFile << "\n";

#ifdef AMREX_USE_CUDA
  jobInfoFile << "CUDA version:  " << buildInfoGetCUDAVersion() << "\n";
  jobInfoFile << "\n";
#endif

  jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Fortran comp:  " << buildInfoGetFName() << "\n";
  jobInfoFile << "Fortran flags: " << buildInfoGetFFlags() << "\n";

  jobInfoFile << "\n";

  jobInfoFile << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  jobInfoFile << "Libraries:     " << buildInfoGetLibraries() << "\n";

  jobInfoFile << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    jobInfoFile << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
  }

  jobInfoFile << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  const char* githash3 = buildInfoGetGitHash(3);
  if (strlen(githash1) > 0) {
    jobInfoFile << "Castro       git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    jobInfoFile << "AMReX        git describe: " << githash2 << "\n";
  }
  if (strlen(githash3) > 0) {
    jobInfoFile << "Microphysics git describe: " << githash3 << "\n";
  }

  const char* buildgithash = buildInfoGetBuildGitHash();
  const char* buildgitname = buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0){
    jobInfoFile << buildgitname << " git describe: " << buildgithash << "\n";
  }

  jobInfoFile << "\n\n";


  // grid information
  jobInfoFile << PrettyLine;
  jobInfoFile << " Grid Information\n";
  jobInfoFile << PrettyLine;

  int f_lev = parent->finestLevel();

  for (int i = 0; i <= f_lev; i++)
    {
      jobInfoFile << " level: " << i << "\n";
      jobInfoFile << "   number of boxes = " << parent->numGrids(i) << "\n";
      jobInfoFile << "   maximum zones   = ";
      for (int n = 0; n < AMREX_SPACEDIM; n++)
        {
          jobInfoFile << parent->Geom(i).Domain().length(n) << " ";
          //jobInfoFile << parent->Geom(i).ProbHi(n) << " ";
        }
      jobInfoFile << "\n\n";
    }

  jobInfoFile << " Boundary conditions\n";
  Vector<int> lo_bc_out(AMREX_SPACEDIM), hi_bc_out(AMREX_SPACEDIM);
  ParmParse pp("castro");
  pp.getarr("lo_bc",lo_bc_out,0,AMREX_SPACEDIM);
  pp.getarr("hi_bc",hi_bc_out,0,AMREX_SPACEDIM);


  // these names correspond to the integer flags setup in the
  // Castro_setup.cpp
  const char* names_bc[] =
    { "interior", "inflow", "outflow",
      "symmetry", "slipwall", "noslipwall" };


  jobInfoFile << "   -x: " << names_bc[lo_bc_out[0]] << "\n";
  jobInfoFile << "   +x: " << names_bc[hi_bc_out[0]] << "\n";
  if (AMREX_SPACEDIM >= 2) {
    jobInfoFile << "   -y: " << names_bc[lo_bc_out[1]] << "\n";
    jobInfoFile << "   +y: " << names_bc[hi_bc_out[1]] << "\n";
  }
  if (AMREX_SPACEDIM == 3) {
    jobInfoFile << "   -z: " << names_bc[lo_bc_out[2]] << "\n";
    jobInfoFile << "   +z: " << names_bc[hi_bc_out[2]] << "\n";
  }

  jobInfoFile << "\n\n";

  jobInfoFile << " Domain geometry info\n";

  jobInfoFile << "     center: " << problem::center[0] << " , " << problem::center[1] << " , " << problem::center[2] << "\n";
  jobInfoFile << "\n";

  jobInfoFile << "     geometry.is_periodic: ";
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    jobInfoFile << geom.isPeriodic(idir) << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     geometry.coord_sys:   " << geom.Coord() << "\n";

  jobInfoFile << "     geometry.prob_lo:     ";
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    jobInfoFile << geom.ProbLo(idir) << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     geometry.prob_hi:     ";
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    jobInfoFile << geom.ProbHi(idir) << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     amr.n_cell:           ";
  const int*  domain_lo = geom.Domain().loVect();
  const int*  domain_hi = geom.Domain().hiVect();
  for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
    jobInfoFile << domain_hi[idir] - domain_lo[idir] + 1 << " ";
  }
  jobInfoFile << "\n";

  int max_level = parent->maxLevel();
  jobInfoFile << "     amr.max_level:        " << max_level << "\n";

  jobInfoFile << "     amr.ref_ratio:        ";
  for (int lev = 1; lev <= max_level; lev++) {
    IntVect ref_ratio = parent->refRatio(lev-1);
    jobInfoFile << ref_ratio[0] << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     amr.n_error_buf:      ";
  for (int lev = 1; lev <= max_level; lev++) {
    int errbuf = parent->nErrorBuf(lev-1);
    jobInfoFile << errbuf << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     amr.regrid_int:       ";
  for (int lev = 1; lev <= max_level; lev++) {
    int regridint = parent->regridInt(lev-1);
    jobInfoFile << regridint << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     amr.blocking_factor:  ";
  for (int lev = 1; lev <= max_level; lev++) {
    IntVect bf = parent->blockingFactor(lev-1);
    jobInfoFile << bf[0] << " ";
  }
  jobInfoFile << "\n";

  jobInfoFile << "     amr.max_grid_size:    ";
  for (int lev = 1; lev <= max_level; lev++) {
    IntVect mgs = parent->maxGridSize(lev-1);
    jobInfoFile << mgs[0] << " ";
  }
  jobInfoFile << "\n\n";

  jobInfoFile << "     amr.subcycling_mode: " << parent->subcyclingMode();

  jobInfoFile << "\n\n";



  // species info
  int mlen = 20;

  jobInfoFile << PrettyLine;
  jobInfoFile << " Species Information\n";
  jobInfoFile << PrettyLine;

  jobInfoFile <<
    std::setw(6) << "index" << SkipSpace <<
    std::setw(mlen+1) << "name" << SkipSpace <<
    std::setw(7) << "A" << SkipSpace <<
    std::setw(7) << "Z" << "\n";
  jobInfoFile << OtherLine;

  for (int i = 0; i < NumSpec; i++)
    {
      jobInfoFile <<
        std::setw(6) << i << SkipSpace <<
        std::setw(mlen+1) << std::setfill(' ') << short_spec_names_cxx[i] << SkipSpace <<
        std::setw(7) << aion[i] << SkipSpace <<
        std::setw(7) << zion[i] << "\n";
    }
  jobInfoFile << "\n\n";


  // runtime parameters
  jobInfoFile << PrettyLine;
  jobInfoFile << " Inputs File Parameters\n";
  jobInfoFile << PrettyLine;

#include <castro_job_info_tests.H>
#ifdef AMREX_PARTICLES
#include <particles_job_info_tests.H>
#endif

#ifdef GRAVITY
  gravity->output_job_info_params(jobInfoFile);
#endif
#ifdef DIFFUSION
  diffusion->output_job_info_params(jobInfoFile);
#endif

#include <prob_job_info_tests.H>

#include <extern_job_info_tests.H>

  jobInfoFile.close();

}


void
Castro::writeBuildInfo ()
{
  std::string PrettyLine = std::string(78, '=') + "\n";

  // build information
  std::cout << PrettyLine;
  std::cout << " Castro Build Information\n";
  std::cout << PrettyLine;

  std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
  std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
  std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
  std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  std::cout << "\n";

  std::cout << "COMP:          " << buildInfoGetComp() << "\n";
  std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  std::cout << "\n";

  std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  std::cout << "\n";

  std::cout << "Fortran comp:  " << buildInfoGetFName() << "\n";
  std::cout << "Fortran flags: " << buildInfoGetFFlags() << "\n";

  std::cout << "\n";

  std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";

  std::cout << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n) << "\n";
  }

  std::cout << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  const char* githash3 = buildInfoGetGitHash(3);
  if (strlen(githash1) > 0) {
    std::cout << "Castro       git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    std::cout << "AMReX        git describe: " << githash2 << "\n";
  }
  if (strlen(githash3) > 0) {
    std::cout << "Microphysics git describe: " << githash3 << "\n";
  }

  const char* buildgithash = buildInfoGetBuildGitHash();
  const char* buildgitname = buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0){
    std::cout << buildgitname << " git describe: " << buildgithash << "\n";
  }

  std::cout << "\n\n";
}

void
Castro::writePlotFile(const std::string& dir,
                      ostream& os,
                      VisMF::How how)
{
  plotFileOutput(dir, os, how, 0);
}


void
Castro::writeSmallPlotFile (const std::string& dir,
                            ostream&       os,
                            VisMF::How     how)
{
  plotFileOutput(dir, os, how, 1);
}


void
Castro::plotFileOutput(const std::string& dir,
                       ostream& os,
                       VisMF::How how,
                       const int is_small)
{
    // We have nothing to do if we're on a level that
    // we are not plotting.

    if (is_small) {
        if (level > parent->smallplotMaxLevel()) {
            return;
        }
    }
    else {
        if (level > parent->plotMaxLevel()) {
            return;
        }
    }

#ifdef AMREX_PARTICLES
  ParticlePlotFile(dir);
#endif

    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++) {
      for (int comp = 0; comp < desc_lst[typ].nComp(); comp++) {
        if (((parent->isStatePlotVar(desc_lst[typ].name(comp)) && is_small == 0) ||
             (parent->isStateSmallPlotVar(desc_lst[typ].name(comp)) && is_small == 1)) &&
            desc_lst[typ].getType() == IndexType::TheCellType()) {
            plot_var_map.emplace_back(typ, comp);
        }
      }
    }

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (const auto & dd : dlist) {

        if ((parent->isDerivePlotVar(dd.name()) && is_small == 0) ||
            (parent->isDeriveSmallPlotVar(dd.name()) && is_small == 1))
        {
#ifdef AMREX_PARTICLES
            if (dd.name() == "particle_count" ||
                dd.name() == "total_particle_count")
            {
                if (Castro::theTracerPC())
                {
                    derive_names.push_back(dd.name());
                    num_derive = num_derive + dd.numDerive();
                }
            } else
#endif
            {
               derive_names.push_back(dd.name());
               num_derive = num_derive + dd.numDerive();
            }
        }
    }

    int n_data_items = static_cast<int>(plot_var_map.size()) + num_derive;

#ifdef RADIATION
    if (Radiation::nplotvar > 0) n_data_items += Radiation::nplotvar;
#endif

#ifdef REACTIONS
#ifndef TRUE_SDC
    if (store_burn_weights) {
        n_data_items += static_cast<int>(Castro::burn_weight_names.size());
    }
#endif
#endif

    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0) {
          amrex::Error("Must specify at least one valid data item to plot");
        }

        os << n_data_items << '\n';

        //
        // Names of variables -- first state, then derived
        //
        for (const auto& [typ, comp] : plot_var_map)
        {
            os << desc_lst[typ].name(comp) << '\n';
        }

        for (auto &name : derive_names)
        {
            const DeriveRec* rec = derive_lst.get(name);
            if (rec->numDerive() > 1) {
                for (int i = 0; i < rec->numDerive(); ++i) {
                    os << rec->variableName(0) + '_' + std::to_string(i) + '\n';
                }
            }
            else {
                os << rec->variableName(0) << '\n';
            }
        }

#ifdef RADIATION
        for (int i=0; i<Radiation::nplotvar; ++i) {
          os << Radiation::plotvar_names[i] << '\n';
        }
#endif

#ifdef REACTIONS
#ifndef TRUE_SDC
        if (store_burn_weights) {
            for (const auto& name: Castro::burn_weight_names) {
                os << name << '\n';
            }
        }
#endif
#endif

        os << AMREX_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        if (is_small) {
            f_lev = std::min(f_lev, parent->smallplotMaxLevel());
        }
        else {
            f_lev = std::min(f_lev, parent->plotMaxLevel());
        }

        os << f_lev << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            os << geom.ProbLo(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            os << geom.ProbHi(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i < f_lev; i++) {
          os << parent->refRatio(i)[0] << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
          os << parent->Geom(i).Domain() << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++) {
          os << parent->levelSteps(i) << ' ';
        }
        os << '\n';
        for (int i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < AMREX_SPACEDIM; k++) {
              os << parent->Geom(i).CellSize()[k] << ' ';
            }
            os << '\n';
        }
        os << (int) geom.Coord() << '\n';
        os << "0\n"; // Write bndry data.

#ifdef RADIATION
        if (do_radiation && Radiation::do_multigroup) {
          std::ofstream groupfile;
          std::string FullPathGroupFile = dir;
          FullPathGroupFile += "/RadiationGroups";
          groupfile.open(FullPathGroupFile.c_str(), std::ios::out);

          radiation->write_groups(groupfile);

          groupfile.close();
        }
#endif
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";
    char buf[64];
    sprintf(buf, "Level_%d", level);
    std::string Level = buf;
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.size()-1] != '/') {
      FullPath += '/';
    }
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor()) {
        if (!amrex::UtilCreateDirectory(FullPath, 0755)) {
            amrex::CreateDirectoryFailed(FullPath);
        }
    }
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (int i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (int n = 0; n < AMREX_SPACEDIM; n++) {
              os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
            }
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    const int nGrow = 0;
    MultiFab  plotMF(grids,dmap,n_data_items,nGrow);
    MultiFab* this_dat = nullptr;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (const auto& [typ, comp] : plot_var_map) {
        this_dat = &state[typ].newData();
        MultiFab::Copy(plotMF, *this_dat, comp, cnt, 1, nGrow);
        cnt++;
    }
    //
    // Cull data from derived variables.
    //
    if (!dlist.empty())
    {
        for (const auto & dd : dlist) {

            if ((parent->isDerivePlotVar(dd.name()) && is_small == 0) ||
                (parent->isDeriveSmallPlotVar(dd.name()) && is_small == 1)) {

                auto derive_dat = derive(dd.variableName(0), cur_time, nGrow);
                MultiFab::Copy(plotMF, *derive_dat, 0, cnt, dd.numDerive(), nGrow);
                cnt = cnt + dd.numDerive();
            }
        }
    }

#ifdef RADIATION
    if (Radiation::nplotvar > 0) {
        MultiFab::Copy(plotMF,*(radiation->plotvar[level]),0,cnt,Radiation::nplotvar,0);
        cnt += Radiation::nplotvar;
    }
#endif

#ifdef REACTIONS
#ifndef TRUE_SDC
    if (store_burn_weights) {
        MultiFab::Copy(plotMF, getLevel(level).burn_weights, 0, cnt, static_cast<int>(Castro::burn_weight_names.size()), 0);
        cnt += static_cast<int>(Castro::burn_weight_names.size());  // NOLINT(clang-analyzer-deadcode.DeadStores)
    }
#endif
#endif

    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;

    const Real io_start_time = ParallelDescriptor::second();

    if (amrex::AsyncOut::UseAsyncOut()) {
        VisMF::AsyncWrite(std::move(plotMF),TheFullPath);
    } else {
        VisMF::Write(plotMF,TheFullPath,how,true);
    }

    const Real io_time = ParallelDescriptor::second() - io_start_time;

    if (level == 0 && ParallelDescriptor::IOProcessor()) {
        writeJobInfo(dir, io_time);
    }

#ifdef GRAVITY
    if (use_point_mass && level == 0) {

        // store current value of the point mass
        std::ofstream PMFile;
        std::string FullPathPMFile = dir;
        FullPathPMFile += "/point_mass";
        PMFile.open(FullPathPMFile.c_str(), std::ios::out);

        PMFile << std::setprecision(17) << point_mass << std::endl;

        PMFile.close();

    }
#endif

}
