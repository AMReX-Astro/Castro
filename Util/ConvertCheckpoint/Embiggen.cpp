
// This is the version that reads a Checkpoint file
// and writes it out again.
// ---------------------------------------------------------------
#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <vector>
#include <string>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_Geometry.H>
#include <AMReX_StateDescriptor.H>
#include <AMReX_StateData.H>
#include <AMReX_BCRec.H>
#include <AMReX_LevelBld.H>
#include <AMReX_AmrLevel.H>

using namespace amrex;

LevelBld *getLevelBld() {
  return 0;
}

#ifdef BL_USE_SETBUF
#define pubsetbuf setbuf
#endif

using std::cout;
using std::cerr;
using std::endl;

using namespace amrex;

std::string CheckFileIn;
std::string CheckFileOut;
int nFiles(64);
bool verbose(true);
int num_new_levels(1);
int      ref_ratio(1);
int   grown_factor(1);
int star_at_center(-1);
int   max_grid_size(4096);
int   coord(-1);
const std::string CheckPointVersion = "CheckPointVersion_1.0";

Vector<int> nsets_save(1);

VisMF::How how = VisMF::OneFilePerCPU;


// ---------------------------------------------------------------
struct FakeStateData {
    struct TimeInterval {
        Real start, stop;
    };
    const StateDescriptor *desc;
    Box domain;
    BoxArray grids;
    TimeInterval new_time;
    TimeInterval old_time;
    MultiFab *new_data;
    MultiFab *old_data;
    Vector< Vector<BCRec> > bc;
};


struct FakeAmrLevel {
    int level;                        // AMR level (0 is coarsest).
    Geometry geom;                    // Geom at this level.
    BoxArray grids;                   // Cell-centered locations of grids.
    IntVect crse_ratio;               // Refinement ratio to coarser level.
    IntVect fine_ratio;               // Refinement ratio to finer level.
    Vector<FakeStateData> state;       // Array of state data.
    Vector<FakeStateData> new_state;   // Array of new state data.
};


struct FakeAmr {
  int                  finest_level;
  Real                 cumtime;
  Vector<Real>          dt_level;
  Vector<int>           level_steps;
  Vector<int>           level_count;
  Vector<int>           n_cycle;
  Vector<Real>          dt_min;
  Vector<IntVect>       ref_ratio;
  Vector<Geometry>      geom;
  Vector< Vector<Real> > dx;
  Vector<FakeAmrLevel> fakeAmrLevels;
};

FakeAmr fakeAmr;

// ---------------------------------------------------------------
static void ScanArguments() {
    ParmParse pp;

    if(pp.contains("checkin")) {
      pp.get("checkin", CheckFileIn);
    }
    if(pp.contains("checkout")) {
      pp.get("checkout", CheckFileOut);
    }
    if(pp.contains("nfiles")) {
      pp.get("nfiles", nFiles);
    }
    if(pp.contains("verbose")) {
      pp.get("verbose", verbose);
    }
    if(pp.contains("ref_ratio")) {
      pp.get("ref_ratio", ref_ratio);
    }

    if(pp.contains("grown_factor")) {
      pp.get("grown_factor", grown_factor);
    }

    pp.get("star_at_center", star_at_center);

    if (star_at_center != 0 && star_at_center != 1)
       amrex::Abort("star_at_center must be 0 or 1");

    if (ref_ratio != 2 && ref_ratio != 4)
       amrex::Abort("ref_ratio must be 2 or 4");

    if (grown_factor <= 1)
        amrex::Abort("must have grown_factor > 1");

    if (star_at_center == 1)
       if (grown_factor != 2 && grown_factor != 3)
          amrex::Abort("must have grown_factor = 2 or 3 for star at center");
}

// ---------------------------------------------------------------
static void PrintUsage (char *progName) {
    cout << "Usage: " << progName << " checkin=filename "
         << "checkout=outfilename "
         << "ref_ratio= 2 or 4 "
         << "grown_factor=integer "
         << "star_at_center =0 or 1  "
         << "[nfiles=nfilesout] "
         << "[verbose=trueorfalse]" << endl;
    exit(1);
}

// ---------------------------------------------------------------

static void ReadCheckpointFile(const std::string& fileName) {
    int i;
    std::string File = fileName;

    File += '/';
    File += "Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream is;
    is.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    is.open(File.c_str(), std::ios::in);
    if( ! is.good()) {
        amrex::FileOpenFailed(File);
    }

    //
    // Read global data.
    //
    // Attempt to differentiate between old and new CheckPointFiles.
    //
    int         spdim;
    bool        new_checkpoint_format = false;
    std::string first_line;

    std::getline(is,first_line);

    if (first_line == CheckPointVersion) {
        new_checkpoint_format = true;
        is >> spdim;
    } else {
        spdim = atoi(first_line.c_str());
    }

    if(spdim != AMREX_SPACEDIM) {
        cerr << "Amr::restart(): bad spacedim = " << spdim << '\n';
        amrex::Abort();
    }

    is >> fakeAmr.cumtime;
    int mx_lev;
    is >> mx_lev;
    is >> fakeAmr.finest_level;

    if(ParallelDescriptor::IOProcessor())
       std::cout << "previous finest_lev is " << fakeAmr.finest_level <<  std::endl;

    // ADDING LEVELS
    int n = num_new_levels;
    mx_lev = mx_lev + n;
    fakeAmr.finest_level = fakeAmr.finest_level + n;

    if(ParallelDescriptor::IOProcessor()) {
       std::cout << "     new finest_lev is " << fakeAmr.finest_level << std::endl;
       std::cout << "previous     mx_lev is " << mx_lev-n << std::endl;
       std::cout << "     new     mx_lev is " << mx_lev << std::endl;
    }

    fakeAmr.geom.resize(mx_lev + 1);
    fakeAmr.ref_ratio.resize(mx_lev);
    fakeAmr.dt_level.resize(mx_lev + 1);
    fakeAmr.dt_min.resize(mx_lev + 1);
    fakeAmr.n_cycle.resize(mx_lev + 1);
    fakeAmr.level_steps.resize(mx_lev + 1);
    fakeAmr.level_count.resize(mx_lev + 1);
    fakeAmr.fakeAmrLevels.resize(mx_lev + 1);

    // READING GEOM, REF_RATIO, DT_LEVEL
    for (i = n; i <= mx_lev; i++) {
      is >> fakeAmr.geom[i];
    }

    if(ParallelDescriptor::IOProcessor()) {
       std::cout << " " << std::endl;
       for (i = 1; i <= mx_lev; i++) {
          std::cout << "Old checkpoint level    " << i-1 << std::endl;
          std::cout << " ... domain is       " << fakeAmr.geom[i].Domain() << std::endl;
          std::cout << " ...     dx is       " << fakeAmr.geom[i].CellSize()[0] << std::endl;
          std::cout << "  " << std::endl;
       }
    }

    // Make sure current domain is divisible by 2*ref_ratio so length of coarsened domain is even
    Box dom0(fakeAmr.geom[0].Domain());
    for (int d = 0; d < AMREX_SPACEDIM; d++)
    {
      int dlen = dom0.size()[d];
      int scaled = dlen / (2*ref_ratio);
      if ( (scaled * 2 * ref_ratio) != dlen )
        amrex::Abort("must have domain divisible by 2*ref_ratio");
    }

    if (grown_factor <= 1)
        amrex::Abort("must have grown_factor > 1");

    for (i = 1; i <  mx_lev; i++) {
      is >> fakeAmr.ref_ratio[i];
    }
    for (i = 1; i <= mx_lev; i++) {
      is >> fakeAmr.dt_level[i];
    }

    Box          domain(fakeAmr.geom[1].Domain());
    RealBox prob_domain(fakeAmr.geom[1].ProbDomain());
    coord = fakeAmr.geom[1].Coord();

    // Define domain for new levels
    domain.coarsen(ref_ratio);
    fakeAmr.geom[0].define(domain,&prob_domain,coord);

    // Define ref_ratio for new levels
    fakeAmr.ref_ratio[0] = ref_ratio * IntVect::TheUnitVector();

    // Define dt_level for new levels
    fakeAmr.dt_level[0] = fakeAmr.dt_level[1] * ref_ratio;

    if (new_checkpoint_format) {
      for (i = 1; i <= mx_lev; i++) is >> fakeAmr.dt_min[i];
      fakeAmr.dt_min[0] = fakeAmr.dt_min[1] * ref_ratio;
    } else {
      for (i = 0; i <= mx_lev; i++) fakeAmr.dt_min[i] = fakeAmr.dt_level[i];
    }

    // READING N_CYCLE, LEVEL_STEPS, LEVEL_COUNT
    for (i = 1; i <= mx_lev; i++) {
      is >> fakeAmr.n_cycle[i];
    }

    for (i = 1; i <= mx_lev; i++) {
      is >> fakeAmr.level_steps[i];
    }
    for (i = 1; i <= mx_lev; i++) {
      is >> fakeAmr.level_count[i];
    }

    // ADDING LEVELS

    // n_cycle is always equal to 1 at the coarsest level
    fakeAmr.n_cycle[0] = 1;

    // At the old coarsest level, which is now level 1, we set n_cycle to ref_ratio
    fakeAmr.n_cycle[1] = ref_ratio;

    fakeAmr.level_steps[0] = fakeAmr.level_steps[1] / ref_ratio;
    if ( (fakeAmr.level_steps[0]*ref_ratio) != fakeAmr.level_steps[1] )
       amrex::Abort("Number of steps in original checkpoint must be divisible by ref_ratio");

    // level_count is how many steps we've taken at this level since the last regrid
    if (fakeAmr.level_count[1] == fakeAmr.level_steps[1])
    {
       fakeAmr.level_count[0] = fakeAmr.level_steps[0];

    // this is actually wrong but should work for now
    } else {
       fakeAmr.level_count[0] = std::min(fakeAmr.level_count[1],fakeAmr.level_steps[0]);;
    }

    int ndesc_save;

    // READ LEVEL DATA
    for(int lev(1); lev <= fakeAmr.finest_level; ++lev) {

      FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[lev];

      is >> falRef.level;
      falRef.level = falRef.level + n;

      is >> falRef.geom;

      falRef.fine_ratio = IntVect::TheUnitVector();
      falRef.fine_ratio.scale(-1);
      falRef.crse_ratio = IntVect::TheUnitVector();
      falRef.crse_ratio.scale(-1);

      if(falRef.level > 0) {
        falRef.crse_ratio = fakeAmr.ref_ratio[falRef.level-1];
      }
      if(falRef.level < mx_lev) {
        falRef.fine_ratio = fakeAmr.ref_ratio[falRef.level];
      }

      falRef.grids.readFrom(is);

      int nstate;
      is >> nstate;
      int ndesc = nstate;

      // This should be the same at all levels
      ndesc_save = ndesc;

      // ndesc depends on which descriptor so we store a value for each
      if (lev == 1) nsets_save.resize(ndesc_save);

      falRef.state.resize(ndesc);
      falRef.new_state.resize(ndesc);

      for(int i = 0; i < ndesc; i++) {
        // ******* StateDescriptor::restart

        is >> falRef.state[i].domain;

        falRef.state[i].grids.readFrom(is);

        is >> falRef.state[i].old_time.start;
        is >> falRef.state[i].old_time.stop;
        is >> falRef.state[i].new_time.start;
        is >> falRef.state[i].new_time.stop;

        int nsets;
        is >> nsets;

        nsets_save[i] = nsets;

        falRef.state[i].old_data = 0;
        falRef.state[i].new_data = 0;

        std::string mf_name;
        std::string FullPathName;

        // This reads the "new" data, if it's there
        if (nsets >= 1) {
           falRef.state[i].new_data = new MultiFab;
           is >> mf_name;
           // Note that mf_name is relative to the Header file.
           // We need to prepend the name of the fileName directory.
           FullPathName = fileName;
           if( ! fileName.empty() && fileName[fileName.length()-1] != '/') {
             FullPathName += '/';
           }
           FullPathName += mf_name;
           VisMF::Read(*(falRef.state[i].new_data), FullPathName);
        }

        // This reads the "old" data, if it's there
        if (nsets == 2) {
          falRef.state[i].old_data = new MultiFab;
          is >> mf_name;
          // Note that mf_name is relative to the Header file.
          // We need to prepend the name of the fileName directory.
          FullPathName = fileName;
          if( ! fileName.empty() && fileName[fileName.length()-1] != '/') {
            FullPathName += '/';
          }
          FullPathName += mf_name;
          VisMF::Read(*(falRef.state[i].old_data), FullPathName);
        }

      }
    }

    FakeAmrLevel &falRef_orig = fakeAmr.fakeAmrLevels[n];

    // Compute the effective max_grid_size
    int max_len = 0;
    BoxArray g(fakeAmr.fakeAmrLevels[n].grids);
    for (int b = 0; b < g.size(); b++)
       for (int d = 0; d < AMREX_SPACEDIM; d++)
         max_len = std::max(max_len, g[b].length(d));
    max_grid_size = max_len;

    // Add new level data
    for(int lev(n-1); lev >= 0; lev--) {
      FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[lev];
      falRef.level = lev;

      // This version breaks up the new coarser domain based on the computed max_grid_size
      BoxArray new_grids(domain);
      new_grids.maxSize(max_grid_size);

      falRef.grids = new_grids;

      falRef.geom.define(domain,&prob_domain,coord);

      if(falRef.level > 0)
        falRef.crse_ratio = ref_ratio * IntVect::TheUnitVector();
      falRef.fine_ratio = ref_ratio * IntVect::TheUnitVector();

      falRef.state.resize(ndesc_save);
      falRef.new_state.resize(ndesc_save);

      for(int i = 0; i < ndesc_save; i++) {

        falRef.state[i].domain = domain;
        falRef.state[i].grids = falRef.grids;
        falRef.state[i].new_time.start = falRef_orig.state[i].new_time.start;
        falRef.state[i].new_time.stop  = falRef_orig.state[i].new_time.stop;
        falRef.state[i].old_time.start = falRef.state[i].new_time.start - fakeAmr.dt_level[lev];
        falRef.state[i].old_time.stop  = falRef.state[i].new_time.stop  - fakeAmr.dt_level[lev];

        falRef.state[i].old_data = 0;
        falRef.state[i].new_data = 0;

        if (nsets_save[i] >= 1) {

           int ncomp = falRef_orig.state[i].new_data->nComp();
           int ngrow = falRef_orig.state[i].new_data->nGrow();

           DistributionMapping dmap {falRef.grids};
           falRef.state[i].new_data = new MultiFab(falRef.grids, dmap, ncomp, ngrow);
           falRef.state[i].new_data->setVal(0.);

           if (nsets_save[i] == 2) {
             falRef.state[i].old_data = new MultiFab(falRef.grids, dmap, ncomp, ngrow);
             falRef.state[i].old_data->setVal(0.);
           }

        }
      }
    }
}

// ---------------------------------------------------------------
static void WriteCheckpointFile(const std::string& inFileName, const std::string &outFileName) {
    VisMF::SetNOutFiles(nFiles);
    // In checkpoint files always write out FABs in NATIVE format.
    FABio::Format thePrevFormat = FArrayBox::getFormat();
    FArrayBox::setFormat(FABio::FAB_NATIVE);

    const std::string ckfile = outFileName;

    // Only the I/O processor makes the directory if it doesn't already exist.
    if(ParallelDescriptor::IOProcessor()) {
      if( ! amrex::UtilCreateDirectory(ckfile, 0755)) {
        amrex::CreateDirectoryFailed(ckfile);
      }
    }
    // Force other processors to wait till directory is built.
    ParallelDescriptor::Barrier();

    // Copy some standard auxiliary files to the new checkpoint.
    // Only one processor (the I/O processor) needs to do this.

    if (ParallelDescriptor::IOProcessor()) {

      std::ifstream oldCastroHeaderFile;

      std::string oldCastroHeaderName = inFileName + "/CastroHeader";
      oldCastroHeaderFile.open(oldCastroHeaderName.c_str(), std::ios::binary);

      if (oldCastroHeaderFile.good()) {

        std::ofstream newCastroHeaderFile;

        std::string newCastroHeaderName = outFileName + "/CastroHeader";
        newCastroHeaderFile.open(newCastroHeaderName.c_str(), std::ios::binary);

        if (newCastroHeaderFile.good()) {
          newCastroHeaderFile << oldCastroHeaderFile.rdbuf();
          newCastroHeaderFile.close();
        }

        oldCastroHeaderFile.close();

      }

    }

    // Write the main header file.

    std::string HeaderFileName = ckfile + "/Header";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    int old_prec(0), i;

    if(ParallelDescriptor::IOProcessor()) {
        // Only the IOProcessor() writes to the header file.
        HeaderFile.open(HeaderFileName.c_str(),
                        std::ios::out|std::ios::trunc|std::ios::binary);

        if( ! HeaderFile.good()) {
          amrex::FileOpenFailed(HeaderFileName);
        }

        old_prec = HeaderFile.precision(15);

        int max_level(fakeAmr.finest_level);
        HeaderFile << CheckPointVersion << '\n'
                   << AMREX_SPACEDIM       << '\n'
                   << fakeAmr.cumtime           << '\n'
                   << max_level                 << '\n'
                   << fakeAmr.finest_level      << '\n';
        //
        // Write out problem domain.
        //
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr.geom[i]        << ' ';
        HeaderFile << '\n';
        for (i = 0; i < max_level; i++)  HeaderFile << fakeAmr.ref_ratio[i]   << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr.dt_level[i]    << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr.dt_min[i]      << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr.n_cycle[i]     << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr.level_steps[i] << ' ';
        HeaderFile << '\n';
        for (i = 0; i <= max_level; i++) HeaderFile << fakeAmr.level_count[i] << ' ';
        HeaderFile << '\n';
    }

    for(int lev(0); lev <= fakeAmr.finest_level; ++lev) {

      std::ostream &os = HeaderFile;
      FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[lev];
      int ndesc = falRef.state.size();

      // Build directory to hold the MultiFabs in the StateData at this level.
      char buf[64];
      sprintf(buf, "Level_%d", lev);
      std::string Level = buf;

      // Now for the full pathname of that directory.
      std::string FullPath = ckfile;
      if( ! FullPath.empty() && FullPath[FullPath.length()-1] != '/') {
        FullPath += '/';
      }
      FullPath += Level;

      // Only the I/O processor makes the directory if it doesn't already exist.
      if(ParallelDescriptor::IOProcessor()) { if( ! amrex::UtilCreateDirectory(FullPath, 0755)) {
          amrex::CreateDirectoryFailed(FullPath);
        }
      }
      // Force other processors to wait till directory is built.
      ParallelDescriptor::Barrier();

      if(ParallelDescriptor::IOProcessor()) {
        os << lev << '\n' << falRef.geom  << '\n';
        falRef.grids.writeOn(os);
        os << ndesc << '\n';
      }
      //
      // Output state data.
      //
      for(int i(0); i < ndesc; ++i) {
        //
        // Now build the full relative pathname of the StateData.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        // There is only one MultiFab written out at each level in HyperCLaw.
        //
        std::string PathNameInHeader = Level;
        sprintf(buf, "/SD_%d", i);
        PathNameInHeader += buf;
        std::string FullPathName = FullPath;
        FullPathName += buf;
        // ++++++++++++ state[i].checkPoint(PathNameInHeader, FullPathName, os, how);
          static const std::string NewSuffix("_New_MF");
          static const std::string OldSuffix("_Old_MF");
          const std::string name(PathNameInHeader);
          const std::string fullpathname(FullPathName);

          bool dump_old(true);
          if(dump_old == true && falRef.state[i].old_data == 0) {
            dump_old = false;
          }

          if(ParallelDescriptor::IOProcessor()) {
            // The relative name gets written to the Header file.
            std::string mf_name_old = name;
            mf_name_old += OldSuffix;
            std::string mf_name_new = name;
            mf_name_new += NewSuffix;

            os << falRef.state[i].domain << '\n';

            falRef.state[i].grids.writeOn(os);

            os << falRef.state[i].old_time.start << '\n'
               << falRef.state[i].old_time.stop  << '\n'
               << falRef.state[i].new_time.start << '\n'
               << falRef.state[i].new_time.stop  << '\n';

            if (nsets_save[i] > 0) {
               if(dump_old) {
                 os << 2 << '\n' << mf_name_new << '\n' << mf_name_old << '\n';
               } else {
                 os << 1 << '\n' << mf_name_new << '\n';
               }
            } else {
              os << 0 << '\n';
            }

          }

          if (nsets_save[i] > 0) {
             BL_ASSERT(falRef.state[i].new_data);
             std::string mf_fullpath_new = fullpathname;
             mf_fullpath_new += NewSuffix;
             VisMF::Write(*(falRef.state[i].new_data),mf_fullpath_new,how);
          }

          if (nsets_save[i] > 1) {
            BL_ASSERT(dump_old);
            BL_ASSERT(falRef.state[i].old_data);
            std::string mf_fullpath_old = fullpathname;
            mf_fullpath_old += OldSuffix;
            VisMF::Write(*(falRef.state[i].old_data),mf_fullpath_old,how);
          }
          // ++++++++++++
      }
      // ========================

       if (ParallelDescriptor::IOProcessor()) {
          if (lev == 0) {
             std::cout << " " << std::endl;
             std::cout << " **************************************** " << std::endl;
             std::cout << " " << std::endl;
          }
          std::cout << "New checkpoint level    " << lev << std::endl;
          std::cout << " ... domain is       " << fakeAmr.geom[lev].Domain() << std::endl;
          std::cout << " ...     dx is       " << fakeAmr.geom[lev].CellSize()[0] << std::endl;
          std::cout << "  " << std::endl;
       }

    }

    if(ParallelDescriptor::IOProcessor()) {
        HeaderFile.precision(old_prec);

        if( ! HeaderFile.good()) {
          amrex::Error("Amr::checkpoint() failed");
        }
    }

    FArrayBox::setFormat(thePrevFormat);
}

// ---------------------------------------------------------------

static void ConvertData() {

   int max_level(fakeAmr.finest_level);

   // Modify the "grids" boxarray *only* at level 0
   // NOTE: must do this before re-defining domain below.

   FakeAmrLevel &falRef0 = fakeAmr.fakeAmrLevels[0];
   Box domain(fakeAmr.geom[0].Domain());
   BoxList(newgrid_list);

   int dlenx = domain.size()[0];
#if (AMREX_SPACEDIM >= 2)
   int dleny = domain.size()[1];
#if (AMREX_SPACEDIM == 3)
   int dlenz = domain.size()[1];
#endif
#endif

   // We treat the r-z case with the star in the middle specially
#if (AMREX_SPACEDIM == 2)
   if (coord == 1 && star_at_center == 1)
   {
     for (int jy = 0; jy < grown_factor; jy++)
      for (int jx = 0; jx < grown_factor; jx++)
        for (int n = 0; n < falRef0.grids.size(); n++)
        {
          int shiftx(jx*dlenx);
          int shifty(jy*dleny);
          IntVect shift_box(AMREX_D_DECL(shiftx,shifty,shiftz));

          Box bx(falRef0.grids[n]);
          bx.shift(shift_box);
          newgrid_list.push_back(bx);
        }

      Box new_domain(fakeAmr.geom[0].Domain());
      new_domain.refine(grown_factor);
      newgrid_list.intersect(new_domain);

      newgrid_list.maxSize(dlenx);

   } else {
#endif

   if (star_at_center == 1 && grown_factor == 2)
   {
   // Here we tile the domain with tiles smaller than the original domain --
   //   we first tile with domain-sized pieces, then intersect with the new domain

#if (AMREX_SPACEDIM == 3)
   for (int jz = 0; jz < 3; jz++)
#endif
#if (AMREX_SPACEDIM >= 2)
   for (int jy = 0; jy < 3; jy++)
#endif
    for (int jx = 0; jx < 3; jx++)
      for (int n = 0; n < falRef0.grids.size(); n++) {
        int shiftx(jx*dlenx - dlenx/2);
#if (AMREX_SPACEDIM >= 2)
        int shifty(jy*dleny - dleny/2);
#if (AMREX_SPACEDIM == 3)
        int shiftz(jz*dlenz - dlenz/2);
#endif
#endif
        IntVect shift_box(AMREX_D_DECL(shiftx,shifty,shiftz));

        Box bx(falRef0.grids[n]);
        bx.shift(shift_box);
        newgrid_list.push_back(bx);
      }

      Box new_domain(fakeAmr.geom[0].Domain());
      new_domain.refine(grown_factor);
      newgrid_list.intersect(new_domain);

   } else {
   // Here we tile the domain with tiles the size of the original domain

#if (AMREX_SPACEDIM == 3)
   for (int jz = 0; jz < grown_factor; jz++)
#endif
#if (AMREX_SPACEDIM >= 2)
   for (int jy = 0; jy < grown_factor; jy++)
#endif
    for (int jx = 0; jx < grown_factor; jx++)
      for (int n = 0; n < falRef0.grids.size(); n++)
      {
        int shiftx(jx*dlenx);
#if (AMREX_SPACEDIM >= 2)
        int shifty(jy*dleny);
#if (AMREX_SPACEDIM == 3)
        int shiftz(jz*dlenz);
#endif
#endif
        IntVect shift_box(AMREX_D_DECL(shiftx,shifty,shiftz));

        Box bx(falRef0.grids[n]);
        bx.shift(shift_box);
        newgrid_list.push_back(bx);
      }

   }  // end of star_at_center test

#if (AMREX_SPACEDIM == 2)
   }  // end of r-z test
#endif

   BoxArray newgrids(newgrid_list);
   falRef0.grids = newgrids;

   int nstatetypes = falRef0.state.size();

   for (int n = 0; n < nstatetypes; n++)
      falRef0.state[n].grids = newgrids;

   // Enlarge the ProbDomain (RealBox) of the geom at each level --
   //   but we only have to do this at level 0 because they are
   //   actually all the same copy
   RealBox rb(fakeAmr.geom[0].ProbDomain());

   // If this is an octant then we always grow only in the high directions
   if (star_at_center == 0)
   {
      // Here we grow only prob_hi, extending the domain in one direction.
      // This works when the star's center is at the origin
      for (int dm = 0; dm < AMREX_SPACEDIM; dm++)
         rb.setHi(dm,grown_factor*rb.hi(dm));
   }

   // We treat the r-z case with the star in the middle specially
#if (AMREX_SPACEDIM == 2)
   else if (coord == 1)
   {
      // Here we grow only prob_hi in the r-direction, but both prob_hi
      //   and prob_lo in the z-direction.
     int dm = 0;
     rb.setHi(dm,grown_factor*rb.hi(dm));

     dm = 1;
     Real dist   = 0.5 * (rb.hi(dm)-rb.lo(dm));
     Real center = 0.5 * (rb.hi(dm)+rb.lo(dm));
     Real newlo = center - grown_factor * dist;
     Real newhi = center + grown_factor * dist;
     rb.setLo(dm,newlo);
     rb.setHi(dm,newhi);
   }
#endif

   // This has star_at_center = 0
   else
   {
      // Here we grow prob_lo and prob_hi, extending the domain in all directions.
      // This works when the star's center is at the center of the domain.
      for (int dm = 0; dm < AMREX_SPACEDIM; dm++)
      {
         Real dist   = 0.5 * (rb.hi(dm)-rb.lo(dm));
         Real center = 0.5 * (rb.hi(dm)+rb.lo(dm));
         Real newlo = center - grown_factor * dist;
         Real newhi = center + grown_factor * dist;
         rb.setLo(dm,newlo);
         rb.setHi(dm,newhi);
      }
   }

   Geometry::ResetDefaultProbDomain(rb);
   for (auto& g : fakeAmr.geom) {
       g.ProbDomain(rb);
   }
   for (auto& l : fakeAmr.fakeAmrLevels) {
       l.geom.ProbDomain(rb);
   }

   const Real* xlo = fakeAmr.geom[0].ProbLo();

   // This sets the CoordSys member "offset" which should be identical to Geometry's problo
   //    but isn't automatically set.
   fakeAmr.geom[0].SetOffset(xlo);

   IntVect shift_iv[max_level+1];

   // Define the shift IntVect for later
   if (star_at_center == 1)
   {
      if (coord == 1) // r-z
      {
         for (int i = 0; i <= max_level; i++)
         {
            Box domain(fakeAmr.geom[i].Domain());
            // We only handle grown_factor = 2
            shift_iv[i][0] = 0;
            shift_iv[i][1] = domain.size()[1] / 2;
         }
      } else if (coord == 0) { // x-y
         for (int i = 0; i <= max_level; i++)
         {
            Box domain(fakeAmr.geom[i].Domain());
            if (grown_factor == 3) {
               shift_iv[i] = domain.size();
            } else if (grown_factor == 2) {
               shift_iv[i] = domain.size() / 2;
            }
         }
      }
   }

   // Enlarge the Domain (Box) of the geom at each level
   for (int i = 0; i <= max_level; i++)
   {
      Box domain(fakeAmr.geom[i].Domain());
      domain.refine(grown_factor);
      fakeAmr.geom[i].Domain(domain);

      FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[i];
      falRef.geom.Domain(domain);
   }

   // Now fix the state data domain
   for (int i = 0; i <= max_level; i++)
   {
      FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[i];
      for (int n = 0; n < nstatetypes; n++)
         falRef.state[n].domain.refine(grown_factor);
   }

   DistributionMapping newdm {newgrids};

   // We need to allocate a MultiFab for new data but don't need to fill it
   for (int n = 0; n < nstatetypes; n++)
   {
      if (falRef0.state[n].new_data != 0) {
         int ncomps = (falRef0.state[n].new_data)->nComp();
         MultiFab * newNewData = new MultiFab(newgrids,newdm,ncomps,1);

         newNewData->setVal(0.);

         if (star_at_center == 1)
            (falRef0.state[n].new_data)->shift(shift_iv[0]);

         falRef0.state[n].new_data = newNewData;

//       newNewData->copy(*(falRef0.state[n].new_data),0,0,ncomps);
      }
   }

   // If we have old_data as well as new_data
   for (int n = 0; n < nstatetypes; n++)
   {
      if (falRef0.state[n].old_data != 0) {
         int ncomps = (falRef0.state[n].old_data)->nComp();
         MultiFab * newOldData = new MultiFab(newgrids,newdm,ncomps,1);
         newOldData->setVal(0.);

         if (star_at_center == 1)
            (falRef0.state[n].old_data)->shift(shift_iv[0]);

         falRef0.state[n].old_data = newOldData;

//       newOldData->copy(*(falRef0.state[n].old_data),0,0,ncomps);
      }
   }

   // Now shift the data at the higher levels
   if (star_at_center == 1) {
      for (int i = 1; i <= max_level; i++)
      {
         FakeAmrLevel &falRef = fakeAmr.fakeAmrLevels[i];

         // Shift the grids associated with each level
         falRef.grids.shift(shift_iv[i]);

         for (int n = 0; n < nstatetypes; n++)
         {
            // Shift the grids associated with each StateData
            falRef.state[n].grids.shift(shift_iv[i]);

            // Shift the grids associated with the MultiFab in each StateData
            if (falRef0.state[n].new_data != 0)
               (falRef.state[n].new_data)->shift(shift_iv[i]);
            if (falRef0.state[n].old_data != 0)
               (falRef.state[n].old_data)->shift(shift_iv[i]);
         }
      }
   }
}


// ---------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc,argv);

    if(argc < 3) {
      PrintUsage(argv[0]);
    }

    ScanArguments();

    if(verbose && ParallelDescriptor::IOProcessor()) {
      cout << " " << std::endl;
      cout << "Reading from old checkpoint file: " <<  CheckFileIn << endl;
      cout << " " << std::endl;
    }

    if(verbose && ParallelDescriptor::IOProcessor()) {
      if (star_at_center == 0) cout << "Star at corner " << endl;
      if (star_at_center == 1) cout << "Star at center " << endl;
      cout << " " << std::endl;
    }

    // Read in the original checkpoint directory and add a coarser level covering the same domain
    ReadCheckpointFile(CheckFileIn);

    // Enlarge the new level 0
    ConvertData();

    // Write out the new checkpoint directory
    WriteCheckpointFile(CheckFileIn, CheckFileOut);

    if(verbose && ParallelDescriptor::IOProcessor()) {
      cout << " " << std::endl;
      cout << "Finished writing to new checkpoint file: " <<  CheckFileOut << endl;
      cout << " " << std::endl;
    }

    amrex::Finalize();
}
// ---------------------------------------------------------------
// ---------------------------------------------------------------
