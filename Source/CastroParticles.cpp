#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "Castro.H"
#include "Castro_F.H"

#ifdef PARTICLES

TracerParticleContainer* Castro::TracerPC =  0;
int Castro::do_tracer_particles           =  0;
int Castro::particle_verbose              =  1;

namespace {
    std::string       particle_init_file;
    std::string       particle_restart_file;
    int               restart_from_nonparticle_chkfile = 0;
    std::string       particle_output_file;
    std::string       timestamp_dir;
    std::vector<int>  timestamp_indices;
    //
    const std::string chk_tracer_particle_file("Tracer");
}

void 
Castro::read_particle_params ()
{
    ParmParse pp("castro");
    
    pp.query("do_tracer_particles",do_tracer_particles);

    //
    // Control the verbosity of the Particle class
    //
    ParmParse ppp("particles");

    ppp.query("v",particle_verbose);
    //
    // Used in initData() on startup to read in a file of particles.
    //
    ppp.query("particle_init_file", particle_init_file);
    //
    // Used in post_restart() to read in a file of particles.
    //
    ppp.query("particle_restart_file", particle_init_file);
    //
    // This must be true the first time you try to restart from a checkpoint
    // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
    // the particle checkpoint stuff (even if there are no active particles).
    // Otherwise the code will fail when trying to read the checkpointed particles.
    //
    ppp.query("restart_from_nonparticle_chkfile", restart_from_nonparticle_chkfile);
    //
    // Used in post_restart() to write out the file of particles.
    //
    ppp.query("particle_output_file", particle_output_file);
    //
    // The directory in which to store timestamp files.
    //
    ppp.query("timestamp_dir", timestamp_dir);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(timestamp_dir, 0755))
            BoxLib::CreateDirectoryFailed(timestamp_dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
}

void
Castro::init_particles ()
{
    BL_PROFILE("Castro::init_particles()");

    if (level > 0)
        return;

    if (do_tracer_particles)
    {
	BL_ASSERT(TracerPC == 0);
	
	TracerPC = new TracerParticleContainer(parent);
	
	TracerPC->SetVerbose(particle_verbose);
	
	if (! particle_init_file.empty())
	{
	    TracerPC->InitFromAsciiFile(particle_init_file,0);
	}
    }
}

void
Castro::ParticleCheckPoint(const std::string& dir)
{
    if (level == 0)
    {
        if (TracerPC)
            TracerPC->Checkpoint(dir, chk_tracer_particle_file);
    }
}

void
Castro::ParticlePostRestart (const std::string& restart_file)
{
    if (level == 0)
    {
        if (do_tracer_particles)
        {
            BL_ASSERT(TracerPC == 0);

            TracerPC = new TracerParticleContainer(parent);

            TracerPC->SetVerbose(particle_verbose);
	    //
	    // We want to be able to add new particles on a restart.
	    // As well as the ability to write the particles out to an ascii file.
	    //
	    if (!restart_from_nonparticle_chkfile)
	    {
		TracerPC->Restart(parent->theRestartFile(), chk_tracer_particle_file);
	    }

	    if (!particle_restart_file.empty())
	    {
		TracerPC->InitFromAsciiFile(particle_restart_file,0);
	    }
	    
	    if (!particle_output_file.empty())
	    {
		TracerPC->WriteAsciiFile(particle_output_file);
	    }
        }
    }
}

MultiFab*
Castro::ParticleDerive(const std::string& name,
                       Real               time,
                       int                ngrow)
{
    BL_PROFILE("Castro::ParticleDerive()");

  if (TracerPC && name == "particle_count")
  {
      MultiFab* derive_dat = new MultiFab(grids,1,0);
      MultiFab    temp_dat(grids,1,0);
      temp_dat.setVal(0);
      TracerPC->Increment(temp_dat,level);
      MultiFab::Copy(*derive_dat,temp_dat,0,0,1,0);
      return derive_dat;
  }
  else if (TracerPC && name == "total_particle_count")
  {
      //
      // We want the total particle count at this level or higher.
      //
      MultiFab* derive_dat = ParticleDerive("particle_count",time,ngrow);

      IntVect trr(D_DECL(1,1,1));

      for (int lev = level+1; lev <= parent->finestLevel(); lev++)
      {
          BoxArray ba = parent->boxArray(lev);

          MultiFab temp_dat(ba,1,0);

          trr *= parent->refRatio(lev-1);

          ba.coarsen(trr);

          MultiFab ctemp_dat(ba,1,0);

          temp_dat.setVal(0);
          ctemp_dat.setVal(0);

          TracerPC->Increment(temp_dat,lev);

          for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
          {
              const FArrayBox& ffab =  temp_dat[mfi];
              FArrayBox&       cfab = ctemp_dat[mfi];
              const Box&       fbx  = ffab.box();

              BL_ASSERT(cfab.box() == BoxLib::coarsen(fbx,trr));

              for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
              {
                  const Real val = ffab(p);
                  if (val > 0)
                      cfab(BoxLib::coarsen(p,trr)) += val;
              }
          }

          temp_dat.clear();

          MultiFab dat(grids,1,0);
          dat.setVal(0);
          dat.copy(ctemp_dat);

          MultiFab::Add(*derive_dat,dat,0,0,1,0);
      }

      return derive_dat;
  }
  else
  {
     return AmrLevel::derive(name,time,ngrow);
  }
}

void
Castro::TimestampParticles (int ngrow)
{
    static bool first = true;
    static int imax = -1;
    if (first)
    {
	first = false;

	ParmParse ppp("particles");

	// have to do it here, not in read_particle_params, because Density, ..., are set after
	// read_particle_params is called.

	int timestamp_density = 1;
	ppp.query("timestamp_density", timestamp_density);
	if (timestamp_density) {
	    timestamp_indices.push_back(Density);
	    std::cout << "Density = " << Density << std::endl;
	}

	int timestamp_temperature = 0;
	ppp.query("timestamp_temperature", timestamp_temperature);
	if (timestamp_temperature) {
	    timestamp_indices.push_back(Temp);
	    std::cout << "Temp = " << Temp << std::endl;
	}
	
	if (!timestamp_indices.empty()) {
	    imax = *(std::max_element(timestamp_indices.begin(), timestamp_indices.end()));
	}
    }

    if ( TracerPC && !timestamp_dir.empty())
    {
	std::string basename = timestamp_dir;
		
	if (basename[basename.length()-1] != '/') basename += '/';
	
	basename += "Timestamp";
	
	int finest_level = parent->finestLevel();
	Real time        = state[State_Type].curTime();

	for (int lev = level; lev <= finest_level; lev++)
	{
	    MultiFab& S_new = parent->getLevel(lev).get_new_data(State_Type);

	    if (imax >= 0) {  // FillPatchIterator will fail otherwise
		int ng = (lev == level) ? ngrow : 1;
		FillPatchIterator fpi(parent->getLevel(lev), S_new, 
				      ng, time, State_Type, 0, imax+1);
		const MultiFab& S = fpi.get_mf();
		TracerPC->Timestamp(basename, S    , lev, time, timestamp_indices);
	    } else {
		TracerPC->Timestamp(basename, S_new, lev, time, timestamp_indices);
	    }
	}
    }	
}

#endif
