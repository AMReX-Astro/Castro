#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include <Castro.H>
#include <Castro_F.H>

#include <particles_params.H>

#include <particles_declares.H>

using namespace amrex;

#ifdef AMREX_PARTICLES

AmrTracerParticleContainer* Castro::TracerPC =  0;

namespace {
    std::vector<int>  timestamp_indices;
    //
    const std::string chk_tracer_particle_file("Tracer");
}

void
Castro::read_particle_params ()
{

  ParmParse pp("particles");

    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(particles::timestamp_dir, 0755))
            amrex::CreateDirectoryFailed(particles::timestamp_dir);
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

        TracerPC = new AmrTracerParticleContainer(parent);

        TracerPC->SetVerbose(particles::particle_verbose);

        if (! particles::particle_init_file.empty())
        {
            TracerPC->InitFromAsciiFile(particles::particle_init_file,0);
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
Castro::ParticlePlotFile(const std::string& dir)
{
    if (level == 0)
    {
      //  We call TracerPC->Checkpoint instead of TracerPC->WritePlotFile
      //  so that the particle ids also get written out.
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

            TracerPC = new AmrTracerParticleContainer(parent);

            TracerPC->SetVerbose(particles::particle_verbose);
            //
            // We want to be able to add new particles on a restart.
            // As well as the ability to write the particles out to an ascii file.
            //
            if (!particles::restart_from_nonparticle_chkfile)
            {
                TracerPC->Restart(parent->theRestartFile(), chk_tracer_particle_file);
            }

            if (!particles::particle_restart_file.empty())
            {
                TracerPC->InitFromAsciiFile(particles::particle_restart_file,0);
            }

            if (!particles::particle_output_file.empty())
            {
                TracerPC->WriteAsciiFile(particles::particle_output_file);
            }
        }
    }
}

std::unique_ptr<MultiFab>
Castro::ParticleDerive(const std::string& name,
                       Real               time,
                       int                ngrow)
{
    BL_PROFILE("Castro::ParticleDerive()");

  if (TracerPC && name == "particle_count")
  {
      auto derive_dat = new MultiFab(grids,dmap,1,0);
      MultiFab    temp_dat(grids,dmap,1,0);
      temp_dat.setVal(0);
      TracerPC->Increment(temp_dat,level);
      MultiFab::Copy(*derive_dat,temp_dat,0,0,1,0);
      return std::unique_ptr<MultiFab>(derive_dat);
  }
  else if (TracerPC && name == "total_particle_count")
  {
      //
      // We want the total particle count at this level or higher.
      //
      auto derive_dat = ParticleDerive("particle_count",time,ngrow);

      IntVect trr(D_DECL(1,1,1));

      for (int lev = level+1; lev <= parent->finestLevel(); lev++)
      {
          BoxArray ba = parent->boxArray(lev);
          const DistributionMapping& dm = parent->DistributionMap(lev);

          MultiFab temp_dat(ba,dm,1,0);

          trr *= parent->refRatio(lev-1);

          ba.coarsen(trr);

          MultiFab ctemp_dat(ba,dm,1,0);

          temp_dat.setVal(0);
          ctemp_dat.setVal(0);

          TracerPC->Increment(temp_dat,lev);

          for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
          {
              const FArrayBox& ffab =  temp_dat[mfi];
              FArrayBox&       cfab = ctemp_dat[mfi];
              const Box&       fbx  = ffab.box();

              BL_ASSERT(cfab.box() == amrex::coarsen(fbx,trr));

              for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
              {
                  const Real val = ffab(p);
                  if (val > 0)
                      cfab(amrex::coarsen(p,trr)) += val;
              }
          }

          temp_dat.clear();

          MultiFab dat(grids,dmap,1,0);
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

        // have to do it here, not in read_particle_params, because Density, ..., are set after
        // read_particle_params is called.
        if (particles::timestamp_density) {
            timestamp_indices.push_back(URHO);
            std::cout << "Density = " << URHO << std::endl;
        }
        if (particles::timestamp_temperature) {
            timestamp_indices.push_back(UTEMP);
            std::cout << "Temp = " << UTEMP << std::endl;
        }

        if (!timestamp_indices.empty()) {
            imax = *(std::max_element(timestamp_indices.begin(), timestamp_indices.end()));
        }
    }

    if ( TracerPC && !particles::timestamp_dir.empty())
    {
        std::string basename = particles::timestamp_dir;

        if (basename[basename.length()-1] != '/') basename += '/';

        basename += "Timestamp";

        int finest_level = parent->finestLevel();
        Real time        = state[State_Type].curTime();

        for (int lev = level; lev <= finest_level; lev++)
        {
            if (TracerPC->NumberOfParticlesAtLevel(lev) <= 0) continue;

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

void
Castro::advance_particles(int iteration, Real time, Real dt)
{
    if (TracerPC)
    {
        int ng = iteration;
        Real t = time + 0.5*dt;

        MultiFab Ucc(grids,dmap,BL_SPACEDIM,ng); // cell centered velocity

        {
            FillPatchIterator fpi(*this, Ucc, ng, t, State_Type, 0, BL_SPACEDIM+1);
            MultiFab& S = fpi.get_mf();

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(Ucc,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox();
                S[mfi].invert(1.0, bx, 0, 1);
                for (int dir=0; dir < BL_SPACEDIM; ++dir) {
                    Ucc[mfi].copy(S[mfi], bx, dir+1, bx, dir, 1);
                    Ucc[mfi].mult(S[mfi], bx, 0, dir);
                }
            }
        }

        TracerPC->AdvectWithUcc(Ucc, level, dt);
    }
}
