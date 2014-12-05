#include <iomanip>
#include "Castro.H"
#ifdef GRAVITY
#include "Gravity.H"
#endif
#include "Castro_F.H"

#ifdef PARTICLES

#include <Particles_F.H>

static std::string ascii_particle_file;

static DarkMatterParticleContainer* DMPC = 0;  // There's really only one of these.
DarkMatterParticleContainer* Castro::theDMPC () { return DMPC; }

bool Castro::do_dm_particles            = false;

std::string Castro::particle_init_type  = "";
std::string Castro::particle_move_type  = "";

bool Castro::particle_initrandom_serialize = false;
Real Castro::particle_initrandom_mass;
long Castro::particle_initrandom_count;
int  Castro::particle_initrandom_iseed;

int Castro::particle_verbose            =  1;
int Castro::write_particles_in_plotfile =  0;

Real Castro::particle_cfl               =  0.5;

static const std::string chk_particle_file("DM");

void 
Castro::read_particle_params ()
{
    ParmParse pp("castro");
    pp.query("do_dm_particles",do_dm_particles);

    pp.get("particle_init_type",particle_init_type);
    pp.get("particle_move_type",particle_move_type);

    if (!do_grav && particle_move_type == "Gravitational") 
    {
        std::cerr << "ERROR:: doesnt make sense to have do_grav=false but move_type = Gravitational" << std::endl;
        BoxLib::Error();
    }

    pp.query("particle_initrandom_serialize",particle_initrandom_serialize);
    pp.query("particle_initrandom_count",particle_initrandom_count);
    pp.query("particle_initrandom_mass",particle_initrandom_mass);
    pp.query("particle_initrandom_iseed",particle_initrandom_iseed);
 
    pp.query("ascii_particle_file",ascii_particle_file);

    if (!ascii_particle_file.empty() && particle_init_type != "AsciiFile") 
    {
        std::cerr << "ERROR::particle_init_type is not AsciiFile but you specified ascii_particle_file" << std::endl;;
        BoxLib::Error();
    }

    //
    // Control the verbosity of the Particle class
    //
    ParmParse ppp("particles");
    ppp.query("v",particle_verbose);
    ppp.query("write_in_plotfile",write_particles_in_plotfile);

    //
    // Set the cfl for particle motion (fraction of cell that a particle can move in a timestep)
    //
    ppp.query("cfl",particle_cfl);
}

void
Castro::init_particles ()
{
    BL_PROFILE("Castro::init_particles()");

     if (level > 0)
        return;
     //
     // Need to initialize particles before defining gravity.
     //
     if (do_dm_particles)
     {
         BL_ASSERT(DMPC == 0);

         DMPC = new DarkMatterParticleContainer(parent);
         //
         // 2 gives more stuff than 1.
         //
         DMPC->SetVerbose(particle_verbose);

         if (particle_init_type == "Random")
         {
             if (particle_initrandom_count <= 0)
             {
                 BoxLib::Abort("Castro::init_particles(): particle_initrandom_count must be > 0");
             }
             if (particle_initrandom_iseed <= 0)
             {
                 BoxLib::Abort("Castro::init_particles(): particle_initrandom_iseed must be > 0");
             }

             if (verbose && ParallelDescriptor::IOProcessor())
             {
                 std::cout << "\nInitializing DM with cloud of " 
                           << particle_initrandom_count
                           << " random particles with initial seed: "
                           << particle_initrandom_iseed << "\n\n";
             }

             DMPC->InitRandom(particle_initrandom_count,particle_initrandom_iseed, 
                              particle_initrandom_mass, particle_initrandom_serialize);
         }
         else if (particle_init_type == "AsciiFile")
         {
             if (verbose && ParallelDescriptor::IOProcessor())
             {
                 std::cout << "\nInitializing DM particles from \""
                           << ascii_particle_file
                           << "\" ...\n\n";
             }
             //
             // The second argument is how many Reals we read into m_data[]
             // after reading in m_pos[].  Here we're reading in the particle
             // mass and velocity.
             //
             DMPC->InitFromAsciiFile(ascii_particle_file,BL_SPACEDIM+1);
         }
         else
         {
             BoxLib::Error("not a valid input for castro.particle_init_type");
         }
     }
     DMPC->SetAllowParticlesNearBoundary(true);
}

void
Castro::ParticleCheckPoint(const std::string& dir)
{
    if (level == 0)
    {
        if (DMPC)
            DMPC->Checkpoint(dir,chk_particle_file);

    }
}

void
Castro::ParticlePlotFile(const std::string& dir)
{
    if (level == 0)
    {
       if (DMPC && write_particles_in_plotfile)
          DMPC->Checkpoint(dir,chk_particle_file);
    }
}

void
Castro::ParticlePostRestart (const std::string& restart_file)
{
    if (level == 0)
    {
        if (do_dm_particles)
        {
            BL_ASSERT(DMPC == 0);

            DMPC = new DarkMatterParticleContainer(parent);
            //
            // 2 gives more stuff than 1.
            //
            DMPC->SetVerbose(particle_verbose);

            DMPC->Restart(restart_file, chk_particle_file);
            //
            // We want the ability to write the particles out to an ascii file.
            //
            ParmParse pp("particles");

            std::string particle_output_file;

            pp.query("particle_output_file", particle_output_file);

            if (!particle_output_file.empty())
            {
                DMPC->WriteAsciiFile(particle_output_file);
            }
        }
    }
}

#ifdef GRAVITY
void
Castro::ParticleEstTimeStep(Real& estdt)
{
    if (DMPC && particle_move_type == "Gravitational")
    {

       const Real      a              = 1.0;
       MultiFab& grav                 = get_new_data(Gravity_Type);
       const Real      estdt_particle = DMPC->estTimestep(grav,a,level,particle_cfl);

       if (estdt_particle > 0.0) 
          estdt = std::min(estdt, estdt_particle);

       if (verbose && ParallelDescriptor::IOProcessor())
       {
           if (estdt_particle > 0.0) 
           {
              std::cout << "...estdt from particles at level " << level << ": " << estdt_particle << '\n';
           } else {
              std::cout << "...there are no particles at level " << level << std::endl;
           }
       }
    }
}
#endif

void
Castro::ParticleRedistribute()
{
    if (DMPC)
        DMPC->Redistribute();
}

void
Castro::ParticleMoveRandom()
{
    if (DMPC && particle_move_type == "Random")
    {
        BL_ASSERT(level == 0);

        DMPC->MoveRandom();
    }
}

MultiFab*
Castro::ParticleDerive(const std::string& name,
                       Real               time,
                       int                ngrow)
{
    BL_PROFILE("Castro::ParticleDerive()");

  if (DMPC && name == "particle_count")
  {
      MultiFab* derive_dat = new MultiFab(grids,1,0);
      MultiFab    temp_dat(grids,1,0);
      temp_dat.setVal(0);
      DMPC->Increment(temp_dat,level);
      MultiFab::Copy(*derive_dat,temp_dat,0,0,1,0);
      return derive_dat;
  }
  else if (DMPC && name == "total_particle_count")
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

          DMPC->Increment(temp_dat,lev);

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
  else if (DMPC && name == "particle_mass_density")
  {
      MultiFab* derive_dat = new MultiFab(grids,1,0);
 
      // We need to do the multilevel AssignDensity even though we're only asking
      //  for one level's worth because otherwise we don't get the coarse-fine
      //  distribution of particles correct.
      PArray<MultiFab> partmf;
      DMPC->AssignDensity(partmf); 

      for (int lev = parent->finestLevel()-1; lev >= 0; lev--)
      {
         const IntVect& ratio = parent->refRatio(lev);
         Gravity::avgDown(partmf[lev],partmf[lev+1],ratio);
      }

      MultiFab::Copy(*derive_dat,partmf[level],0,0,1,0);
      return derive_dat;
  }
  else
  {
     return AmrLevel::derive(name,time,ngrow);
  }
}
#endif
