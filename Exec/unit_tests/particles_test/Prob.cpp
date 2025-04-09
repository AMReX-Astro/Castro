/* Implementations of functions in Problem.H go here */

#include <Castro.H>
#include <AMReX_TracerParticles.H>
#include <AMReX_AmrParticles.H>

using namespace amrex;
using ParticleType = Particle<1, 1>;

#ifdef DO_PROBLEM_POST_SIMULATION
void Castro::problem_post_simulation(Vector<std::unique_ptr<AmrLevel> >& amr_level) {

    int nlevels = amr_level.size();

    Real err = -1.e30;

    auto lev = 0;
    auto ipass = 0;

    Castro& castro = dynamic_cast<Castro&>(*amr_level[lev]);

    const Real strttime = amrex::second();
    const Real*     dx       = castro.geom.CellSize();
    const Real*     plo      = castro.geom.ProbLo();
    const Real*     phi      = castro.geom.ProbHi();

    // read initial positions from Ascii file

    const int MyProc   = ParallelDescriptor::MyProc();

    Vector<ParticleType> nparticles;

    std::string file = "particle_file";

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

    std::ifstream ifs;

    ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    ifs.open(file.c_str(), std::ios::in);

    if (!ifs.good())
        amrex::FileOpenFailed(file);

    int num_particles = 0;

    ifs >> num_particles >> std::ws;

    ParticleLocData pld;
    ParticleType p, p_rep;

    for (int i = 0; i < num_particles; i++)
    {
        AMREX_D_TERM(ifs >> p.m_pos[0]; ,
                     ifs >> p.m_pos[1]; ,
                     ifs >> p.m_pos[2]; );

        auto extradata = 0;

        for (int n = 0; n < extradata; n++)
        {
            ifs >> p.m_rdata[AMREX_SPACEDIM+n];
        }

        p.id()  = ParticleType::NextID();
        p.cpu() = MyProc;

        nparticles.push_back(p);
    }

    // compute the change in position of the particles wrt to their initial positions

    auto& pmap = TracerPC->GetParticles(lev);

    double total_err = 0;

    for (auto& kv : pmap) {
        int grid = kv.first.first;
        auto& pbox = kv.second.GetArrayOfStructs();
        const int n = pbox.size();

        for (int i = 0; i < n; i++)
        {
            auto& p = pbox[i];

            bool match = false;

            auto it = nparticles.begin();

            // find the original particle and calculate change in position
            for (; !match && it != nparticles.end(); ++it) {
                if (it->id() == p.id())
                    match = true;
            }

            if (!match) Print() << "haven't found a match :(" << std::endl;

            // calculate change in position
            auto deltax = it->m_pos[0] - p.m_pos[0];
            auto deltay = it->m_pos[1] - p.m_pos[1];

            auto delta = std::sqrt(deltax*deltax + deltay*deltay);

            // quick nan check here
            if (delta == delta)
                total_err += delta;

            // Print() << "delta r = " << delta << std::endl;

        }
    }

    const std::string stars(78,'*');
    amrex::Print() << stars << "\n"
                   << " particles problem post_simulation() \n"
                   << " Average change in position: " << total_err / num_particles << "\n"
                   << stars << "\n\n";
}
#endif
