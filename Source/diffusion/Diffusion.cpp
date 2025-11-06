#include <AMReX_ParmParse.H>
#include <Diffusion.H>
#include <Castro.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <MGutils.H>

#include <castro_limits.H>

using namespace amrex;

Diffusion::Diffusion(Amr* Parent, BCRec* _phys_bc)
  :
    parent(Parent),
    LevelData(MAX_LEV),
    grids(MAX_LEV),
    volume(MAX_LEV),
    area(MAX_LEV),
    phys_bc(_phys_bc)
{
    AMREX_ALWAYS_ASSERT(parent->maxLevel() < MAX_LEV);

    make_mg_bc();
}


void
Diffusion::output_job_info_params(std::ostream& jobInfoFile)
{
#include <diffusion_job_info_tests.H>
}


void
Diffusion::install_level (int                   level,
                          AmrLevel*             level_data,
                          MultiFab&             _volume,
                          MultiFab*             _area)
{
    if (verbose) {
        amrex::Print() << "Installing Diffusion level " << level << '\n';
    }

    LevelData[level] = level_data;

    volume[level] = &_volume;

    area[level] = _area;

    BoxArray ba(LevelData[level]->boxArray());
    grids[level] = ba;
}

void
Diffusion::applyop (int level, MultiFab& Temperature,
                    MultiFab& CrseTemp, MultiFab& DiffTerm,
                    Vector<std::unique_ptr<MultiFab> >& temp_cond_coef)
{
    applyop_mlmg(level, Temperature, CrseTemp, DiffTerm, temp_cond_coef);
}


void
Diffusion::make_mg_bc ()
{
    const Geometry& geom = parent->Geom(0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            mlmg_lobc[idim] = MLLinOp::BCType::Periodic;
            mlmg_hibc[idim] = MLLinOp::BCType::Periodic;
        } else {
            if (phys_bc->lo(idim) == amrex::PhysBCType::inflow) {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            }
            if (phys_bc->hi(idim) == amrex::PhysBCType::symmetry) {
                mlmg_hibc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_hibc[idim] = MLLinOp::BCType::Neumann;
            }
        }
    }

    // Set Neumann bc at r=0.
    const auto* problo = geom.ProbLo();
    if ((geom.IsSPHERICAL() || geom.IsRZ()) && problo[0] == 0.0_rt) {
        mlmg_lobc[0] = MLLinOp::BCType::Neumann;
    }

}

void
Diffusion::applyop_mlmg (int level, MultiFab& Temperature,
                         MultiFab& CrseTemp, MultiFab& DiffTerm,
                         Vector<std::unique_ptr<MultiFab> >& temp_cond_coef)
{
    BL_PROFILE("Diffusion::applyop_mlmg()");

    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... compute diffusive term at level " << level << '\n';
    }

    const Geometry& geom = parent->Geom(level);
    const BoxArray& ba = Temperature.boxArray();
    const DistributionMapping& dm = Temperature.DistributionMap();

    LPInfo info;
    info.setMetricTerm(true);
    info.setMaxCoarseningLevel(0);
    info.setAgglomeration(false);
    info.setConsolidation(false);

    MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
    mlabec.setMaxOrder(diffusion::mlmg_maxorder);

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);

    if (level > 0) {
        const auto& rr = parent->refRatio(level-1);
        mlabec.setCoarseFineBC(&CrseTemp, rr[0]);
    }
    mlabec.setLevelBC(0, &Temperature);

    mlabec.setScalars(0.0, -1.0);
    mlabec.setBCoeffs(0, Array<MultiFab const*, AMREX_SPACEDIM>{AMREX_D_DECL(temp_cond_coef[0].get(),
                                                                             temp_cond_coef[1].get(),
                                                                             temp_cond_coef[2].get())});

    MLMG mlmg(mlabec);
    mlmg.setVerbose(verbose);
    mlmg.apply({&DiffTerm}, {&Temperature});
}
