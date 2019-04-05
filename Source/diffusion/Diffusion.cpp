#include <AMReX_ParmParse.H>
#include "Diffusion.H"
#include "Castro.H"
#include "Castro_F.H"
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>

#define MAX_LEV 15

#include "diffusion_defaults.H"

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
    read_params();
    make_mg_bc();
}

Diffusion::~Diffusion() {}

void
Diffusion::read_params ()
{
    static bool done = false;

    if (!done)
    {
        ParmParse pp("diffusion");

#include "diffusion_queries.H"

        done = true;
    }
}


void 
Diffusion::output_job_info_params(std::ostream& jobInfoFile)
{
#include "diffusion_job_info_tests.H"
}


void
Diffusion::install_level (int                   level,
                          AmrLevel*             level_data,
                          MultiFab&             _volume,
                          MultiFab*             _area)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Installing Diffusion level " << level << '\n';

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

#if (BL_SPACEDIM < 3)
void
Diffusion::weight_cc(int level, MultiFab& cc)
{
    const Real* dx = parent->Geom(level).CellSize();
    const int coord_type = Geometry::Coord();
#ifdef _OPENMP
#pragma omp parallel	  
#endif
    for (MFIter mfi(cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_weight_cc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		     BL_TO_FORTRAN_ANYD(cc[mfi]),
                     AMREX_REAL_ANYD(dx), coord_type);
    }
}

void
Diffusion::unweight_cc(int level, MultiFab& cc)
{
    const Real* dx = parent->Geom(level).CellSize();
    const int coord_type = Geometry::Coord();
#ifdef _OPENMP
#pragma omp parallel	  
#endif
    for (MFIter mfi(cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

#pragma gpu
        ca_unweight_cc(AMREX_INT_ANYD(bx.loVect()), AMREX_INT_ANYD(bx.hiVect()),
		       BL_TO_FORTRAN_ANYD(cc[mfi]),
                       AMREX_REAL_ANYD(dx), coord_type);
    }
}
#endif

void
Diffusion::make_mg_bc ()
{
    const Geometry& geom = parent->Geom(0);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom.isPeriodic(idim)) {
            mlmg_lobc[idim] = MLLinOp::BCType::Periodic;
            mlmg_hibc[idim] = MLLinOp::BCType::Periodic;
        } else {
            if (phys_bc->lo(idim) == Inflow) {
                mlmg_lobc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_lobc[idim] = MLLinOp::BCType::Neumann;
            }
            if (phys_bc->hi(idim) == Symmetry) {
                mlmg_hibc[idim] = MLLinOp::BCType::Dirichlet;
            } else {
                mlmg_hibc[idim] = MLLinOp::BCType::Neumann;
            }
        }
    }

    // Set Neumann bc at r=0.
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() ) {
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
    
    MLABecLaplacian mlabec({geom}, {ba}, {dm},
                           LPInfo().setMetricTerm(true).setMaxCoarseningLevel(0));
    mlabec.setMaxOrder(mlmg_maxorder);

    mlabec.setDomainBC(mlmg_lobc, mlmg_hibc);

    if (level > 0) {
        const auto& rr = parent->refRatio(level-1);
        mlabec.setCoarseFineBC(&CrseTemp, rr[0]);
    }
    mlabec.setLevelBC(0, &Temperature);

    mlabec.setScalars(0.0, -1.0);
    mlabec.setBCoeffs(0, {AMREX_D_DECL(temp_cond_coef[0].get(),
                                       temp_cond_coef[1].get(),
                                       temp_cond_coef[2].get())});

    MLMG mlmg(mlabec);
    mlmg.setVerbose(verbose);
    mlmg.apply({&DiffTerm}, {&Temperature});
}
