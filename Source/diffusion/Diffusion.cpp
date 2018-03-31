#include <AMReX_ParmParse.H>
#include "Diffusion.H"
#include "Castro.H"
#include "Castro_F.H"
#include <AMReX_FMultiGrid.H>

#define MAX_LEV 15

#include "diffusion_defaults.H"

using namespace amrex;

int  Diffusion::stencil_type = CC_CROSS_STENCIL;
 
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
    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... compute diffusive term at level " << level << '\n';
    }

    Vector<std::unique_ptr<MultiFab> > coeffs_curv;
#if (BL_SPACEDIM < 3)
    // NOTE: we just pass DiffTerm here to use in the MFIter loop...
    if (Geometry::IsRZ() || Geometry::IsSPHERICAL())
    {
	coeffs_curv.resize(BL_SPACEDIM);

	for (int i = 0; i< BL_SPACEDIM; ++i) {
	    coeffs_curv[i].reset(new MultiFab(temp_cond_coef[i]->boxArray(),
					      temp_cond_coef[i]->DistributionMap(),
					      1, 0));
	    MultiFab::Copy(*coeffs_curv[i], *temp_cond_coef[i], 0, 0, 1, 0);
	}

	applyMetricTerms(0, DiffTerm, coeffs_curv);
    }
#endif

    auto & coeffs = (coeffs_curv.size() > 0) ? coeffs_curv : temp_cond_coef;

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                   : IntVect::TheZeroVector();

    FMultiGrid fmg(parent->Geom(level), level, crse_ratio);

    if (level == 0) {
	fmg.set_bc(mg_bc, Temperature);
    } else {
	fmg.set_bc(mg_bc, CrseTemp, Temperature);
    }

    fmg.set_diffusion_coeffs(amrex::GetVecOfPtrs(coeffs));

    fmg.applyop(Temperature, DiffTerm);

#if (BL_SPACEDIM < 3)
    // Do this to unweight Res
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
	unweight_cc(level, DiffTerm);
#endif
}

#if (BL_SPACEDIM == 1)
void
Diffusion::applyViscOp (int level, MultiFab& Vel, 
                        MultiFab& CrseVel, MultiFab& ViscTerm, 
                        Vector<std::unique_ptr<MultiFab> >& visc_coeff)
{
    if (verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "   " << '\n';
        std::cout << "... compute second part of viscous term at level " << level << '\n';
    }

    IntVect crse_ratio = level > 0 ? parent->refRatio(level-1)
                                   : IntVect::TheZeroVector();

    FMultiGrid fmg(parent->Geom(level), level, crse_ratio);

    if (level == 0) {
	fmg.set_bc(mg_bc, Vel);
    } else {
	fmg.set_bc(mg_bc, CrseVel, Vel);
    }

    // Here we DO NOT multiply the coefficients by (1/r^2) for spherical coefficients
    // because we are computing (1/r^2) d/dr (const * d/dr(r^2 u))
    fmg.set_diffusion_coeffs(amrex::GetVecOfPtrs(visc_coeff));

#if (BL_SPACEDIM < 3)
    // Here we weight the Vel going into the FMG applyop
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
	weight_cc(level, Vel);
#endif

    fmg.applyop(Vel, ViscTerm);

#if (BL_SPACEDIM < 3)
    // Do this to unweight Res
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
	unweight_cc(level, ViscTerm);
#endif
}
#endif

#if (BL_SPACEDIM < 3)
void
Diffusion::applyMetricTerms(int level, MultiFab& Rhs, Vector<std::unique_ptr<MultiFab> >& coeffs)
{
    const Real* dx = parent->Geom(level).CellSize();
    const int coord_type = Geometry::Coord();
#ifdef _OPENMP
#pragma omp parallel	  
#endif
    for (MFIter mfi(Rhs,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
	D_TERM(const Box xbx = mfi.nodaltilebox(0);,
	       const Box ybx = mfi.nodaltilebox(1);,
	       const Box zbx = mfi.nodaltilebox(2);)
        // Modify Rhs and coeffs with the appropriate metric terms.
        ca_apply_metric(bx.loVect(), bx.hiVect(),
			D_DECL(xbx.loVect(),
			       ybx.loVect(),
			       zbx.loVect()),
			D_DECL(xbx.hiVect(),
			       ybx.hiVect(),
			       zbx.hiVect()),
			BL_TO_FORTRAN(Rhs[mfi]),
			D_DECL(BL_TO_FORTRAN((*coeffs[0])[mfi]),
			       BL_TO_FORTRAN((*coeffs[1])[mfi]),
			       BL_TO_FORTRAN((*coeffs[2])[mfi])),
			dx,&coord_type);
    }
}
#endif

#if (BL_SPACEDIM < 3)
void
Diffusion::weight_cc(int level, MultiFab& cc)
{
    const Real* dx = parent->Geom(level).CellSize();
    const int coord_type = Geometry::Coord();
#ifdef _OPENMP
#pragma omp parallel	  
#endif
    for (MFIter mfi(cc,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        ca_weight_cc(bx.loVect(), bx.hiVect(),
		     BL_TO_FORTRAN(cc[mfi]),dx,&coord_type);
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
    for (MFIter mfi(cc,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        ca_unweight_cc(bx.loVect(), bx.hiVect(),
		       BL_TO_FORTRAN(cc[mfi]),dx,&coord_type);
    }
}
#endif

void
Diffusion::make_mg_bc ()
{
    const Geometry& geom = parent->Geom(0);
    for ( int dir = 0; dir < BL_SPACEDIM; ++dir )
    {
        if ( geom.isPeriodic(dir) )
        {
            mg_bc[2*dir + 0] = 0;
            mg_bc[2*dir + 1] = 0;
        }
        else
        {
            if (phys_bc->lo(dir) == Symmetry   || 
                phys_bc->lo(dir) == SlipWall   || 
                phys_bc->lo(dir) == NoSlipWall || 
                phys_bc->lo(dir) == Outflow)   
            { 
              mg_bc[2*dir + 0] = MGT_BC_NEU;
            }
            else if (phys_bc->lo(dir) ==  Inflow)
            { 
              mg_bc[2*dir + 0] = MGT_BC_DIR;
            } else {
              amrex::Error("Failed to set lo mg_bc in Diffusion::make_mg_bc" );
            }

            if (phys_bc->hi(dir) == Symmetry   || 
                phys_bc->hi(dir) == SlipWall   || 
                phys_bc->hi(dir) == NoSlipWall || 
                phys_bc->hi(dir) == Outflow)   
            {
              mg_bc[2*dir + 1] = MGT_BC_NEU;
            } 
            else if (phys_bc->hi(dir) ==  Inflow)
            { 
              mg_bc[2*dir + 0] = MGT_BC_DIR;
            } else {
              amrex::Error("Failed to set hi mg_bc in Diffusion::make_mg_bc" );
            }
        }
    }

    // Set Neumann bc at r=0.
    if (Geometry::IsSPHERICAL() || Geometry::IsRZ() )
        mg_bc[0] = MGT_BC_NEU;
}

