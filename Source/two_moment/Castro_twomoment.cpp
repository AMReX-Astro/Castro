#include "Castro.H"
#include "AMReX_ParmParse.H"
#include "TwoMoment_F.H"

using std::string;
using namespace amrex;

void
Castro::read_thornado_params ()
{
  ParmParse pp("thornado");

  pp.query("eL"    , thornado_eL);
  pp.query("eR"    , thornado_eR);
  pp.query("ndimse", thornado_ndimse);
  pp.query("zoome" , thornado_zoome);
}

int
Castro::init_thornado()
{
    int nDimsX   = BL_SPACEDIM;
    int nDimsE   = thornado_ndimse; // number of energy groups in thornado
    int nSpecies = THORNADO_NSPECIES;

    Real zoomE    = thornado_zoome;  // geometric zoom factor for energy groups

    amrex::Print() << "*****Calling InitThornado " << std::endl; 
    InitThornado(&nDimsX, &nDimsE, &zoomE, &nSpecies);

    int ncomp_thornado;

    ca_get_rad_ncomp(&ncomp_thornado);

    return ncomp_thornado;
}

void
Castro::init_thornado_data()
{
    MultiFab& Fluid_new = get_new_data(State_Type);
    MultiFab&  Thor_new = get_new_data(Thornado_Type);

    int my_ngrow = 2;  // two fluid ghost cells

    const Real* dx = geom.CellSize();

// *************************************************************
    int swE = 0; // stencil for energy array; no energy ghost cells needed
    Vector<Real> grid_lo(3);
    Vector<Real> grid_hi(3);

    int * boxlen = new int[3];
    const Real* prob_lo   = geom.ProbLo();

    amrex::Print() << "***THORNADO *** Using eL and eR = " << thornado_eL << " " << thornado_eR << std::endl;

    int swX[3];
    swX[0] = my_ngrow;
    swX[1] = my_ngrow;
#if (BL_SPACEDIM > 2)
    swX[2] = my_ngrow;
#else
    swX[2] = 0;
#endif

    int nr =  Thor_new.nComp();
    const Real  cur_time = state[Thornado_Type].curTime();

    amrex::Print() << "*****Calling InitThornado_Patch and init_thornado on each patch " << std::endl; 

    // For now we will not allowing logical tiling
    for (MFIter mfi(Fluid_new, false); mfi.isValid(); ++mfi) 
    {
        Box bx = mfi.validbox();

        grid_lo[0] = prob_lo[0] +  bx.smallEnd(0)  * dx[0] / 100.0; // Factor of 100 because Thornado uses m, not cm
        grid_hi[0] = prob_lo[0] + (bx.bigEnd(0)+1) * dx[0] / 100.0;
        boxlen[0] = bx.length(0);

        grid_lo[1] = prob_lo[1] +  bx.smallEnd(1)  * dx[1] / 100.0;
        grid_hi[1] = prob_lo[1] + (bx.bigEnd(1)+1) * dx[1] / 100.0;
        boxlen[1] = bx.length(1);

#if (BL_SPACEDIM > 2)
        grid_lo[2] = prob_lo[2] +  bx.smallEnd(2)  * dx[2] / 100.0;
        grid_hi[2] = prob_lo[2] + (bx.bigEnd(2)+1) * dx[2] / 100.0;
        boxlen[2] = bx.length(2);
#else
        grid_lo[2] = 0.;
        grid_hi[2] = 1.;
        boxlen[2]  = 1;
#endif

       // Note thornado_eL and thornado_eR in units of MeV; we will convert 
       //      to thornado units inside InitThornado_Patch
        InitThornado_Patch(boxlen, swX,
            grid_lo.dataPtr(), grid_hi.dataPtr(),
            &swE, &thornado_eL, &thornado_eR);

       RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
       const int* lo      = bx.loVect();
       const int* hi      = bx.hiVect();
  
        ca_init_thornado_data
	  (level, cur_time, lo, hi,
           nr, BL_TO_FORTRAN(Thor_new[mfi]), 
               BL_TO_FORTRAN(Fluid_new[mfi]), 
           dx, gridloc.lo(), gridloc.hi());

        FreeThornado_Patch();
    }
    delete boxlen;
}

void
Castro::create_thornado_source(Real dt)
{
    MultiFab& S_new = get_new_data(State_Type);

    MultiFab& U_R_old = get_old_data(Thornado_Type);
    MultiFab& U_R_new = get_new_data(Thornado_Type);

    int my_ngrow = 2;  // two fluid ghost cells

    // This fills the ghost cells of the fluid MultiFab which we will pass into Thornado
    MultiFab S_border(grids, dmap, NUM_STATE, my_ngrow);
    const Real  prev_time = state[State_Type].prevTime();
    AmrLevel::FillPatch(*this, S_border, my_ngrow, prev_time, State_Type, 0, NUM_STATE);

    // This fills the ghost cells of the radiation MultiFab which we will pass into Thornado
    MultiFab R_border(grids, dmap, U_R_old.nComp(), my_ngrow);
    AmrLevel::FillPatch(*this, R_border, my_ngrow, prev_time, Thornado_Type, 0, U_R_old.nComp());

    // int n_sub = GetNSteps(dt); // From thornado
    int n_sub = 1; // THIS IS JUST A HACK TO MAKE IT COMPILE 

    const Real* dx = geom.CellSize();

    int n_fluid_dof = THORNADO_FLUID_NDOF;
    int n_moments   = THORNADO_NMOMENTS;

    // For right now create a temporary holder for the source term -- we'll incorporate it 
    //    more permanently later.  
    MultiFab dS(grids, dmap, S_new.nComp(), S_new.nGrow());
    const Real* prob_lo   = geom.ProbLo();

    Real dt_sub = dt / n_sub;
    int * boxlen = new int[3];

    int swE = 0; // stencil for energy array; no energy ghost cells needed
    Vector<Real> grid_lo(3);
    Vector<Real> grid_hi(3);

    // Note these are in units of MeV; we will convert to thornado units inside
    //      InitThornado_Patch
    amrex::Print() << "***THORNADO *** Using eL and eR = " << thornado_eL << " " << thornado_eR << std::endl;

    int swX[3];
    swX[0] = my_ngrow;
    swX[1] = my_ngrow;
#if (BL_SPACEDIM > 2)
    swX[2] = my_ngrow;
#else
    swX[2] = 0;
#endif

    amrex::Print() << "*****Calling InitThornado_Patch " << std::endl; 

    for (int i = 0; i < n_sub; i++)
    {

      // Make sure to zero dS here since not all the 
      //    components will be filled in call_to_thornado
      //    and we don't want to re-add terms from the last iteration
      dS.setVal(0.);

      // For now we will not allowing logical tiling
      for (MFIter mfi(S_border, false); mfi.isValid(); ++mfi) 
      {
        Box bx = mfi.validbox();

        if (i == 0)  
        {
           grid_lo[0] = prob_lo[0] +  bx.smallEnd(0)  * dx[0] / 100.0; // Factor of 100 because Thornado uses m, not cm
           grid_hi[0] = prob_lo[0] + (bx.bigEnd(0)+1) * dx[0] / 100.0;
           boxlen[0] = bx.length(0);

           grid_lo[1] = prob_lo[1] +  bx.smallEnd(1)  * dx[1] / 100.0;
           grid_hi[1] = prob_lo[1] + (bx.bigEnd(1)+1) * dx[1] / 100.0;
           boxlen[1] = bx.length(1);

#if (BL_SPACEDIM > 2)
           grid_lo[2] = prob_lo[2] +  bx.smallEnd(2)  * dx[2] / 100.0;
           grid_hi[2] = prob_lo[2] + (bx.bigEnd(2)+1) * dx[2] / 100.0;
           boxlen[2] = bx.length(2);
#else
           grid_lo[2] = 0.;
           grid_hi[2] = 1.;
           boxlen[2]  = 1;
#endif
           std::cout << "CALLING PATCH ON BOX " << bx << std::endl;

           // Note thornado_eL and thornado_eR in units of MeV; we will convert 
           //      to thornado units inside InitThornado_Patch
           InitThornado_Patch(boxlen, swX,
               grid_lo.dataPtr(), grid_hi.dataPtr(),
               &swE, &thornado_eL, &thornado_eR);
        }

        std::cout << "CALLING SOURCE ON BOX " << bx << std::endl;
        call_to_thornado(BL_TO_FORTRAN_BOX(bx), &dt_sub,
                         BL_TO_FORTRAN_FAB(S_border[mfi]),
                         BL_TO_FORTRAN_FAB(dS[mfi]),
                         BL_TO_FORTRAN_FAB(R_border[mfi]),
                         BL_TO_FORTRAN_FAB(U_R_new[mfi]), 
                         &n_fluid_dof, &n_moments, &my_ngrow);
        std::cout << "DONE CALLING SOURCE ON BOX " << bx << std::endl;

        // Add the source term to all components even though there should
        //     only be non-zero source terms for (Rho, Xmom, Ymom, Zmom, RhoE, UFX)
        MultiFab::Add(S_new, dS, Density, 0, S_new.nComp(), 0);

//      Note we can't call FreeThornado_Patch here because it is only initialized 
//        at the begining of the run, not every timestep
        if (i == (n_sub-1)) FreeThornado_Patch();
      }
      S_border.FillBoundary();
    }
    delete boxlen;
}
