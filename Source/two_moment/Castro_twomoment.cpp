
#include "Castro.H"
#include "TwoMoment_F.H"

using std::string;
using namespace amrex;

int
Castro::init_thornado()
{
    int nDimsX   = BL_SPACEDIM;
    int nDimsE   = 1;
    int nSpecies = THORNADO_NSPECIES;

    InitThornado(&nDimsX, &nDimsE, &nSpecies);
   
    int ncomp_thornado;

    ca_get_rad_ncomp(&ncomp_thornado);

    return ncomp_thornado;
}

void
Castro::init_thornado_data()
{
    MultiFab& Thor_new = get_new_data(Thornado_Type);

    int nc = Thor_new.nComp();
    const Real* dx = geom.CellSize();
    const Real  cur_time = state[Thornado_Type].curTime();

    for (MFIter mfi(Thor_new); mfi.isValid(); ++mfi)
    {
       RealBox gridloc = RealBox(grids[mfi.index()],geom.CellSize(),geom.ProbLo());
       const Box& box     = mfi.validbox();
       const int* lo      = box.loVect();
       const int* hi      = box.hiVect();
  
        ca_init_thornado_data
	  (level, cur_time, lo, hi, nc,
           BL_TO_FORTRAN(Thor_new[mfi]), dx,
           gridloc.lo(), gridloc.hi());
    }
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

    int swE = 0; // stencil for energy array; no energy ghost cells needed
    Vector<Real> grid_lo(3);
    Vector<Real> grid_hi(3);

    const Real* dx = geom.CellSize();

    int n_fluid_dof = THORNADO_FLUID_NDOF;
    int n_moments   = THORNADO_NMOMENTS;

    // For right now create a temporary holder for the source term -- we'll incorporate it 
    //    more permanently later.  
    MultiFab dS(grids, dmap, S_new.nComp(), S_new.nGrow());

    Real eL = 0.;
    Real eR = 1.;

    int swX[3];
    swX[0] = my_ngrow;
    swX[1] = my_ngrow;
    swX[2] = my_ngrow;

    int * boxlen = new int[3];

    for (int i = 0; i < n_sub; i++)
    {
      Real dt_sub = dt / n_sub;

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
           grid_lo[0] =  bx.smallEnd(0)  * dx[0];
           grid_lo[1] =  bx.smallEnd(1)  * dx[1];
           grid_lo[2] =  bx.smallEnd(2)  * dx[2];
           grid_hi[0] = (bx.bigEnd(0)+1) * dx[0];
           grid_hi[1] = (bx.bigEnd(1)+1) * dx[1];
           grid_hi[2] = (bx.bigEnd(2)+1) * dx[2];

           boxlen[0] = bx.length(0);
           boxlen[1] = bx.length(1);
           boxlen[2] = bx.length(2);
  
           InitThornado_Patch(boxlen, swX,
               grid_lo.dataPtr(), grid_hi.dataPtr(),
               &swE, &eL, &eR);
        }

        call_to_thornado(BL_TO_FORTRAN_BOX(bx), &dt_sub,
                         BL_TO_FORTRAN_FAB(S_border[mfi]),
                         BL_TO_FORTRAN_FAB(dS[mfi]),
                         BL_TO_FORTRAN_FAB(R_border[mfi]),
                         BL_TO_FORTRAN_FAB(U_R_new[mfi]), 
                         &n_fluid_dof, &n_moments, &my_ngrow);

        // Add the source term to all components even though there should
        //     only be non-zero source terms for (Rho, Xmom, Ymom, Zmom, RhoE, UFX)
        MultiFab::Add(S_new, dS, Density, 0, S_new.nComp(), 0);

        if (i == (n_sub-1)) FreeThornado_Patch();
      }
      S_border.FillBoundary();
    }
    delete boxlen;
}
