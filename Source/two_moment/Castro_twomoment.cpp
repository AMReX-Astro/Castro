
#include "Castro.H"
#include "TwoMoment_F.H"

using std::string;
using namespace amrex;

void
Castro::init_thornado()
{
    int nDimsX = THORNADO_NDIMS_X;
    int nDimsE = THORNADO_NDIMS_E;
    InitThornado(&nDimsX, &nDimsE);
}

void
Castro::create_thornado_source(Real dt)
{
    MultiFab& S_new = get_new_data(State_Type);

    int my_ncomp = BL_SPACEDIM+3;  // rho, rho*u, rho*v, rho*w, rho*E, Y_e
    int my_ngrow = 2;  // two fluid ghost cells

    // This fills the ghost cells of the MultiFab which we will then copy into U_F
    MultiFab S_with_ghost_cells(grids, dmap, NUM_STATE, my_ngrow);
    const Real  cur_time = state[State_Type].curTime();
    AmrLevel::FillPatch(*this, S_with_ghost_cells, my_ngrow, cur_time, State_Type, 0, NUM_STATE);

    // Create U_F to hold the right order of the right variables
    MultiFab U_F(grids, dmap, my_ncomp, my_ngrow);

    // Copy into U_F just these variables: rho, rho*u, rho*v, rho*w, rho*E, rho*ne
    int      cnt = 0;
    MultiFab::Copy(U_F, S_with_ghost_cells, Density , cnt, 1, my_ngrow); cnt++;
    MultiFab::Copy(U_F, S_with_ghost_cells, Xmom    , cnt, 1, my_ngrow); cnt++;
    MultiFab::Copy(U_F, S_with_ghost_cells, Ymom    , cnt, 1, my_ngrow); cnt++;
#if (BL_SPACEDIM == 3)
    MultiFab::Copy(U_F, S_with_ghost_cells, Zmom    , cnt, 1, my_ngrow); cnt++;
#endif
    MultiFab::Copy(U_F, S_with_ghost_cells, Eden    , cnt, 1, my_ngrow); cnt++;
    MultiFab::Copy(U_F, S_with_ghost_cells, FirstAux, cnt, 1, my_ngrow); 

    MultiFab& U_R_old = get_old_data(Thornado_Type);
    MultiFab& U_R_new = get_new_data(Thornado_Type);

    // int n_sub = GetNSteps(dt); // From thornado
    int n_sub = 1; // THIS IS JUST A HACK TO MAKE IT COMPILE 

    int swE = 0; // stencil for energy array; no energy ghost cells needed
    Vector<Real> grid_lo(3);
    Vector<Real> grid_hi(3);

    const Real* dx = geom.CellSize();

    int n_fluid_dof = THORNADO_FLUID_NDOF;
    int n_rad_dof   = THORNADO_RAD_NDOF;
    int n_energy    = THORNADO_NENERGY;
    int n_species   = THORNADO_NSPECIES;
    int n_moments   = THORNADO_NMOMENTS;

    // For right now create a temporary holder for the source term -- we'll incorporate it 
    //    more permanently later.  
    MultiFab dS(grids, dmap, S_new.nComp(), S_new.nGrow());

    // ASA -- WE NEED TO SET THESE FOR REAL
    Real eL = 0.;
    Real eR = 0.;

    int swX[3];
    swX[0] = 1;
    swX[1] = 1;
    swX[2] = 1;

    int * boxlen = new int[3];
    int nr_comp = U_R_new.nComp();

    for (int i = 0; i < n_sub; i++)
    {
      Real dt_sub = dt / n_sub;

      // Make sure to zero dS here since not all the 
      //    components will be filled in call_to_thornado
      //    and we don't want to re-add terms from the last iteration
      dS.setVal(0.);

      // For now we will not allowing logical tiling
      for (MFIter mfi(U_F, false); mfi.isValid(); ++mfi) 
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
               &n_energy, &swE, &eL, &eR, &n_species);
        }

        // n_fluid_dof = THORNADO_FLUID_NDOF;
        // n_rad_dof   = THORNADO_RAD_NDOF;
        // n_energy    = THORNADO_NENERGY;
        // n_species   = THORNADO_NSPECIES;
        // n_moments   = THORNADO_NMOMENTS;

        call_to_thornado(BL_TO_FORTRAN_BOX(bx), &dt_sub,
                         S_new[mfi].dataPtr(),
                         BL_TO_FORTRAN_FAB(dS[mfi]),
                         U_R_old[mfi].dataPtr(),
                         BL_TO_FORTRAN_FAB(U_R_new[mfi]), 
                         &n_fluid_dof, &n_energy, &n_species, 
                         &n_rad_dof, &n_moments);

        // Add the source term to all components even though there should
        //     only be non-zero source terms for (Rho, Xmom, Ymom, Zmom, RhoE, UFX)
        MultiFab::Add(S_new, dS, Density, 0, S_new.nComp(), 0);

        if (i == (n_sub-1)) FreeThornado_Patch();
      }
      U_F.FillBoundary();
    }
    delete boxlen;
}
