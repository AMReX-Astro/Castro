#ifdef THORNDAO
#include "Castro.H"
#include "Castro_F.H"

using std::string;
using namespace amrex;

void
Castro::init_thornado()
{
    int nDimsX = BL_SPACEDIM;
    int nDimsE = 1;  // This is the number of energy groups -- for now we'll test with just 1
    InitThornado(nDimsX, nDimsE);
}
void
Castro::init_thornado_on_patch()
{
    int nDimsX = BL_SPACEDIM;
    int nDimsE = 1;  // This is the number of energy groups -- for now we'll test with just 1
    InitThornado(nDimsX, nDimsE);
}
void
Castro::create_thornado_source()
{
    int my_ncomp = BL_SPACEDIM+3;  // rho, rho*u, rho*v, rho*w, rho*E, Y_e
    int my_ngrow = 1;  // one fluid ghost cell
    // Create a temporary so it has the right order of the right variables and no ghost cells
    MultiFab U_F(grids, dmap, my_ncomp, my_ngrow);

    // For right now create a temporary holder for the source term -- we'll incorporate it 
    // more permanently later
    MultiFab dU_F(grids, dmap, my_ncomp, my_ngrow);

    // Copy the current state S_new into U_F
    // Note that the first five components get copied as is
    int src_comp = 0;
    int dst_comp = 0;
    int   n_comp = BL_SPACEDIM+2;
    MultiFab::Copy(U_F, S_new, src_comp, dst_comp, n_comp, my_ngrow); // rho, rho*u, rho*v, rho*w, rho*E

    // Now copy Y_e, which is first "aux" variable, into the 6th spot of U_F
    src_comp = FirstAux;
    dst_comp = BL_SPACEDIM+2;
    MultiFab::Copy(U_F, S_new, src_comp, dst_comp, 1, my_ngrow);

    MultiFab& U_R_old = get_old_data(Thornado_Type);
    MultiFab& U_R_new = get_new_data(Thornado_Type);

    // int n_sub = GetNSteps(dt); // From thornado
    int n_sub = 1; // THIS IS JUST A HACK TO MAKE IT COMPILE 

    for (int i = 0; i < n_sub; i++)
    {
      Real dt_sub = dt / n_sub;

      // For now we will not allowing logical tiling
      for (MFIter mfi(U_F, false); mfi.isValid(); ++mfi) 
      {
        if (i == 0)  init_thornado_on_patch();

        const Box& bx = mfi.validbox();
        call_to_thornado(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_FAB(U_F[mfi]),
                         U_R_old[mfi].dataPtr(),
                         BL_TO_FORTRAN_FAB(U_R_new[mfi]),
                         BL_TO_FORTRAN_FAB(dU_F[mfi]),&dt_sub);
      }
      U_F.FillBoundary();
    }
}
#endif
