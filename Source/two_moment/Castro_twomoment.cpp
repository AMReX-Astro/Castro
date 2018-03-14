#ifdef THORNDAO
#include "Castro.H"
#include "Castro_F.H"

using std::string;
using namespace amrex;

Real
Castro::init_thornado()
{
    int my_ncomp = BL_SPACEDIM+3;  \\ rho, rho*u, rho*v, rho*w, rho*E, Y_e
    int my_ngrow = 1;  \\ one fluid ghost cell
    // Create a temporary so it has the right order of the right variables and no ghost cells
    MultiFab U_F(grids, dmap, my_ncomp, my_ngrow);

    // Copy the current state S_new into U_F
    // Note that the first five components get copied as is
    int src_comp = 0;
    int dst_comp = 0;
    int   n_comp = BL_SPACEDIM+2;
    MultiFab::Copy(U_F, S_new, src_comp, dst_comp, n_comp, my_ngrow); \\ rho, rho*u, rho*v, rho*w, rho*E

    // Now copy Y_e, which is first "aux" variable, into the 6th spot of U_F
    int src_comp = FirstAux;
    int dst_comp = BL_SPACEDIM+2;
    MultiFab::Copy(U_F, S_new, src_comp, dst_comp, 1, my_grow);

    MultiFab& U_R_old = get_old_data(Thornado_Type);
    MultiFab& U_R_new = get_new_data(Thornado_Type);

    // For now we will not allowing logical tiling
    for (MFIter mfi(U_F, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        init_thornado(BL_TO_FORTRAN_BOX(bx),
                      BL_TO_FORTRAN_FAB(U_F[mfi]),
                      BL_TO_FORTRAN_FAB(U_R_old[mfi]));
    }
    MultiFab::Copy(U_R_new, U_R_old, 0, 0, U_R_old.nComp(), U_R_old.nGrow());
}
#endif
