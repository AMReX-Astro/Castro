#include "Castro.H"
#include "Castro_F.H"

using namespace amrex;

void
Castro::do_sdc_update(int m_start, int m_end, Real dt_m) {

  // this routine needs to do the update from time node m to m+1
  //
  // We come in with:
  //   A_new[m_start] : this is the advective update at node m_start
  //   A_old[:] : this is the advective source for all nodes at the old iterate
  //
  //   R_old[:] : this is the reaction source for all nodes at the old iterate

  // If we do advection only, then the update is explicit.  If we do
  // reactions, then the update is implicit within a zone.

#ifdef REACTIONS
  amrex::Abort("SDC update not implemented with reactions");
#endif

  for (MFIter mfi(*k_new[0]); mfi.isValid(); ++mfi) {

    const Box& bx = mfi.tilebox();

#ifdef REACTIONS
    // advection + reactions
    if (mol_order == 2) {
      ca_sdc_update_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                       BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                       BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                       BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                       BL_TO_FORTRAN_3D((*A_old[0])[mfi]),
                       BL_TO_FORTRAN_3D((*A_old[1])[mfi]),
                       &m_start);
    } else {
      amrex::Abort("mol_order != 2 not implemented");
    }
#else
    // pure advection
    if (mol_order == 2) {
      ca_sdc_update_advection_o2(BL_TO_FORTRAN_BOX(bx), &dt_m,
                                 BL_TO_FORTRAN_3D((*k_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*k_new[m_end])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_new[m_start])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[0])[mfi]),
                                 BL_TO_FORTRAN_3D((*A_old[1])[mfi]),
                                 &m_start);
    } else {
      amrex::Abort("mol_order != 2 not implemented");
    }
#endif

  }
}


void
Castro::construct_old_react_source() {

  // this routine simply fills R_old with the old-iteration reactive
  // source

}
