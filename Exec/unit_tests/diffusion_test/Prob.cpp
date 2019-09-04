/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

using namespace amrex;

#ifdef DO_PROBLEM_POST_SIMULATION
void Castro::problem_post_simulation(Vector<std::unique_ptr<AmrLevel> >& amr_level) {

  // compute the norm of the solution vs. the analytic solution

  int nlevels = amr_level.size();

  Real err = -1.e30;


  for (int n = 0; n < nlevels; ++n) {

    // the Castro object for this level
    Castro& castro = dynamic_cast<Castro&>(*amr_level[n]);
    Real time = castro.get_state_data(State_Type).curTime();

    const int* domain_lo = castro.geom.Domain().loVect();
    const int* domain_hi = castro.geom.Domain().hiVect();

    // the state data
    MultiFab& S = castro.get_new_data(State_Type);

    // derive the analytic solution
    auto analytic = castro.derive("analytic", time, 1);

    // if we are fourth-order, we need to convert to averages
    if (mol_order == 4 || sdc_order == 4) {
      for (MFIter mfi(*analytic); mfi.isValid(); ++mfi) {

        const Box& gbx = mfi.growntilebox(1);
        ca_make_fourth_in_place(AMREX_INT_ANYD(gbx.loVect()), AMREX_INT_ANYD(gbx.hiVect()),
                                BL_TO_FORTRAN_FAB((*analytic)[mfi]),
                                AMREX_INT_ANYD(domain_lo), AMREX_INT_ANYD(domain_hi));

      }
    }

    // compute the norm of the error
    MultiFab::Subtract(*analytic, S, Temp, 0, 1, 0);

    err = std::max(err, analytic->norm0());
    

  }

  const std::string stars(78,'*');
  amrex::Print() << stars << "\n"
                 << " diffusion problem post_simulation() \n"
                 << " L-inf error against analytic solution: " << err << "\n"
                 << stars << "\n";
}
#endif
