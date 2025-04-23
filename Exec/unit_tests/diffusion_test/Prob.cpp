/* Implementations of functions in Problem.H go here */

#include <Castro.H>

using namespace amrex;

#ifdef DO_PROBLEM_POST_SIMULATION
void Castro::problem_post_simulation(Vector<std::unique_ptr<AmrLevel> >& amr_level) {

  // compute the norm of the solution vs. the analytic solution

  int nlevels = static_cast<int>(amr_level.size());

  Real L0 = -1.e30;
  Real L2 = -1.e30;

  for (int n = 0; n < nlevels; ++n) {

    // the Castro object for this level
    auto& castro = dynamic_cast<Castro&>(*amr_level[n]);
    Real time = castro.get_state_data(State_Type).curTime();

    // the state data
    MultiFab& S = castro.get_new_data(State_Type);

    // derive the analytic solution
    const int ngrow = 1;
    auto analytic = castro.derive("analytic", time, ngrow);

#ifdef TRUE_SDC
    // if we are fourth-order, we need to convert to averages
    if (sdc_order == 4) {
        FArrayBox tmp;

        for (MFIter mfi(*analytic); mfi.isValid(); ++mfi) {

            const Box& bx = mfi.tilebox();
            tmp.resize(bx, 1);
            auto tmp_arr = tmp.array();

            castro.make_fourth_in_place(bx,
                                        analytic->array(mfi),
                                        tmp_arr,
                                        domain_lo, domain_hi);

      }
    }
#endif

    // compute the norm of the error
    MultiFab::Subtract(*analytic, S, UTEMP, 0, 1, 0);

    L0 = std::max(L0, analytic->norm0());
    L2 = std::max(L2, analytic->norm2());

  }

  const std::string stars(78,'*');
  amrex::Print() << stars << "\n"
                 << " diffusion problem post_simulation() \n"
                 << " L-inf error against analytic solution: " << L0 << "\n"
                 << " L-2 error against analytic solution: " << L2 << "\n"
                 << stars << "\n";
}
#endif
