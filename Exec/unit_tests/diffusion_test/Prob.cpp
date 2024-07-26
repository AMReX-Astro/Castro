/* Implementations of functions in Problem.H go here */

#include <Castro.H>

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

    auto domain_lo = castro.geom.Domain().loVect3d();
    auto domain_hi = castro.geom.Domain().hiVect3d();

    // the state data
    MultiFab& S = castro.get_new_data(State_Type);

    // derive the analytic solution
    auto analytic = castro.derive("analytic", time, 1);

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

    err = std::max(err, analytic->norm0());


  }

  const std::string stars(78,'*');
  amrex::Print() << stars << "\n"
                 << " diffusion problem post_simulation() \n"
                 << " L-inf error against analytic solution: " << err << "\n"
                 << stars << "\n";
}
#endif
