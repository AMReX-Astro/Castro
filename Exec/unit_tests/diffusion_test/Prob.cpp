/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

void Castro::post_simulation(PArray<AmrLevel>& amr_level) {

  // compute the norm of the solution vs. the analytic solution

  int nlevels = amr_level.size();

  for (int n = 0; n < nlevels; ++n) {

    // the Castro object for this level
    Castro* castro = dynamic_cast<Castro*>(&amr_level[n]);
    Real time = castro->get_state_data(State_Type).curTime();

    // the state data
    MultiFab& S = castro->get_new_data(State_Type);

    // derive the analytic solution
    MultiFab *analytic = castro->derive("analytic", time, 0);
    
    // compute the norm of the error
    std::cout << "norms: " << S.norm0(Temp, 0) << " " << analytic->norm0() << std::endl;
    
    // cleanup
    delete analytic;

  }
   

}
