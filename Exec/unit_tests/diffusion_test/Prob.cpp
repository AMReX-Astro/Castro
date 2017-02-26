/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

void Castro::post_simulation(PArray<AmrLevel>& amr_level) {

  // compute the norm of the solution vs. the analytic solution

  int nlevels = amr_level.size();

  Real err = -1.e30;

  for (int n = 0; n < nlevels; ++n) {

    // the Castro object for this level
    Castro* castro = dynamic_cast<Castro*>(&amr_level[n]);
    Real time = castro->get_state_data(State_Type).curTime();

    // the state data
    MultiFab& S = castro->get_new_data(State_Type);

    // derive the analytic solution
    std::unique_ptr<MultiFab> analytic(castro->derive("analytic", time, 0));
    
    // compute the norm of the error
    MultiFab::Subtract(*analytic, S, Temp, 0, 1, 0);

    err = std::max(err, analytic->norm0());
    

  }
   
  std::cout << std::string(78, '*') << std::endl;
  std::cout << " diffusion problem post_simulation() " << std::endl;
  std::cout << " L-inf error against analytic solution: " << err << std::endl;
  std::cout << std::string(78, '*') << std::endl;


}
