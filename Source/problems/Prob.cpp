/* Implementations of functions in Problem.H go here */

#include <Castro.H>

using namespace amrex;

#ifdef DO_PROBLEM_POST_SIMULATION
void Castro::problem_post_simulation(Vector<std::unique_ptr<AmrLevel> >& amr_level) {

  // this is a stub post_simulation() routine

  // you can put this in your problem directory to create a custom
  // diagnostic to run at the end of a simulation

  // all of the needed data comes in though the amr_level array,
  // which is a PArray of AmrLevel objects.  The number of levels
  // is simply the size of this array:

  // int nlevels = amr_level.size();

  // To access data, cast a level to a Castro object, e.g. for
  // level 0:

  // Castro* castro = dynamic_cast<Castro*>(&amr_level[0]);

  // then you can get the data, e.g. for state data as:

  // MultiFab& S = castro->get_new_data(State_Type);

  // and if needed, the state descriptor:

  // const StateDescriptor* desc = &castro->desc_lst[State_Type];

  // and then you can get the names of the state data as desc->name(comp)

}
#endif

