/* Implementations of functions in Problem.H go here */

#include "Castro.H"
#include "Castro_F.H"
#include "Problem_F.H"

void Castro::post_simulation(PArray<AmrLevel>& amr_level) {

  // here we have all of the amr_level data -- this is essentially a
  // Castro object for each level

  // get the variable info
  Castro* castro = dynamic_cast<Castro*>(&amr_level[0]);
  MultiFab& S = castro->get_new_data(State_Type);

  const StateDescriptor* desc = &castro->desc_lst[State_Type];

  int nvar = S.nComp();

  int nlevels = amr_level.size();

  // storage of min / max of each state variable over levels
  Array<Real> vmin(nvar);
  Array<Real> vmax(nvar);

  for (int n = 0; n < nlevels; ++n) {

    Castro* castro = dynamic_cast<Castro*>(&amr_level[n]);
    MultiFab& S = castro->get_new_data(State_Type);

    for (int k = 0; k < S.nComp(); ++k) {
      if (n == 0) {
	vmin[k] = S.min(k);
	vmax[k] = S.max(k);
      } else {
	vmin[k] = std::min(vmin[k], S.min(k));
	vmax[k] = std::max(vmax[k], S.max(k));
      }
    }
  }
   
  if (ParallelDescriptor::IOProcessor())  {
    for (int k = 0; k < S.nComp(); ++k) {
      std::cout << desc->name(k) << " " << vmin[k] << " " << vmax[k] << std::endl;
    }
  }

}
