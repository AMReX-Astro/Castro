
#include "Castro.H"
#include "Castro_error_F.H"
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif
#include <Castro_prob_err_F.H>

using std::string;
using namespace amrex;

typedef StateDescriptor::BndryFunc BndryFunc;

void
Castro::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //

    // This is the only uncommented line which means we will only tag on Ye / grad(Ye)
    err_list.add("Ye",1,ErrorRec::Special,ca_ye_error);

    // err_list.add("density",1,ErrorRec::Special,ca_denerror);

#if 0
    // laplac_error uses the special evaluation of the second derivative
    // and can be called with any variable.  Note that two ghost cells are needed.
    err_list.add("density",2,ErrorRec::Special,ca_laplac_error);
    err_list.add("pressure",2,ErrorRec::Special,ca_laplac_error);
    err_list.add("Temp",1,ErrorRec::Special,ca_temperror);
    err_list.add("pressure",1,ErrorRec::Special,ca_presserror);
    err_list.add("x_velocity",1,ErrorRec::Special,ca_velerror);
    err_list.add("y_velocity",1,ErrorRec::Special,ca_velerror);
#if (BL_SPACEDIM == 3)
    err_list.add("z_velocity",1,ErrorRec::Special,ca_velerror);
#endif
#endif

    // Save the number of built-in functions; this will help us
    // distinguish between those, and the ones the user is about to add.

    num_err_list_default = err_list.size();

#include <Castro_prob_err_list.H>

}
