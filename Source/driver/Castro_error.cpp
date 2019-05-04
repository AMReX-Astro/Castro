
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

    BL_PROFILE("Castro::ErrorSetUp()");

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //

    // This routine uses the special evaluation of the second derivative
    //   and can be called with any variable.  Note that two ghost cells are needed.
//  err_list.add("density",2,ErrorRec::Special,ca_laplac_error);
//  err_list.add("pressure",2,ErrorRec::Special,ca_laplac_error);

    err_list_names.push_back("density");
    err_list_ng.push_back(1);

    err_list_names.push_back("Temp");
    err_list_ng.push_back(1);

    err_list_names.push_back("pressure");
    err_list_ng.push_back(1);

    err_list_names.push_back("x_velocity");
    err_list_ng.push_back(1);

#if (BL_SPACEDIM >= 2)
    err_list_names.push_back("y_velocity");
    err_list_ng.push_back(1);
#endif

#if (BL_SPACEDIM == 3)
    err_list_names.push_back("z_velocity");
    err_list_ng.push_back(1);
#endif

#ifdef REACTIONS
    err_list_names.push_back("t_sound_t_enuc");
    err_list_ng.push_back(0);

    err_list_names.push_back("enuc");
    err_list_ng.push_back(0);
#endif

#ifdef RADIATION
    if (do_radiation && !Radiation::do_multigroup) {
      err_list_names.push_back("rad");
      err_list_ng.push_back(1);
    }
#endif

    // Save the number of built-in functions; this will help us
    // distinguish between those, and the ones the user is about to add.

    num_err_list_default = err_list_names.size();

#include <Castro_prob_err_list.H>

}
