#include <winstd.H>

#include "Castro.H"
#include "Castro_error_F.H"
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif

using std::string;

typedef StateDescriptor::BndryFunc BndryFunc;

void
Castro::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //

    // This routine uses the special evaluation of the second derivative
    //   and can be called with any variable.  Note that two ghost cells are needed.
//  err_list.add("density",2,ErrorRec::Special,ca_laplac_error);
//  err_list.add("pressure",2,ErrorRec::Special,ca_laplac_error);

    err_list.add("density",1,ErrorRec::Special,ca_denerror);
//  err_list.add("Temp",1,ErrorRec::Special,ca_temperror);
    err_list.add("pressure",1,ErrorRec::Special,ca_presserror);
    err_list.add("x_velocity",1,ErrorRec::Special,ca_velerror);

    err_list.add("entropy",1,ErrorRec::Special,ca_enterror);

#ifdef RADIATION
    if (do_radiation && !Radiation::do_multigroup) {
      err_list.add("rad",1,ErrorRec::Special,ca_raderror);
    }
#ifdef NEUTRINO
    err_list.add("Ye",1,ErrorRec::Special,ca_yeerror);
#if (BL_SPACEDIM == 1)
    err_list.add("grav_x",1,ErrorRec::Special,ca_graverror);
#else
    err_list.add("maggrav",1,ErrorRec::Special,ca_graverror);
#endif
#endif

#endif
}
