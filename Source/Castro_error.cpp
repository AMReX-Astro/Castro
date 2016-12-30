
#include "Castro.H"
#include "Castro_error_F.H"
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif

using std::string;
using namespace amrex;

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
    err_list.add("Temp",1,ErrorRec::Special,ca_temperror);
    err_list.add("pressure",1,ErrorRec::Special,ca_presserror);
    err_list.add("x_velocity",1,ErrorRec::Special,ca_velerror);
#if (BL_SPACEDIM >= 2)
    err_list.add("y_velocity",1,ErrorRec::Special,ca_velerror);
#endif
#if (BL_SPACEDIM == 3)
    err_list.add("z_velocity",1,ErrorRec::Special,ca_velerror);
#endif

//   err_list.add("entropy",1,ErrorRec::Special,ca_enterror);

#ifdef REACTIONS
    err_list.add("t_sound_t_enuc",0,ErrorRec::Special,ca_nucerror);
#endif

#ifdef RADIATION
    if (do_radiation && !Radiation::do_multigroup) {
      err_list.add("rad",1,ErrorRec::Special,ca_raderror);
    }
#endif

}
