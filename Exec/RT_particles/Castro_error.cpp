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
    err_list.add("density",2,ErrorRec::Special,
  	       BL_FORT_PROC_CALL(CA_LAPLAC_ERROR,ca_laplac_error));
//  err_list.add("pressure",2,ErrorRec::Special,
//	       BL_FORT_PROC_CALL(CA_LAPLAC_ERROR,ca_laplac_error));

    err_list.add("density",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(CA_DENERROR,ca_denerror));
//  err_list.add("Temp",1,ErrorRec::Special,
//		 BL_FORT_PROC_CALL(CA_TEMPERROR,ca_temperror));
    err_list.add("pressure",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(CA_PRESSERROR,ca_presserror));
    err_list.add("x_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(CA_VELERROR,ca_velerror));
#if (BL_SPACEDIM >= 2)
    err_list.add("y_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(CA_VELERROR,ca_velerror));
#endif
#if (BL_SPACEDIM == 3)
    err_list.add("z_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(CA_VELERROR,ca_velerror));
#endif

//   err_list.add("entropy",1,ErrorRec::Special,
//		 BL_FORT_PROC_CALL(CA_ENTERROR,ca_enterror));

#ifdef RADIATION
    if (do_radiation && !Radiation::do_multigroup) {
      err_list.add("rad",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(CA_RADERROR,ca_raderror));
    }
#endif

#ifdef LEVELSET
    err_list.add("LSphi",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(CA_LSERROR,ca_lserror));
#endif
}
