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
    err_list.add("x_velocity",2,ErrorRec::Special,
  	       BL_FORT_PROC_CALL(CA_LAPLAC_ERROR,ca_laplac_error));

    // These routines access the variable indicated and test on value or gradient.
    //   and can be called with any variable.  Note that two ghost cells are needed.
    err_list.add("density",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(CA_DENERROR,ca_denerror));
    err_list.add("x_velocity",1,ErrorRec::Special,
                   BL_FORT_PROC_CALL(CA_VELERROR,ca_velerror));
}
