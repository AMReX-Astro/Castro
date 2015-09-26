#include <winstd.H>

#include "Castro.H"
#include "Castro_error_F.H"
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif

using std::string;

typedef StateDescriptor::BndryFunc BndryFunc;

// note that "StateErr" here is defined as a derived variable in
// Castro_setup.cpp and constructed to contain 3 quantities: density,
// temperature, and the first species

void
Castro::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    err_list.add("StateErr",1,ErrorRec::Special,
		 BL_FORT_PROC_CALL(CA_STATE_ERROR,ca_state_error));
}
