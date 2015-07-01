#include <winstd.H>

#include "Castro.H"
#include "Castro_F.H"
#ifdef RADIATION
# include "Radiation.H"
# include "RAD_F.H"
#endif

using std::string;

static Box the_same_box (const Box& b) { return b; }
static Box grow_box_by_one (const Box& b) { return BoxLib::grow(b,1); }
static Box grow_box_by_two (const Box& b) { return BoxLib::grow(b,2); }

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
