#include <winstd.H>

#include "Castro.H"
#include "Castro_error_F.H"

using std::string;

static Box the_same_box (const Box& b) { return b; }
static Box grow_box_by_one (const Box& b) { return BoxLib::grow(b,1); }
static Box grow_box_by_two (const Box& b) { return BoxLib::grow(b,2); }

typedef StateDescriptor::BndryFunc BndryFunc;

void
Castro::ErrorSetUp ()
{
    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //

  //  err_list.add("Temp",2,ErrorRec::Special,
  //  	       BL_FORT_PROC_CALL(CA_LAPLAC_ERROR,ca_laplac_error));
  err_list.add("Temp",1,ErrorRec::Special,
	       BL_FORT_PROC_CALL(CA_TEMPERROR,ca_temperror));
}
