#include <winstd.H>

#include <new>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>

#ifndef WIN32
#include <unistd.h>
#endif

#include <CArena.H>
#include <REAL.H>
#include <Utility.H>
#include <IntVect.H>
#include <Box.H>
#include <Amr.H>
#include <ParmParse.H>
#include <ParallelDescriptor.H>
#include <AmrLevel.H>

#include <time.h>

#ifdef HAS_DUMPMODEL
#include <DumpModel1d.H>
#endif

#ifdef HAS_XGRAPH
#include <XGraph1d.H>
#endif

#include "Castro_io.H"


extern "C"
{
   void do_burn();
}


std::string inputs_name = "";

int
main (int   argc,
      char* argv[])
{

    //
    // Make sure to catch new failures.
    //
    BoxLib::Initialize(argc,argv);

    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
	inputs_name = argv[1];
      }
    }

    do_burn();

    BoxLib::Finalize();

    return 0;
}
