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

    BL_PROFILE_VAR("main()", pmain);

    Real dRunTime1 = ParallelDescriptor::second();

    std::cout << std::setprecision(10);

    int  max_step;
    Real strt_time;
    Real stop_time;
    ParmParse pp; 

    max_step  = -1;
    strt_time =  0.0;
    stop_time = -1.0;

    pp.query("max_step",max_step);
    pp.query("strt_time",strt_time);
    pp.query("stop_time",stop_time);

    if (strt_time < 0.0)
    {
        BoxLib::Abort("MUST SPECIFY a non-negative strt_time"); 
    }

    if (max_step < 0 && stop_time < 0.0) {
      BoxLib::Abort(
       "Exiting because neither max_step nor stop_time is non-negative.");
    }

    // Print the current date and time.

    time_t time_type;

    struct tm* time_pointer;

    time(&time_type);

    time_pointer = gmtime(&time_type);

    if (ParallelDescriptor::IOProcessor()) 
      std::cout << std::setfill('0') << "\nStarting run at "
		<< std::setw(2) << time_pointer->tm_hour << ":"
		<< std::setw(2) << time_pointer->tm_min << ":"
		<< std::setw(2) << time_pointer->tm_sec << " UTC on "
		<< std::setw(2) << time_pointer->tm_mday << "/"
		<< std::setw(2) << time_pointer->tm_mon + 1 << "/"
		<< time_pointer->tm_year + 1900 << "." << std::endl;
    
    //
    // Initialize random seed after we're running in parallel.
    //

    Amr* amrptr = new Amr;

    amrptr->init(strt_time,stop_time);

#ifdef HAS_DUMPMODEL
    DumpModel *dumpmodelptr = new DumpModel();
#endif

#ifdef HAS_XGRAPH
    XGraph1d *xgraphptr = new XGraph1d(*amrptr);
#endif

#ifdef HAS_XGRAPH
    xgraphptr->draw(amrptr->levelSteps(0),amrptr->cumTime());
#endif

    // If we set the regrid_on_restart flag and if we are *not* going to take
    //    a time step then we want to go ahead and regrid here.
    if ( amrptr->RegridOnRestart() && 
         ( (amrptr->levelSteps(0) >= max_step) ||
           (amrptr->cumTime() >= stop_time) ) )
           {
           //
           // Regrid only!
           //
           amrptr->RegridOnly(amrptr->cumTime());
           }

    Real dRunTime2 = ParallelDescriptor::second();

    while ( amrptr->okToContinue()                            &&
           (amrptr->levelSteps(0) < max_step || max_step < 0) &&
           (amrptr->cumTime() < stop_time || stop_time < 0.0) )

    {
        //
        // Do a timestep.
        //
        amrptr->coarseTimeStep(stop_time);

#ifdef HAS_DUMPMODEL
	dumpmodelptr->dump(amrptr);
#endif

#ifdef HAS_XGRAPH
	xgraphptr->draw(amrptr->levelSteps(0),amrptr->cumTime());
#endif

    }

    if (!amrptr->okToContinue() && ParallelDescriptor::IOProcessor()) {
      std::cout << "Stopping simulation because we are not OK to continue." << std::endl;
    }

#ifdef HAS_DUMPMODEL
    dumpmodelptr->dump(amrptr, 1);
    delete dumpmodelptr;
#endif

#ifdef HAS_XGRAPH
    xgraphptr->draw(amrptr->levelSteps(0),amrptr->cumTime(), 1);
#endif

#ifdef HAS_XGRAPH
    delete xgraphptr;
#endif

    // Write final checkpoint and plotfile
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
        amrptr->checkPoint();
    }

    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
        amrptr->writePlotFile();
    }

    time(&time_type);

    time_pointer = gmtime(&time_type);

    if (ParallelDescriptor::IOProcessor())
      std::cout << std::setfill('0') << "\nStarting run at "
		<< std::setw(2) << time_pointer->tm_hour << ":"
		<< std::setw(2) << time_pointer->tm_min << ":"
		<< std::setw(2) << time_pointer->tm_sec << " UTC on "
		<< std::setw(2) << time_pointer->tm_mday << "/"
		<< std::setw(2) << time_pointer->tm_mon + 1 << "/"
		<< time_pointer->tm_year + 1900 << "." << std::endl;
    
    delete amrptr;
    //
    // This MUST follow the above delete as ~Amr() may dump files to disk.
    //
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    Real dRunTime3 = ParallelDescriptor::second();

    Real runtime_total = dRunTime3 - dRunTime1;
    Real runtime_timestep = dRunTime3 - dRunTime2;

    ParallelDescriptor::ReduceRealMax(runtime_total,IOProc);
    ParallelDescriptor::ReduceRealMax(runtime_timestep,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "Run time = " << runtime_total << std::endl;
        std::cout << "Run time w/o init = " << runtime_timestep << std::endl;
    }

    if (CArena* arena = dynamic_cast<CArena*>(BoxLib::The_Arena()))
    {
        //
        // A barrier to make sure our output follows that of RunStats.
        //
        ParallelDescriptor::Barrier();
        //
        // We're using a CArena -- output some FAB memory stats.
        //
        // This'll output total # of bytes of heap space in the Arena.
        //
        // It's actually the high water mark of heap space required by FABs.
        //
        char buf[256];

        sprintf(buf,
                "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %ld",
                ParallelDescriptor::MyProc(),
                arena->heap_space_used());

        std::cout << buf << std::endl;
    }

    BL_PROFILE_VAR_STOP(pmain);
    BL_PROFILE_SET_RUN_TIME(dRunTime2);

    BoxLib::Finalize();

    return 0;
}
