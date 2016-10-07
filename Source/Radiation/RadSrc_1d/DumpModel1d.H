#ifndef  _DumpModel1d_H_
#define  _DumpModel1d_H_

// This class is currently defined only for the 1D neutrino code, and
// Make.Castro only turns on HAS_DUMPMODEL for this case.  What is does
// is dump a file called "modelDump" containing the fluid and neutrino
// state at all exposed cells of the 1D problem.  This file can then be
// used to initialize 1D, 2D, or 3D versions of the calculation.

// Further comments are in DumpModel1d.cpp

#ifdef HAS_DUMPMODEL

#if (BL_SPACEDIM == 1)

class DumpModel {

 public:

  DumpModel();

  void dump(class Amr* parent, int force_dump = 0);

 private:

  int interval;
  int verbose;

};

#endif

#endif

#endif
