#ifndef  _XGraph1d_H_
#define  _XGraph1d_H_   "@(#)XGraph1d.H	1.1\t10/2/92"

#if (BL_SPACEDIM == 1)

#include <iostream>
using std::istream;
using std::ostream;

//#include <Box.H>
//#include <PArray.H>
#include <MultiFab.H>
//#include <FArrayBox.H>
//#include <BaseFab.H>

#ifdef _MEM_CHECK_
#include "mem.H"
#endif

class XGraph1d;

class XGFrame {
   public:
      friend class XGraph1d;
      friend ostream& operator << (ostream &os, XGraph1d &cont);

      XGFrame(std::string& fn, std::string& var_nm, int freq, int lev);
      XGFrame(std::string& fn, std::string& var_nm, std::string& plot_nm, int freq, int lev);

     ~XGFrame();
     std::string* varPtr() { return &var_name; }

   private:
      std::string    file_name;   // root name of variable
      std::string    var_name;    // name of quantity to draw
      std::string    plot_name;   // name of plotfile quantity to draw
      int        interval;    // how often to plot
      int        level;       // maximum level of plot
};

// forward ref
class Amr;

static const int XGPtMXGY=256; // Man is weak...
class XGPt {
public:
   Real x;
   Real y[XGPtMXGY];
};
class XGGrid {
public:
   std::list<XGPt> point;
};
inline bool operator<(const XGGrid& lhs,const XGGrid& rhs) {
   return lhs.point.begin()->x<rhs.point.begin()->x;
}


class XGraph1d {
   friend ostream& operator << (ostream &os, XGraph1d &cont);
   public:
      XGraph1d(Amr& amrsys);

      ~XGraph1d();

      // manually clear XGraph1d space
      void clear();

      void addVar(std::string& file_nm, std::string& var_nm, int freq, int lev);

      void draw(int nstep, Real time, int force_draw = 0);

      void dumpXGraph(MultiFab& q, MultiFab& qf, int comp,
         std::list<XGGrid>& xgg, int& lev, int& maxlev);

   private:

      // each file contains a single variable
      void draw_single(int nstep, Real time, int force_draw);
      // all variables at the time are stored in file
      void draw_all(int nstep, Real time, int force_draw);

      std::list<XGFrame> frames;
      Amr *amrptr;

      static int verbose;

      std::string format;

      int use_xmgrace_legend;
      int use_xmgrace_title;
};


#endif

#endif
