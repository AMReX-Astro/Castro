#if    (BL_SPACEDIM == 1)

#include "XGraph1d.H"

#include <ParallelDescriptor.H>
#include <vector>
#include <map>

#include <Amr.H>
#include <AmrLevel.H>
#include <ParmParse.H>
#include <Array.H>
#include <Misc.H>

#ifdef _NO_INLINE_
#       define  INLINE
#       include "XGraph1d.I"
#endif

#include <iomanip>
#include <sstream>
#include <cstdio>
#include <cctype>
#include <algorithm>
//using std::setprecision;
//using namespace std;
int XGraph1d::verbose = false;


XGFrame::XGFrame(std::string& file_nm, std::string& var_nm, int freq, int lev)
{
   file_name = file_nm;
   var_name = var_nm;
   plot_name = var_nm;
   interval   = freq;
   level      = lev;
}

XGFrame::XGFrame(std::string& file_nm, std::string& var_nm, std::string& plot_nm, int freq, int lev)
{
   file_name = file_nm;
   var_name = var_nm;
   plot_name = plot_nm;
   interval   = freq;
   level      = lev;
}

XGFrame::~XGFrame()
{
}

//   ############################################################
//   ##### XGraph1d members
//   ############################################################

XGraph1d::XGraph1d(Amr& amrsys ){
	amrptr = &amrsys;

    //parse input file
    ParmParse pp("xgraph");
    pp.query("v",verbose); pp.query("verbose",verbose);
    int n = pp.countname("graph");
    if (n>XGPtMXGY) BoxLib::Error("xgraph hardwired max number of vars");
    int k;
    int freq0 = -1, lev0 = -1, newinput = 0;
    std::string fname0;

    use_xmgrace_legend = 1;
    pp.query("use_xmgrace_legend",use_xmgrace_legend);
    use_xmgrace_title = 1;
    pp.query("use_xmgrace_title",use_xmgrace_title);

    if (n>0) {
       std::string full_format="ult";
       pp.query("format",full_format);
       format=full_format.substr(0,3);
       std::transform(format.begin(),format.end(),format.begin(),(int(*)(int)) std::tolower);
       if(ParallelDescriptor::IOProcessor()) {
          std::cout << "###### XGraph (xgraph:*) parameters ######" << std::endl;
          std::cout << "format    = " << format << std::endl;
          if (newinput) {
            std::cout << "filename  = " << fname0 << std::endl;
            std::cout << "frequency = " << freq0 << std::endl;
            std::cout << "level     = " << lev0 << std::endl;
          }
          std::cout << "verbose   = " << verbose << std::endl;
          std::cout << "##########################################" << std::endl;
          std::cout << std::endl;
       }
    }
    for (k = 0; k < n; k++) {
      int freq = -1, lev=-1;
      std::string fname;
      std::string vname;
      if (!newinput) { // old xgraph input format
        pp.getkth("graph",k,freq,2);
        if (freq > 0) {
            pp.getkth("graph",k,fname,0);
            pp.getkth("graph",k,vname,1);
            pp.getkth("graph",k,freq,2);
            pp.getkth("graph",k,lev,3);
            if (lev>-1&&freq>-1&&ParallelDescriptor::IOProcessor()) std::cout << "XGraph::plots at the same time will have the same refinement level!" << std::endl;
        } else continue;
      } else {
        pp.getkth("graph",k,vname);
        fname=fname0;
        freq=freq0;
        lev=lev0;
      }
      if (vname!="ALL") {
	     XGFrame cf = XGFrame(fname, vname, freq, lev);
	     frames.push_back(cf);
      } else {

        const DescriptorList& desc_lst = AmrLevel::get_desc_lst();
        for (int typ = 0; typ < desc_lst.size(); typ++) {
          for (int comp = 0; comp < desc_lst[typ].nComp();comp++) {
            std::string var_name = desc_lst[typ].name(comp);
            if (Amr::isStatePlotVar(var_name)) {
               if(var_name.substr(0,4)!="+par") {
                 // remove invalid characters from plot/file names
                 std::string this_name=var_name;
                 std::string alphanum("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_");
                 std::string::size_type pos = this_name.find_first_not_of(alphanum);
                 while (pos != std::string::npos) {
                   this_name.erase(pos,1);
                   pos = this_name.find_first_not_of(alphanum);
                 }
                 std::string plot_name = var_name, file_name=fname;
                 if (format=="gnu") plot_name = this_name;
                 else if (format=="ult") file_name += std::string("_") + this_name;
                 XGFrame cf = XGFrame(file_name, var_name, plot_name, freq, lev);
                 frames.push_back(cf);
               }
            }
          }
        }

        DeriveList& derive_lst = AmrLevel::get_derive_lst();
        const std::list<DeriveRec>& dlist = derive_lst.dlist();

        for (std::list<DeriveRec>::const_iterator it = dlist.begin();
             it != dlist.end();
             ++it)
        {
            std::string var_name = it->name();

            if (Amr::isDerivePlotVar(var_name)) {
               if (var_name.substr(0,4)!="+par") {
                 // remove invalid characters from plot/file names
                 std::string this_name=var_name;
                 std::string alphanum("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_");
                 std::string::size_type pos = this_name.find_first_not_of(alphanum);
                 while (pos != std::string::npos) {
                   this_name.erase(pos,1);
                   pos = this_name.find_first_not_of(alphanum);
                 }
                 std::string plot_name = var_name, file_name=fname;
                 if (format=="gnu") plot_name = this_name;
                 else if (format=="ult") file_name += std::string("_") + this_name;
                 XGFrame cf = XGFrame(file_name, var_name, plot_name, freq, lev);
                 frames.push_back(cf);
               }
            }
        }

      }
    }
}

XGraph1d::~XGraph1d()
{
   clear();
}

void XGraph1d::clear()
{
//   ListIterator<XGFrame*> li(frames);

	frames.std::list<XGFrame>::clear();
}

std::ostream& operator << (std::ostream &os, XGraph1d &xgr)
{
   for (std::list<XGFrame>::iterator li = xgr.frames.begin();
         li != xgr.frames.end(); li++) {

        os << "XGraph1d  " << li->var_name
           << "  " << li->interval
           << "  " << li->level << std::endl;
   }
   return os;
}

void XGraph1d::addVar(std::string& file_nm, std::string& var_nm, int freq, int lev)
{
   XGFrame c = XGFrame(file_nm,var_nm,freq,lev);
   frames.push_back(c);
}

void XGraph1d::draw(int nstep, Real time, int force_draw)
{
   if (frames.size()==0) return;
   std::list<XGFrame>::iterator lib=frames.begin();
   // determine whether any variable will be plotted this step
   // (all variables are derived whenever any one is plotted)
   // create vector of filenames and mappings from variable names
   int draw = 0;
   int maxlev = 0;
   std::vector<std::string> fname_vec;
   std::map<std::string,int> fname_map;
   for(std::list<XGFrame>::iterator li = lib; li != frames.end(); li++) {
      int myinterval = li->interval;
      std::string myvarname = li->var_name, myfilename = li->file_name;
      if ((myinterval > 0) && (nstep%myinterval == 0)||force_draw) {
         bool fname_used=false;
         int findex = 0;
         for (std::vector<std::string>::iterator vi=fname_vec.begin(); vi!=fname_vec.end(); ++vi,++findex) {
            if (myfilename==*vi) {
               fname_used=true;
               break;
            }
         }
         if (fname_map.count(myvarname)==0) {
           fname_map[myvarname]=findex;
         } else {
           if(ParallelDescriptor::IOProcessor()) std::cerr << "XGraph::draw::duplicate specification of variable: " << myvarname << std::endl;
           abort();
         }
         draw=1;
         if (!fname_used) fname_vec.push_back(myfilename);
         if (maxlev>-1){ 
           int mylev = li->level;
           if (mylev>maxlev||mylev<0) maxlev=mylev;
         }
      }
   }
   // return if no variables will be plotted
   if(!draw) return;
   int amrlev = amrptr->finestLevel();
   if(maxlev < 0||maxlev>amrlev) maxlev = amrlev;
   Array<MultiFab*> soln(amrlev+1);
   for(int lev = 0; lev <= amrlev; lev++)
      soln[lev]=new MultiFab(amrptr->getLevel(lev).boxArray(),frames.size(),0);
   // set up arrays of output streams/filenames
   BL_ASSERT(fname_vec.size()<=XGPtMXGY);
   std::ofstream os[XGPtMXGY];
   std::vector<std::string> fname(fname_vec.size());
   // write headers for each file
   std::string xvar = (amrptr->Geom(0).Coord()==0) ? "x" : "r";
   for(int i=0; i<fname_vec.size(); ++i) {  
     std::ostringstream ofname;
     std::string xgfname = fname_vec[i];
     if (format=="ult") {
       int i1 = fname_vec[i].find("ALL_");
       if (i1==0) {
         xgfname.erase(0,4);
       }
     }
     ofname << xgfname << "_" << std::setfill('0') << std::setw(4) << nstep;
     if (format=="ult") {
       ofname << ".ultra";
     } else if (format=="xmg") {
       ofname << ".xmgr";
     } else {
       ofname << "." << format;
     }
     fname[i]=ofname.str();
     if(verbose&&ParallelDescriptor::IOProcessor()) std::cout << "Writing data " << fname[i] << std::endl;
     if (ParallelDescriptor::IOProcessor()) {
       os[i].open(fname[i].c_str(),std::ios::out);
       if (format=="gnu") {
         os[i] << "; set title \"step = " << nstep << ", t = " << time;
         if (maxlev<amrlev) os[i] << " level=" << maxlev;
         os[i] << "\"" << std::endl;
         os[i] << "; set datafile comment ';'" << std::endl;
         os[i] << "; _(x)=column(x)" << std::endl;
         os[i] << std::setprecision(12);
         os[i] << "; set xlabel \"" << xvar << "\"" << std::endl;
       } else if (format == "pla") {
        // Dont write anything -- this is used for reading in by another code
       } else {
         os[i].setf(std::ios::scientific);
         if (format=="xmg") {
           if (use_xmgrace_legend) {
             os[i] << std::setprecision(6);
             if (use_xmgrace_title) {
               os[i] << "@TITLE \"step = " << nstep << ", t = " << time;
	     } else {
               os[i] << "# step = " << nstep << ", t = " << time;
	     }
             os[i] << std::setprecision(12);
             if (maxlev<amrlev) os[i] << " level=" << maxlev;
             os[i] << "\"" << std::endl;
           }
           os[i] << "@xaxis label \"" << xvar << "\"" << std::endl;
         } else {
           os[i] << " " << std::endl;
           os[i] << std::setprecision(3);
           int i1 = fname_vec[i].find("ALL_");
           std::string varname = fname_vec[i];
           if (i1==0) {
             varname.erase(0,4);
           }
           os[i] << "%# " << varname << ",step=" << nstep << ",t=" << time << std::endl;
         }
         os[i] << std::setprecision(12);
         os[i] << " " << std::endl;
       }
     }
   }
   int cntall=0;
   std::vector<int> cnt(fname_vec.size(),0);
   std::vector<std::string> plot_string(fname_vec.size());
   for(std::list<XGFrame>::iterator li = lib; li != frames.end(); li++,cntall++) {
      // determine whether variable is plotted and appropriate output stream
      int findex = fname_map.count(li->var_name) ? fname_map[li->var_name] : -1;
      if (findex>-1&&ParallelDescriptor::IOProcessor()) {
        if (format=="gnu") {
          std::string li_plot_name = li->plot_name;
          if (li_plot_name[0]=='+') li_plot_name.erase(0,1);
          os[findex] << "; " << li_plot_name << '=' << cnt[findex]+2 << std::endl;
          // set up plot variables in gnuplot file
          if (cnt[findex]==0) plot_string[findex] = "; plot ";
          else plot_string[findex] += ", ";
          plot_string[findex] += "'" + fname[findex] + "'u 1:(_(" + li_plot_name + "))";
        } else {
          if (format=="xmg") {
            // write out legend for xmgrace files
            os[findex] << "@S" << cnt[findex] << " LINESTYLE 1" << std::endl; 
            if(use_xmgrace_legend) {
              os[findex] << "@LEGEND STRING " << cnt[findex] << " \"" << li->plot_name << "\"" << std::endl;
            } else {
              int oldprecision = os[findex].precision();
              os[findex] << std::setprecision(6);
              if (use_xmgrace_title) {
                os[findex] << "@TITLE \"" << li->plot_name << " at step = " << nstep << ", t = " << time << "\"" << std::endl;
	      } else {
                os[findex] << "# " << li->plot_name << " at step = " << nstep << ", t = " << time << "\"" << std::endl;
	      }
              os[findex] << std::setprecision(oldprecision);
              os[findex] << "@yaxis label \" " << li->plot_name << "\"" << std::endl;
            }
          }
        }
        cnt[findex]++;
      }

      // derive data and compute min and max for all variables
      double v_min = 1.e30;
      double v_max = -1.e30;

      for(int lev = 0; lev <= maxlev; lev++) {
         // void derive() doesn't work correctly for derived variables
         // amrptr->getLevel(lev).derive(li->var_name,time,*soln[lev],cntall);
         MultiFab* temp = amrptr->derive(li->var_name,time,lev,0);
         MultiFab::Copy(*soln[lev],*temp,0,cntall,1,0);
         Real g_min = soln[lev]->min(cntall);
         Real g_max = soln[lev]->max(cntall);
         v_min = Min(v_min,g_min);
         v_max = Max(v_max,g_max);
	 delete temp;
      }

      ParallelDescriptor::ReduceRealMin(v_min);
      ParallelDescriptor::ReduceRealMax(v_max);
      // only display min/max if variable is plotted
      if (verbose&&findex>-1&&ParallelDescriptor::IOProcessor())
         std::cout << li->var_name << " min,max= " << v_min << ", " << v_max << std::endl;
   }
   std::list<XGGrid> XGGlist;
   for(int lev = 0; lev <= maxlev; lev++) {
      dumpXGraph(*soln[lev],*soln[Min(lev+1,maxlev)],frames.size(),XGGlist,lev,maxlev);
   }

   if (ParallelDescriptor::IOProcessor()) {
     for (int findex=0;findex<fname_vec.size();++findex) {
       if (format=="gnu") {
         // write out plot string to gnuplot file
         os[findex] << plot_string[findex] << std::endl;
         os[findex] << "; exit" << std::endl;
       }
     }
     for(std::list<XGGrid>::iterator it=XGGlist.begin();it!=XGGlist.end();) {
       if(it->point.empty()) {
         it=XGGlist.erase(it);
       } else {
         ++it;
       }
     }
     XGGlist.sort();
     for(std::list<XGGrid>::iterator it=XGGlist.begin();it!=XGGlist.end();++it) {
        for(std::list<XGPt>::iterator itt=it->point.begin();itt!=it->point.end();++itt) {
           // write out x-variable to each output stream
           for (int findex=0;findex<fname_vec.size();++findex) os[findex] << itt->x;
           int i=0;
           for(std::list<XGFrame>::iterator li = lib; li != frames.end(); li++,i++) {
             // write out requested variables to appropriate output streams
             int findex = fname_map.count(li->var_name) ? fname_map[li->var_name] : -1;
             if (findex>-1) os[findex] << " " << itt->y[i];
           }
           // newline for each output stream
           for (int findex=0;findex<fname_vec.size();++findex) os[findex] << std::endl;
        }
     }
   }
   for(int lev = 0; lev <= amrlev; lev++)
      delete soln[lev];
}
static const Real blackout = 4.444444444e+44;
void XGraph1d::dumpXGraph(MultiFab& q, MultiFab& qf,int comp,
                          std::list<XGGrid>& xgg, int& lev, int& maxlev)
{
   // first black out regions that are covered by fine grids
   if(lev < maxlev) {
      int ratio = amrptr->refRatio(lev)[0];
      for (MFIter mfi(q); mfi.isValid(); ++mfi) {
         const int i = mfi.index();
         Box cbox(q.box(i));
         for(int j=0;j<qf.boxArray().size();j++) {
            Box fbox(qf.box(j));
            fbox.coarsen(ratio);
            fbox &= cbox;
            if(fbox.ok()) {
               q[i].setVal(blackout,fbox,0);
            }
         }
      }
   }
   // now write out portions that are not blacked out
   Real xl,dx = (amrptr->Geom(lev)).CellSize(0);
   
   std::vector< std::list<XGGrid> > xgg_lev(q.boxArray().size());
   for (MFIter mfi(q); mfi.isValid(); ++mfi) {
      const int i = mfi.index();
      const int* nx = q[i].length();
      amrptr->Geom(lev).CellCenter(q[i].box().smallEnd(),&xl);
      bool is_data=false;
      XGGrid xggtmp;
      xgg_lev[i].push_back(xggtmp);
      for(int j = 0; j < *nx; j++) {
         if(*(q[i].dataPtr(0)+j) != blackout) {
            XGPt ptmp;
            ptmp.x=xl+dx*j;
            for(int c=0;c<comp;++c) ptmp.y[c]=*(q[i].dataPtr(c)+j);
            xgg_lev[i].back().point.push_back(ptmp);
            is_data=true;
         } else if(is_data) {
            xgg_lev[i].push_back(xggtmp);
            is_data=false;
         }
      }
      is_data=false;
   }
   const int IOProc = ParallelDescriptor::IOProcessorNumber();
   const DistributionMapping dm=q.DistributionMap();
   for(int i=0;i<q.boxArray().size();i++) {
      int q_proc=dm[i];
      int grids = xgg_lev[i].size();
      ParallelDescriptor::ReduceIntMax(grids);
      if (IOProc!=q_proc) {
        int pts;
        if (ParallelDescriptor::IOProcessor()) {
          xgg_lev[i].resize(grids);
          for (std::list<XGGrid>::iterator it=xgg_lev[i].begin();it!=xgg_lev[i].end();++it) {
            ParallelDescriptor::Recv(&pts,1,q_proc,0);
            std::list<XGPt> xggtmp(pts);
            it->point=xggtmp;
            for (std::list<XGPt>::iterator itt=it->point.begin();itt!=it->point.end();++itt) {
              ParallelDescriptor::Recv(&(itt->x),1,q_proc,1);
              ParallelDescriptor::Recv(itt->y,comp,q_proc,2);
            }
          }
        } else if (ParallelDescriptor::MyProc()==q_proc) {
          for (std::list<XGGrid>::iterator it=xgg_lev[i].begin();it!=xgg_lev[i].end();++it) {
            pts = it->point.size();
            ParallelDescriptor::Send(&pts,1,IOProc,0);
            for (std::list<XGPt>::iterator itt=it->point.begin();itt!=it->point.end();++itt) {
              ParallelDescriptor::Send(&(itt->x),1,IOProc,1);
              ParallelDescriptor::Send(itt->y,comp,IOProc,2);
            }
          }
        }
      }
      for (std::list<XGGrid>::iterator it=xgg_lev[i].begin();it!=xgg_lev[i].end();++it) {
        xgg.push_back(*it);
      }
   }
}
#endif
