
// This class is currently defined only for the 1D neutrino code, and
// Make.Castro only turns on HAS_DUMPMODEL for this case.  What is does
// is dump a file called "modelDump" containing the fluid and neutrino
// state at all exposed cells of the 1D problem.  This file can then be
// used to initialize 1D, 2D, or 3D versions of the calculation.  The
// underlying assumption is that a collapse calculation can be run for
// a while in 1D, then later switched to a higher-dimensional run once
// it is time for the flow to become unstable and complicated.

// The dump works when the 1D calculation uses parallel AMR.  Output
// to the file is serialized through the I/O processor.  Exposed
// portions of all levels are sorted so that cells are written out in
// radial order.  The format of the fluid portion of the state is
// based on that of the modelInput file used to initialize neutrino
// calculations from another code.

// There are assumptions that make sense only for neutrino calculations.
// In particular, Ye is printed, and also the entire radiation state
// with Radiation::nGroups components.  With minor changes, though, this
// class could be modified to do similar dumps for other types of Castro
// calculations.

#ifdef HAS_DUMPMODEL

#if (BL_SPACEDIM == 1)

#include "Radiation.H"
#include "DumpModel1d.H"

#include <AMReX_ParmParse.H>

#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

DumpModel::DumpModel() : interval(1000), verbose(1)
{
  ParmParse pp("dumpmodel");
  pp.query("interval", interval);
  pp.query("v", verbose);
  pp.query("verbose", verbose);
}

void DumpModel::dump(Amr* parent, int force_dump)
{

  // Exit if it isn't the right time to do the dump

  if (force_dump == 0 &&
      parent->levelSteps(0) % interval != 0) {
    return;
  }

  // Create flattened list of ordered, disjoint boxes from the grid hierarchy

  list<Box> bxlist;
  list<int> lnlist;

  vector<int> rat(parent->finestLevel(), 1);

  for (int ln = 0; ln <= parent->finestLevel(); ln++) {

    // rat[lnp] will be ratio between levels lnp and ln

    if (ln > 0) {
      for (int lnp = 0; lnp < ln; lnp++) {
        rat[lnp] *= parent->refRatio(ln - 1)[0];
      }
    }

    // Insert grids of level ln into list, trimming coarser grids as needed

    const BoxArray& grids = parent->boxArray(ln);
    for (int igr = 0; igr < grids.size(); igr++) {
      const Box& reg = grids[igr];
      list<Box>::iterator bi = bxlist.begin();
      list<int>::iterator li = lnlist.begin();
      for ( ; bi != bxlist.end(); ) {
        if (*li == ln) {
          if (bi->smallEnd(0) > reg.bigEnd(0)) {
            // Insert reg before bi and break
            bxlist.insert(bi, reg);
            lnlist.insert(li, ln);
            break;
          }
          // Increment and continue loop as reg not inserted yet
        }
        else { // *li < ln
          Box bx(amrex::refine(*bi,rat[*li]));
          if (bx.bigEnd(0) >= reg.smallEnd(0)) {
            if (bx.smallEnd(0) < reg.smallEnd(0) &&
                bx.bigEnd(0) <= reg.bigEnd(0)) {
              // Trim overlap with reg from *bi leaving low end of *bi
              bx.setBig(0, reg.smallEnd(0) - 1);
              bx.coarsen(rat[*li]);
              *bi = bx;
              // Increment and continue loop as reg not inserted yet
            }
            else if (bx.smallEnd(0) >= reg.smallEnd(0) &&
                     bx.bigEnd(0) <= reg.bigEnd(0)) {
              // Remove *bi from list as reg covers it completely
              list<Box>::iterator bt = bi;
              list<int>::iterator lt = li;
              ++bi;
              ++li;
              bxlist.erase(bt);
              lnlist.erase(lt);
              // Continue loop as reg not inserted yet; but bypass
              // increment at end of loop since we've done it already.
              continue;
            }
            else {
              // To get this far (bx.bigEnd(0) > reg.bigEnd(0)) must be true
              if (bx.smallEnd(0) < reg.smallEnd(0)) {
                // Trim overlap with reg from middle of *bi, extract low frag
                Box bxx(bx);
                bxx.setBig(0, reg.smallEnd(0) - 1);
                bxx.coarsen(rat[*li]);
                // Insert this low fragment of *bi before remainder of *bi
                bxlist.insert(bi, bxx);
                lnlist.insert(li, *li);
              }
              if (bx.smallEnd(0) <= reg.bigEnd(0)) {
                // Trim overlap with reg leaving end high end of *bi
                bx.setSmall(0, reg.bigEnd(0) + 1);
                bx.coarsen(rat[*li]);
                *bi = bx;
              }
              // Insert reg before high end of *bi and break
              bxlist.insert(bi, reg);
              lnlist.insert(li, ln);
              break;
            }
          }
        }

        // Increment and continue loop as reg not inserted yet
        ++bi;
        ++li;

      } // end of loop over bi, li

      if (bi == bxlist.end()) {
        // We made it to the end without finding another place to insert reg
        bxlist.push_back(reg);
        lnlist.push_back(ln);
      }
    }
  }

  // bxlist now contains an ordered list of disjoint boxes representing
  // the exposed portions of all levels of the hierarchy.  Each box is
  // at the resolution of its respective level, and the corresponding
  // entry of lnlist contains the level number.

  if (ParallelDescriptor::IOProcessor() && verbose >= 2) {
    cout << "Printing disjoint box list" << endl;

    list<Box>::iterator bi = bxlist.begin();
    list<int>::iterator li = lnlist.begin();
    for ( ; bi != bxlist.end(); ++bi, ++li) {
      cout << "Level " << *li << ", box " << *bi << endl;
    }
  }

  if (ParallelDescriptor::IOProcessor() && verbose >= 1) {
    cout << "Creating modelDump" << endl;
  }

  ofstream dumpfile;
  if (ParallelDescriptor::IOProcessor()) {
    dumpfile.open("modelDump", std::ios::out);
  }

  const int bufsiz = 200; // must hold 9 fields of width 16 plus a little more
  char buf[bufsiz];

  // The following are not parallel loops.
  // We are doing the Fab copy portion on all processors because
  // all must participate, but the write is done only on IOProcessor.

  int k = 1;
  list<Box>::iterator bi = bxlist.begin();
  list<int>::iterator li = lnlist.begin();
  for ( ; bi != bxlist.end(); ++bi, ++li) {
    const Box& reg = *bi;
    int ln         = *li;
    Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(ln));

    MultiFab& S_new = castro->get_new_data(State_Type);

    Fab stmp(reg, S_new.nComp());

    S_new.copy(stmp);

    if (ParallelDescriptor::IOProcessor()) {
      for (int i = reg.smallEnd(0); i <= reg.bigEnd(0); i++) {
        IntVect p(i);

        // Quantities written as 0.0 are not currently used, so we
        // don't bother to derive them.

        sprintf(buf,
                "%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E%16.8E %d\n",
                parent->Geom(ln).CellCenter(i, 0), // radius
                stmp(p, URHO),                      // density
                stmp(p, UMX) / stmp(p, URHO),      // velocity
                0.0,                               // pressure
                stmp(p, UTEMP),                    // temperature (MeV ?)
                stmp(p, UEINT) / stmp(p, URHO),      // internal energy e
                0.0,                               // entropy
                0.0,                               // cumulative mass
                stmp(p, UFX / stmp(p, URHO),       // Ye
                k++);
        std::string Buf = buf;
        dumpfile << Buf;
      }
    }
  }

  if (ParallelDescriptor::IOProcessor()) {
    dumpfile << Radiation::nGroups << '\n';
  }

  k = 1;
  bi = bxlist.begin();
  li = lnlist.begin();
  for ( ; bi != bxlist.end(); ++bi, ++li) {
    const Box& reg = *bi;
    int ln         = *li;
    Castro *castro = dynamic_cast<Castro*>(&parent->getLevel(ln));

    MultiFab& R_new = castro->get_new_data(Rad_Type);

    Fab rtmp(reg, R_new.nComp());

    R_new.copy(rtmp);

    if (ParallelDescriptor::IOProcessor()) {
      for (int i = reg.smallEnd(0); i <= reg.bigEnd(0); i++) {
        IntVect p(i);
        char* bufp = buf;
        for (int j = 0; ; j++) {
          sprintf(bufp, "%16.8E", rtmp(p, j));
          bufp += 16;
          if (j + 1 == Radiation::nGroups) {
            sprintf(bufp, " %d\n", k++);
            std::string Buf = buf;
            dumpfile << Buf;
            break;
          }
          else if (j % 4 == 3) {
            sprintf(bufp, "\n");
            std::string Buf = buf;
            dumpfile << Buf;
            bufp = buf;
          }
        }
      }
    }
  }

  if (ParallelDescriptor::IOProcessor()) {
    dumpfile.close();
  }
}

#endif

#endif
