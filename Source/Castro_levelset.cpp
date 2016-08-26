#include "Castro.H"
#include "Castro_F.H"

void
Castro::advance_levelset(Real time, Real dt)
{
    BL_PROFILE("Castro::advance_levelset()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... update levelset\n";

    MultiFab&  LS_new = get_new_data(LS_State_Type);
    LS_new.copy(get_old_data(LS_State_Type));

    for (MFIter mfi(LStype); mfi.isValid(); ++mfi)
    {
        IntFab& type = LStype[mfi];
        IntFab& nband = LSnband[mfi];
        IntFab& mine = LSmine[mfi];
        const Box& box = mfi.validbox();
        int nbandsize = nband.box().numPts();
        int minesize = mine.box().numPts();
        
        // Set nband data based on type
	BL_FORT_PROC_CALL(LS_NARROWBAND,ls_narrowband)
	     (BL_TO_FORTRAN(type),
              nband.dataPtr(), &nbandsize,
              mine.dataPtr(), &minesize,
              box.loVect(), box.hiVect());
    }
    
    Real phidt = dt;
    Real phit = dt;
    
    while(phit > 0) {
        
        int nGrowLS = 1;
        int nCompLS = 1;
        const Real* dx = geom.CellSize();
        
        phidt = phit;
        
	LStype.FillBoundary(0, 1, geom.periodicity());
        
        for (FillPatchIterator fpi(*this,LS_new,nGrowLS,state[LS_State_Type].curTime(),LS_State_Type,0,nCompLS);
             fpi.isValid(); ++fpi)
        {
            
            const FArrayBox& lsfab = fpi();  // Note that this is a temporary fab with nCompLS components filled by FPI
            const Box& box = fpi.validbox();
            IntFab& type = LStype[fpi];
            
            IntFab& nband = LSnband[fpi];
            IntFab& mine = LSmine[fpi];
            int nbandsize = nband.box().numPts();
            int minesize = mine.box().numPts();
            
            const FArrayBox& uadv = u_gdnv[0][fpi];
            const FArrayBox& vadv = u_gdnv[1][fpi];
#if (BL_SPACEDIM == 3)
            const FArrayBox& wadv = u_gdnv[2][fpi];
#endif

            Real lscfl;
            BL_FORT_PROC_CALL(LS_CFL,ls_cfl)
                  (&lscfl, 
                   BL_TO_FORTRAN(lsfab), BL_TO_FORTRAN(uadv),  BL_TO_FORTRAN(vadv),
#if (BL_SPACEDIM == 3)
                   BL_TO_FORTRAN(wadv), 
#endif
                   nband.dataPtr(), &nbandsize, mine.dataPtr(), &minesize,
                   box.loVect(), box.hiVect(), &phit, dx,
                   BL_TO_FORTRAN(type));
            phidt = std::min(phidt,lscfl);
        }
        ParallelDescriptor::ReduceRealMin(phidt);
        
        phit = phit - phidt;
	if (verbose && ParallelDescriptor::IOProcessor())
	  {
	    std::cout<<"phit"<<phit<<"\n";
	    std::cout<<"phidt"<<phidt<<"\n";
	  }
        bool reinit = false;
        
        for (FillPatchIterator fpi(*this,LS_new,nGrowLS,state[LS_State_Type].curTime(),LS_State_Type,0,nCompLS);
             fpi.isValid(); ++fpi)
        {
            const FArrayBox& lsfab = fpi();  // Note that this is a temporary fab with nCompLS components filled by FPI
            FArrayBox& phinew = LS_new[fpi];  // This is a fab of the actual state, with phi in the 0 component
            const Box& box = fpi.validbox();
            const FArrayBox& uadv = u_gdnv[0][fpi];
            const FArrayBox& vadv = u_gdnv[1][fpi];
#if (BL_SPACEDIM == 3)
            const FArrayBox& wadv = u_gdnv[2][fpi];
#endif

            IntFab& type = LStype[fpi];
            
            IntFab& nband = LSnband[fpi];
            IntFab& mine = LSmine[fpi];
            int nbandsize = nband.box().numPts();
            int minesize = mine.box().numPts();
            
            // Advance level set for this fab, result directly into state
            int reinit_flag;
	    BL_FORT_PROC_CALL(LS_PHIUPD,ls_phiupd)
                  (&reinit_flag,BL_TO_FORTRAN(lsfab), BL_TO_FORTRAN(phinew), 
                   BL_TO_FORTRAN(uadv), BL_TO_FORTRAN(vadv), 
#if (BL_SPACEDIM == 3)
                   BL_TO_FORTRAN(wadv), 
#endif
                   nband.dataPtr(), &nbandsize, mine.dataPtr(), &minesize,
                   box.loVect(), box.hiVect(), &phidt, dx,
                   BL_TO_FORTRAN(type)); 
            reinit |= reinit_flag;
        }
        
        ParallelDescriptor::ReduceBoolOr(reinit);
        
        if (reinit)
        {
            reinit_phi(time);
        }
    }
    
    delete [] u_gdnv;

}

void
Castro::reinit_phi(Real time)
{
    BL_PROFILE("Castro::reinit_phi()");

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...... reinitialzing levelset\n";

    MultiFab& LS_new = get_new_data(LS_State_Type);
    const Real* dx   = geom.CellSize();
    int nGrowLS = 2;
    int nCompLS = 1;

    LStype.FillBoundary(0, 1, geom.periodicity());
        
    // Load valid region of phi
    Array<int> intfacep, intfacen, heap, heaploc;
    for (FillPatchIterator fpi(*this,LS_new,nGrowLS,state[LS_State_Type].curTime(),LS_State_Type,0,nCompLS);
         fpi.isValid(); ++fpi)
    {
        const FArrayBox& lsfab = fpi();  // Note that this is a temporary fab with nCompLS components filled by FPI
        FArrayBox& phinew = LS_new[fpi];  // This is a fab of the actual state, with phi in the 0 component
        IntFab& type = LStype[fpi];
        const Box& box = fpi.validbox();
        
        int intfacesize = box.numPts();
        int intfacenump = 0;
        int intfacenumn = 0;
        intfacep.resize(BL_SPACEDIM*intfacesize);
        intfacen.resize(BL_SPACEDIM*intfacesize);
        
        IntFab& nband = LSnband[fpi];
        IntFab& mine = LSmine[fpi];
        int nbandsize = nband.box().numPts();

        int nbandnum = 0;
        heap.resize(BL_SPACEDIM*nbandsize);
        int  typesize = type.box().numPts();
        heaploc.resize(typesize);

        // Set type for all cells in band to "outside band"
	BL_FORT_PROC_CALL(LS_RETYPIFY,ls_retypify)
            (BL_TO_FORTRAN(type), nband.dataPtr(), &nbandsize);

        // Set list of positive and negative narrowband points to be filled by FASTMARCH
	BL_FORT_PROC_CALL(LS_FINDINTERFACE,ls_findinterface)
            (BL_TO_FORTRAN(lsfab),  BL_TO_FORTRAN(phinew),  
             BL_TO_FORTRAN(type),  
             box.loVect(), box.hiVect(), dx, &intfacenump, &intfacenumn,
             intfacep.dataPtr(), intfacen.dataPtr(),
             nband.dataPtr(), &nbandsize, &intfacesize);

        // Fill in narrow band (+ve/-ve sides), set mines
        if(intfacenump > 0)
        {
            int positive = 1;
	    BL_FORT_PROC_CALL(LS_FASTMARCH,ls_fastmarch)
                 (BL_TO_FORTRAN(phinew),  BL_TO_FORTRAN(type),  
                  box.loVect(), box.hiVect(), dx, &intfacenump, intfacep.dataPtr(),
                  nband.dataPtr(), &nbandsize, &nbandnum,
                  mine.dataPtr(), &positive,&intfacesize,
                  heap.dataPtr(), heaploc.dataPtr());
        }
        
        if(intfacenumn > 0)
        {
            int negative = -1;
	    BL_FORT_PROC_CALL(LS_FASTMARCH,ls_fastmarch)
                 (BL_TO_FORTRAN(phinew),  BL_TO_FORTRAN(type),  
                  box.loVect(), box.hiVect(), dx,  &intfacenumn,  intfacen.dataPtr(),
            	  nband.dataPtr(), &nbandsize, &nbandnum,
            	  mine.dataPtr(), &negative,&intfacesize,
            	  heap.dataPtr(), heaploc.dataPtr());
        }
    }

    bool notdone = true;
    bool redo = false;
    
    // Check grow region and see if anything changes due to neighboring grids
    while (notdone)
    {
	LS_new.FillBoundary(geom.periodicity());
	LStype.FillBoundary(0, 1, geom.periodicity());
        
        for (MFIter mfi(LS_new); mfi.isValid(); ++mfi)
        {
            FArrayBox& phinew = LS_new[mfi];  // This is a fab of the actual state, with phi in the 0 component
            IntFab& type = LStype[mfi];
            const Box& box = mfi.validbox();
            
            IntFab& nband = LSnband[mfi];
            int nbandsize = nband.box().numPts();
            
            int nbandnum = 0;
            
            int  numtemptype = type.box().numPts();            
            heaploc.resize(numtemptype);
            
	    BL_FORT_PROC_CALL(LS_NBANDNUMIFY,ls_nbandnumify)
               (nband.dataPtr(), &nbandsize,&nbandnum);
            
            int positive = 1;
            int redo_flag;
	    BL_FORT_PROC_CALL(LS_FASTMARCH2,ls_fastmarch2)
                 (&redo_flag, BL_TO_FORTRAN(phinew),BL_TO_FORTRAN(type),
                  box.loVect(), box.hiVect(), dx,
                  nband.dataPtr(), &nbandsize, &nbandnum, &positive, heaploc.dataPtr());
            redo |= redo_flag;
            
            int negative = -1;
	    BL_FORT_PROC_CALL(LS_FASTMARCH2,ls_fastmarch2)
                 (&redo_flag, BL_TO_FORTRAN(phinew),BL_TO_FORTRAN(type),
                  box.loVect(), box.hiVect(), dx,
                  nband.dataPtr(), &nbandsize, &nbandnum, &negative, heaploc.dataPtr());
            redo |= redo_flag;
        }
        ParallelDescriptor::ReduceBoolOr(redo);
        notdone = redo;
        redo = false;
    }

    for (MFIter mfi(LStype); mfi.isValid(); ++mfi)
    {
        IntFab& type = LStype[mfi];
        IntFab& nband = LSnband[mfi];
        IntFab& mine = LSmine[mfi];
        const Box& box = mfi.validbox();
        
        int nbandsize = nband.box().numPts();
        int minesize = mine.box().numPts();
        
	BL_FORT_PROC_CALL(LS_MINE,ls_mine)
             (BL_TO_FORTRAN(type),
              nband.dataPtr(), &nbandsize, mine.dataPtr(), &minesize,
              box.loVect(), box.hiVect());
    }
}
