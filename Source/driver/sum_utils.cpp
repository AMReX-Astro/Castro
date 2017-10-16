#include <iomanip>

#include <Castro.H>
#include <Castro_F.H>

#ifdef SELF_GRAVITY
#include <Gravity.H>
#include <Gravity_F.H>
#endif

using namespace amrex;

Real
Castro::sumDerive (const std::string& name,
                   Real               time,
		   bool               local)
{
    Real sum     = 0.0;
    auto mf = derive(name, time, 0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
    {
	for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
	{
	    sum += (*mf)[mfi].sum(mfi.tilebox(),0);
	}
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSum (const std::string& name,
                   Real               time,
		   bool               local,
		   bool               finemask)
{
    BL_PROFILE("Castro::volWgtSum()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto mf = derive(name,time,0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel() && finemask)
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

	Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	ca_summass(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
		   ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSquaredSum (const std::string& name,
                          Real               time,
			  bool               local)
{
    BL_PROFILE("Castro::volWgtSquaredSum()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto mf = derive(name,time,0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];
    
        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	ca_sumsquared(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
		      ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSum (const std::string& name,
                   Real               time,
                   int                idir,
		   bool               local)
{
    BL_PROFILE("Castro::locWgtSum()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto mf = derive(name,time,0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& fab = (*mf)[mfi];
    
        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	ca_sumlocmass(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
		      ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s,idir);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSum2D (const std::string& name,
                     Real               time,
                     int                idir1,
                     int                idir2,
		     bool               local)
{
    BL_PROFILE("Castro::locWgtSum2D()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto mf = derive(name,time,0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& fab = (*mf)[mfi];
    
        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	ca_sumlocmass2d(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
			ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s,idir1,idir2);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSumMF (const MultiFab& mf, int comp, bool local) 
{
    BL_PROFILE("Castro::volWgtSumMF()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& fab = mf[mfi];

        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //

	ca_summass(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_N_3D(fab,comp),
		   ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::volWgtSumOneSide (const std::string& name,
                          Real               time, 
                          int                side,
                          int                bdir,
			  bool               local)
{
    BL_PROFILE("Castro::volWgtSumOneSide()");

    // This function is a clone of volWgtSum except it computes the result only on half of the domain.
    // The lower half corresponds to side == 0 and the upper half corresponds to side == 1.
    // The argument bdir gives the direction along which to bisect.
    // Only designed to work in Cartesian geometries.

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto        mf      = derive(name,time,0);
    const int* domhi    = geom.Domain().hiVect();

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];
    
        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        int hiLeft[3]       = { *hi, *(hi+1), *(hi+2) };
        int loRight[3]      = { *lo, *(lo+1), *(lo+2) };

        hiLeft[bdir]        = *(domhi+bdir) / 2;
        loRight[bdir]       = *(domhi+bdir) / 2 + 1;

        const int* hiLeftPtr  = hiLeft;
        const int* loRightPtr = loRight;
        const int* loFinal;
        const int* hiFinal;

        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        //
        
        bool doSum = false;
        if ( side == 0 && *(lo + bdir) <= *(hiLeftPtr + bdir) ) {
          doSum = true;
          if ( *(hi + bdir) <= *(hiLeftPtr + bdir) ) {
            loFinal   = lo;
            hiFinal   = hi;
	  }
	  else {
            loFinal   = lo;
            hiFinal   = hiLeftPtr;
	  }
	}  
        else if ( side == 1 && *(hi + bdir) >= *(loRightPtr + bdir) ) {
          doSum = true;
          if ( *(lo + bdir) >= *(loRightPtr + bdir) ) {
            loFinal   = lo;
            hiFinal   = hi;
          }
          else {
            loFinal   = loRightPtr;
            hiFinal   = hi;
          }
	}

        if ( doSum ) {

          ca_summass(ARLIM_3D(loFinal),ARLIM_3D(hiFinal),BL_TO_FORTRAN_3D(fab),
		     ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);

        }
        
        sum += s;
		
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locWgtSumOneSide (const std::string& name,
                          Real               time,
                          int                idir, 
                          int                side,
                          int                bdir,
			  bool               local)
{
    BL_PROFILE("Castro::locWgtSumOneSide()");

  // This function is a clone of locWgtSum except that it only sums over one half of the domain.
  // The lower half corresponds to side == 0, and the upper half corresponds to side == 1.
  // The argument idir (x == 0, y == 1, z == 2) gives the direction to location weight by,
  // and the argument bdir gives the direction along which to bisect.

    Real sum            = 0.0;
    const Real* dx      = geom.CellSize();
    auto        mf      = derive(name,time,0); 
    const int* domhi    = geom.Domain().hiVect(); 

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif        
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        FArrayBox& fab = (*mf)[mfi];

        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

        int hiLeft[3]       = { *hi, *(hi+1), *(hi+2) };
        int loRight[3]      = { *lo, *(lo+1), *(lo+2) };

        hiLeft[bdir]        = *(domhi+bdir) / 2;
        loRight[bdir]       = (*(domhi+bdir) / 2) + 1;

        const int* hiLeftPtr  = hiLeft;
        const int* loRightPtr = loRight;
        const int* loFinal;
        const int* hiFinal;        
        //
        // Note that this routine will do a volume weighted sum of
        // whatever quantity is passed in, not strictly the "mass".
        // 

        bool doSum = false;
        if ( side == 0 && *(lo + bdir) <= *(hiLeftPtr + bdir) ) {
          doSum = true;
          if ( *(hi + bdir) <= *(hiLeftPtr + bdir) ) {
            loFinal   = box.loVect();
            hiFinal   = box.hiVect();
	  }
	  else {
            loFinal   = box.loVect();
            hiFinal   = hiLeftPtr;
	  }
	}  
        else if ( side == 1 && *(hi + bdir) >= *(loRightPtr + bdir) ) {
          doSum = true;
          if ( *(lo + bdir) >= *(loRightPtr + bdir) ) {
            loFinal   = box.loVect();
            hiFinal   = box.hiVect();
          }
          else {
            loFinal   = loRightPtr;
            hiFinal   = box.hiVect();
          }
	}

        if ( doSum ) {

          ca_sumlocmass(ARLIM_3D(loFinal),ARLIM_3D(hiFinal),BL_TO_FORTRAN_3D(fab),
			ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s,idir);

        }
     
        sum += s;
        
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;

}

Real
Castro::volProductSum (const std::string& name1, 
                       const std::string& name2,
                       Real time, bool local)
{
    BL_PROFILE("Castro::volProductSum()");

    Real        sum = 0.0;
    const Real* dx  = geom.CellSize();
    auto        mf1 = derive(name1,time,0);
    auto        mf2 = derive(name2,time,0);

    BL_ASSERT(mf1);
    BL_ASSERT(mf2);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf1, mask, 0, 0, 1, 0);
	MultiFab::Multiply(*mf2, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf1,true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& fab1 = (*mf1)[mfi];
        const FArrayBox& fab2 = (*mf2)[mfi];
    
        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	ca_sumproduct(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab1),
		      BL_TO_FORTRAN_3D(fab2),ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s);
        
        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

Real
Castro::locSquaredSum (const std::string& name,
                       Real               time,
                       int                idir,
		       bool               local)
{
    BL_PROFILE("Castro::locSquaredSum()");

    Real        sum     = 0.0;
    const Real* dx      = geom.CellSize();
    auto        mf      = derive(name,time,0);

    BL_ASSERT(mf);

    if (level < parent->finestLevel())
    {
	const MultiFab& mask = getLevel(level+1).build_fine_mask();
	MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
    }

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif    
    for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
    {
        const FArrayBox& fab = (*mf)[mfi];
    
        Real s = 0.0;
        const Box& box  = mfi.tilebox();
        const int* lo   = box.loVect();
        const int* hi   = box.hiVect();

	ca_sumlocsquaredmass(ARLIM_3D(lo),ARLIM_3D(hi),BL_TO_FORTRAN_3D(fab),
			     ZFILL(dx),BL_TO_FORTRAN_3D(volume[mfi]),&s,idir);

        sum += s;
    }

    if (!local)
	ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

