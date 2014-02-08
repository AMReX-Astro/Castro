
#include <EdgeVar.H>
#include <EDGEVAR_F.H>

#include <Using.H>

// Constructors destructor

EdgeVar::EdgeVar()
{
}

EdgeVar::EdgeVar(const BoxArray& boxes,
                           const int num_var,
                           const int tran_grow ,
                           const int norm_grow,
                           const Real initval,
                           const FabAlloc mem_mode)
{

  BL_ASSERT(num_var > 0);
  ba.define(boxes);
  ncomp = num_var;
  ngrow_norm = norm_grow;
  ngrow_tran = tran_grow;
  allocate = mem_mode ;
  Edge = 0 ;
  if (allocate) {
    Edge = new MultiFab[BL_SPACEDIM] ;

    for( int dim=0; dim<BL_SPACEDIM; dim++ ) {
      BoxArray bedge(ba);
      bedge.grow(dim, norm_grow);
      for( int tran=0; tran<BL_SPACEDIM; tran++ ) {
	if( dim != tran ) bedge.grow(tran, tran_grow);
      }
      bedge.surroundingNodes(dim);

      Edge[dim].define(bedge, num_var, 0, Fab_allocate );
      Edge[dim].setVal(initval) ;
    }
  }
}

EdgeVar::EdgeVar (const EdgeVar & ev )   // Copy constructor
  : ba(ev.ba),
    ncomp(ev.ncomp),
    ngrow_norm(ev.ngrow_norm),
    ngrow_tran(ev.ngrow_tran),
    allocate(ev.allocate),
    Edge(0)
{
  if (allocate) {
    Edge = new MultiFab[BL_SPACEDIM] ;

    for( int dim=0; dim<BL_SPACEDIM; dim++ ) {
      BoxArray bedge(ba);
      bedge.grow(dim, ngrow_norm);
      for( int tran=0; tran<BL_SPACEDIM; tran++ ) {
	if( dim != tran ) bedge.grow(tran, ngrow_tran);
      }
      bedge.surroundingNodes(dim);

      Edge[dim].define(bedge, ncomp, 0, Fab_allocate );
      for (MFIter ei(Edge[dim]); ei.isValid(); ++ei) {
	int i = ei.index();
	Edge[dim][i].copy(ev.Edge[dim][i]);
      }
    }
  }
}

EdgeVar::~EdgeVar() 
{
   if (allocate) {
      delete [] Edge ;
   }
}

#if 0

// Removed from pBoxlib version

// Insertation overload

ostream& operator<< (ostream& os, const EdgeVar& l)
{
#if BL_SPACEDIM == 2
   os << "-------------- EdgeVar -------------" << endl 
      << "Number of variables = " << l.ncomp  << endl
      << "Normal grow factor = " << l.ngrow_norm << endl
      << "Transverse grow factor = " << l.ngrow_tran << endl ;
      if (  (&(l.Edge)[0]) && (&(l.Edge)[1])  ) {
         os  << "X-normal MultiFab --- " << endl << (l.Edge)[0] 
             << "Y-normal MultiFab --- " << endl << (l.Edge)[1]  ;
      }
      else {
         os << "No MultiFab's Allocated" << endl ;
      }
      os << "-------------------------------" << endl ;
#endif
#if BL_SPACEDIM == 3
   os << "-------------- EdgeVar -------------" << endl 
      << "Number of variables = " << ncomp  << endl
      << "Normal grow factor = " << ngrow_norm << endl
      << "Transverse grow factor = " << ngrow_tran << endl ;
      if ((&(l.Edge[0]) && (&(l.Edge)[1]) && (&(l.Edge)[2]) ) {
         os << "X-normal MultiFab --- " << endl << (l.Edge)[0] 
            << "Y-normal MultiFab --- " << endl << (l.Edge)[1]
            << "Z-normal MultiFab --- " << endl << (l.Edge)[2] ;
      }
      else {
         os << "No MultiFab's Allocated" << endl ;
      }
      os << "-------------------------------" << endl ;
#endif
   return os ;
}
#endif

void EdgeVar::setVal(Real val)
{
   for (int i = 0 ; i < BL_SPACEDIM ; i++) {
     for (MFIter ei(Edge[i]); ei.isValid(); ++ei) {
       int j = ei.index();
       Edge[i][j].setVal(val) ;
     }
   }
}

void EdgeVar::Extensive(const Real* dx, Real factor)
{
   Real volume = 1.0 ;

   for( int n = 0 ; n < BL_SPACEDIM ; n++ ) {
      volume *= dx[n] ;
   }
   factor *= volume ;
   for( int n = 0 ; n < BL_SPACEDIM ; n++ ) {
      for (MFIter ei(Edge[n]); ei.isValid(); ++ei) {
	 int i = ei.index();
         Edge[n][i] *= (factor/dx[n]) ;
      }
   }
}

void EdgeVar::Intensive(const Real* dx, Real factor)
{
   Real volume = 1.0 ;

   for( int n = 0 ; n < BL_SPACEDIM ; n++ ) {
      volume *= dx[n] ;
   }
   for( int n = 0 ; n < BL_SPACEDIM ; n++ ) {
      for (MFIter ei(Edge[n]); ei.isValid(); ++ei) {
	 int i = ei.index();
         Edge[n][i] *= ((factor*dx[n])/volume) ;
      }
   }
}


EdgeVar& EdgeVar::operator*=(Real fact)
{
   int nghost = 0 ;
   for(int i = 0 ; i < BL_SPACEDIM ; i++) {
     for (MFIter ei(Edge[i]); ei.isValid(); ++ei) {
       int j = ei.index();
       Edge[i][j].mult(fact,nghost) ;
     }
   }
   return *this ;
}

EdgeVar& EdgeVar::operator*=(const EdgeVar& factor)
{
   BL_ASSERT( this->ba == factor.ba ) ;
   BL_ASSERT( this->ncomp == factor.ncomp ) ;

   for(int n = 0 ; n < BL_SPACEDIM ; n++) {
      for (MFIter ei(Edge[n]); ei.isValid(); ++ei) {
	 int i = ei.index();
         Edge[n][i] *= factor[n][i] ;
      }
   }
   return *this ;
}

EdgeVar& EdgeVar::operator/=(const EdgeVar& factor)
{
   BL_ASSERT( this->ba == factor.ba ) ;
   BL_ASSERT( this->ncomp == factor.ncomp ) ;

   for(int n = 0 ; n < BL_SPACEDIM ; n++) {
      for (MFIter ei(Edge[n]); ei.isValid(); ++ei) {
	 int i = ei.index();
         Edge[n][i] /= factor[n][i] ;
      }
   }
   return *this ;
}


EdgeVar& EdgeVar::operator+=(const EdgeVar& val)
{
   for(int i = 0 ; i < BL_SPACEDIM ; i++) {
     for (MFIter ei(Edge[i]); ei.isValid(); ++ei) {
       int j = ei.index();
       Edge[i][j].plus(val[i][j], 0, 0, ncomp) ;
     }
   }
   return *this ;
}

EdgeVar& EdgeVar::operator=(const EdgeVar& val)
{

   // Note : only works for EdgeVar's with same size

   BL_ASSERT( val.ba == ba ) ;
   BL_ASSERT( val.ncomp == ncomp ) ;

   int ngrid = val.ba.size() ;

   for(int n = 0 ; n < BL_SPACEDIM ; n++) {
      for (MFIter ei(Edge[n]); ei.isValid(); ++ei) {
         int j = ei.index();
         Edge[n][j].copy(val[n][j]) ;
      }
   }

   return *this ;
}


void EdgeVar::AvgDown(const EdgeVar & Fine, IntVect & nref)
{

   int nfine = Fine.ba.size() ;
   int ncoarse = (*this).ba.size() ;

   BoxArray crse_box( Fine.Boxes() ) ;
   crse_box.coarsen(nref) ;
   EdgeVar Coarse(crse_box,ncomp) ;

   cout << "EdgeVar::AvgDown not parallel yet" << endl;
   exit(1);

   for ( int i=0 ; i < nfine ; i++ ) {
#if (BL_SPACEDIM == 2)
      FORT_EDGE_AVG_DOWN(
         Fine[0][i].dataPtr(),   dimlist(Fine[0][i].box()),
         Fine[1][i].dataPtr(),   dimlist(Fine[1][i].box()),
         nref.getVect(),
         Coarse[0][i].dataPtr(), dimlist(Coarse[0][i].box()),
         Coarse[1][i].dataPtr(), dimlist(Coarse[1][i].box()) );
#elif (BL_SPACEDIM == 3)
      FORT_EDGE_AVG_DOWN(
         Fine[0][i].dataPtr(),   dimlist(Fine[0][i].box()),
         Fine[1][i].dataPtr(),   dimlist(Fine[1][i].box()),
         Fine[2][i].dataPtr(),   dimlist(Fine[2][i].box()),
         nref.getVect(),
         Coarse[0][i].dataPtr(), dimlist(Coarse[0][i].box()),
         Coarse[1][i].dataPtr(), dimlist(Coarse[1][i].box()),
         Coarse[2][i].dataPtr(), dimlist(Coarse[2][i].box()) );
#endif
      for (int j = 0 ; j < ncoarse ; j++ ) {
         for (int n = 0 ; n < BL_SPACEDIM ; n++ ) {
            ((*this)[n][j]).copy(Coarse[n][i]) ;
         }
      }
   }
}

void EdgeVar::Interp(const EdgeVar & Crse, IntVect nref)
{

   int ncoarse = Crse.ba.size() ;
   int nfine = (*this).ba.size() ;

   BoxArray fine_box( Crse.Boxes() ) ;
   fine_box.refine(nref) ;
   EdgeVar Fine(fine_box,ncomp) ;

   for (MFIter mfi(Crse[0]); mfi.isValid(); ++mfi) {
      int i = mfi.index();
#if (BL_SPACEDIM == 1)
      FORT_EDGE_INTERP(
         Crse[0][i].dataPtr(), dimlist(Crse[0][i].box()),
         nref.getVect(),
         Fine[0][i].dataPtr(), dimlist(Fine[0][i].box()));
#elif (BL_SPACEDIM == 2)
      FORT_EDGE_INTERP(
         Crse[0][i].dataPtr(), dimlist(Crse[0][i].box()),
         Crse[1][i].dataPtr(), dimlist(Crse[1][i].box()),
         nref.getVect(),
         Fine[0][i].dataPtr(), dimlist(Fine[0][i].box()),
         Fine[1][i].dataPtr(), dimlist(Fine[1][i].box()));
#elif (BL_SPACEDIM == 3)
      FORT_EDGE_INTERP(
         Crse[0][i].dataPtr(), dimlist(Crse[0][i].box()),
         Crse[1][i].dataPtr(), dimlist(Crse[1][i].box()),
         Crse[2][i].dataPtr(), dimlist(Crse[2][i].box()),
         nref.getVect(),
         Fine[0][i].dataPtr(), dimlist(Fine[0][i].box()),
         Fine[1][i].dataPtr(), dimlist(Fine[1][i].box()),
         Fine[2][i].dataPtr(), dimlist(Fine[2][i].box()));
#endif
   }

   for (int n = 0 ; n < BL_SPACEDIM ; n++ ) {
      (*this)[n].copy(Fine[n]);
   }
#if 0
   for ( int i=0 ; i < ncoarse ; i++ ) {
      FORT_EDGE_INTERP(
         Crse[0][i].dataPtr(), dimlist(Crse[0][i].box()),
         Crse[1][i].dataPtr(), dimlist(Crse[1][i].box()),
         nref.getVect(),
         Fine[0][i].dataPtr(), dimlist(Fine[0][i].box()),
         Fine[1][i].dataPtr(), dimlist(Fine[1][i].box()));

      for (int j = 0 ; j < nfine ; j++ ) {
         for (int n = 0 ; n < BL_SPACEDIM ; n++ ) {
            ((*this)[n][j]).copy(Fine[n][i]) ;
         }
      }
   }
#endif
}

void 
EdgeVar::SetMultiFabPtr(MultiFab* data) 
{
   for (int i = 0 ; i < BL_SPACEDIM ; i++) {
      BL_ASSERT( data[i].boxArray() == ba ) ;
   }
   Edge = data ;
}

void
EdgeVar::Gradient(const MultiFab & scalar, const Real* dx) 
{

   BL_ASSERT (ba == scalar.boxArray()) ;
   BL_ASSERT (ncomp == scalar.nComp()) ;

   for (MFIter ei(scalar); ei.isValid(); ++ei) {
      int i = ei.index();
#if (BL_SPACEDIM == 1)
      FORT_SET_GRADIENT( scalar[i].dataPtr(),  dimlist(scalar[i].box()),
                         dx,
                        &ncomp,
                         Edge[0][i].dataPtr(), dimlist(Edge[0][i].box()) ) ;
#elif (BL_SPACEDIM == 2)
      FORT_SET_GRADIENT( scalar[i].dataPtr(),  dimlist(scalar[i].box()),
                         dx,
                        &ncomp,
                         Edge[0][i].dataPtr(), dimlist(Edge[0][i].box()),
                         Edge[1][i].dataPtr(), dimlist(Edge[1][i].box()) ) ;
#elif (BL_SPACEDIM == 3)
      FORT_SET_GRADIENT( scalar[i].dataPtr(),  dimlist(scalar[i].box()),
                         dx,
                        &ncomp,
                         Edge[0][i].dataPtr(), dimlist(Edge[0][i].box()),
                         Edge[1][i].dataPtr(), dimlist(Edge[1][i].box()),
                         Edge[2][i].dataPtr(), dimlist(Edge[2][i].box()) ) ;
#endif
   }
}

void
EdgeVar::AvgDivergence(MultiFab & div, const Real* dx) const
{

   BL_ASSERT (ba == div.boxArray()) ;
   BL_ASSERT (ncomp == div.nComp()) ;

   for (MFIter ei(div); ei.isValid(); ++ei) {
      int i = ei.index();
#if (BL_SPACEDIM == 1)
      FORT_AVG_DIV( Edge[0][i].dataPtr(), dimlist(Edge[0][i].box()),
                    dx,
                   &ncomp,
                    div[i].dataPtr(),     dimlist(div[i].box()) ) ;
#elif (BL_SPACEDIM == 2)
      FORT_AVG_DIV( Edge[0][i].dataPtr(), dimlist(Edge[0][i].box()),
		    Edge[1][i].dataPtr(), dimlist(Edge[1][i].box()),
                    dx,
                   &ncomp,
                    div[i].dataPtr(),     dimlist(div[i].box()) ) ;
#elif (BL_SPACEDIM == 3)
      FORT_AVG_DIV( Edge[0][i].dataPtr(), dimlist(Edge[0][i].box()),
		    Edge[1][i].dataPtr(), dimlist(Edge[1][i].box()),
		    Edge[2][i].dataPtr(), dimlist(Edge[2][i].box()),
                    dx,
                   &ncomp,
                    div[i].dataPtr(),     dimlist(div[i].box()) ) ;
#endif
   }
}
