#ifndef _EDGE_VAR_H_
#define _EDGE_VAR_H_

#include <iostream>

#include <BoxArray.H>
#include <FabArray.H>
#include <MultiFab.H>

//
// EdgeVar - 
//	Contains edge centered MultiFabs in all
//	coordinate directions based on a cell centered BoxArray.
//
//      Growth in normal and transverse directions allowed (default 0)
//
//      Initialization to a value allowed (default 0.0) 
//

class EdgeVar
{

public:

// constructors , destructor
  EdgeVar() ;
  EdgeVar(
     const BoxArray& boxes,   // Input Box Array around which fluxes
                              // are to be defined.
     const int num_var,       // number of components
     const int tran_grow =  0,  // Amount to grow in transverse
                                // direction.
     const int norm_grow =  0,  // Amount to grow in normal
                                // direction.
     const Real initval = 0.0,  // Value to initialize MultiFabs with.
     FabAlloc mem_mode = Fab_allocate // Allocate by default
   );
  EdgeVar (const EdgeVar& );   // Copy constructor
  ~EdgeVar();

   void Extensive(const Real* dx, Real factor = 1.0 ) ;
   void Intensive(const Real* dx, Real factor = 1.0 ) ;
   void AvgDown(const EdgeVar & Fine, IntVect & nref) ;
   void Interp(const EdgeVar & Crse, IntVect nref) ;
   void Gradient(const MultiFab & data, const Real* dx) ;
   void AvgDivergence( MultiFab & div, const Real* dx) const ;
   void setVal(Real value) ;

// operators
EdgeVar& operator*=(Real fact) ;
EdgeVar& operator*=(const EdgeVar& factor) ;
EdgeVar& operator/=(const EdgeVar& factor) ;
EdgeVar& operator+=(const EdgeVar& val) ;
EdgeVar& operator=(const EdgeVar& val) ;

// access functions

    inline MultiFab& operator[](int dir)  {return Edge[dir];};
    inline const MultiFab& operator[](int dir)  const {return Edge[dir];};
    inline const BoxArray& Boxes() const {return ba;} ;
    inline BoxArray& Boxes()  {return ba;} ;
    inline int num_components() {return ncomp;} ;
    inline int normal_growfactor() {return ngrow_norm;} ;
    inline int transverse_growfactor() {return ngrow_tran;} ;
    inline MultiFab* FabPtr() {return Edge;} ;
    void SetMultiFabPtr( MultiFab* data) ;

// i/o

    //friend ostream& operator << (ostream&, const EdgeVar&) ;

protected :
    BoxArray ba;
    int ncomp;
    int ngrow_norm ;
    int ngrow_tran ;
    FabAlloc allocate ;
    MultiFab* Edge;

};

#endif
