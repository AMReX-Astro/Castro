#ifndef _Diffusion_H_
#define _Diffusion_H_

#include <AmrLevel.H>

class Diffusion {

public:

  Diffusion(Amr* Parent, BCRec* _phys_bc);
  ~Diffusion ();
  void read_params ();

  void install_level (int                   level,
                      AmrLevel*             level_data,
                      MultiFab&             volume,
                      MultiFab*             area);

  void applyop(int level,MultiFab& Temperature,MultiFab& CrseTemp,
               MultiFab& DiffTerm, Array<std::unique_ptr<MultiFab> >& temp_cond_coef);

  void applyViscOp(int level,MultiFab& Vel, MultiFab& CrseVel,
                   MultiFab& ViscTerm, Array<std::unique_ptr<MultiFab> >& visc_coeff);

  void make_mg_bc();

  void GetCrsePhi(int level, 
                  MultiFab& phi_crse,
                  Real time);

  void zeroPhiFluxReg (int level);


protected:
  //
  // Pointers to amr,amrlevel.
  //
  Amr*             parent;
  Array<AmrLevel*> LevelData;

  //
  // BoxArray at each level
  //
  Array<BoxArray> grids;
  //
  // Volume and area fractions.
  //
  Array<MultiFab*> volume;
  Array<MultiFab*> area;

  int mg_bc[2*BL_SPACEDIM];

  BCRec*       phys_bc;

#include "diffusion_params.H"

  static int   stencil_type;

#if (BL_SPACEDIM < 3)
  void applyMetricTerms(int level,MultiFab& Rhs, Array<std::unique_ptr<MultiFab> >& coeffs);
  void   weight_cc(int level,MultiFab& cc);
  void unweight_cc(int level,MultiFab& cc);
#endif
};
#endif

