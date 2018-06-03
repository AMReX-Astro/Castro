#ifndef _Diffusion_H_
#define _Diffusion_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_MLLinOp.H>

class Diffusion {

public:

  Diffusion(amrex::Amr* Parent, amrex::BCRec* _phys_bc);
  ~Diffusion ();
  void read_params ();

  void output_job_info_params(std::ostream& jobInfoFile);

  void install_level (int                   level,
                      amrex::AmrLevel*             level_data,
                      amrex::MultiFab&             volume,
                      amrex::MultiFab*             area);

  void applyop(int level,amrex::MultiFab& Temperature,amrex::MultiFab& CrseTemp,
               amrex::MultiFab& DiffTerm, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& temp_cond_coef);

  void applyViscOp(int level,amrex::MultiFab& Vel, amrex::MultiFab& CrseVel,
                   amrex::MultiFab& ViscTerm, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& visc_coeff);

  void make_mg_bc();

protected:
  //
  // Pointers to amr,amrlevel.
  //
  amrex::Amr*             parent;
  amrex::Vector<amrex::AmrLevel*> LevelData;

  //
  // BoxArray at each level
  //
  amrex::Vector<amrex::BoxArray> grids;
  //
  // Volume and area fractions.
  //
  amrex::Vector<amrex::MultiFab*> volume;
  amrex::Vector<amrex::MultiFab*> area;

  int mg_bc[2*BL_SPACEDIM];

  amrex::BCRec*       phys_bc;

  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;

#include "diffusion_params.H"

  static int   stencil_type;

#if (BL_SPACEDIM < 3)
  void applyMetricTerms(int level,amrex::MultiFab& Rhs, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& coeffs);
  void   weight_cc(int level,amrex::MultiFab& cc);
  void unweight_cc(int level,amrex::MultiFab& cc);
#endif

  void applyop_mlmg(int level,amrex::MultiFab& Temperature,amrex::MultiFab& CrseTemp,
                    amrex::MultiFab& DiffTerm, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& temp_cond_coef);

  void applyViscOp_mlmg(int level,amrex::MultiFab& Vel, amrex::MultiFab& CrseVel,
                        amrex::MultiFab& ViscTerm, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& visc_coeff);

};
#endif

