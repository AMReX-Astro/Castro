#ifndef DIFFUSION_H
#define DIFFUSION_H

#include <AMReX_AmrLevel.H>
#include <AMReX_MLLinOp.H>

#include <diffusion_params.H>

///
/// @class Diffusion
/// @brief
///
class Diffusion {

public:


///
/// @param Parent
/// @param _phys_bc
///
  Diffusion(amrex::Amr* Parent, amrex::BCRec* _phys_bc);
  ~Diffusion () = default;

///
/// Remove copy/move constructors/assignment operators.
///
  Diffusion (const Diffusion&) = delete;
  Diffusion (Diffusion&&) = delete;
  Diffusion& operator= (const Diffusion&) = delete;
  Diffusion& operator= (Diffusion&&) = delete;

  void read_params ();


///
/// @param jobInfoFile
///
  static void output_job_info_params(std::ostream& jobInfoFile);


///
/// @param level
/// @param level_data
/// @param volume
/// @param area
///
  void install_level (int                   level,
                      amrex::AmrLevel*             level_data,
                      amrex::MultiFab&             volume,
                      amrex::MultiFab*             area);


///
/// @param level
/// @param Temperature
/// @param CrseTemp
/// @param DiffTerm
/// @param temp_cond_coef
///
  void applyop(int level,amrex::MultiFab& Temperature,amrex::MultiFab& CrseTemp,
               amrex::MultiFab& DiffTerm, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& temp_cond_coef);

  void make_mg_bc();

protected:

///
/// Pointers to amr,amrlevel.
///
  amrex::Amr*             parent;
  amrex::Vector<amrex::AmrLevel*> LevelData;


///
/// BoxArray at each level
///
  amrex::Vector<amrex::BoxArray> grids;

///
/// Volume and area fractions.
///
  amrex::Vector<amrex::MultiFab*> volume;
  amrex::Vector<amrex::MultiFab*> area;

  amrex::BCRec*       phys_bc;

  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_lobc;
  std::array<amrex::MLLinOp::BCType,AMREX_SPACEDIM> mlmg_hibc;


///
/// @param level
/// @param Temperature
/// @param CrseTemp
/// @param DiffTerm
/// @param temp_cond_coef
///
  void applyop_mlmg(int level,amrex::MultiFab& Temperature,amrex::MultiFab& CrseTemp,
                    amrex::MultiFab& DiffTerm, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& temp_cond_coef);

};
#endif
