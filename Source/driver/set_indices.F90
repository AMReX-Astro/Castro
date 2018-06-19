
! DO NOT EDIT!!!

! This file is automatically created by set_variables.py.  To update
! or add variable indices, please edit _variables and then rerun the
! script.


subroutine check_equal(index1, index2)

  use amrex_error_module

  implicit none

  integer, intent(in) :: index1, index2

#ifndef AMREX_USE_CUDA
  if (index1 /= index2) then
    call amrex_error("ERROR: mismatch of indices")
  endif
#endif

end subroutine check_equal


subroutine ca_set_auxiliary_indices()


  use meth_params_module
  use network, only: naux, nspec
#ifdef RADIATION
  use rad_params_module, only : ngroups
#endif
  implicit none

  QGAMC = 1

  QC = 2

  QDPDR = 3

  QDPDE = 4

#ifdef RADIATION
  QGAMCG = 5
#endif

#ifdef RADIATION
  QCG = 6
#endif

#ifdef RADIATION
  QLAMS = 7
#endif

end subroutine ca_set_auxiliary_indices

subroutine ca_set_conserved_indices( &
#ifdef HYBRID_MOMENTUM
                                    Rmom, &
#endif
#ifdef HYBRID_MOMENTUM
                                    Lmom, &
#endif
#ifdef HYBRID_MOMENTUM
                                    Pmom, &
#endif
#ifdef SHOCK_VAR
                                    Shock, &
#endif
                                    Density, &
                                    Xmom, &
                                    Ymom, &
                                    Zmom, &
                                    Eden, &
                                    Eint, &
                                    Temp, &
                                    FirstAdv, &
                                    FirstSpec, &
                                    FirstAux &
                                   )


  use meth_params_module
  use network, only: naux, nspec
#ifdef RADIATION
  use rad_params_module, only : ngroups
#endif
  implicit none
  integer, intent(in) :: Density
  integer, intent(in) :: Xmom
  integer, intent(in) :: Ymom
  integer, intent(in) :: Zmom
#ifdef HYBRID_MOMENTUM
  integer, intent(in) :: Rmom
#endif
#ifdef HYBRID_MOMENTUM
  integer, intent(in) :: Lmom
#endif
#ifdef HYBRID_MOMENTUM
  integer, intent(in) :: Pmom
#endif
  integer, intent(in) :: Eden
  integer, intent(in) :: Eint
  integer, intent(in) :: Temp
  integer, intent(in) :: FirstAdv
  integer, intent(in) :: FirstSpec
  integer, intent(in) :: FirstAux
#ifdef SHOCK_VAR
  integer, intent(in) :: Shock
#endif

  URHO = 1
  call check_equal(URHO,Density+1)

  UMX = 2
  call check_equal(UMX,Xmom+1)

  UMY = 3
  call check_equal(UMY,Ymom+1)

  UMZ = 4
  call check_equal(UMZ,Zmom+1)

#ifdef HYBRID_MOMENTUM
  UMR = 5
  call check_equal(UMR,Rmom+1)
#endif

#ifdef HYBRID_MOMENTUM
  UML = 6
  call check_equal(UML,Lmom+1)
#endif

#ifdef HYBRID_MOMENTUM
  UMP = 7
  call check_equal(UMP,Pmom+1)
#endif

  UEDEN = 8
  call check_equal(UEDEN,Eden+1)

  UEINT = 9
  call check_equal(UEINT,Eint+1)

  UTEMP = 10
  call check_equal(UTEMP,Temp+1)

  if (nadv > 0) then
    UFA = 11
  else
    UFA = 0
  endif
  call check_equal(UFA,FirstAdv+1)

  if (nspec > 0) then
    UFS = 11 + nadv
  else
    UFS = 0
  endif
  call check_equal(UFS,FirstSpec+1)

  if (naux > 0) then
    UFX = 11 + nadv + nspec
  else
    UFX = 0
  endif
  call check_equal(UFX,FirstAux+1)

#ifdef SHOCK_VAR
  USHK = 11 + nadv + nspec + naux
  call check_equal(USHK,Shock+1)
#endif

end subroutine ca_set_conserved_indices

subroutine ca_set_godunov_indices()


  use meth_params_module
  use network, only: naux, nspec
#ifdef RADIATION
  use rad_params_module, only : ngroups
#endif
  implicit none

  GDRHO = 1

  GDU = 2

  GDV = 3

  GDW = 4

  GDPRES = 5

  GDGAME = 6

#ifdef RADIATION
  GDLAMS = 7
#endif

#ifdef RADIATION
  GDERADS = 7 + ngroups
#endif

end subroutine ca_set_godunov_indices

subroutine ca_set_primitive_indices()


  use meth_params_module
  use network, only: naux, nspec
#ifdef RADIATION
  use rad_params_module, only : ngroups
#endif
  implicit none

  QRHO = 1

  QU = 2

  QV = 3

  QW = 4

  QGAME = 5

  QPRES = 6

  QREINT = 7

#ifdef MHD
  QMAGX = 8
#endif

#ifdef MHD
  QMAGY = 9
#endif

#ifdef MHD
  QMAGZ = 10
#endif

  QTEMP = 11

  QFA = 12

  QFS = 12 + nadv

  QFX = 12 + nadv + nspec

#ifdef RADIATION
  QPTOT = 12 + nadv + nspec + naux
#endif

#ifdef RADIATION
  QREITOT = 13 + nadv + nspec + naux
#endif

#ifdef RADIATION
  QRAD = 14 + nadv + nspec + naux
#endif

end subroutine ca_set_primitive_indices

