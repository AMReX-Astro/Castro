subroutine ca_set_godunov_indices( &
                                )

  use meth_params_module

  NGDNV = 1

  GDRHO = NGDNV
  NGDNV = NGDNV + 1

  GDU = NGDNV
  NGDNV = NGDNV + 1

  GDV = NGDNV
  NGDNV = NGDNV + 1

  GDW = NGDNV
  NGDNV = NGDNV + 1

  GDPRES = NGDNV
  NGDNV = NGDNV + 1

  GDGAME = NGDNV
  NGDNV = NGDNV + 1

#ifdef RADIATION
  GDLAMS = NGDNV
  NGDNV = NGDNV + ngroups

  GDERADS = NGDNV
  NGDNV = NGDNV + ngroups

#endif
  NGDNV = NGDNV - 1
end subroutine ca_set_godunov_indices

subroutine ca_set_conserved_indices( &
                                   Eint, &
                                   Temp, &
                                   Zmom, &
                                   Density, &
                                   Eden, &
                                   FirstAdv, &
                                   Shock, &
                                   FirstSpec, &
                                   Ymom, &
                                   FirstAux, &
                                   Xmom, &
                                  )

  use meth_params_module
  integer, intent(in) :: Eint
  integer, intent(in) :: Temp
  integer, intent(in) :: Zmom
  integer, intent(in) :: Density
  integer, intent(in) :: Eden
  integer, intent(in) :: FirstAdv
  integer, intent(in) :: Shock
  integer, intent(in) :: FirstSpec
  integer, intent(in) :: Ymom
  integer, intent(in) :: FirstAux
  integer, intent(in) :: Xmom

  NVAR = 1

  URHO = NVAR
  NVAR = NVAR + 1
  call check_equal(URHO,Density+1)

  UMX = NVAR
  NVAR = NVAR + 1
  call check_equal(UMX,Xmom+1)

  UMY = NVAR
  NVAR = NVAR + 1
  call check_equal(UMY,Ymom+1)

  UMZ = NVAR
  NVAR = NVAR + 1
  call check_equal(UMZ,Zmom+1)

  UEDEN = NVAR
  NVAR = NVAR + 1
  call check_equal(UEDEN,Eden+1)

  UEINT = NVAR
  NVAR = NVAR + 1
  call check_equal(UEINT,Eint+1)

  UTEMP = NVAR
  NVAR = NVAR + 1
  call check_equal(UTEMP,Temp+1)

  UFA = NVAR
  NVAR = NVAR + numadv
  call check_equal(UFA,FirstAdv+1)

  UFS = NVAR
  NVAR = NVAR + nspec
  call check_equal(UFS,FirstSpec+1)

  UFX = NVAR
  NVAR = NVAR + naux
  call check_equal(UFX,FirstAux+1)

  USHK = NVAR
  NVAR = NVAR + 1
  call check_equal(USHK,Shock+1)

#ifdef HYBRID_MOMENTUM
  UMR = NVAR
  NVAR = NVAR + 1

  UML = NVAR
  NVAR = NVAR + 1

  UMP = NVAR
  NVAR = NVAR + 1

#endif
  NVAR = NVAR - 1
end subroutine ca_set_conserved_indices

subroutine ca_set_auxillary_indices( &
                                  )

  use meth_params_module

  NQAUX = 1

  QGAMC = NQAUX
  NQAUX = NQAUX + 1

  QC = NQAUX
  NQAUX = NQAUX + 1

  QDPDR = NQAUX
  NQAUX = NQAUX + 1

  QDPDE = NQAUX
  NQAUX = NQAUX + 1

#ifdef RADIATION
  QGAMCC = NQAUX
  NQAUX = NQAUX + 1

  QCG = NQAUX
  NQAUX = NQAUX + 1

  QLAMS = NQAUX
  NQAUX = NQAUX + ngroups

#endif
  NQAUX = NQAUX - 1
end subroutine ca_set_auxillary_indices

subroutine ca_set_primitive_indices( &
                                  )

  use meth_params_module

  NQ = 1

  QVAR = 1

  QRHO = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QU = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QV = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QW = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QGAME = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QPRES = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QREINT = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QTEMP = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QFA = NQ
  NQ = NQ + numadv
  QVAR = QVAR + numadv

  QFS = NQ
  NQ = NQ + nspec
  QVAR = QVAR + nspec

  QFX = NQ
  NQ = NQ + naux
  QVAR = QVAR + naux

#ifdef RADIATION
  QPTOT = NQ
  NQ = NQ + 1

  QREITOT = NQ
  NQ = NQ + 1

  QRAD = NQ
  NQ = NQ + ngroups

#endif
#ifdef MHD
  QMAGX = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QMAGY = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

  QMAGZ = NQ
  NQ = NQ + 1
  QVAR = QVAR + 1

#endif
  NQ = NQ - 1
  QVAR = QVAR - 1
end subroutine ca_set_primitive_indices

