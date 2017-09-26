
subroutine ca_advance_mhd() bind(C, name="ca_advance_mhd")

  use meth_params_module, only : MAXADV, NMAG, QMAGX

  print *, "MAXADV ", MAXADV
  print *, "NMAG ", NMAG 

end subroutine
