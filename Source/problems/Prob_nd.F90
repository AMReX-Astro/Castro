subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer          :: init, namlen
  integer          :: name(namlen)
  real(rt)         :: problo(3), probhi(3)

end subroutine amrex_probinit
