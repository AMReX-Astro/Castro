
subroutine filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

  use prob_params_module, only: dim
  use amrex_fort_module, only: rt => amrex_real
  use amrex_filcc_module, only: filccn

  implicit none
  
  integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
  integer,  intent(in   ) :: bc(dim,2)
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: delta(3), xlo(3)
  real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

  call filccn(adv_lo, adv_hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

end subroutine filcc_nd
