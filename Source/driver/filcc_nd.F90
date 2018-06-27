AMREX_LAUNCH subroutine filcc_nd_wrapper(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)
  use prob_params_module, only: dim
  use amrex_fort_module, only: rt => amrex_real, get_loop_bounds
  use amrex_filcc_module, only: amrex_filccn

  implicit none

  integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
  integer,  intent(in   ) :: bc(dim,2)
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: delta(3), xlo(3)
  real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

  integer :: blo(3), bhi(3)

  call get_loop_bounds(blo, bhi, adv_lo, adv_hi)

  call amrex_filccn(blo, bhi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

end subroutine filcc_nd_wrapper


subroutine filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

  use prob_params_module, only: dim
  use amrex_fort_module, only: rt => amrex_real
#ifdef CUDA
  use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

  implicit none
  
  integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
  integer,  intent(in   ) :: bc(dim,2)
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(in   ) :: delta(3), xlo(3)
  real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

#ifdef CUDA
  attributes(device) :: adv, adv_lo, adv_hi, bc, domlo, domhi, delta, xlo
#endif

  call filcc_nd_wrapper &
#ifdef CUDA
       <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
       (adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, bc)

end subroutine filcc_nd
