! Given a dimension agnostic array, call the relevant copy of filcc
! depending on the dimensionality of the problem. Note that we rely on
! the assumption that the unused indices equal zero for higher dimensions.

subroutine filcc_nd(adv,adv_lo,adv_hi,domlo,domhi,delta,xlo,bc)

  use prob_params_module, only: dim
  
  use amrex_fort_module, only : rt => amrex_real
  implicit none
  
  include 'AMReX_bc_types.fi'  
  
  integer          :: adv_lo(3),adv_hi(3)
  integer          :: bc(dim,2)
  integer          :: domlo(3), domhi(3)
  real(rt)         :: delta(3), xlo(3)
  real(rt)         :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))

#if (BL_SPACEDIM == 1)

     call filcc(adv(:,0,0),adv_lo(1),adv_hi(1),domlo(1:dim),domhi(1:dim),delta(1:dim),xlo(1:dim),bc)

#elif (BL_SPACEDIM == 2)

     call filcc(adv(:,:,0),adv_lo(1),adv_lo(2),adv_hi(1),adv_hi(2),domlo(1:dim),domhi(1:dim),delta(1:dim),xlo(1:dim),bc)

#else

     call filcc(adv(:,:,:),adv_lo(1),adv_lo(2),adv_lo(3),adv_hi(1),adv_hi(2),adv_hi(3),domlo,domhi,delta,xlo,bc)

#endif


end subroutine filcc_nd
