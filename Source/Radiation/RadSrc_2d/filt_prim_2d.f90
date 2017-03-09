subroutine ca_filt_prim(lo, hi, &
     Stmp, Stmp_l1, Stmp_l2, Stmp_h1, Stmp_h2, &
     Snew, Snew_l1, Snew_l2, Snew_h1, Snew_h2, &
     mask, mask_l1, mask_l2, mask_h1, mask_h2, &
     filt_T, S, domlo,domhi, &
     delta,xlo,problo,time,level) bind(C, name="ca_filt_prim")

  use network, only : naux, nspec
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
       small_temp, small_dens, nadv
  use filter_module

  use amrex_fort_module, only : rt => c_real
  implicit none

  integer, intent(in) :: lo(2), hi(2), domlo(2), domhi(2), level
  integer, intent(in) :: filt_T, S
  real(rt)        , intent(in) :: delta(2), xlo(2), problo(2), time
  integer, intent(in) ::   Stmp_l1,Stmp_h1,Stmp_l2,Stmp_h2
  integer, intent(in) ::   Snew_l1,Snew_h1,Snew_l2,Snew_h2
  integer, intent(in) ::   mask_l1,mask_h1,mask_l2,mask_h2
  real(rt)         :: Stmp(Stmp_l1:Stmp_h1,Stmp_l2:Stmp_h2,NVAR)
  real(rt)         :: Snew(Snew_l1:Snew_h1,Snew_l2:Snew_h2,NVAR)
  real(rt)         :: mask(mask_l1:mask_h1,mask_l2:mask_h2)
  ! mask has three possible values: -1.e0_rt, 0.e0_rt, and 1.e0_rt.
  ! -1.e0_rt appears only in cells that are covered by neither this level nor the finer level.
  !       It can only appear in ghost cells. 
  !  0.e0_rt appears only in cells that are covered by only this level, but not the finer level.
  !  1.e0_rt appears only in cells that are covered by the finer level.
  !       It can appear in either valid cells or ghost cells. 

  integer :: i,j
  real(rt)         :: p, e, X(nspec+naux)

  ! this is a stub -- a problem can override this in its own directory
  ! to implement filtering

end subroutine ca_filt_prim

