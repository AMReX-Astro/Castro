subroutine ca_filt_prim(lo, hi, &
     Stmp, Stmp_l1, Stmp_h1, &
     Snew, Snew_l1, Snew_h1, &
     mask, mask_l1, mask_h1, &
     filt_T, S, domlo,domhi, &
     delta,xlo,problo,time,level) bind(C, name="ca_filt_prim")

  use network, only : naux, nspec
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
       small_temp, small_dens, nadv
  use filter_module

  implicit none

  integer, intent(in) :: lo(1), hi(1), domlo(1), domhi(1), level
  integer, intent(in) :: filt_T, S
  double precision, intent(in) :: delta(1), xlo(1), problo(1), time
  integer, intent(in) ::   Stmp_l1,Stmp_h1
  integer, intent(in) ::   Snew_l1,Snew_h1
  integer, intent(in) ::   mask_l1,mask_h1
  double precision :: Stmp(Stmp_l1:Stmp_h1,NVAR)
  double precision :: Snew(Snew_l1:Snew_h1,NVAR)
  double precision :: mask(mask_l1:mask_h1)
  ! mask has three possible values: -1.d0, 0.d0, and 1.d0.
  ! -1.d0 appears only in cells that are covered by neither this level nor the finer level.
  !       It can only appear in ghost cells. 
  !  0.d0 appears only in cells that are covered by only this level, but not the finer level.
  !  1.d0 appears only in cells that are covered by the finer level.
  !       It can appear in either valid cells or ghost cells. 


  ! this is a stub -- a problem can override this in its own directory
  ! to implement filtering

end subroutine ca_filt_prim

