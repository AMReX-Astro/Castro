! no tiling
subroutine ca_correct_dterm(dfx, dfx_l1, dfx_h1, &
     re, rc) bind(C, name="ca_correct_dterm")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: dfx_l1, dfx_h1
  real(rt)        , intent(inout) :: dfx(dfx_l1:dfx_h1)
  real(rt)        , intent(in) :: re(dfx_l1:dfx_h1), rc(1)

  integer :: i

  do i=dfx_l1, dfx_h1
     dfx(i) = dfx(i) / (re(i) + 1.e-50_rt)
  end do

end subroutine ca_correct_dterm
