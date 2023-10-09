! no tiling
subroutine ca_correct_dterm(  &
                            dfx, dfx_l1, dfx_l2, dfx_h1, dfx_h2, &
                            dfy, dfy_l1, dfy_l2, dfy_h1, dfy_h2, &
                            re, rc) bind(C, name="ca_correct_dterm")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer, intent(in) :: dfx_l1, dfx_l2, dfx_h1, dfx_h2
  integer, intent(in) :: dfy_l1, dfy_l2, dfy_h1, dfy_h2
  real(rt)        , intent(inout) :: dfx(dfx_l1:dfx_h1,dfx_l2:dfx_h2)
  real(rt)        , intent(inout) :: dfy(dfy_l1:dfy_h1,dfy_l2:dfy_h2)
  real(rt)        , intent(in) :: re(dfx_l1:dfx_h1), rc(dfy_l1:dfy_h1)

  integer :: i, j

  do j=dfx_l2, dfx_h2
     do i=dfx_l1, dfx_h1
        dfx(i,j) = dfx(i,j) / (re(i) + 1.e-50_rt)
     end do
  end do

  do j=dfy_l2, dfy_h2
     do i=dfy_l1, dfy_h1
        dfy(i,j) = dfy(i,j) / rc(i)
     end do
  end do

end subroutine ca_correct_dterm
