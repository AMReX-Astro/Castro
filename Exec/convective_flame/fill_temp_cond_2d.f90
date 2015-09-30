subroutine ca_fill_temp_cond(lo,hi, &
                             state,s_l1,s_l2,s_h1,s_h2, &
                             coefx,cx_l1,cx_l2,cx_h1,cx_h2, &
                             coefy,cy_l1,cy_l2,cy_h1,cy_h2, dx)
  use network
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX
  use probdata_module   , only : thermal_conductivity

  implicit none

  integer         , intent(in   ) :: lo(2), hi(2)
  integer         , intent(in   ) :: s_l1, s_l2, s_h1, s_h2
  integer         , intent(in   ) :: cx_l1, cx_l2, cx_h1, cx_h2
  integer         , intent(in   ) :: cy_l1, cy_l2, cy_h1, cy_h2
  real (kind=dp_t), intent(in   ) :: state(s_l1:s_h1,s_l2:s_h2,NVAR)
  real (kind=dp_t), intent(inout) :: coefx(cx_l1:cx_h1,cx_l2:cx_h2)
  real (kind=dp_t), intent(inout) :: coefy(cy_l1:cy_h1,cy_l2:cy_h2)
  real (kind=dp_t), intent(in   ) :: dx(2)

  ! local variables
  integer          :: i,j,n
  double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

  ! fill the cell-centered conductivity
  do j = lo(2)-1,hi(2)+1
     do i = lo(1)-1,hi(1)+1
        coef_cc(i,j) = thermal_conductivity
     enddo
  enddo

  ! average to the interfaces
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        coefx(i,j) = 0.5d0 * (coef_cc(i,j) + coef_cc(i-1,j))
     end do
  end do
  
  do j = lo(2),hi(2)+1
     do i = lo(1),hi(1)
        coefy(i,j) = 0.5d0 * (coef_cc(i,j) + coef_cc(i,j-1))
     end do
  end do
  
end subroutine ca_fill_temp_cond
