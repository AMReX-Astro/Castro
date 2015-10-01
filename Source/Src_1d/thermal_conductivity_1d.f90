! This routine computes the thermal conductivity at cell-centers

subroutine ca_therm_cond_cc(lo, hi, &
                            state,s_l1,s_h1, &
                            cond,c_l1,c_h1)
  use network
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX
  use probdata_module   , only : thermal_conductivity

  implicit none

  integer         , intent(in   ) :: lo(1), hi(1)
  integer         , intent(in   ) :: s_l1, s_h1
  integer         , intent(in   ) :: c_l1, c_h1
  real (kind=dp_t), intent(in   ) :: state(s_l1:s_h1,NVAR)
  real (kind=dp_t), intent(inout) :: cond(c_l1:c_h1)

  ! local variables
  integer          :: i

  ! fill the cell-centered conductivity

  do i = lo(1), hi(1)
     cond(i) = thermal_conductivity
  enddo

end subroutine ca_therm_cond_cc


! This routine fills the thermal conductivity on the edges of a zone
! by calling the cell-centered conductivity routine and averaging to
! the interfaces

subroutine ca_fill_temp_cond(lo,hi, &
                             state,s_l1,s_h1, &
                             coefx,cx_l1,cx_h1, &
                             coefy,cy_l1,cy_h1, dx)
  use network
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX
  use probdata_module   , only : thermal_conductivity

  implicit none

  integer         , intent(in   ) :: lo(1), hi(1)
  integer         , intent(in   ) :: s_l1, s_h1
  integer         , intent(in   ) :: cx_l1, cx_h1
  integer         , intent(in   ) :: cy_l1, cy_h1
  real (kind=dp_t), intent(in   ) :: state(s_l1:s_h1,NVAR)
  real (kind=dp_t), intent(inout) :: coefx(cx_l1:cx_h1)
  real (kind=dp_t), intent(inout) :: coefy(cy_l1:cy_h1)
  real (kind=dp_t), intent(in   ) :: dx(1)

  ! local variables
  integer          :: i
  double precision :: coef_cc(lo(1)-1:hi(1)+1)

  ! fill the cell-centered conductivity

  do i = lo(1)-1,hi(1)+1
     coef_cc(i) = thermal_conductivity
  enddo

  ! average to the interfaces
  do i = lo(1),hi(1)+1
     coefx(i) = 0.5d0 * (coef_cc(i) + coef_cc(i-1))
  end do

end subroutine ca_fill_temp_cond
