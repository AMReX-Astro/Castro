! This routine fills the thermal conductivity on the edges of a zone
! by calling the cell-centered conductivity routine and averaging to
! the interfaces

subroutine ca_fill_temp_cond(lo,hi, &
                             state,s_l1,s_l2,s_h1,s_h2, &
                             coefx,cx_l1,cx_l2,cx_h1,cx_h2, &
                             coefy,cy_l1,cy_l2,cy_h1,cy_h2, dx)

  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX
  use conductivity_module
  use eos_type_module

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
  integer          :: i, j, n
  double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

  type (eos_t) :: eos_state
  double precision :: cond

  ! fill the cell-centered conductivity
  do j = lo(2)-1,hi(2)+1
     do i = lo(1)-1,hi(1)+1
        eos_state%rho    = state(i,j,URHO)
        eos_state%T      = state(i,j,UTEMP)
        eos_state%xn(:)  = state(i,j,UFS:UFS-1+nspec)
        eos_state%aux(:) = state(i,j,UFX:UFX-1+naux)

        call thermal_conductivity(eos_state, cond)

        coef_cc(i,j) = cond
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
