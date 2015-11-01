! This routine fills the thermal conductivity on the edges of a zone
! by calling the cell-centered conductivity routine and averaging to
! the interfaces

subroutine ca_fill_temp_cond(lo,hi, &
                             state,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                             coefx,cx_l1,cx_l2,cx_l3,cx_h1,cx_h2,cx_h3, &
                             coefy,cy_l1,cy_l2,cy_l3,cy_h1,cy_h2,cy_h3, &
                             coefz,cz_l1,cz_l2,cz_l3,cz_h1,cz_h2,cz_h3, dx)

  use bl_constants_module
  use network, only: nspec, naux
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX, diffuse_cutoff_density
  use conductivity_module
  use eos_type_module

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3)
  integer         , intent(in   ) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
  integer         , intent(in   ) :: cx_l1, cx_l2, cx_l3, cx_h1, cx_h2, cx_h3
  integer         , intent(in   ) :: cy_l1, cy_l2, cy_l3, cy_h1, cy_h2, cy_h3
  integer         , intent(in   ) :: cz_l1, cz_l2, cz_l3, cz_h1, cz_h2, cz_h3
  real (kind=dp_t), intent(in   ) :: state(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NVAR)
  real (kind=dp_t), intent(inout) :: coefx(cx_l1:cx_h1,cx_l2:cx_h2,cx_l3:cx_h3)
  real (kind=dp_t), intent(inout) :: coefy(cy_l1:cy_h1,cy_l2:cy_h2,cy_l3:cy_h3)
  real (kind=dp_t), intent(inout) :: coefz(cz_l1:cz_h1,cz_l2:cz_h2,cz_l3:cz_h3)
  real (kind=dp_t), intent(in   ) :: dx(3)

  ! local variables
  integer          :: i, j, k
  double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

  type (eos_t) :: eos_state
  double precision :: cond

  ! fill the cell-centered conductivity

  do k = lo(3)-1,hi(3)+1
     do j = lo(2)-1,hi(2)+1
        do i = lo(1)-1,hi(1)+1
           eos_state%rho    = state(i,j,k,URHO)
           eos_state%T      = state(i,j,k,UTEMP)
           eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
           eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

           if (eos_state%rho > diffuse_cutoff_density) then
              call thermal_conductivity(eos_state, cond)
           else
              cond = ZERO
           endif

           coef_cc(i,j,k) = cond
        enddo
     enddo
  enddo

  ! average to the interfaces
  do k = lo(3),hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
           coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1,j,k))
        end do
     end do
  enddo

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)+1
        do i = lo(1),hi(1)
           coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1,k))
        end do
     end do
  enddo

  do k = lo(3),hi(3)+1
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1))
        end do
     end do
  enddo

end subroutine ca_fill_temp_cond
