! This routine fills the species conductivity on the edges of a zone
! by calling the cell-centered conductivity routine and averaging to
! the interfaces

subroutine ca_fill_spec_cond(lo,hi, &
                             state,s_lo,s_hi, &
                             coefx,cx_lo,cx_hi, &
                             coefy,cy_lo,cy_hi, &
                             coefz,cz_lo,cz_hi, dx)

  use bl_constants_module
  use network, only: nspec, naux
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX, diffuse_cutoff_density
  use prob_params_module, only : dg
  use conductivity_module
  use eos_type_module

  implicit none

  integer         , intent(in   ) :: lo(3), hi(3)
  integer         , intent(in   ) :: s_lo(3), s_hi(3)
  integer         , intent(in   ) :: cx_lo(3), cx_hi(3), cy_lo(3), cy_hi(3), cz_lo(3), cz_hi(3)
  real (kind=dp_t), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
  real (kind=dp_t), intent(inout) :: coefx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
  real (kind=dp_t), intent(inout) :: coefy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
  real (kind=dp_t), intent(inout) :: coefz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))
  real (kind=dp_t), intent(in   ) :: dx(3)

  ! local variables
  integer          :: i, j, k
  double precision :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

  type (eos_t) :: eos_state
  double precision :: cond

  ! fill the cell-centered conductivity

  do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
     do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
        do i = lo(1)-1*dg(1),hi(1)+1*dg(1)
           eos_state%rho    = state(i,j,k,URHO)
           eos_state%T      = state(i,j,k,UTEMP)
           eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
           eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

           if (eos_state%rho > diffuse_cutoff_density) then
              call thermal_conductivity(eos_state, cond)
              cond = cond / eos_state%cp
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
        do i = lo(1),hi(1)+1*dg(1)
           coefx(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i-1*dg(1),j,k))
        end do
     end do
  enddo

  do k = lo(3),hi(3)
     do j = lo(2),hi(2)+1*dg(2)
        do i = lo(1),hi(1)
           coefy(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j-1*dg(2),k))
        end do
     end do
  enddo

  do k = lo(3),hi(3)+1*dg(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)
           coefz(i,j,k) = 0.5d0 * (coef_cc(i,j,k) + coef_cc(i,j,k-1*dg(3)))
        end do
     end do
  enddo

end subroutine ca_fill_spec_cond
