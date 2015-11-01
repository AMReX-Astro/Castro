! This routine fills the thermal conductivity on the edges of a zone
! by calling the cell-centered conductivity routine and averaging to
! the interfaces

subroutine ca_fill_temp_cond(lo,hi, &
                             state,s_l1,s_h1, &
                             coefx,cx_l1,cx_h1, &
                             coefy,cy_l1,cy_h1, dx)

  use bl_constants_module
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UEDEN, UTEMP, UFS, UFX, diffuse_cutoff_density
  use conductivity_module
  use eos_type_module

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

  type (eos_t) :: eos_state
  double precision :: cond

  ! fill the cell-centered conductivity

  do i = lo(1)-1,hi(1)+1
     eos_state%rho    = state(i,URHO)
     eos_state%T      = state(i,UTEMP)
     eos_state%xn(:)  = state(i,UFS:UFS-1+nspec)
     eos_state%aux(:) = state(i,UFX:UFX-1+naux)

     if (eos_state%rho > diffuse_cutoff_density) then
        call thermal_conductivity(eos_state, cond)
     else
        cond = ZERO
     endif

     coef_cc(i) = cond
  enddo

  ! average to the interfaces
  do i = lo(1),hi(1)+1
     coefx(i) = 0.5d0 * (coef_cc(i) + coef_cc(i-1))
  end do

end subroutine ca_fill_temp_cond
