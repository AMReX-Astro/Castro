module enthalpy_module

  implicit none

  public

contains

  ! This routine defines enthalpy from the EOS

  subroutine make_enthalpy(lo,hi, &
       state,s_lo,s_hi, &
       enth ,e_lo,e_hi) &
       bind(C, name="make_enthalpy")

    use bl_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UFS, UFX, diffuse_cutoff_density
    use conductivity_module
    use eos_type_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: e_lo(3), e_hi(3)
    real (kind=dp_t), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real (kind=dp_t), intent(inout) ::  enth(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))

    ! local variables
    integer          :: i, j, k

    type (eos_t) :: eos_state

    ! fill the cell-centered diffusion coefficient

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)

             if (eos_state%rho > diffuse_cutoff_density) then
                enth(i,j,k) = eos_state%h
             else
                enth(i,j,k) = ZERO
             endif
          enddo
       enddo
    enddo

  end subroutine make_enthalpy
  
end module enthalpy_module
