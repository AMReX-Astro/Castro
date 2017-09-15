module enthalpy_module

  use amrex_fort_module, only : rt => amrex_real
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
    use meth_params_module, only : NVAR, URHO, UTEMP, UEINT, UFS, UFX, diffuse_cutoff_density
    use eos_type_module
    use eos_module, only : eos
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: e_lo(3), e_hi(3)
    real(rt)        , intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        , intent(inout) ::  enth(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))

    ! local variables
    integer          :: i, j, k
    integer          :: i_lo, i_hi
    integer          :: j_lo, j_hi
    integer          :: k_lo, k_hi

    type (eos_t) :: eos_state

    i_lo = max(e_lo(1),s_lo(1))
    i_hi = min(e_hi(1),s_hi(1))
    j_lo = max(e_lo(2),s_lo(2))
    j_hi = min(e_hi(2),s_hi(2))
    k_lo = max(e_lo(3),s_lo(3))
    k_hi = min(e_hi(3),s_hi(3))

    ! Fill the cell-centered diffusion coefficient
    do k = k_lo, k_hi
       do j = j_lo, j_hi
          do i = i_lo, i_hi
             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)
             eos_state%e      = state(i,j,k,UEINT)/state(i,j,k,URHO)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)/ state(i,j,k,URHO)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)/ state(i,j,k,URHO)
             call eos(eos_input_re,eos_state)

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
