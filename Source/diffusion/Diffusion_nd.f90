module diffusion_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains


  !> @brief This routine fills the thermal conductivity on the edges of a zone
  !! by calling the cell-centered conductivity routine and averaging to
  !! the interfaces
  !!
  !! @note Binds to C function ``ca_fill_temp_cond``
  !!
  !! @param[in] lo integer
  !! @param[in] s_lo integer
  !! @param[in] cx_lo integer
  !! @param[in] state real(rt)
  !! @param[inout] coefx real(rt)
  !! @param[inout] coefy real(rt)
  !! @param[inout] coefz real(rt)
  !!
  subroutine ca_fill_temp_cond(lo,hi, &
       state,s_lo,s_hi, &
       coefx,cx_lo,cx_hi, &
       coefy,cy_lo,cy_hi, &
       coefz,cz_lo,cz_hi) &
       bind(C, name="ca_fill_temp_cond")



    use amrex_constants_module
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UEINt, UFS, UFX, &
         diffuse_cutoff_density, diffuse_cutoff_density_hi, diffuse_cond_scale_fac, &
         small_temp
    use prob_params_module, only : dg
    use conductivity_module
    use eos_type_module
    use eos_module, only : eos
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: s_lo(3), s_hi(3)
    integer         , intent(in   ) :: cx_lo(3), cx_hi(3), cy_lo(3), cy_hi(3), cz_lo(3), cz_hi(3)
    real(rt)        , intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt)        , intent(inout) :: coefx(cx_lo(1):cx_hi(1),cx_lo(2):cx_hi(2),cx_lo(3):cx_hi(3))
    real(rt)        , intent(inout) :: coefy(cy_lo(1):cy_hi(1),cy_lo(2):cy_hi(2),cy_lo(3):cy_hi(3))
    real(rt)        , intent(inout) :: coefz(cz_lo(1):cz_hi(1),cz_lo(2):cz_hi(2),cz_lo(3):cz_hi(3))

    ! local variables
    integer          :: i, j, k
    real(rt)         :: coef_cc(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    type (eos_t) :: eos_state
    real(rt) :: multiplier

    ! fill the cell-centered conductivity

    do k = lo(3)-1*dg(3),hi(3)+1*dg(3)
       do j = lo(2)-1*dg(2),hi(2)+1*dg(2)
          do i = lo(1)-1*dg(1),hi(1)+1*dg(1)

             eos_state%rho    = state(i,j,k,URHO)
             eos_state%T      = state(i,j,k,UTEMP)   ! needed as an initial guess
             eos_state%e      = state(i,j,k,UEINT)/state(i,j,k,URHO)
             eos_state%xn(:)  = state(i,j,k,UFS:UFS-1+nspec)/ state(i,j,k,URHO)
             eos_state%aux(:) = state(i,j,k,UFX:UFX-1+naux)/ state(i,j,k,URHO)

             if (eos_state%e < ZERO) then
                eos_state%T = small_temp
                call eos(eos_input_rt,eos_state)
             else
                call eos(eos_input_re,eos_state)
             endif

             if (eos_state%rho > diffuse_cutoff_density) then
                call conductivity(eos_state)

                if (eos_state%rho < diffuse_cutoff_density_hi) then
                    multiplier = (eos_state%rho - diffuse_cutoff_density) / &
                            (diffuse_cutoff_density_hi - diffuse_cutoff_density)
                    eos_state % conductivity = eos_state % conductivity * multiplier
                endif
             else
                eos_state % conductivity = ZERO
             endif
             coef_cc(i,j,k) = diffuse_cond_scale_fac * eos_state % conductivity
          enddo
       enddo
    enddo

    ! average to the interfaces
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1*dg(1)
             coefx(i,j,k) = 0.5e0_rt * (coef_cc(i,j,k) + coef_cc(i-1*dg(1),j,k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1*dg(2)
          do i = lo(1),hi(1)
             coefy(i,j,k) = 0.5e0_rt * (coef_cc(i,j,k) + coef_cc(i,j-1*dg(2),k))
          end do
       end do
    enddo

    do k = lo(3),hi(3)+1*dg(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             coefz(i,j,k) = 0.5e0_rt * (coef_cc(i,j,k) + coef_cc(i,j,k-1*dg(3)))
          end do
       end do
    enddo

  end subroutine ca_fill_temp_cond

end module diffusion_module
