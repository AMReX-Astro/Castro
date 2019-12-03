module diffusion_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine ca_fill_temp_cond(lo, hi, &
                               state, s_lo, s_hi, &
                               coef, c_lo, c_hi) &
                               bind(C, name="ca_fill_temp_cond")
    ! This routine fills the thermal conductivity at zone centers
    ! by calling the cell-centered conductivity routine
    !
    ! .. note::
    !    Binds to C function ``ca_fill_temp_cond``

    use amrex_constants_module, only: ZERO, ONE
    use network, only: nspec, naux
    use meth_params_module, only : NVAR, URHO, UTEMP, UEINt, UFS, UFX, &
                                   diffuse_cutoff_density, diffuse_cutoff_density_hi, diffuse_cond_scale_fac, &
                                   small_temp
    use conductivity_module, only: conductivity
    use eos_type_module, only: eos_input_rt, eos_input_re, eos_t
    use eos_module, only: eos

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3) ! lower limit
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: coef(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))

    ! local variables
    integer  :: i, j, k

    type (eos_t) :: eos_state
    real(rt) :: multiplier, rhoinv

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             eos_state % rho    = state(i,j,k,URHO)
             rhoinv = ONE/eos_state % rho
             eos_state % T      = state(i,j,k,UTEMP)   ! needed as an initial guess
             eos_state % e      = state(i,j,k,UEINT) * rhoinv
             eos_state % xn(:)  = state(i,j,k,UFS:UFS-1+nspec) * rhoinv
             eos_state % aux(:) = state(i,j,k,UFX:UFX-1+naux) * rhoinv

             if (eos_state%e < ZERO) then
                eos_state%T = small_temp
                call eos(eos_input_rt,eos_state)
             else
                call eos(eos_input_re,eos_state)
             endif

             if (eos_state%rho > diffuse_cutoff_density) then
                call conductivity(eos_state)

                if (eos_state%rho < diffuse_cutoff_density_hi) then
                    multiplier = (eos_state % rho - diffuse_cutoff_density) / &
                            (diffuse_cutoff_density_hi - diffuse_cutoff_density)
                    eos_state % conductivity = eos_state % conductivity * multiplier
                endif
             else
                eos_state % conductivity = ZERO
             endif
             coef(i,j,k) = diffuse_cond_scale_fac * eos_state % conductivity
          end do
       end do
    end do

  end subroutine ca_fill_temp_cond

  subroutine ca_average_coef_cc_to_ec(lo, hi, &
                                      coef_c, c_lo, c_hi, &
                                      coef_e, e_lo, e_hi, &
                                      dir) bind(c, name="ca_average_coef_cc_to_ec")
    ! This routine averages cell-centered conductivity coefficients to zone edges
    !
    ! .. note::
    !    Binds to C function ``ca_average_coef_cc_to_ec``

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: e_lo(3), e_hi(3)
    real(rt), intent(in   ) :: coef_c(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(inout) :: coef_e(e_lo(1):e_hi(1),e_lo(2):e_hi(2),e_lo(3):e_hi(3))
    integer,  intent(in   ), value :: dir

    integer :: i, j, k

    !$gpu

    if (dir == 1) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                coef_e(i,j,k) = 0.5e0_rt * (coef_c(i,j,k) + coef_c(i-1,j,k))
             end do
          end do
       end do

    else if (dir == 2) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                coef_e(i,j,k) = 0.5e0_rt * (coef_c(i,j,k) + coef_c(i,j-1,k))
             end do
          end do
       end do

    else

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                coef_e(i,j,k) = 0.5e0_rt * (coef_c(i,j,k) + coef_c(i,j,k-1))
             end do
          end do
       end do

    end if

  end subroutine ca_average_coef_cc_to_ec

end module diffusion_module
