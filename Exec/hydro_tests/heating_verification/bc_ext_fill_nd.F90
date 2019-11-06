module bc_ext_fill_module
    ! this module contains different routines for filling the
    ! hydrodynamics boundary conditions

    ! .. note::
    !    the hydrostatic boundary conditions here rely on
    !    constant gravity

  use amrex_constants_module, only: ZERO, HALF, ONE, TWO
#ifndef AMREX_USE_CUDA
  use castro_error_module, only: castro_error
#endif
  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
                                 UEDEN, UEINT, UFS, UTEMP, const_grav, T_guess
  use prob_params_module, only: dim

  implicit none

  include 'AMReX_bc_types.fi'

contains


  subroutine ext_fill(lo, hi, adv, adv_lo, adv_hi, &
                      domlo, domhi, delta, xlo, time, bc) &
                      bind(C, name="ext_fill")

    use prob_params_module, only : problo, dim, center
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only: nspec
    use probdata_module
    use actual_eos_module, only : gamma_const

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2,NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k, n
    integer :: jmin, jmax, kmin, kmax

    real(rt) :: xn(nspec)
    real(rt) :: z, H, const
    integer :: npts_1d

    type (eos_t) :: eos_state

    real(rt), allocatable :: pressure(:), density(:), temp(:), eint(:)

    if (bc(1,1,URHO) == EXT_DIR .and. lo(1) < domlo(1)) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = domlo(1)-1, adv_lo(1), -1

                adv(i,j,k,URHO) = rho_in
                adv(i,j,k,UFS:UFS-1+nspec) = 0.0_rt
                adv(i,j,k,UFS) = rho_in

                adv(i,j,k,UMX) = rho_in * u_in
                adv(i,j,k,UMY) = 0.0_rt
                adv(i,j,k,UMZ) = 0.0_rt

                adv(i,j,k,UEINT) = rho_in * e_in
                adv(i,j,k,UEDEN) = rho_in * e_in + 0.5_rt * rho_in * u_in**2
                adv(i,j,k,UTEMP) = T_in

             end do
          end do
       end do

    endif


#ifndef AMREX_USE_CUDA
    if (bc(1,2,URHO) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for +x BC, d")
    end if

#if AMREX_SPACEDIM >= 2
    if (bc(2,1,URHO) == EXT_DIR .and. lo(2) < domlo(2)) then
       call castro_error("ERROR: HSE boundaries not implemented for -y BC")
    endif

    if (bc(2,2,URHO) == EXT_DIR .and. hi(2) > domhi(2)) then
       call castro_error("ERROR: HSE boundaries not implemented for +y BC, d")
    end if
#endif

#if AMREX_SPACEDIM == 3
    if (bc(3,1,URHO) == EXT_DIR .and. lo(3) < domlo(3)) then
       call castro_error("ERROR: HSE boundaries not implemented for -z BC")
    endif

    if (bc(3,2,URHO) == EXT_DIR .and. hi(3) > domhi(3)) then
       call castro_error("ERROR: HSE boundaries not implemented for +z BC, d")
    end if
#endif
#endif

  end subroutine ext_fill


  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) &
                         bind(C, name="ext_denfill")

    use prob_params_module, only: problo, center
#ifndef AMREX_USE_CUDA
    use castro_error_module, only: castro_error
#endif
    use probdata_module
    use actual_eos_module, only : gamma_const

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(dim,2)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k
    integer :: jmin, jmax, kmin, kmax

    real(rt) :: z, H, const
    integer :: npts_1d
    real(rt), allocatable :: density(:)

    if (bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = domlo(1)-1, adv_lo(1), -1
                adv(i,j,k) = rho_in
             end do
          end do
       end do

    endif


#ifndef AMREX_USE_CUDA
    if (bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for +x BC, d")
    end if

    if (bc(2,1) == EXT_DIR .and. lo(2) < domlo(2)) then
       call castro_error("ERROR: HSE boundaries not implemented for -y BC")
    endif

    if (bc(2,2) == EXT_DIR .and. hi(2) > domhi(2)) then
       call castro_error("ERROR: HSE boundaries not implemented for +y BC, d")
    end if

    if (bc(3,1) == EXT_DIR .and. lo(3) < domlo(3)) then
       call castro_error("ERROR: HSE boundaries not implemented for -z BC")
    endif

    if (bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then
       call castro_error("ERROR: HSE boundaries not implemented for +z BC, d")
    end if
#endif


  end subroutine ext_denfill

end module bc_ext_fill_module
