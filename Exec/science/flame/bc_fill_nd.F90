module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real
  use meth_params_module, only: NVAR

  implicit none

  include 'AMReX_bc_types.fi'

  public

contains

  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only : amrex_filccn
    use bc_ext_fill_module, only : ext_fill
    use meth_params_module, only : URHO, UTEMP, UMX, UMY, UMZ, UEINT, UEDEN, UFS
    use probdata_module
    use amrex_constants_module, only : ZERO, HALF

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

    if (AMREX_SPACEDIM == 1) then

       if (hi(1) > domhi(1)) then

          ! +x boundary
          if (bc(1,2,UMX) == FOEXTRAP) then
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   do i = domhi(1)+1, hi(1)

                      !adv(i,UMX) = adv(domhi(1),UMX)
                      adv(i,j,k,URHO)  = rho_fuel
                      adv(i,j,k,UMX)   = rho_fuel * v_inflow
                      adv(i,j,k,UMY)   = ZERO
                      adv(i,j,k,UMZ)   = ZERO
                      adv(i,j,k,UTEMP) = T_fuel
                      adv(i,j,k,UEINT) = rho_fuel * e_fuel
                      adv(i,j,k,UEDEN) = rho_fuel * e_fuel + HALF * sum(adv(i,j,k,UMX:UMZ)**2)/rho_fuel
                      adv(i,j,k,UFS:UFS+nspec-1) = rho_fuel * xn_fuel(:)
                   end do
                end do
             end do

          end if
       end if

    end if

  end subroutine hypfill


  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    use amrex_filcc_module, only: amrex_filccn
#ifndef AMREX_USE_CUDA
    use bc_ext_fill_module, only: ext_denfill
#endif

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    !$gpu

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, 1, domlo, domhi, delta, xlo, bc)

  end subroutine denfill


#ifdef REACTIONS
  subroutine reactfill(lo, hi, react, react_lo, react_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="reactfill")

    use amrex_filcc_module, only: amrex_filccn

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: react_lo(3), react_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: react(react_lo(1):react_hi(1),react_lo(2):react_hi(2),react_lo(3):react_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: d
    integer :: bc_temp(AMREX_SPACEDIM,2)

    !$gpu

    ! handle an external BC via extrapolation here
    bc_temp(:,:) = bc(:,:)

    do d = 1, AMREX_SPACEDIM
       if (bc(d,1) == EXT_DIR .and. react_lo(d) < domlo(d)) then
          bc_temp(d,1) = FOEXTRAP
       end if

       if (bc(d,2) == EXT_DIR .and. react_hi(d) > domhi(d)) then
          bc_temp(d,2) = FOEXTRAP
       end if
    end do

    call amrex_filccn(lo, hi, react, react_lo, react_hi, 1, domlo, domhi, delta, xlo, bc_temp)

  end subroutine reactfill
#endif


end module bc_fill_module
