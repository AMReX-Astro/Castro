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

    ! compute background state.  Here we assume that the background is
    ! in the vertical direction.  These BCs only make sense for the
    ! boundary along the background direction (dimension
    ! AMREX_SPACEDIM).

    ! first make a 1D initial model for the entire domain + 4 ghost cells on either end
    npts_1d = (TWO*center(AMREX_SPACEDIM) + 1.e-8_rt) / delta(AMREX_SPACEDIM)

    allocate(pressure(-5:npts_1d+4))
    allocate(density (-5:npts_1d+4))
    allocate(temp    (-5:npts_1d+4))
    allocate(eint    (-5:npts_1d+4))

    const = pres_base/dens_base**gamma_const

    pressure(0) = pres_base
    density(0)  = dens_base

    ! only initialize the first species
    xn(1) = ONE

    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(const_grav)

    do j = 0, npts_1d+4

       ! initial guess
       temp(j) = T_guess

       if (do_isentropic) then
          z = dble(j) * delta(AMREX_SPACEDIM)
          density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
               (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
       else
          z = (dble(j)+HALF) * delta(AMREX_SPACEDIM)
          density(j) = dens_base * exp(-z/H)
       end if

       if (j > 0) then
          pressure(j) = pressure(j-1) - &
               delta(AMREX_SPACEDIM) * HALF * (density(j)+density(j-1)) * abs(const_grav)
       end if

       eos_state%rho = density(j)
       eos_state%T = temp(j)
       eos_state%p = pressure(j)
       eos_state%xn(:) = xn(:)

       call eos(eos_input_rp, eos_state)

       eint(j) = eos_state%e
       temp(j) = eos_state%T

    end do

    do j=-1,-5,-1

       ! initial guess
       temp(j) = T_guess

       if (do_isentropic) then
          z = dble(j) * delta(AMREX_SPACEDIM)
          density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
               (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
       else
          z = (dble(j)+HALF) * delta(AMREX_SPACEDIM)
          density(j) = dens_base * exp(-z/H)
       end if

       pressure(j) = pressure(j+1) + &
            delta(AMREX_SPACEDIM) * HALF * (density(j)+density(j+1)) * abs(const_grav)

       eos_state%rho = density(j)
       eos_state%T = temp(j)
       eos_state%p = pressure(j)
       eos_state%xn(:) = xn(:)

       call eos(eos_input_rp, eos_state)

       eint(j) = eos_state%e
       temp(j) = eos_state%T

    end do

    do n = 1, NVAR

#ifndef AMREX_USE_CUDA
       if (bc(1,1,n) == EXT_DIR .and. lo(1) < domlo(1)) then
          call castro_error("ERROR: HSE boundaries not implemented for -x BC")
       endif

       if (bc(1,2,n) == EXT_DIR .and. hi(1) > domhi(1)) then
          call castro_error("ERROR: HSE boundaries not implemented for +x BC, d")
       end if
#endif

#if AMREX_SPACEDIM >= 2
       !-------------------------------------------------------------------------
       ! y boundaries
       !-------------------------------------------------------------------------

       ! YLO
       if (bc(2,1,n) == EXT_DIR .and. lo(2) < domlo(2)) then

          ! we will fill all the variables when we consider URHO
          if (n == URHO) then
             do k = lo(3), hi(3)
                do i = lo(1), hi(1)

                   ! integrate downward
                   jmin = adv_lo(2)
                   jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
                   ! For CUDA, this should only be one thread doing the work:
                   ! we'll arbitrary choose the zone with index domlo(2) - 1.
                   if (hi(2) /= jmax) then
                      jmax = jmin - 1
                   end if
#endif
                   do j = jmax, jmin, -1

                      ! zero transverse momentum
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMZ) = ZERO

                      if (boundary_type .eq. 1) then
                         ! extrapolate normal momentum
                         ! enforces pi=0 at boundary
                         adv(i,j,k,UMY) = adv(i,domlo(2),k,UMY)
                      else
                         ! zero normal momentum
                         ! permits pi to pass through boundary
                         adv(i,j,k,UMY) = ZERO
                      end if

                      adv(i,j,k,URHO) = density(j)
                      adv(i,j,k,UFS) = adv(i,j,k,URHO)
                      adv(i,j,k,UEINT) = eint(j)*adv(i,j,k,URHO)
                      adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                           + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UTEMP) = temp(j)

                   end do
                end do
             end do
          endif  ! n == URHO

       endif


       ! YHI
       if (bc(2,2,n) == EXT_DIR .and. hi(2) > domhi(2)) then

          ! we will fill all the variables when we consider URHO
          if (n == URHO) then

             do k = lo(3), hi(3)
                do i = lo(1), hi(1)

                   jmin = domhi(2)+1
                   jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
                   if (lo(2) /= jmin) then
                      jmin = jmax + 1
                   end if
#endif
                   do j = jmin, jmax

                      ! zero transverse momentum
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMZ) = ZERO

                      if (boundary_type .eq. 1) then
                         ! extrapolate normal momentum
                         ! enforces pi=0 at boundary
                         adv(i,j,k,UMY) = adv(i,domhi(2),k,UMY)
                      else
                         ! zero normal momentum
                         ! permits pi to pass through boundary
                         adv(i,j,k,UMY) = ZERO
                      end if

                      adv(i,j,k,URHO) = density(j)
                      adv(i,j,k,UFS) = adv(i,j,k,URHO)
                      adv(i,j,k,UEINT) = eint(j)*adv(i,j,k,URHO)
                      adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                           + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UTEMP) = temp(j)

                   end do

                end do
             end do
          end if

       endif
#endif

#if AMREX_SPACEDIM == 3

       !-------------------------------------------------------------------------
       ! z boundaries
       !-------------------------------------------------------------------------

       ! ZLO
       if (bc(3,1,n) == EXT_DIR .and. lo(3) < domlo(3)) then

          ! we will fill all the variables when we consider URHO
          if (n == URHO) then
             do j= lo(2), hi(2)
                do i = lo(1), hi(1)

                   ! integrate downward
                   kmin = adv_lo(3)
                   kmax = domlo(3) - 1
#ifdef AMREX_USE_CUDA
                   if (hi(3) /= kmax) then
                      kmax = kmin - 1
                   end if
#endif
                   do k = kmax, kmin, -1

                      ! zero transverse momentum
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMY) = ZERO

                      if (boundary_type .eq. 1) then
                         ! extrapolate normal momentum
                         ! enforces pi=0 at boundary
                         adv(i,j,k,UMZ) = adv(i,j,domlo(3),UMZ)
                      else
                         ! zero normal momentum
                         ! permits pi to pass through boundary
                         adv(i,j,k,UMZ) = ZERO
                      end if

                      adv(i,j,k,URHO) = density(k)
                      adv(i,j,k,UFS) = adv(i,j,k,URHO)
                      adv(i,j,k,UEINT) = eint(k)*adv(i,j,k,URHO)
                      adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                           + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UTEMP) = temp(k)

                   end do
                end do
             end do

          endif  ! n == URHO

       endif

       ! ZHI
       if (bc(3,2,n) == EXT_DIR .and. hi(3) > domhi(3)) then

          ! we will fill all the variables when we consider URHO
          if (n == URHO) then
             do j= lo(2), hi(2)
                do i = lo(1), hi(1)

                   kmin = domhi(3) + 1
                   kmax = adv_hi(3)
#ifdef AMREX_USE_CUDA
                   if (lo(3) /= kmin) then
                      kmin = kmax + 1
                   end if
#endif
                   do k = kmin, kmax

                      ! zero transverse momentum
                      adv(i,j,k,UMX) = ZERO
                      adv(i,j,k,UMY) = ZERO

                      if (boundary_type .eq. 1) then
                         ! extrapolate normal momentum
                         ! enforces pi=0 at boundary
                         adv(i,j,k,UMZ) = adv(i,j,domhi(3),UMZ)
                      else
                         ! zero normal momentum
                         ! permits pi to pass through boundary
                         adv(i,j,k,UMZ) = ZERO
                      end if

                      adv(i,j,k,URHO) = density(k)
                      adv(i,j,k,UFS) = adv(i,j,k,URHO)
                      adv(i,j,k,UEINT) = eint(k)*adv(i,j,k,URHO)
                      adv(i,j,k,UEDEN) = adv(i,j,k,UEINT) &
                           + HALF*sum(adv(i,j,k,UMX:UMZ)**2)/adv(i,j,k,URHO)
                      adv(i,j,k,UTEMP) = temp(k)

                   end do

                end do
             end do
          end if

       end if
#endif

    end do

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

    ! Note: this function should not be needed, technically, but is
    ! provided to filpatch because there are many times in the algorithm
    ! when just the density is needed.  We try to rig up the filling so
    ! that the same function is called here and in hypfill where all the
    ! states are filled.

    ! compute background state.  Here we assume that the background is
    ! in the vertical direction.  These BCs only make sense for the
    ! boundary along the background direction (dimension
    ! AMREX_SPACEDIM).

    ! first make a 1D initial model for the entire domain + 4 ghost cells on either end
    npts_1d = (TWO*center(AMREX_SPACEDIM) + 1.e-8_rt) / delta(AMREX_SPACEDIM)

    allocate(density (-5:npts_1d+4))

    const = pres_base/dens_base**gamma_const

    density(0)  = dens_base


    ! compute the pressure scale height (for an isothermal, ideal-gas
    ! atmosphere)
    H = pres_base / dens_base / abs(const_grav)

    do j = 0, npts_1d+4

       if (do_isentropic) then
          z = dble(j) * delta(AMREX_SPACEDIM)
          density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
               (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
       else
          z = (dble(j)+HALF) * delta(AMREX_SPACEDIM)
          density(j) = dens_base * exp(-z/H)
       end if

    end do

    do j=-1,-5,-1

       if (do_isentropic) then
          z = dble(j) * delta(AMREX_SPACEDIM)
          density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
               (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
       else
          z = (dble(j)+HALF) * delta(AMREX_SPACEDIM)
          density(j) = dens_base * exp(-z/H)
       end if
    end do

#ifndef AMREX_USE_CUDA
    if (bc(1,1) == EXT_DIR .and. lo(1) < domlo(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for -x BC")
    endif

    if (bc(1,2) == EXT_DIR .and. hi(1) > domhi(1)) then
       call castro_error("ERROR: HSE boundaries not implemented for +x BC, d")
    end if
#endif

#if AMREX_SPACEDIM >= 2
    !-------------------------------------------------------------------------
    ! y boundaries
    !-------------------------------------------------------------------------

    ! YLO
    if (bc(2,1) == EXT_DIR .and. lo(2) < domlo(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             ! integrate downward
             jmin = adv_lo(2)
             jmax = domlo(2)-1
#ifdef AMREX_USE_CUDA
             ! For CUDA, this should only be one thread doing the work:
             ! we'll arbitrary choose the zone with index domlo(2) - 1.
             if (hi(2) /= jmax) then
                jmax = jmin - 1
             end if
#endif
             do j = jmax, jmin, -1
                adv(i,j,k) = density(j)
             end do
          end do
       end do

    endif


    ! YHI
    if (bc(2,2) == EXT_DIR .and. hi(2) > domhi(2)) then

       do k = lo(3), hi(3)
          do i = lo(1), hi(1)

             jmin = domhi(2)+1
             jmax = adv_hi(2)
#ifdef AMREX_USE_CUDA
             if (lo(2) /= jmin) then
                jmin = jmax + 1
             end if
#endif
             do j = jmin, jmax
                adv(i,j,k) = density(j)
             end do

          end do
       end do
    end if
#endif

#if AMREX_SPACEDIM == 3

    !-------------------------------------------------------------------------
    ! z boundaries
    !-------------------------------------------------------------------------

    ! ZLO
    if (bc(3,1) == EXT_DIR .and. lo(3) < domlo(3)) then

       do j= lo(2), hi(2)
          do i = lo(1), hi(1)

             ! integrate downward
             kmin = adv_lo(3)
             kmax = domlo(3) - 1
#ifdef AMREX_USE_CUDA
             if (hi(3) /= kmax) then
                kmax = kmin - 1
             end if
#endif
             do k = kmax, kmin, -1
                adv(i,j,k) = density(k)
             end do
          end do
       end do

    endif

    ! ZHI
    if (bc(3,2) == EXT_DIR .and. hi(3) > domhi(3)) then

       do j= lo(2), hi(2)
          do i = lo(1), hi(1)

             kmin = domhi(3) + 1
             kmax = adv_hi(3)
#ifdef AMREX_USE_CUDA
             if (lo(3) /= kmin) then
                kmin = kmax + 1
             end if
#endif
             do k = kmin, kmax
                adv(i,j,k) = density(k)
             end do

          end do
       end do
    end if
#endif

  end subroutine ext_denfill

end module bc_ext_fill_module
