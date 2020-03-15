module bc_fill_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  ! Given a zone state, fill it with ambient material.

  subroutine fill_ambient(state, s_lo, s_hi, i, j, k, x, time)

    use amrex_constants_module, only: ZERO, HALF
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UTEMP, UEINT, UEDEN, UFS
#ifdef GRAVITY
    use meth_params_module, only: gravity_type_int, const_grav
#endif
    use network, only: nspec
    use prob_params_module, only: problo, probhi
    use probdata_module

    implicit none

    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    integer,  intent(in   ) :: i, j, k
    real(rt), intent(in   ) :: x, time

    real(rt) :: c_T

    real(rt) :: vel_g

    !$gpu

    c_T = problo(1) + center_T * (probhi(1) - problo(1))

    state(i,j,k,URHO) = ambient_dens
    state(i,j,k,UFS:UFS-1+nspec) = ambient_dens * ambient_comp

    if (x < c_T) then
       state(i,j,k,UTEMP) = T_l
       state(i,j,k,UEINT) = state(i,j,k,URHO) * ambient_e_l
       state(i,j,k,UMX) = state(i,j,k,URHO) * vel
#ifdef GRAVITY
       if (gravity_type_int == 0) then
          state(i,j,k,UMX) = state(i,j,k,UMX) + state(i,j,k,URHO) * const_grav * time
       end if
#endif
    else
       state(i,j,k,UTEMP) = T_r
       state(i,j,k,UEINT) = state(i,j,k,URHO) * ambient_e_r
       state(i,j,k,UMX) = -state(i,j,k,URHO) * vel
#ifdef GRAVITY
       if (gravity_type_int == 0) then
          state(i,j,k,UMX) = state(i,j,k,UMX) + state(i,j,k,URHO) * const_grav * time
       end if
#endif
    end if

    state(i,j,k,UMY:UMZ) = ZERO
    state(i,j,k,UEDEN) = state(i,j,k,UEINT) + HALF * sum(state(i,j,k,UMX:UMZ)**2) / state(i,j,k,URHO)

  end subroutine fill_ambient



  subroutine hypfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="hypfill")

    use amrex_filcc_module, only: amrex_filccn
    use amrex_constants_module, only: HALF
    use meth_params_module, only: NVAR
    use probdata_module, only: fill_ambient_bc
    use prob_params_module, only : problo

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k, n
    real(rt) :: x

    !$gpu

    ! First, use the generic filling routines to make sure we have valid data everywhere
    ! on physical domain ghost cells.

    call amrex_filccn(lo, hi, adv, adv_lo, adv_hi, NVAR, domlo, domhi, delta, xlo, bc)

    ! Override the generic routine at the physical boundaries by
    ! setting the material to the ambient state. Note that we
    ! don't want to do this for interior/periodic boundaries,
    ! which have bc == 0, or for reflecting boundaries, which have
    ! bc == -1 or bc == 1.

    if (fill_ambient_bc) then

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                x = problo(1) + (dble(i) + HALF)*delta(1)

                if (AMREX_SPACEDIM .ge. 1) then
                   if (i .lt. domlo(1) .and. (bc(1,1,1) .ne. -1 .and. bc(1,1,1) .ne. 0 .and. bc(1,1,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   else if (i .gt. domhi(1) .and. (bc(1,2,1) .ne. -1 .and. bc(1,2,1) .ne. 0 .and. bc(1,2,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   end if
                end if

                if (AMREX_SPACEDIM .ge. 2) then
                   if (j .lt. domlo(2) .and. (bc(2,1,1) .ne. -1 .and. bc(2,1,1) .ne. 0 .and. bc(2,1,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   else if (j .gt. domhi(2) .and. (bc(2,2,1) .ne. -1 .and. bc(2,2,1) .ne. 0 .and. bc(2,2,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   end if
                end if
                
                if (AMREX_SPACEDIM .eq. 3) then
                   if (k .lt. domlo(3) .and. (bc(3,1,1) .ne. -1 .and. bc(3,1,1) .ne. 0 .and. bc(3,1,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   else if (k .gt. domhi(3) .and. (bc(3,2,1) .ne. -1 .and. bc(3,2,1) .ne. 0 .and. bc(3,2,1) .ne. 1)) then
                      call fill_ambient(adv, adv_lo, adv_hi, i, j, k, x, time)
                   endif
                endif

             enddo
          enddo
       enddo

    endif

  end subroutine hypfill



  subroutine denfill(lo, hi, adv, adv_lo, adv_hi, domlo, domhi, delta, xlo, time, bc) bind(C, name="denfill")

    use amrex_filcc_module, only: amrex_filccn

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

end module bc_fill_module
