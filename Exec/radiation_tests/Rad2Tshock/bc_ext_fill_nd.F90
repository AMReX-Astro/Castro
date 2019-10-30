module bc_ext_fill_module

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  ! ::: -----------------------------------------------------------

  subroutine ext_fill(lo, hi, adv, adv_lo, adv_hi, &
                      domlo, domhi, delta, xlo, time, bc) bind(C, name="ext_fill")

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UFX, UTEMP
    use network, only : nspec, naux
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_rt
    use probdata_module

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2, NVAR)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3),NVAR)
    real(rt), intent(in   ), value :: time

    integer i, j, k, n


    ! The strategy here is to set Dirichlet condition for inflow and
    ! outflow boundaries, and let the Riemann solver sort out the
    ! proper upwinding.  However, this decision makes this routine
    ! look somewhat non-orthodox, in that we need to set external
    ! values in either case....how do we know it's Outflow?  We have
    ! to assume that the setup routines converted Outflow to FOEXTRAP.

    !$gpu

    if (idir == 1) then

       !XLO
       if ( bc(1,1,1) .eq. EXT_DIR .and. adv_lo(1) .lt. domlo(1)) then
          do k = adv_lo(3), adv_hi(3)
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), domlo(1)-1
                   adv(i,j,k,URHO) = rho0
                   adv(i,j,k,UMX) = rho0*v0
                   adv(i,j,k,UMY) = 0.e0_rt
                   adv(i,j,k,UMZ) = 0.e0_rt
                   adv(i,j,k,UFS) = adv(i,j,k,URHO)
                   if (naux > 0) then
                      adv(i,j,k,UFX) = adv(i,j,k,URHO)
                   end if
                   adv(i,j,k,UEINT) = eint0
                   adv(i,j,k,UEDEN) = etot0
                   adv(i,j,k,UTEMP) = T0
                end do
             end do
          end do
       end if

       !     XHI
       if ( bc(1,2,1) .eq. EXT_DIR .and. adv_hi(1) .gt. domhi(1)) then
          do k = adv_lo(3), adv_hi(3)
             do j = adv_lo(2), adv_hi(2)
                do i = domhi(1)+1, adv_hi(1)
                   adv(i,j,k,URHO) = rho1
                   adv(i,j,k,UMX) = rho1*v1
                   adv(i,j,k,UMY) = 0.e0_rt
                   adv(i,j,k,UMZ) = 0.e0_rt
                   adv(i,j,k,UFS) = adv(i,j,k,URHO)
                   if (naux > 0) then
                      adv(i,j,k,UFX) = adv(i,j,k,URHO)
                   end if
                   adv(i,j,k,UEINT) = eint1
                   adv(i,j,k,UEDEN) = etot1
                   adv(i,j,k,UTEMP) = T1
                end do
             end do
          end do
       end if

    else if (idir == 2) then
       !YLO
       if ( bc(2,1,1) .eq. EXT_DIR .and. adv_lo(2) .lt. domlo(2)) then
          do k = adv_lo(3), adv_hi(3)
             do j = adv_lo(2), domlo(2)-1
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k,URHO) = rho0
                   adv(i,j,k,UMX) = 0.e0_rt
                   adv(i,j,k,UMY) = rho0*v0
                   adv(i,j,k,UMZ) = 0.e0_rt
                   adv(i,j,k,UFS) = adv(i,j,k,URHO)
                   if (naux > 0) then
                      adv(i,j,k,UFX) = adv(i,j,k,URHO)
                   end if
                   adv(i,j,k,UEINT) = eint0
                   adv(i,j,k,UEDEN) = etot0
                   adv(i,j,k,UTEMP) = T0
                end do
             end do
          end do
       end if

       !     YHI
       if ( bc(2,2,1) .eq. EXT_DIR .and. adv_hi(2) .gt. domhi(2)) then
          do k = adv_lo(3), adv_hi(3)
             do j = domhi(2)+1, adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k,URHO) = rho1
                   adv(i,j,k,UMX) = 0.e0_rt
                   adv(i,j,k,UMY) = rho1*v1
                   adv(i,j,k,UMZ) = 0.e0_rt
                   adv(i,j,k,UFS) = adv(i,j,k,URHO)
                   if (naux > 0) then
                      adv(i,j,k,UFX) = adv(i,j,k,URHO)
                   end if
                   adv(i,j,k,UEINT) = eint1
                   adv(i,j,k,UEDEN) = etot1
                   adv(i,j,k,UTEMP) = T1
                end do
             end do
          end do
       end if

    else if (idir == 3) then
       !ZLO
       if ( bc(3,1,1) .eq. EXT_DIR .and. adv_lo(3) .lt. domlo(3)) then
          do k = adv_lo(3), domlo(3)-1
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k,URHO) = rho0
                   adv(i,j,k,UMX) = 0.e0_rt
                   adv(i,j,k,UMY) = 0.e0_rt
                   adv(i,j,k,UMZ) = rho0*v0
                   adv(i,j,k,UFS) = adv(i,j,k,URHO)
                   if (naux > 0) then
                      adv(i,j,k,UFX) = adv(i,j,k,URHO)
                   end if
                   adv(i,j,k,UEINT) = eint0
                   adv(i,j,k,UEDEN) = etot0
                   adv(i,j,k,UTEMP) = T0
                enddo
             end do
          end do
       end if

       !     ZHI
       if ( bc(3,2,1) .eq. EXT_DIR .and. adv_hi(3) .gt. domhi(3)) then
          do k = domhi(3)+1, adv_hi(3)
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k,URHO) = rho1
                   adv(i,j,k,UMX) = 0.e0_rt
                   adv(i,j,k,UMY) = 0.e0_rt
                   adv(i,j,k,UMZ) = rho1*v1
                   adv(i,j,k,UFS) = adv(i,j,k,URHO)
                   if (naux > 0) then
                      adv(i,j,k,UFX) = adv(i,j,k,URHO)
                   end if
                   adv(i,j,k,UEINT) = eint1
                   adv(i,j,k,UEDEN) = etot1
                   adv(i,j,k,UTEMP) = T1
                end do
             end do
          end do
       end if

    end if

  end subroutine ext_fill

  ! :::
  ! ::: -----------------------------------------------------------
  ! :::

  subroutine ext_denfill(lo, hi, adv, adv_lo, adv_hi, &
                         domlo, domhi, delta, xlo, time, bc) bind(C, name="ext_denfill")

    use probdata_module

    implicit none

    include 'AMReX_bc_types.fi'

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: adv_lo(3), adv_hi(3)
    integer,  intent(in   ) :: domlo(3), domhi(3)
    integer,  intent(in   ) :: bc(AMREX_SPACEDIM, 2)
    real(rt), intent(in   ) :: delta(3), xlo(3)
    real(rt), intent(inout) :: adv(adv_lo(1):adv_hi(1),adv_lo(2):adv_hi(2),adv_lo(3):adv_hi(3))
    real(rt), intent(in   ), value :: time

    integer :: i, j, k

    !$gpu

    !     Note: this function should not be needed, technically, but is provided
    !     to filpatch because there are many times in the algorithm when just
    !     the density is needed.  We try to rig up the filling so that the same
    !     function is called here and in hypfill where all the states are filled.

    if (idir == 1) then
       !     XLO
       if ( bc(1,1) .eq. EXT_DIR .and. adv_lo(1) .lt. domlo(1)) then
          do k = adv_lo(3), adv_hi(3)
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), domlo(1)-1
                   adv(i,j,k) = rho0
                end do
             end do
          end do
       end if

       !     XHI
       if ( bc(1,2) .eq. EXT_DIR .and. adv_hi(1) .gt. domhi(1)) then
          do k = adv_lo(3), adv_hi(3)
             do j = adv_lo(2), adv_hi(2)
                do i = domhi(1)+1, adv_hi(1)
                   adv(i,j,k) = rho1
                end do
             end do
          end do
       end if

    else if (idir == 2) then
       !     YLO
       if ( bc(2,1) .eq. EXT_DIR .and. adv_lo(2) .lt. domlo(2)) then
          do k = adv_lo(3), adv_hi(3)
             do j = adv_lo(2), domlo(2)-1
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k) = rho0
                end do
             end do
          end do
       end if

       !     YHI
       if ( bc(2,2) .eq. EXT_DIR .and. adv_hi(2) .gt. domhi(2)) then
          do k = adv_lo(3), adv_hi(3)
             do j = domhi(2)+1, adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k) = rho1
                end do
             end do
          end do
       end if

    else if (idir == 3) then
       !     ZLO
       if ( bc(3,1) .eq. EXT_DIR .and. adv_lo(3) .lt. domlo(3)) then
          do k = adv_lo(3), domlo(3)-1
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k) = rho0
                end do
             end do
          end do
       end if

       !     ZHI
       if ( bc(3,2) .eq. EXT_DIR .and. adv_hi(3) .gt. domhi(3)) then
          do k = domhi(3)+1, adv_hi(3)
             do j = adv_lo(2), adv_hi(2)
                do i = adv_lo(1), adv_hi(1)
                   adv(i,j,k) = rho1
                end do
             end do
          end do
       end if

    end if

  end subroutine ext_denfill

end module bc_ext_fill_module
