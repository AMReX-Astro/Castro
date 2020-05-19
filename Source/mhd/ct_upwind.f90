module ct_upwind
  use amrex_mempool_module, only : bl_allocate, bl_deallocate
  use amrex_fort_module, only : rt => amrex_real
  use hlld_solver, only : hlld
  use meth_params_module
  use mhd_state_module

  implicit none

  private primtocons
  public corner_transport

  ! note: in this module, we use left and right to mean with respect
  ! to the interface.  So qleft, uleft, ul, ... are the left state on
  ! an interface and qright, uright, ur, ... are the right state on an
  ! interface.

  ! Miniati and Martin use + and - with respect to a zone center, so
  ! u+ would be on the right edge of a zone, which is the left state
  ! to an interface.

contains


  subroutine corner_transport(lo, hi, &
                              q, q_lo, q_hi, &
                              qleft, ql_lo, ql_hi, &
                              qright, qr_lo, qr_hi, &
                              flxx, flxx_lo , flxx_hi, &
                              flxy, flxy_lo , flxy_hi, &
                              flxz, flxz_lo , flxz_hi, &
                              Ex, ex_lo, ex_hi, &
                              Ey, ey_lo, ey_hi, &
                              Ez, ez_lo, ez_hi, &
                              dx, dt) bind(C, name="corner_transport")

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NQ
    use electric_field
    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: qr_lo(3), qr_hi(3)
    integer, intent(in) :: ql_lo(3), ql_hi(3)

    integer, intent(in) :: ex_lo(3), ex_hi(3)
    integer, intent(in) :: ey_lo(3), ey_hi(3)
    integer, intent(in) :: ez_lo(3), ez_hi(3)

    integer, intent(in) :: flxx_lo(3), flxx_hi(3)
    integer, intent(in) :: flxy_lo(3), flxy_hi(3)
    integer, intent(in) :: flxz_lo(3), flxz_hi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(in) :: qright(qr_lo(1):qr_hi(1), qr_lo(2):qr_hi(2), qr_lo(3):qr_hi(3), NQ, 3)
    real(rt), intent(in) :: qleft(ql_lo(1):ql_hi(1), ql_lo(2):ql_hi(2), ql_lo(3):ql_hi(3), NQ, 3)

    ! fluxes should be NVAR+3
    real(rt), intent(out) :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3), NVAR+3)   !Half Step Fluxes
    real(rt), intent(out) :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3), NVAR+3)   !Half Step Fluxes
    real(rt), intent(out) :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3), NVAR+3)   !Half Step Fluxes

    real(rt), intent(out)  :: Ex(ex_lo(1):ex_hi(1), ex_lo(2):ex_hi(2), ex_lo(3):ex_hi(3))
    real(rt), intent(out)  :: Ey(ey_lo(1):ey_hi(1), ey_lo(2):ey_hi(2), ey_lo(3):ey_hi(3))
    real(rt), intent(out)  :: Ez(ez_lo(1):ez_hi(1), ez_lo(2):ez_hi(2), ez_lo(3):ez_hi(3))

    ! these are conserved + magnetic field (cell centered)

    real(rt), pointer :: ux_left(:,:,:,:)
    real(rt), pointer :: ux_right(:,:,:,:)
    real(rt), pointer :: uy_left(:,:,:,:)
    real(rt), pointer :: uy_right(:,:,:,:)
    real(rt), pointer :: uz_left(:,:,:,:)
    real(rt), pointer :: uz_right(:,:,:,:)

    ! temporary conserved variables
    real(rt), pointer :: utmp_left(:,:,:,:)
    real(rt), pointer :: utmp_right(:,:,:,:)

    real(rt), pointer :: qtmp_left(:,:,:,:)
    real(rt), pointer :: qtmp_right(:,:,:,:)

    integer :: fx_lo(3), fx_hi(3)
    integer :: fy_lo(3), fy_hi(3)
    integer :: fz_lo(3), fz_hi(3)

    integer :: fxy_lo(3), fxy_hi(3)
    integer :: fyx_lo(3), fyx_hi(3)
    integer :: fxz_lo(3), fxz_hi(3)
    integer :: fzx_lo(3), fzx_hi(3)
    integer :: fyz_lo(3), fyz_hi(3)
    integer :: fzy_lo(3), fzy_hi(3)

    integer :: u_lo(3), u_hi(3)
    integer :: ut_lo(3), ut_hi(3)

    real(rt), pointer :: flxx1D(:,:,:,:)
    real(rt), pointer :: flxy1D(:,:,:,:)
    real(rt), pointer :: flxz1D(:,:,:,:)

    real(rt), pointer :: flx_xy(:,:,:,:)
    real(rt), pointer :: flx_xz(:,:,:,:)
    real(rt), pointer :: flx_yx(:,:,:,:)
    real(rt), pointer :: flx_yz(:,:,:,:)
    real(rt), pointer :: flx_zx(:,:,:,:)
    real(rt), pointer :: flx_zy(:,:,:,:)

    real(rt)  :: q2D(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt)  :: dx(3)
    real(rt), intent(in), value :: dt

    integer  :: i, j, k, work_lo(3), work_hi(3)

    ! MM CTU Step 1
    ! Calculate Flux 1D, eq.35

    ! x-dir
    ![lo(1)-2, lo(2)-3, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
    fx_lo = (/ lo(1)-2, lo(2)-3, lo(3)-3 /)
    fx_hi = (/ hi(1)+3, hi(2)+3, hi(3)+3 /)

    call bl_allocate(flxx1D, fx_lo, fx_hi, NVAR+3)

    call hlld(fx_lo, fx_hi, &
              qleft(:,:,:,:,1), ql_lo, ql_hi, &
              qright(:,:,:,:,1), qr_lo, qr_hi, &
              flxx1D, fx_lo, fx_hi, 1)

    !y-dir
    ![lo(1)-3, lo(2)-2, lo(3)-3] [hi(1)+3, hi(2)+3, hi(3)+3]
    fy_lo = (/ lo(1)-3, lo(2)-2, lo(3)-3 /)
    fy_hi = (/ hi(1)+3, hi(2)+3, hi(3)+3 /)

    call bl_allocate(flxy1D, fy_lo, fy_hi, NVAR+3)

    call hlld(fy_lo, fy_hi, &
              qleft(:,:,:,:,2), ql_lo, ql_hi, &
              qright(:,:,:,:,2), qr_lo, qr_hi, &
              flxy1D, fy_lo, fy_hi, 2)

    !z-dir
    ![lo(1)-3, lo(2)-3, lo(3)-2] [hi(1)+3, hi(2)+3, hi(3)+3]
    fz_lo = (/ lo(1)-3, lo(2)-3, lo(3)-2 /)
    fz_hi = (/ hi(1)+3, hi(2)+3, hi(3)+3 /)

    call bl_allocate(flxz1D, fz_lo, fz_hi, NVAR+3)

    call hlld(fz_lo, fz_hi, &
              qleft(:,:,:,:,3), ql_lo, ql_hi, &
              qright(:,:,:,:,3), qr_lo, qr_hi, &
              flxz1D, fz_lo, fz_hi, 3)

    !Prim to Cons
    u_lo = [lo(1)-2, lo(2)-2, lo(3)-2]
    u_hi = [hi(1)+2, hi(2)+2, hi(3)+2]

    call bl_allocate(ux_left, u_lo, u_hi, NVAR+3)
    call bl_allocate(ux_right, u_lo, u_hi, NVAR+3)

    call bl_allocate(uy_left, u_lo, u_hi, NVAR+3)
    call bl_allocate(uy_right, u_lo, u_hi, NVAR+3)

    call bl_allocate(uz_left, u_lo, u_hi, NVAR+3)
    call bl_allocate(uz_right, u_lo, u_hi, NVAR+3)

    call PrimToCons(u_lo, u_hi, qleft(:,:,:,:,1), ql_lo, ql_hi, ux_left, u_lo, u_hi)
    call PrimToCons(u_lo, u_hi, qright(:,:,:,:,1), qr_lo, qr_hi, ux_right, u_lo, u_hi)

    call PrimToCons(u_lo, u_hi, qleft(:,:,:,:,2), ql_lo, ql_hi, uy_left, u_lo, u_hi)
    call PrimToCons(u_lo, u_hi, qright(:,:,:,:,2), qr_lo, qr_hi, uy_right, u_lo, u_hi)

    call PrimToCons(u_lo, u_hi, qleft(:,:,:,:,3), ql_lo, ql_hi, uz_left, u_lo, u_hi)
    call PrimToCons(u_lo, u_hi, qright(:,:,:,:,3), qr_lo, qr_hi, uz_right, u_lo, u_hi)


    ! MM CTU Step 2
    ! Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields, eq.36

    ![lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+2, hi(2)+3, hi(3)+3]
    work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
    work_hi = (/ hi(1)+2, hi(2)+3, hi(3)+3 /)
    call electric_edge_x(work_lo, work_hi, &
                         q, q_lo, q_hi, &
                         Ex, ex_lo, ex_hi, &
                         flxy1D, fy_lo, fy_hi, &
                         flxz1D, fz_lo, fz_hi)

    ![lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+2, hi(3)+3]
    work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
    work_hi = (/ hi(1)+3, hi(2)+2, hi(3)+3 /)
    call electric_edge_y(work_lo, work_hi, &
                         q, q_lo, q_hi, &
                         Ey, ey_lo, ey_hi, &
                         flxx1D, fx_lo, fx_hi, &
                         flxz1D, fz_lo, fz_hi)

    ![lo(1)-2, lo(2)-2, lo(3)-2][hi(1)+3, hi(2)+3, hi(3)+2]
    work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-2 /)
    work_hi = (/ hi(1)+3, hi(2)+3, hi(3)+2 /)
    call electric_edge_z(work_lo, work_hi, &
                         q, q_lo, q_hi, &
                         Ez, ez_lo, ez_hi, &
                         flxx1D, fx_lo, fx_hi, &
                         flxy1D, fy_lo, fy_hi)

    ! MM CTU Steps 3, 4, and 5
    ! Corner Couple, eq. 37, 38 and 39 Correct Conservative vars using Transverse Fluxes

    ! X direction

    ! affected by Y Flux
    ![lo(1)-2, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(2)+2]
    work_lo = (/ lo(1)-1 , lo(2)-2, lo(3)-2 /)
    work_hi = (/ hi(1)+2 , hi(2)+2, hi(3)+2 /)

    ut_lo = [lo(1)-2, lo(2)-2, lo(3)-2]
    ut_hi = [hi(1)+2, hi(2)+2, hi(3)+2]

    call bl_allocate(utmp_left, ut_lo, ut_hi, NVAR+3)
    call bl_allocate(utmp_right, ut_lo, ut_hi, NVAR+3)

    call corner_couple(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       ux_right, u_lo, u_hi, &
                       ux_left, u_lo, u_hi, &
                       flxy1D, fy_lo, fy_hi, &
                       !d1 = x, d2 = y
                       1, 2, &
                       dx(1), dt) !qmpxy

    ! magnetic field components on x interface affected by y fluxes
    call corner_couple_mag(work_lo, work_hi, &
                           utmp_right, ut_lo, ut_hi, &
                           utmp_left, ut_lo, ut_hi, &
                           ux_right, u_lo, u_hi, &
                           ux_left, u_lo, u_hi, &
                           Ex, ex_lo, ex_hi, &
                           Ez, ez_lo, ez_hi, &
                           !x,y,z
                           1, 2, 3, dx(1), dt)

    call bl_allocate(qtmp_left, ut_lo, ut_hi, NQ)
    call bl_allocate(qtmp_right, ut_lo, ut_hi, NQ)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)

    ! Calculate Flux 2D eq. 40
    ! [lo(1)-1, lo(2)-2, lo(3)-2][hi(1)+2,hi(2)+2,hi(3)+2]
    fxy_lo = (/ lo(1)-1, lo(2)-2, lo(3)-2 /)
    fxy_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call bl_allocate(flx_xy, fxy_lo, fxy_hi, NVAR+3)

    print *, "calling hlld"
    call hlld(fxy_lo, fxy_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flx_xy, fxy_lo, fxy_hi, 1) !F^{x|y}

    print *, "done"

    ! affected by Z Flux
    call corner_couple(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       ux_right, u_lo, u_hi, &
                       ux_left, u_lo, u_hi, &
                       flxz1D, fz_lo, fz_hi, &
                       1, 3, &
                       dx(1), dt) !qrpxz

    call corner_couple_mag(work_lo, work_hi, &
                           utmp_right, ut_lo, ut_hi, &
                           utmp_left, ut_lo, ut_hi, &
                           ux_right, u_lo, u_hi, &
                           ux_left, u_lo, u_hi, &
                           Ex, ex_lo, ex_hi, &
                           Ey, ey_lo, ey_hi, &
                           1, 3, 2, dx(1), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)

    fxz_lo = (/ lo(1)-1, lo(2)-2, lo(3)-2 /)
    fxz_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call bl_allocate(flx_xz, fxz_lo, fxz_hi, NVAR+3)

    call hlld(fxz_lo, fxz_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flx_xz, fxz_lo, fxz_hi, 1) !F^{x|z}


    !Y direction

    print *, "here"

    ! affected by X Flux
    ![lo(1)-2, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(3)+2]
    work_lo = (/ lo(1)-2, lo(2)-1, lo(3)-2 /)
    work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call corner_couple(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       uy_right, u_lo, u_hi, &
                       uy_left, u_lo, u_hi, &
                       flxx1D, fx_lo, fx_hi, &
                       2, 1, &
                       dx(2), dt) !qrpyx

    call corner_couple_mag(work_lo, work_hi, &
                           utmp_right, ut_lo, ut_hi, &
                           utmp_left, ut_lo, ut_hi, &
                           uy_right, u_lo, u_hi, &
                           uy_left, u_lo, u_hi, &
                           Ey, ey_lo, ey_hi, &
                           Ez, ez_lo, ez_hi, &
                           2, 1, 3, dx(2), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)


    ![lo(1)-2, lo(2)-1, lo(3)-2][hi(1)+2,hi(2)+2,hi(3)+2]
    fyx_lo = (/ lo(1)-2, lo(2)-1, lo(3)-2 /)
    fyx_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call bl_allocate(flx_yx, fyx_lo, fyx_hi, NVAR+3)

    call hlld(fyx_lo, fyx_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flx_yx, fyx_lo, fyx_hi, 2) !F^{y|x}

    print *, "z direction"

    ! affected by Z Flux
    call corner_couple(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       uy_right, u_lo, u_hi, &
                       uy_left, u_lo, u_hi, &
                       flxz1D, fz_lo, fz_hi, &
                       2, 3, &
                       dx(2), dt) !qrpyz

    call corner_couple_mag(work_lo, work_hi, &
                           utmp_right, ut_lo, ut_hi, &
                           utmp_left, ut_lo, ut_hi, &
                           uy_right, u_lo, u_hi, &
                           uy_left, u_lo, u_hi, &
                           Ey, ey_lo, ey_hi, &
                           Ex, ex_lo, ex_hi, &
                           2, 3, 1, dx(2), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)

    print *, "here2"

    fyz_lo = (/ lo(1)-2, lo(2)-1, lo(3)-2 /)
    fyz_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call bl_allocate(flx_yz, fyz_lo, fyz_hi, NVAR+3)

    call hlld(fyz_lo, fyz_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flx_yz, fyz_lo, fyz_hi, 2) !F^{y|z}

    !Z direction

    ! affected by X Flux
    ![lo(1)-2, lo(2)-2, lo(3)-2] [hi(1)+2, hi(2)+2, hi(3)+2]
    work_lo = (/ lo(1)-2, lo(2)-2, lo(3)-1 /)
    work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2/)
    call corner_couple(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       uz_right, u_lo, u_hi, &
                       uz_left, u_lo, u_hi, &
                       flxx1D, fx_lo, fx_hi, &
                       3, 1, &
                       dx(3), dt) !qrpzx

    call corner_couple_mag(work_lo, work_hi, &
                           utmp_right, ut_lo, ut_hi, &
                           utmp_left, ut_lo, ut_hi, &
                           uz_right, u_lo, u_hi, &
                           uz_left, u_lo, u_hi, &
                           Ez, ez_lo, ez_hi, &
                           Ey, ey_lo, ey_hi, &
                           3, 1, 2, dx(3), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)

    print *, "fxz"

    ![lo(1)-2,lo(2)-2,lo(3)-1][h1(1)+2, h1(2)+2, h1(3)+2]
    fzx_lo = (/ lo(1)-2, lo(2)-2, lo(3)-1 /)
    fzx_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call bl_allocate(flx_zx, fzx_lo, fzx_hi, NVAR+3)

    call hlld(fzx_lo, fzx_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flx_zx, fzx_lo, fzx_hi, 3) !F^{z|x}

    ! affected by Y Flux
    call corner_couple(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       uz_right, u_lo, u_hi, &
                       uz_left, u_lo, u_hi, &
                       flxy1D, fy_lo, fy_hi, &
                       3, 2, &
                       dx(3), dt) !qrpzy

    call corner_couple_mag(work_lo, work_hi, &
                           utmp_right, ut_lo, ut_hi, &
                           utmp_left, ut_lo, ut_hi, &
                           uz_right, u_lo, u_hi, &
                           uz_left, u_lo, u_hi, &
                           Ez, ez_lo, ez_hi, &
                           Ex, ex_lo, ex_hi, &
                           3, 2, 1, dx(3), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)

    fzy_lo = (/ lo(1)-2, lo(2)-2, lo(3)-1 /)
    fzy_hi = (/ hi(1)+2, hi(2)+2, hi(3)+2 /)

    call bl_allocate(flx_zy, fzy_lo, fzy_hi, NVAR+3)

    call hlld(fzy_lo, fzy_hi, &
         qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
         flx_zy, fzy_lo, fzy_hi, 3) !F^{z|y}


    ! MM CTU Step 6
    ! Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
    ! eq. 42 and 43

    do k = lo(3)-2, hi(3)+2
       do j = lo(2)-2, hi(2)+2
          do i = lo(1)-1, hi(1)+2
             flxx1D(i,j,k,:) = 0.5d0*(flx_xy(i,j,k,:) + flx_xz(i,j,k,:))
          end do
       end do
    end do

    do k = lo(3)-2, hi(3)+2
       do j = lo(2)-1, hi(2)+2
          do i = lo(1)-2, hi(1)+2
             flxy1D(i,j,k,:) = 0.5d0*(flx_yx(i,j,k,:) + flx_yz(i,j,k,:))
          end do
       end do
    end do

    do k = lo(3)-1, hi(3)+2
       do j = lo(2)-2, hi(2)+2
          do i = lo(1)-2, hi(1)+2
             flxz1D(i,j,k,:) = 0.5d0*(flx_zx(i,j,k,:) + flx_zy(i,j,k,:))
          end do
       end do
    end do


    print *, "after flux averaging"

    ! eq. 41
    ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+2, hi(3)+2]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+2, hi(3)+2 /)
    call electric_edge_x(work_lo, work_hi, &
                         q, q_lo, q_hi, &
                         Ex, ex_lo, ex_hi, &
                         flxy1D, fy_lo, fy_hi, &
                         flxz1D, fz_lo, fz_hi)

    ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+1, hi(3)+2]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+2, hi(2)+1, hi(3)+2 /)
    call electric_edge_y(work_lo, work_hi, &
                         q, q_lo, q_hi, &
                         Ey, ey_lo, ey_hi, &
                         flxx1D, fx_lo, fx_hi, &
                         flxz1D, fz_lo, fz_hi)

    ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+2, hi(2)+2, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+2, hi(2)+2, hi(3)+1 /)
    call electric_edge_z(work_lo, work_hi, &
                         q, q_lo, q_hi, &
                         Ez, ez_lo, ez_hi, &
                         flxx1D, fx_lo, fx_hi, &
                         flxy1D, fy_lo, fy_hi)

    ! MM CTU Step 7, 8, and 9
    ! Half Step conservative vars eq.44, eq.45, eq.46
    ! Here we reuse utmp_left/right to denote the half-time conservative state

    print *, "x half step"

    !for x direction
    ![lo(1)-1,lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call half_step(work_lo, work_hi, &
                   utmp_right, ut_lo, ut_hi, &
                   utmp_left, ut_lo, ut_hi, &
                   ux_right, u_lo, u_hi, &
                   ux_left, u_lo, u_hi, &
                   flx_yz, fyz_lo, fyz_hi, &
                   flx_zy, fzy_lo, fzy_hi, &
                   !dir = x, d1 =y, d2 =z
                   1, 2, 3, dx(1), dt)

    call half_step_mag(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       ux_right, u_lo, u_hi, &
                       ux_left, u_lo, u_hi, &
                       Ex, ex_lo, ex_hi, &
                       Ey, ey_lo, ey_hi, &
                       Ez, ez_lo, ez_hi, &
                       !d=x, d1=y, d2=z,
                       1, 2, 3, &
                       dx(1), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)
    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    ! Final Fluxes eq.47

    ! We need to compute these on a box 1 larger in the transverse directions
    ! than we'd need for hydro alone due to the electric update

    ![lo(1), lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call hlld(work_lo, work_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flxx, flxx_lo, flxx_hi, 1)


    print *, "y half step"

    !for y direction

    ![lo(1)-1,lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call half_step(work_lo, work_hi, &
                   utmp_right, ut_lo, ut_hi, &
                   utmp_left, ut_lo, ut_hi, &
                   uy_right, u_lo, u_hi, &
                   uy_left, u_lo, u_hi, &
                   flx_xz, fxz_lo, fxz_hi, &
                   flx_zx, fzx_lo, fzx_hi, &
                   2, 1, 3, dx(2), dt)

    call half_step_mag(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       uy_right, u_lo, u_hi, &
                       uy_left, u_lo, u_hi, &
                       Ey, ey_lo, ey_hi, &
                       Ex, ex_lo, ex_hi, &
                       Ez, ez_lo, ez_hi, &
                       !d, d1, d2, UMAGD UMAGD1, UMAGD2, sgn,
                       2, 1, 3, &
                       dx(2), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)
    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    ![lo(1)-1, lo(2), lo(3)-1][hi(1)+1,hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call hlld(work_lo, work_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flxy, flxy_lo, flxy_hi, 2)

    print *, "z half step"

    !for z direction

    ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call half_step(work_lo, work_hi, &
                   utmp_right, ut_lo, ut_hi, &
                   utmp_left, ut_lo, ut_hi, &
                   uz_right, u_lo, u_hi, &
                   uz_left, u_lo, u_hi, &
                   flx_xy, fxy_lo, fxy_hi, &
                   flx_yx, fyx_lo, fyx_hi, &
                   3, 1, 2, dx(3), dt)

    call half_step_mag(work_lo, work_hi, &
                       utmp_right, ut_lo, ut_hi, &
                       utmp_left, ut_lo, ut_hi, &
                       uz_right, u_lo, u_hi, &
                       uz_left, u_lo, u_hi, &
                       Ez, ez_lo, ez_hi, &
                       Ex, ex_lo, ex_hi, &
                       Ey, ey_lo, ey_hi, &
                       !d, d1, d2, UMAGD UMAGD1, UMAGD2, sgn,
                       3, 1, 2, &
                       dx(3), dt)

    call ConsToPrim(work_lo, work_hi, &
                    qtmp_right, ut_lo, ut_hi, utmp_right, ut_lo, ut_hi)
    call ConsToPrim(work_lo, work_hi, &
                    qtmp_left, ut_lo, ut_hi, utmp_left, ut_lo, ut_hi)

    ![lo(1)-1,lo(2)-1,lo(3)][hi(1)+1, hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call hlld(work_lo, work_hi, &
              qtmp_left, ut_lo, ut_hi, qtmp_right, ut_lo, ut_hi, &
              flxz, flxz_lo, flxz_hi, 3)


    ! MM CTU Step 10
    ! Primitive update eq. 48
    ![lo(1)-1, lo(2)-1, lo(3)-1][hi(1)+1, hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1)-1, lo(2)-1, lo(3)-1 /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3)+1 /)
    call prim_half(work_lo, work_hi, &
                   q2D, q_lo, q_hi, &
                   q, q_lo, q_hi, &
                   flxx1D, fx_lo, fx_hi, &
                   flxy1D, fy_lo, fy_hi, &
                   flxz1D, fz_lo, fz_hi, &
                   dx(1), dx(2), dx(3), dt)

    !Final Electric Field Update eq.48
    ![lo(1), lo(2), lo(3)][hi(1), hi(2)+1, hi(3)+1]
    work_lo = (/ lo(1), lo(2), lo(3) /)
    work_hi = (/ hi(1), hi(2)+1, hi(3)+1 /)
    call electric_edge_x(work_lo, work_hi, &
                         q2D, q_lo, q_hi, &
                         Ex, ex_lo, ex_hi, &
                         flxy, flxy_lo, flxy_hi, &
                         flxz, flxz_lo, flxz_hi)

    ![lo(1), lo(2), lo(3)][hi(1)+1, hi(2), hi(3)+1]
    work_lo = (/ lo(1), lo(2), lo(3) /)
    work_hi = (/ hi(1)+1, hi(2), hi(3)+1 /)
    call electric_edge_y(work_lo, work_hi, &
                         q2D, q_lo, q_hi, &
                         Ey, ey_lo, ey_hi, &
                         flxx, flxx_lo, flxx_hi, &
                         flxz, flxz_lo, flxz_hi)

    ![lo(1), lo(2), lo(3)][hi(1)+1, hi(2)+1 ,hi(3)]
    work_lo = (/ lo(1), lo(2), lo(3) /)
    work_hi = (/ hi(1)+1, hi(2)+1, hi(3) /)
    call electric_edge_z(work_lo, work_hi, &
                         q2D, q_lo, q_hi, &
                         Ez, ez_lo, ez_hi, &
                         flxx, flxx_lo, flxx_hi, &
                         flxy, flxy_lo, flxy_hi)


    call bl_deallocate(flxx1D)
    call bl_deallocate(flxy1D)
    call bl_deallocate(flxz1D)

    call bl_deallocate(flx_xy)
    call bl_deallocate(flx_xz)
    call bl_deallocate(flx_yx)
    call bl_deallocate(flx_yz)
    call bl_deallocate(flx_zx)
    call bl_deallocate(flx_zy)

    call bl_deallocate(ux_left)
    call bl_deallocate(ux_right)

    call bl_deallocate(uy_left)
    call bl_deallocate(uy_right)

    call bl_deallocate(uz_left)
    call bl_deallocate(uz_right)

    call bl_deallocate(utmp_left)
    call bl_deallocate(utmp_right)

    call bl_deallocate(qtmp_left)
    call bl_deallocate(qtmp_right)

  end subroutine corner_transport

  !================================================= Calculate the Conservative Variables ===============================================

  subroutine PrimToCons(lo, hi, q, q_lo, q_hi, u, u_lo, u_hi)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module
    use eos_type_module, only : eos_t, eos_input_rp
    use eos_module, only: eos
    use network, only: nspec

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)

    real(rt), intent(in) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3), NQ)
    real(rt), intent(out) :: u(u_lo(1):u_hi(1), u_lo(2):u_hi(2), u_lo(3):u_hi(3), NVAR+3)

    integer :: i ,j ,k

    type(eos_t) :: eos_state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             u(i,j,k,URHO)  = q(i,j,k,QRHO)
             u(i,j,k,UMX)    = q(i,j,k,QRHO)*q(i,j,k,QU)
             u(i,j,k,UMY)    = q(i,j,k,QRHO)*q(i,j,k,QV)
             u(i,j,k,UMZ)    = q(i,j,k,QRHO)*q(i,j,k,QW)

             eos_state % rho = q(i, j, k, QRHO)
             eos_state % p   = q(i, j, k, QPRES)
             eos_state % T   = 100.d0 !dummy initial g.
             eos_state % xn  = q(i, j, k, QFS:QFS+nspec-1)

             call eos(eos_input_rp, eos_state)

             u(i,j,k,UEDEN) = eos_state % rho * eos_state % e &
                  + 0.5d0*q(i,j,k,QRHO)*dot_product(q(i,j,k,QU:QW),q(i,j,k,QU:QW)) &
                  + 0.5d0*(dot_product(q(i,j,k,QMAGX:QMAGZ),q(i,j,k,QMAGX:QMAGZ)))

             u(i,j,k,UEINT) = eos_state % rho * eos_state % e
             u(i,j,k,UTEMP) = eos_state % T

             u(i,j,k,UMAGX:UMAGZ) = q(i,j,k,QMAGX:QMAGZ)

             ! species
             u(i,j,k,UFS:UFS-1+nspec) = q(i,j,k,QRHO) * q(i,j,k,QFS:QFS-1+nspec)

          enddo
       enddo
    enddo
  end subroutine PrimToCons

  !================================================= Calculate the Primitve Variables ===============================================

  subroutine ConsToPrim(w_lo, w_hi, q, q_lo, q_hi, u, u_lo, u_hi)

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use network, only : nspec

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: u_lo(3), u_hi(3)
    real(rt), intent(in) :: u(u_lo(1):u_hi(1), u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR+3)
    real(rt), intent(out) :: q(q_lo(1):q_hi(1), q_lo(2):q_hi(2), q_lo(3):q_hi(3),NQ)

    integer :: i ,j ,k

    type (eos_t) :: eos_state

    !q = u
    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)
             q(i,j,k,QRHO)  = u(i,j,k,URHO)
             q(i,j,k,QU)    = u(i,j,k,UMX)/q(i,j,k,QRHO)
             q(i,j,k,QV)    = u(i,j,k,UMY)/q(i,j,k,QRHO)
             q(i,j,k,QW)    = u(i,j,k,UMZ)/q(i,j,k,QRHO)
             q(i,j,k,QREINT) = u(i,j,k,UEDEN) - 0.5d0*q(i,j,k,QRHO)*dot_product(q(i,j,k,QU:QW),q(i,j,k,QU:QW)) &
                  - 0.5d0*dot_product(u(i,j,k,UMAGX:UMAGZ), u(i,j,k,UMAGX:UMAGZ)) !NVAR+1->UMAGX

             ! species
             q(i,j,k,QFS:QFS-1+nspec) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)

             eos_state % rho = q(i, j, k, QRHO)
             eos_state % e   = q(i, j, k, QREINT)/eos_state % rho
             eos_state % xn  = q(i, j, k, QFS:QFS+nspec-1)
             eos_state % T   = 100.d0  ! initial guess

             call eos(eos_input_re, eos_state)

             q(i,j,k,QTEMP) = eos_state % T
             q(i,j,k,QPRES) = eos_state % p
             q(i,j,k,QMAGX:QMAGZ) = u(i,j,k,UMAGX:UMAGZ)
          enddo
       enddo
    enddo
  end subroutine ConsToPrim

  subroutine corner_couple(w_lo, w_hi, &
                           ur_out, uro_lo, uro_hi, &
                           ul_out, ulo_lo, ulo_hi, &
                           ur, ur_lo, ur_hi, &
                           ul, ul_lo, ul_hi, &
                           flxd2, flxd2_lo, flxd2_hi, &
                           d1, d2, &
                           dx, dt)

    ! take conservative interface states ul and ur and update them
    ! with with the transverse flux difference (corner coupling) to
    ! produce ul_out and ur_out
    !
    ! This implements MM step 3 of the CTU algorithm.
    ! the normal direction (for the interface states) is d1
    ! the transverse direction (for the flux difference) is d2

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module
    use network, only : nspec

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: ulo_lo(3), ulo_hi(3)
    integer, intent(in) :: uro_lo(3), uro_hi(3)
    integer, intent(in) :: ur_lo(3), ur_hi(3)
    integer, intent(in) :: ul_lo(3), ul_hi(3)
    integer, intent(in) :: flxd2_lo(3), flxd2_hi(3)
    integer, intent(in) :: d1, d2
    real(rt), intent(in) :: dx, dt

    real(rt), intent(in) :: ur(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3)
    real(rt), intent(in) :: ul(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3)

    real(rt), intent(out) :: flxd2(flxd2_lo(1):flxd2_hi(1),flxd2_lo(2):flxd2_hi(2),flxd2_lo(3):flxd2_hi(3),NVAR+3)

    real(rt), intent(out) :: ur_out(uro_lo(1):uro_hi(1),uro_lo(2):uro_hi(2),uro_lo(3):uro_hi(3),NVAR+3)
    real(rt), intent(out) :: ul_out(ulo_lo(1):ulo_hi(1),ulo_lo(2):ulo_hi(2),ulo_lo(3):ulo_hi(3),NVAR+3)

    real(rt) :: u, v, w, cdtdx
    integer  :: i, j, k, n
    integer  :: dl(3), dr(3) !for the addition of +1 to either i,j,k depending on d2

    ! update the state on interface direction d1 with the input flux in direction d2

    ! for the flux difference, F_r - F_l, we need to shift the indices in the first flux (F_r)
    ! in d2 to get a difference across the interface.  We also need to shift by a zone in d1
    ! for the left interface.  dr(:) and dl(:) will hold these shifts.
    dr(:) = 0
    dl(:) = 0

    ! the first term of the flxd2 substraction is shifted by 1 on the direction d2
    dr(d2) = 1

    cdtdx = dt/(3.d0*dx)

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
             ! MM Eq. 37

             do n = 1, NVAR
                if (n == UTEMP) cycle

                ur_out(i,j,k,n) = ur(i,j,k,n) - &
                     cdtdx*(flxd2(i+dr(1),j+dr(2),k+dr(3),n) - flxd2(i,j,k,n))
             end do

             ! now fix up the internal energy

             u = ur_out(i,j,k,UMX)/ur_out(i,j,k,URHO)
             v = ur_out(i,j,k,UMY)/ur_out(i,j,k,URHO)
             w = ur_out(i,j,k,UMZ)/ur_out(i,j,k,URHO)
             ur_out(i,j,k,UEINT) = ur_out(i,j,k,UEDEN) - 0.5d0*ur_out(i,j,k,URHO)*(u**2 + v**2 + w**2)

          end do
       end do
    end do


    ! left interface (e.g., U_{i-1/2,j,k,L} or the "+" state in MM notation)
    ! note: this uses information one zone to the left in d1
    dl(d1) = -1
    dr(d1) = -1

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             do n = 1, NVAR
                if (n == UTEMP) cycle

                ul_out(i,j,k,n) = ul(i,j,k,n) - &
                     cdtdx*(flxd2(i+dr(1),j+dr(2),k+dr(3),n) - &
                            flxd2(i+dl(1),j+dl(2),k+dl(3),n))
             end do

             ! fix up the internal energy

             u = ul_out(i,j,k,UMX)/ul_out(i,j,k,URHO)
             v = ul_out(i,j,k,UMY)/ul_out(i,j,k,URHO)
             w = ul_out(i,j,k,UMZ)/ul_out(i,j,k,URHO)
             ul_out(i,j,k,UEINT) = ul_out(i,j,k,UEDEN) - 0.5d0*ul_out(i,j,k,URHO)*(u**2 + v**2 + w**2)

          enddo
       enddo
    enddo
  end subroutine corner_couple


  subroutine corner_couple_mag(w_lo, w_hi, &
                               ur_out, uro_lo, uro_hi, &
                               ul_out, ulo_lo, ulo_hi, &
                               ur, ur_lo, ur_hi, &
                               ul, ul_lo, ul_hi, &
                               Ed1, ed1_lo, ed1_hi, &
                               Ed3, ed3_lo, ed3_hi, &
                               d1, d2, d3, &
                               dx, dt)
    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, UEINT

    ! Use 1D Electric Fields to Transverse correct the Temporary Magnetic Fields
    ! Correction using Faraday's Law

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: uro_lo(3), uro_hi(3)
    integer, intent(in) :: ulo_lo(3), ulo_hi(3)
    integer, intent(in) :: ur_lo(3), ur_hi(3)
    integer, intent(in) :: ul_lo(3), ul_hi(3)
    integer, intent(in) :: ed1_lo(3), ed1_hi(3)
    integer, intent(in) :: ed3_lo(3), ed3_hi(3)
    integer, intent(in) :: d1, d2, d3

    real(rt), intent(inout) :: ur_out(uro_lo(1):uro_hi(1), uro_lo(2):uro_hi(2), uro_lo(3):uro_hi(3), NVAR+3)
    real(rt), intent(inout) :: ul_out(ulo_lo(1):ulo_hi(1), ulo_lo(2):ulo_hi(2), ulo_lo(3):ulo_hi(3), NVAR+3)
    real(rt), intent(in)    :: ur(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3)
    real(rt), intent(in)    :: ul(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3)

    real(rt), intent(in) :: Ed1(ed1_lo(1):ed1_hi(1),ed1_lo(2):ed1_hi(2),ed1_lo(3):ed1_hi(3))
    real(rt), intent(in) :: Ed3(ed3_lo(1):ed3_hi(1),ed3_lo(2):ed3_hi(2),ed3_lo(3):ed3_hi(3))

    real(rt) :: dx, dt
    integer :: i ,j ,k
    integer :: d(3), a1(3), a2(3), a3(3), a4(3) !for the additions of +1 to i,j,k

    integer :: UMAGD1, UMAGD2, UMAGD3   !UMAGD1 corresponds to d1, and UMAGD2 to d2, UMAGD3 to d3
    integer :: sgn
    real(rt) :: cdtdx

    ! update the magnetic field components on the d1 face with the
    ! flux in the d2 direction.
    !
    ! There are 3 components here:
    !
    !  * B_d1 is perpendicular to the face and it is updated by the
    !    transverse flux difference d2 on the face.
    !  * B_d2 is unchanged, since it points in the direction d2
    !  * B_d3 is in the plane of the face, but perpendicular to d2,
    !    and is updated by a d2 flux difference inside the cell

    sgn = epsilon_ijk(d1, d2, d3)
    cdtdx = dt/(3.d0*dx)

    UMAGD1 = UMAGX - 1 + d1
    UMAGD2 = UMAGX - 1 + d2
    UMAGD3 = UMAGX - 1 + d3

    d(:) = 0
    a1(:) = 0
    a2(:) = 0
    a3(:) = 0
    a4(:) = 0

    d(d2)  = 1   !for example if d2 = y, j+1 on Ed3

    a1(d2) = 1   !j+1 on first and third term of addition Ed1
    a1(d3) = 1

    a2(d3) = 1   !the second term of addition Ed1 increments by 1 the i,j,k

    a3(d2) = 1

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

             ! d1 -- this is perpendicular to the face, MM Eq. 38
             ! (note MM Eq. 38 has a sign error) e.g., for d1 = x and
             ! d2 = y, this gets updated as
             !
             ! Bx|y_{i-1/2,j,k,R} = Bx_{i-1/2,j,k) -
             !     dt/3dx (Ez_{i-1/2,j+1/2,k) - Ez_{i-1/2,j-1/2,k})
             !
             ! we use d(:) to captured the j+1/2 indexing into Ez

             ur_out(i,j,k,UMAGD1) = ur(i,j,k,UMAGD1) - sgn * cdtdx * &
                  (Ed3(i+d(1),j+d(2),k+d(3)) - Ed3(i,j,k))

             ! d3 -- this is in the plane of the face, MM Eq. 39
             ! e.g.,g for d1 = x, and d2 = y, this gets updated as
             !
             ! Bz|y_{i-1/2,j,k,R} = Bz_{i-1/2,j,k} + 1/2 dt/3dx
             !     (Ex_{i,j+1/2,k+1/2} - Ex_{i,j-1/2,k+1/2} +
             !      Ex_{i,j+1/2,k-1/2} - Ex_{i,j-1/2,k-1/2})
             !
             ! we use a1(:) for the first E term, a2(:) for the
             ! second, and a3(:) for the third

             ur_out(i,j,k,UMAGD3) = ur(i,j,k,UMAGD3) + sgn * 0.5_rt * cdtdx * &
                  ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i+a2(1),j+a2(2),k+a2(3))) + &
                   (Ed1(i+a3(1),j+a3(2),k+a3(3)) - Ed1(i,j,k)))

             ! the component pointing in the transverse update direction, d2, is unchanged
             ur_out(i,j,k,UMAGD2) = ur(i,j,k,UMAGD2)

             ur_out(i,j,k,UEINT) = ur_out(i,j,k,UEINT) - &
                  0.5d0*dot_product(ur_out(i,j,k,UMAGX:UMAGZ), ur_out(i,j,k,UMAGX:UMAGZ))

          enddo
       enddo
    enddo

    ! The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
    ! in direction d1
    a1(d1) = -1
    a2(d1) = -1
    a3(d1) = -1
    a4(d1) = -1

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! left state on the interface (e.g. B_{i-1/2,j,k,L} or `+` in MM notation)

             ul_out(i,j,k,UMAGD1) = ul(i,j,k,UMAGD1) - sgn * cdtdx * &
                  (Ed3(i+d(1),j+d(2),k+d(3)) - Ed3(i,j,k))

             ul_out(i,j,k,UMAGD3) = ul(i,j,k,UMAGD3) + sgn * 0.5_rt * cdtdx * &
                  ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i+a2(1),j+a2(2),k+a2(3))) + &
                   (Ed1(i+a3(1),j+a3(2),k+a3(3)) - Ed1(i+a4(1),j+a4(2),k+a4(3))))

             ul_out(i,j,k,UMAGD2) = ul(i,j,k,UMAGD2)

             ul_out(i,j,k,UEINT) = ul_out(i,j,k,UEINT) - &
                  0.5d0*dot_product(ul_out(i,j,k,UMAGX:UMAGZ), ul_out(i,j,k,UMAGX:UMAGZ))

          enddo
       enddo
    enddo

  end subroutine corner_couple_mag


  subroutine half_step(w_lo, w_hi, &
                       ur_out, uro_lo, uro_hi, &
                       ul_out, ulo_lo, ulo_hi, &
                       ur, ur_lo, ur_hi, &
                       ul, ul_lo, ul_hi, &
                       flxd1, flxd1_lo, flxd1_hi, &
                       flxd2, flxd2_lo, flxd2_hi, &
                       dir, d1, d2, dx, dt)

    ! Final transverse flux corrections to the conservative state

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, URHO, UEDEN, UMX, UMY, UMZ, URHO, UEINT, UFS
    use network, only: nspec

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: uro_lo(3), uro_hi(3)
    integer, intent(in) :: ulo_lo(3), ulo_hi(3)
    integer, intent(in) :: ur_lo(3), ur_hi(3)
    integer, intent(in) :: ul_lo(3), ul_hi(3)

    integer, intent(in) :: flxd1_lo(3), flxd1_hi(3)
    integer, intent(in) :: flxd2_lo(3), flxd2_hi(3)

    real(rt), intent(inout) :: ur_out(uro_lo(1):uro_hi(1), uro_lo(2):uro_hi(2), uro_lo(3):uro_hi(3), NVAR+3)
    real(rt), intent(inout) :: ul_out(ulo_lo(1):ulo_hi(1), ulo_lo(2):ulo_hi(2), ulo_lo(3):ulo_hi(3), NVAR+3)
    real(rt), intent(in)    :: ur(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3)
    real(rt), intent(in)    :: ul(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3)

    real(rt), intent(in)  :: dx, dt !dx will be dx, dy or dz
    integer, intent(in)   :: dir, d1, d2 ! following notation of eq. 44

    real(rt), intent(in)  :: flxd1(flxd1_lo(1):flxd1_hi(1),flxd1_lo(2):flxd1_hi(2),flxd1_lo(3):flxd1_hi(3),NVAR+3)
    real(rt), intent(in)  :: flxd2(flxd2_lo(1):flxd2_hi(1),flxd2_lo(2):flxd2_hi(2),flxd2_lo(3):flxd2_hi(3),NVAR+3)

    real(rt) :: u, v, w
    integer  :: i ,j ,k, n
    integer  :: dl_1(3), dr_1(3), dl_2(3), dr_2(3) !for the shift in i,j,k
    real(rt) :: hdtdx

    hdtdx = 0.5_rt * dt/dx

    dl_1(:) = 0
    dr_1(:) = 0
    dl_2(:) = 0
    dr_2(:) = 0

    dr_1(d1) = 1  ! add +1 to the d1 direction in the first flxd1 term of the subtraction
    dr_2(d2) = 1  ! add +1 to the d2 direction in the first flxd2 term of the subtraction


    ! right interface (e.g. U_{i-1/2,j,k,R} or the "-" state in MM notation)
    ! MM Eq. 44

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             do n = 1, NVAR
                if (n == UTEMP) cycle

                ur_out(i,j,k,n) = ur(i,j,k,n) - &
                     hdtdx * (flxd1(i+dr_1(1),j+dr_1(2),k+dr_1(3),n) - flxd1(i,j,k,n)) - &
                     hdtdx * (flxd2(i+dr_2(1),j+dr_2(2),k+dr_2(3),n) - flxd2(i,j,k,n))

             end do

             u = ur_out(i,j,k,UMX)/ur_out(i,j,k,URHO)
             v = ur_out(i,j,k,UMY)/ur_out(i,j,k,URHO)
             w = ur_out(i,j,k,UMZ)/ur_out(i,j,k,URHO)

             ur_out(i,j,k,UEINT) = ur_out(i,j,k,UEDEN) - &
                  0.5d0*ur_out(i,j,k,URHO)*(u**2 + v**2 + w**2)

          end do
       end do
    end do

    ! for the left state B components on the face dir, the flux
    ! difference is in the zone to the left
    dr_1(dir) = -1
    dl_1(dir) = -1
    dr_2(dir) = -1
    dl_2(dir) = -1

    ! left interface (e.g., U_{i+1/2,j,k,L} or the "+" state in MM notation)
    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             do n = 1, NVAR
                if (n == UTEMP) cycle

                ul_out(i,j,k,n) = ul(i,j,k,n) - &
                     hdtdx * (flxd1(i+dr_1(1),j+dr_1(2),k+dr_1(3),n) - &
                              flxd1(i+dl_1(1),j+dl_1(2),k+dl_1(3),n)) - &
                     hdtdx * (flxd2(i+dr_2(1),j+dr_2(2),k+dr_2(3),n) - &
                              flxd2(i+dl_2(1),j+dl_2(2),k+dl_2(3),n))

             end do

             u = ul_out(i,j,k,UMX)/ul_out(i,j,k,URHO)
             v = ul_out(i,j,k,UMY)/ul_out(i,j,k,URHO)
             w = ul_out(i,j,k,UMZ)/ul_out(i,j,k,URHO)
             ul_out(i,j,k,UEINT) = ul_out(i,j,k,UEDEN) - &
                  0.5d0*ul_out(i,j,k,URHO)*(u**2 + v**2 + w**2)

          end do
       end do
    end do

  end subroutine


  subroutine half_step_mag(w_lo, w_hi, &
                           ur_out, uro_lo, uro_hi, &
                           ul_out, ulo_lo, ulo_hi, &
                           ur, ur_lo, ur_hi, &
                           ul, ul_lo, ul_hi, &
                           Ed, ed_lo, ed_hi, &
                           Ed1, ed1_lo, ed1_hi, &
                           Ed2, ed2_lo, ed2_hi, &
                           d, d1, d2, &
                           dx, dt)

    ! Final Magnetic Corrections

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR, UEINT

    !Correction using Faraday's Law
    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: uro_lo(3), uro_hi(3)
    integer, intent(in) :: ulo_lo(3), ulo_hi(3)
    integer, intent(in) :: ur_lo(3), ur_hi(3)
    integer, intent(in) :: ul_lo(3), ul_hi(3)
    integer, intent(in) :: ed_lo(3), ed_hi(3)
    integer, intent(in) :: ed1_lo(3), ed1_hi(3)
    integer, intent(in) :: ed2_lo(3), ed2_hi(3)

    integer, intent(in)   :: d, d1, d2

    real(rt), intent(inout) :: ur_out(uro_lo(1):uro_hi(1), uro_lo(2):uro_hi(2), uro_lo(3):uro_hi(3), NVAR+3)
    real(rt), intent(inout) :: ul_out(ulo_lo(1):ulo_hi(1), ulo_lo(2):ulo_hi(2), ulo_lo(3):ulo_hi(3), NVAR+3)
    real(rt), intent(in)    :: ur(ur_lo(1):ur_hi(1), ur_lo(2):ur_hi(2), ur_lo(3):ur_hi(3), NVAR+3)
    real(rt), intent(in)    :: ul(ul_lo(1):ul_hi(1), ul_lo(2):ul_hi(2), ul_lo(3):ul_hi(3), NVAR+3)

    real(rt), intent(in) :: Ed(ed_lo(1):ed_hi(1),ed_lo(2):ed_hi(2),ed_lo(3):ed_hi(3))
    real(rt), intent(in) :: Ed1(ed1_lo(1):ed1_hi(1),ed1_lo(2):ed1_hi(2),ed1_lo(3):ed1_hi(3))
    real(rt), intent(in) :: Ed2(ed2_lo(1):ed2_hi(1),ed2_lo(2):ed2_hi(2),ed2_lo(3):ed2_hi(3))

    real(rt), intent(in)  :: dx, dt

    integer :: i ,j ,k

    !to manage the +1 shifts on  i,j,k
    integer :: a1(3), a2(3), b1(3), b2(3), b3(3), b4(3), b5(3), b6(3), b7(3)

    integer :: UMAGD, UMAGD1, UMAGD2
    integer :: sgn

    real(rt) :: hdtdx

    hdtdx = 0.5_rt * dt/dx

    sgn = -1 * epsilon_ijk(d, d1, d2)

    UMAGD = UMAGX - 1 + d
    UMAGD1 = UMAGX - 1 + d1
    UMAGD2 = UMAGX - 1 + d2

    a1(:) = 0
    a2(:) = 0
    b1(:) = 0
    b2(:) = 0
    b3(:) = 0
    b4(:) = 0
    b5(:) = 0
    b6(:) = 0
    b7(:) = 0

    ! for Bd
    a1(d2) = 1 ! shift on first term of Ed1 substraction, in d2 direction
    a2(d1) = 1 ! shift on first term of Ed2 substraction, in d1 direction

    ! for Bd1 and Bd2
    b1(d1) = 1 ! shift on 1st term for Bd1 and Bd2
    b1(d2) = 1 ! in d1 and d2 components

    b2(d1) = 1 !shift on 2nd and 6th term for Bd1, and 3rd term for Bd2

    b3(d2) = 1 !shift on 3rd term for Bd1, and 2nd and 6th term for Bd2

    b4(d)  = 1 !shift on 5th term for Bd1, on d and d1 components
    b4(d1) = 1

    b5(d)  = 1 !shift on 7th term for Bd1 and Bd2

    b6(d)  = 1 !shift on 5th term for Bd2, on d and d2 components
    b6(d2) = 1


    ! d-faces

    ! right state on the interface (e.g. B_{i-1/2,j,k,R} or `-` in MM notation)

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! Bd -- this is perpendicular to the face. Note MM eq.45
             ! in Miniati has a sign error in the epsilon term

             ur_out(i,j,k,UMAGD) = ur(i,j,k,UMAGD) - sgn * hdtdx * &
                  ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i,j,k)) - &
                   (Ed2(i+a2(1),j+a2(2),k+a2(3)) - Ed2(i,j,k)))

             ! Bd1 -- this is one of the components of B in the plane of the face d
             ! Eq.46 in Miniati

             ur_out(i,j,k,UMAGD1) = ur(i,j,k,UMAGD1) + sgn * hdtdx * &
                  ((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b2(1),j+b2(2),k+b2(3))) + &
                   (Ed(i+b3(1),j+b3(2),k+b3(3)) - Ed(i,j,k)) - &
                   (Ed2(i+b4(1),j+b4(2),k+b4(3)) - Ed2(i+b2(1),j+b2(2),k+b2(3))) - &
                   (Ed2(i+b5(1),j+b5(2),k+b5(3)) - Ed2(i,j,k)))

             ! Bd2 -- this is the other component of B in the plane of the face d
             ! Eq. 46 in Miniati

             ur_out(i,j,k,UMAGD2) = ur(i,j,k,UMAGD2) - sgn * hdtdx * &
                  ((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b3(1),j+b3(2),k+b3(3))) + &
                   (Ed(i+b2(1),j+b2(2),k+b2(3)) - Ed(i,j,k)) - &
                   (Ed1(i+b6(1),j+b6(2),k+b6(3)) - Ed1(i+b3(1),j+b3(2),k+b3(3))) - &
                   (Ed1(i+b5(1),j+b5(2),k+b5(3)) - Ed1(i,j,k)))

             ! correct the energy

             ur_out(i,j,k,UEINT) = ur_out(i,j,k,UEINT) - &
                  0.5d0*dot_product(ur_out(i,j,k,UMAGX:UMAGZ), ur_out(i,j,k,UMAGX:UMAGZ))


          end do
       end do
    end do

    ! left state on the interface (e.g., B_{i-1/2,j,k,L} or `+` in MM notation)

    ! The in-plane B component at B_{i-1/2,j,k,L} uses the information one zone to the left
    ! in direction d1

    b1(d) = b1(d) - 1
    b2(d) = b2(d) - 1
    b3(d) = b3(d) - 1
    b4(d) = b4(d) - 1
    b5(d) = b5(d) - 1
    b6(d) = b6(d) - 1
    b7(d) = b7(d) - 1

    do k = w_lo(3), w_hi(3)
       do j = w_lo(2), w_hi(2)
          do i = w_lo(1), w_hi(1)

             ! for the + case, the shifts mentioned above in b6, b5, and b4
             ! also correspond to the 1st, 2nd and 4th, and 3rd term respectevely

             ! Bd -- this is perpendicular to the face (MM Eq. 45 with sign fix)

             ! this is the same face as the right state, so the update the identical
             ul_out(i,j,k,UMAGD) = ul(i,j,k,UMAGD) - sgn * hdtdx * &
                  ((Ed1(i+a1(1),j+a1(2),k+a1(3)) - Ed1(i,j,k)) - &
                   (Ed2(i+a2(1),j+a2(2),k+a2(3)) - Ed2(i,j,k)))

             ! Bd1 -- first component on face d, eq. 46 in Miniati
             ul_out(i,j,k,UMAGD1) = ul(i,j,k,UMAGD1) + sgn * hdtdx * &
                  ((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b2(1),j+b2(2),k+b2(3))) + &
                   (Ed(i+b3(1),j+b3(2),k+b3(3)) - Ed(i+b7(1),j+b7(2),k+b7(3))) - &
                   (Ed2(i+b4(1),j+b4(2),k+b4(3)) - Ed2(i+b2(1),j+b2(2),k+b2(3))) - &
                  ( Ed2(i+b5(1),j+b5(2),k+b5(3)) - Ed2(i+b7(1),j+b7(2),k+b7(3))))

             ! Bd2 -- second component on face d, eq. 46 in Miniati
             ul_out(i,j,k,UMAGD2) = ul(i,j,k,UMAGD2) - sgn * hdtdx * &
                  ((Ed(i+b1(1),j+b1(2),k+b1(3)) - Ed(i+b3(1),j+b3(2),k+b3(3))) + &
                   (Ed(i+b2(1),j+b2(2),k+b2(3)) - Ed(i+b7(1),j+b7(2),k+b7(3))) - &
                   (Ed1(i+b6(1),j+b6(2),k+b6(3)) - Ed1(i+b3(1),j+b3(2),k+b3(3))) - &
                   (Ed1(i+b5(1),j+b5(2),k+b5(3)) - Ed1(i+b7(1),j+b7(2),k+b7(3))))

             ul_out(i,j,k,UEINT) = ul_out(i,j,k,UEINT) - &
                  0.5d0*dot_product(ul_out(i,j,k,UMAGX:UMAGZ), ul_out(i,j,k,UMAGX:UMAGZ))

          end do
       end do
    end do

  end subroutine half_step_mag

  subroutine prim_half(w_lo, w_hi, &
                       q2D, q2_lo, q2_hi, &
                       q, q_lo, q_hi, &
                       flxx, flxx_lo, flxx_hi, &
                       flxy, flxy_lo, flxy_hi, &
                       flxz, flxz_lo, flxz_hi, &
                       dx, dy, dz, dt)

    ! Find the 2D corrected primitive variables

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module, only : NVAR

    implicit none

    integer, intent(in) :: w_lo(3), w_hi(3)
    integer, intent(in) :: q_lo(3), q_hi(3)
    integer, intent(in) :: q2_lo(3), q2_hi(3)
    integer, intent(in) :: flxx_lo(3), flxx_hi(3)
    integer, intent(in) :: flxy_lo(3), flxy_hi(3)
    integer, intent(in) :: flxz_lo(3), flxz_hi(3)

    real(rt), intent(in)  :: q(q_lo(1):q_hi(1),q_lo(2):q_hi(2),q_lo(3):q_hi(3),NQ)
    real(rt), intent(in)  :: flxx(flxx_lo(1):flxx_hi(1),flxx_lo(2):flxx_hi(2),flxx_lo(3):flxx_hi(3),NVAR+3)
    real(rt), intent(in)  :: flxy(flxy_lo(1):flxy_hi(1),flxy_lo(2):flxy_hi(2),flxy_lo(3):flxy_hi(3),NVAR+3)
    real(rt), intent(in)  :: flxz(flxz_lo(1):flxz_hi(1),flxz_lo(2):flxz_hi(2),flxz_lo(3):flxz_hi(3),NVAR+3)

    real(rt), intent(out) :: q2D(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NQ)

    real(rt)  :: divF(NVAR+3)
    real(rt)  :: divF_q(NQ)
    real(rt)  :: dx, dy, dz, dt
    integer   :: i, j, k

    do k = w_lo(3),w_hi(3)
       do j = w_lo(2),w_hi(2)
          do i = w_lo(1),w_hi(1)

             divF(:) = (flxx(i+1,j,k,:) - flxx(i,j,k,:)) / dx + &
                       (flxy(i,j+1,k,:) - flxy(i,j,k,:)) / dy + &
                       (flxz(i,j,k+1,:) - flxz(i,j,k,:)) / dz

             ! that is a flux of conserved variables -- transform it to primitive
             call qflux(divF_q, divF, q(i,j,k,:))

             ! Right below eq. 48
             q2D(i,j,k,:) = q(i,j,k,:) - 0.5d0*dt * divF_q

          enddo
       enddo
    enddo
  end subroutine prim_half


  subroutine qflux(qflx,flx,q)

    ! Calculate the C to P Jacobian applied to the fluxes

    use amrex_fort_module, only : rt => amrex_real
    use meth_params_module !,only : QRHO, QU, QV, QW, QPRES, QMAGX, QMAGY, QMAGZ, QVAR, NVAR
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    ! this is step 10 in the paper, just after Eq. 48

    ! this implements dW/dU . qflux, where dW/dU is the Jacobian of
    ! the primitive quantities (W) with respect to conserved quantities (U)

    real(rt), intent(in) :: flx(NVAR+3), q(NQ)
    real(rt), intent(out) :: qflx(NQ)
    real(rt) :: dedp, dedrho, totalE

    type (eos_t) :: eos_state

    qflx = 0.d0
    qflx(QRHO) = flx(URHO)
    qflx(QU) = ( flx(UMX) - flx(URHO) * q(QU) )/q(QRHO)
    qflx(QV) = ( flx(UMY) - flx(URHO) * q(QV) )/q(QRHO)
    qflx(QW) = ( flx(UMZ) - flx(URHO) * q(QW) )/q(QRHO)

    qflx(QFS:QFS+nspec-1) = ( flx(UFS:UFS+nspec-1) - flx(URHO) * q(QFS:QFS+nspec-1) )/q(QRHO)

    eos_state % rho = q(QRHO)
    eos_state % p   = q(QPRES)
    eos_state % T   = 100.d0 !dummy initial guess
    eos_state % xn  = q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    dedrho = eos_state % dedr - eos_state % dedT * eos_state % dPdr * 1.0d0/eos_state % dPdT
    dedp = eos_state % dedT * 1.0d0/eos_state % dPdT

    qflx(QPRES) = ( -q(QMAGX)*flx(UMAGX) - q(QMAGY)*flx(UMAGY) - q(QMAGZ)*flx(UMAGZ) + &
         flx(UEDEN) - flx(UMX)*q(QU) - flx(UMY)*q(QV) - &
         flx(UMZ)*q(QW) + flx(URHO)*(0.5*(q(QU)**2+q(QV)**2+q(QW)**2) - &
         eos_state % e -q(QRHO)*dedrho) ) / ( dedp * q(QRHO))

    qflx(QMAGX) = flx(UMAGX)
    qflx(QMAGY) = flx(UMAGY)
    qflx(QMAGZ) = flx(UMAGZ)

    qflx(QTEMP) = 0.0_rt

  end subroutine qflux

end module ct_upwind
