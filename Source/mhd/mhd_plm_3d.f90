module mhd_plm_module

  ! Module that gives a piecewise linear interpolation for the primitive variables
  ! They are projected onto the characteristic variables for tracing.

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module
  implicit none

  private centerdif,vanleer, lvecx, lvecy, lvecz, rvecx, rvecy, rvecz, evals, slope

  public plm

  integer, parameter :: NEIGN = 7

  integer, parameter :: IEIGN_RHO = 1
  integer, parameter :: IEIGN_U = 2
  integer, parameter :: IEIGN_V = 3
  integer, parameter :: IEIGN_W = 4
  integer, parameter :: IEIGN_P = 5

  ! perpendicular magnetic field components
  integer, parameter :: IEIGN_BT = 6
  integer, parameter :: IEIGN_BTT = 7

contains

  subroutine plm(lo, hi, &
                 idir, &
                 s, s_lo, s_hi, &
                 flatn, f_lo, f_hi, &
                 bx, bxlo, bxhi, &
                 by, bylo, byhi, &
                 bz, bzlo, bzhi, &
                 Ip, Ip_lo, Ip_hi, &
                 Im, Im_lo, Im_hi, &
                 srcQ, srcq_lo, srcq_hi, &
                 dx, dt) bind(C, name="plm")

    use network, only: nspec
    use eos_module
    use eos_type_module, only: eos_t, eos_input_rp
    use meth_params_module, only: small_dens, small_pres

    implicit none

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in), value :: idir
    integer, intent(in) :: s_lo(3), s_hi(3)
    integer, intent(in) :: f_lo(3), f_hi(3)
    integer, intent(in) :: bxlo(3), bxhi(3)
    integer, intent(in) :: bylo(3), byhi(3)
    integer, intent(in) :: bzlo(3), bzhi(3)
    integer, intent(in) :: Ip_lo(3), Ip_hi(3)
    integer, intent(in) :: Im_lo(3), Im_hi(3)
    integer, intent(in) :: srcq_lo(3), srcq_hi(3)

    real(rt), intent(in) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)

   ! face-centered magnetic fields
    real(rt), intent(in   ) ::  bx(bxlo(1):bxhi(1), bxlo(2):bxhi(2), bxlo(3):bxhi(3))
    real(rt), intent(in   ) ::  by(bylo(1):byhi(1), bylo(2):byhi(2), bylo(3):byhi(3))
    real(rt), intent(in   ) ::  bz(bzlo(1):bzhi(1), bzlo(2):bzhi(2), bzlo(3):bzhi(3))

    real(rt), intent(in) :: srcQ(srcq_lo(1):srcq_hi(1), srcq_lo(2):srcq_hi(2), srcq_lo(3):srcq_hi(3), NQSRC)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))

    real(rt), intent(out) :: Ip(Ip_lo(1):Ip_hi(1), Ip_lo(2):Ip_hi(2), Ip_lo(3):Ip_hi(3), NQ, 3)
    real(rt), intent(out) :: Im(Ip_lo(1):Ip_hi(1), Ip_lo(2):Ip_hi(2), Ip_lo(3):Ip_hi(3), NQ, 3)

    real(rt), intent(in   ) :: dx(3)
    real(rt), intent(in), value :: dt

    real(rt) :: dQL(NEIGN), dQR(NEIGN), dW, dL, dR
    real(rt) :: leig(NEIGN,NEIGN), reig(NEIGN,NEIGN), lam(NEIGN), summ_p(NEIGN), summ_m(NEIGN)
    real(rt) :: smhd(NEIGN)

    integer  :: ii, i , j, k, n
    real(rt) :: un, dtdx

    type(eos_t) :: eos_state

    ! Ip and Im are the interface states in each dimension (the last index '3' is the
    ! direction

    Ip(Ip_lo(1):Ip_hi(1),Ip_lo(2):Ip_hi(2),Ip_lo(3):Ip_hi(3),:,:) = 0.d0
    Im(Im_lo(1):Im_hi(1),Im_lo(2):Im_hi(2),Im_lo(3):Im_hi(3),:,:) = 0.d0

    ! these loops are over cell-centers and for each cell-center, we find the left and
    ! right interface states

    dtdx = dt/dx(idir)

    !=========================== PLM =========================================
    do k = s_lo(3)+1, s_hi(3)-1
       do j = s_lo(2)+1, s_hi(2)-1
          do i = s_lo(1)+1, s_hi(1)-1

             !=========================== X Direction ========================
             summ_p = 0.d0
             summ_m = 0.d0

             dW = 0.d0
             reig = 0.d0
             leig = 0.d0
             lam = 0.d0


             if (idir == 1) then
                ! we use a reduced eigensystem, the normal B field
                ! component (Bx) is omitted

                dQL(IEIGN_RHO) = s(i,j,k,QRHO) - s(i-1,j,k,QRHO)
                dQL(IEIGN_U) = s(i,j,k,QU) - s(i-1,j,k,QU)
                dQL(IEIGN_V) = s(i,j,k,QV) - s(i-1,j,k,QV)
                dQL(IEIGN_W) = s(i,j,k,QW) - s(i-1,j,k,QW)
                dQL(IEIGN_P) = s(i,j,k,QPRES) - s(i-1,j,k,QPRES)
                dQL(IEIGN_BT) = s(i,j,k,QMAGY) - s(i-1,j,k,QMAGY)
                dQL(IEIGN_BTT) = s(i,j,k,QMAGZ) - s(i-1,j,k,QMAGZ)

                dQR(IEIGN_RHO) = s(i+1,j,k,QRHO) - s(i,j,k,QRHO)
                dQR(IEIGN_U) = s(i+1,j,k,QU) - s(i,j,k,QU)
                dQR(IEIGN_V) = s(i+1,j,k,QV) - s(i,j,k,QV)
                dQR(IEIGN_W) = s(i+1,j,k,QW) - s(i,j,k,QW)
                dQR(IEIGN_P) = s(i+1,j,k,QPRES) - s(i,j,k,QPRES)
                dQR(IEIGN_BT) = s(i+1,j,k,QMAGY) - s(i,j,k,QMAGY)
                dQR(IEIGN_BTT) = s(i+1,j,k,QMAGZ) - s(i,j,k,QMAGZ)

             else if (idir == 2) then
                ! we use a reduced eigensystem, the normal B field
                ! component (By) is omitted

                dQL(IEIGN_RHO) = s(i,j,k,QRHO) - s(i,j-1,k,QRHO)
                dQL(IEIGN_U) = s(i,j,k,QU) - s(i,j-1,k,QU)
                dQL(IEIGN_V) = s(i,j,k,QV) - s(i,j-1,k,QV)
                dQL(IEIGN_W) = s(i,j,k,QW) - s(i,j-1,k,QW)
                dQL(IEIGN_P) = s(i,j,k,QPRES) - s(i,j-1,k,QPRES)
                dQL(IEIGN_BT) = s(i,j,k,QMAGX) - s(i,j-1,k,QMAGX)
                dQL(IEIGN_BTT) = s(i,j,k,QMAGZ) - s(i,j-1,k,QMAGZ)

                dQR(IEIGN_RHO) = s(i,j+1,k,QRHO) - s(i,j,k,QRHO)
                dQR(IEIGN_U) = s(i,j+1,k,QU) - s(i,j,k,QU)
                dQR(IEIGN_V) = s(i,j+1,k,QV) - s(i,j,k,QV)
                dQR(IEIGN_W) = s(i,j+1,k,QW) - s(i,j,k,QW)
                dQR(IEIGN_P) = s(i,j+1,k,QPRES) - s(i,j,k,QPRES)
                dQR(IEIGN_BT) = s(i,j+1,k,QMAGX) - s(i,j,k,QMAGX)
                dQR(IEIGN_BTT) = s(i,j+1,k,QMAGZ) - s(i,j,k,QMAGZ)

             else
                ! we use a reduced eigensystem, the normal B field
                ! component (Bz) is omitted

                dQL(IEIGN_RHO) = s(i,j,k,QRHO) - s(i,j,k-1,QRHO)
                dQL(IEIGN_U) = s(i,j,k,QU) - s(i,j,k-1,QU)
                dQL(IEIGN_V) = s(i,j,k,QV) - s(i,j,k-1,QV)
                dQL(IEIGN_W) = s(i,j,k,QW) - s(i,j,k-1,QW)
                dQL(IEIGN_P) = s(i,j,k,QPRES) - s(i,j,k-1,QPRES)
                dQL(IEIGN_BT) = s(i,j,k,QMAGX) - s(i,j,k-1,QMAGX)
                dQL(IEIGN_BTT) = s(i,j,k,QMAGY) - s(i,j,k-1,QMAGY)

                dQR(IEIGN_RHO) = s(i,j,k+1,QRHO) - s(i,j,k,QRHO)
                dQR(IEIGN_U) = s(i,j,k+1,QU) - s(i,j,k,QU)
                dQR(IEIGN_V) = s(i,j,k+1,QV) - s(i,j,k,QV)
                dQR(IEIGN_W) = s(i,j,k+1,QW) - s(i,j,k,QW)
                dQR(IEIGN_P) = s(i,j,k+1,QPRES) - s(i,j,k,QPRES)
                dQR(IEIGN_BT) = s(i,j,k+1,QMAGX) - s(i,j,k,QMAGX)
                dQR(IEIGN_BTT) = s(i,j,k+1,QMAGY) - s(i,j,k,QMAGY)

             end if

             ! compute the eigenvectors and eigenvalues for this coordinate direction

             call evals(lam,  s(i,j,k,:), idir)

             if (idir == 1) then
                call lvecx(leig, s(i,j,k,:))
                call rvecx(reig, s(i,j,k,:))

             else if (idir == 2) then
                call lvecy(leig, s(i,j,k,:))
                call rvecy(reig, s(i,j,k,:))

             else
                call lvecz(leig, s(i,j,k,:))
                call rvecz(reig, s(i,j,k,:))
             end if

             ! MHD Source Terms -- from the Miniati paper, Eq. 32 and 33
             smhd(IEIGN_RHO) = 0.0d0
             smhd(IEIGN_U) = s(i,j,k,QMAGX)/s(i,j,k,QRHO)
             smhd(IEIGN_V) = s(i,j,k,QMAGY)/s(i,j,k,QRHO)
             smhd(IEIGN_W) = s(i,j,k,QMAGZ)/s(i,j,k,QRHO)
             smhd(IEIGN_P) = s(i,j,k,QMAGX)*s(i,j,k,QU) + &
                             s(i,j,k,QMAGY)*s(i,j,k,QV) + &
                             s(i,j,k,QMAGZ)*s(i,j,k,QW)

             if (idir == 1) then
                smhd(IEIGN_BT) = s(i,j,k,QV)
                smhd(IEIGN_BTT) = s(i,j,k,QW)

                ! cross-talk of normal magnetic field direction
                smhd(:) = smhd(:)*(bx(i+1,j,k) - bx(i,j,k))/dx(idir)

             else if (idir == 2) then
                smhd(IEIGN_BT) = s(i,j,k,QU)
                smhd(IEIGN_BTT) = s(i,j,k,QW)

                ! cross-talk of normal magnetic field direction
                smhd(:) = smhd(:)*(by(i,j+1,k) - by(i,j,k))/dx(idir)

             else
                smhd(IEIGN_BT) = s(i,j,k,QU)
                smhd(IEIGN_BTT) = s(i,j,k,QV)

                ! cross-talk of normal magnetic field direction
                smhd(:) = smhd(:)*(bz(i,j,k+1) - bz(i,j,k))/dx(idir)

             end if


             ! Perform the characteristic projection.  Since we are using
             ! Using HLLD, we sum over all eigenvalues -- see the discussion after Eq. 31
             do ii = 1, NEIGN
                dL = dot_product(leig(ii,:), dQL)
                dR = dot_product(leig(ii,:), dQR)
                call slope(dW, dL, dR, flatn(i,j,k))

                summ_p(:) = summ_p(:) + (1.0d0 - dtdx*lam(ii)) * dW * reig(:,ii)
                summ_m(:) = summ_m(:) - (1.0d0 + dtdx*lam(ii)) * dW * reig(:,ii)
             enddo

             ! left state at i+1/2
             Ip(i,j,k,QRHO,idir) = max(small_dens, s(i,j,k,QRHO) + 0.5d0*summ_p(IEIGN_RHO) + 0.5d0*dt*smhd(IEIGN_RHO))
             Ip(i,j,k,QU,idir) = s(i,j,k,QU) + 0.5d0*summ_p(IEIGN_U) + 0.5d0*dt*smhd(IEIGN_U)
             Ip(i,j,k,QV,idir) = s(i,j,k,QV) + 0.5d0*summ_p(IEIGN_V) + 0.5d0*dt*smhd(IEIGN_V)
             Ip(i,j,k,QW,idir) = s(i,j,k,QW) + 0.5d0*summ_p(IEIGN_W) + 0.5d0*dt*smhd(IEIGN_W)
             Ip(i,j,k,QPRES,idir) = max(small_pres, s(i,j,k,QPRES) + 0.5d0*summ_p(IEIGN_P) + 0.5d0*dt*smhd(5))

             if (idir == 1) then
                Ip(i,j,k,QMAGX,idir) = bx(i+1,j,k) !! Bx stuff
                Ip(i,j,k,QMAGY,idir) = s(i,j,k,QMAGY) + 0.5d0*summ_p(IEIGN_BT) + 0.5d0*dt*smhd(IEIGN_BT)
                Ip(i,j,k,QMAGZ,idir) = s(i,j,k,QMAGZ) + 0.5d0*summ_p(IEIGN_BTT) + 0.5d0*dt*smhd(IEIGN_BTT)

             else if (idir == 2) then
                Ip(i,j,k,QMAGX,idir) = s(i,j,k,QMAGX) + 0.5d0*summ_p(IEIGN_BT) + 0.5d0*dt*smhd(IEIGN_BT)
                Ip(i,j,k,QMAGY,idir) = by(i,j+1,k) !! By stuff
                Ip(i,j,k,QMAGZ,idir) = s(i,j,k,QMAGZ) + 0.5d0*summ_p(IEIGN_BTT) + 0.5d0*dt*smhd(IEIGN_BTT)

             else
                Ip(i,j,k,QMAGX,idir) = s(i,j,k,QMAGX) + 0.5d0*summ_p(IEIGN_BT) + 0.5d0*dt*smhd(IEIGN_BT)
                Ip(i,j,k,QMAGY,idir) = s(i,j,k,QMAGY) + 0.5d0*summ_p(IEIGN_BTT) + 0.5d0*dt*smhd(IEIGN_BTT)
                Ip(i,j,k,QMAGZ,idir) = bz(i,j,k+1) !! Bz stuff
             end if

             ! right state at i-1/2
             Im(i,j,k,QRHO,idir) = max(small_dens, s(i,j,k,QRHO) + 0.5d0*summ_m(IEIGN_RHO) + 0.5d0*dt*smhd(IEIGN_RHO))
             Im(i,j,k,QU,idir) = s(i,j,k,QU) + 0.5d0*summ_m(IEIGN_U) + 0.5d0*dt*smhd(IEIGN_U)
             Im(i,j,k,QV,idir) = s(i,j,k,QV) + 0.5d0*summ_m(IEIGN_V) + 0.5d0*dt*smhd(IEIGN_V)
             Im(i,j,k,QW,idir) = s(i,j,k,QW) + 0.5d0*summ_m(IEIGN_W) + 0.5d0*dt*smhd(IEIGN_W)
             Im(i,j,k,QPRES,idir) = max(small_pres, s(i,j,k,QPRES) + 0.5d0*summ_m(IEIGN_P) + 0.5d0*dt*smhd(IEIGN_P))

             if (idir == 1) then
                Im(i,j,k,QMAGX,idir) = bx(i,j,k) !! Bx stuff
                Im(i,j,k,QMAGY,idir) = s(i,j,k,QMAGY) + 0.5d0*summ_m(IEIGN_BT) + 0.5d0*dt*smhd(IEIGN_BT)
                Im(i,j,k,QMAGZ,idir) = s(i,j,k,QMAGZ) + 0.5d0*summ_m(IEIGN_BTT) + 0.5d0*dt*smhd(IEIGN_BTT)

             else if (idir == 2) then
                Im(i,j,k,QMAGX,idir) = s(i,j,k,QMAGX) + 0.5d0*summ_m(IEIGN_BT) + 0.5d0*dt*smhd(IEIGN_BT)
                Im(i,j,k,QMAGY,idir) = by(i,j,k) !! By stuff
                Im(i,j,k,QMAGZ,idir) = s(i,j,k,QMAGZ) + 0.5d0*summ_m(IEIGN_BTT) + 0.5d0*dt*smhd(IEIGN_BTT)

             else
                Im(i,j,k,QMAGX,idir) = s(i,j,k,QMAGX) + 0.5d0*summ_m(IEIGN_BT) + 0.5d0*dt*smhd(IEIGN_BT)
                Im(i,j,k,QMAGY,idir) = s(i,j,k,QMAGY) + 0.5d0*summ_m(IEIGN_BTT) + 0.5d0*dt*smhd(IEIGN_BTT)
                Im(i,j,k,QMAGZ,idir) = bz(i,j,k) !! Bz stuff
             endif

             ! species
             do ii = QFS, QFS+nspec-1
                if (idir == 1) then
                   dL = s(i,j,k,ii) - s(i-1,j,k,ii)
                   dR = s(i+1,j,k,ii) - s(i,j,k,ii)
                   un = s(i,j,k,QU)

                else if (idir == 2) then
                   dL = s(i,j,k,ii) - s(i,j-1,k,ii)
                   dR = s(i,j+1,k,ii) - s(i,j,k,ii)
                   un = s(i,j,k,QV)

                else
                   dL = s(i,j,k,ii) - s(i,j,k-1,ii)
                   dR = s(i,j,k+1,ii) - s(i,j,k,ii)
                   un = s(i,j,k,QW)
                end if

               call slope(dW, dL, dR, flatn(i,j,k))

               Ip(i,j,k,ii,idir) = s(i,j,k,ii) + 0.5d0*(1.0d0 - dtdx*un) * dW
               Im(i,j,k,ii,idir) = s(i,j,k,ii) - 0.5d0*(1.0d0 + dtdx*un) * dW
             enddo

             ! rho e
             eos_state % rho = Ip(i,j,k,QRHO,idir)
             eos_state % p   = Ip(i,j,k,QPRES,idir)
             eos_state % T   = s(i,j,k,QTEMP) !some initial guess?
             eos_state % xn  = Ip(i,j,k,QFS:QFS+nspec-1,idir)

             call eos(eos_input_rp, eos_state)
             Ip(i,j,k,QREINT,idir) = eos_state % e * eos_state % rho

             eos_state % rho = Im(i,j,k,QRHO,idir)
             eos_state % p   = Im(i,j,k,QPRES,idir)
             eos_state % xn  = Im(i,j,k,QFS:QFS+nspec-1,idir)

             call eos(eos_input_rp, eos_state)
             Im(i,j,k,QREINT,idir) = eos_state % e * eos_state % rho

             ! add source terms
             Ip(i,j,k,QRHO,idir) = max(small_dens, Ip(i,j,k,QRHO,idir) + 0.5d0*dt*srcQ(i,j,k,QRHO))
             Ip(i,j,k,QU,idir) = Ip(i,j,k,QU,idir) + 0.5d0*dt*srcQ(i,j,k,QU)
             Ip(i,j,k,QV,idir) = Ip(i,j,k,QV,idir) + 0.5d0*dt*srcQ(i,j,k,QV)
             Ip(i,j,k,QW,idir) = Ip(i,j,k,QW,idir) + 0.5d0*dt*srcQ(i,j,k,QW)
             Ip(i,j,k,QPRES,idir) = Ip(i,j,k,QPRES,idir) + 0.5d0*dt*srcQ(i,j,k,QPRES)
             Ip(i,j,k,QREINT,idir) = Ip(i,j,k,QREINT,idir) + 0.5d0*dt*srcQ(i,j,k,QREINT)

             Im(i,j,k,QRHO,idir) = max(small_dens, Im(i,j,k,QRHO,idir) + 0.5d0*dt*srcQ(i,j,k,QRHO))
             Im(i,j,k,QU,idir) = Im(i,j,k,QU,idir) + 0.5d0*dt*srcQ(i,j,k,QU)
             Im(i,j,k,QV,idir) = Im(i,j,k,QV,idir) + 0.5d0*dt*srcQ(i,j,k,QV)
             Im(i,j,k,QW,idir) = Im(i,j,k,QW,idir) + 0.5d0*dt*srcQ(i,j,k,QW)
             Im(i,j,k,QPRES,idir) = Im(i,j,k,QPRES,idir) + 0.5d0*dt*srcQ(i,j,k,QPRES)
             Im(i,j,k,QREINT,idir) = Im(i,j,k,QREINT,idir) + 0.5d0*dt*srcQ(i,j,k,QREINT)

          enddo
       enddo
    enddo

  end subroutine plm

  !======================================== Minmod TVD slope limiter =========================================
  subroutine minmod(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in) :: WR, WL
    real(rt), intent(out) :: dW
    dW = 0.d0

    if (abs(WR) < abs(WL) .and. WR*WL > 0.d0) then
       dW = WR
    else if (abs(WL) < abs(WR) .and. WR*WL > 0.d0) then
       dW = WL
    endif

  end subroutine minmod


  !========================================= VanLeer TVD slope limiter =======================================
  subroutine vanleer(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in )       ::  WR, WL
    real(rt), intent(out)       ::  dW
    dW = 0.0d0

    if( WR*WL .gt. 0.0d0 ) then
       dW = 2.0d0*WR*WL/(WR + WL)
    endif

  end subroutine vanleer

  !========================================== centered difference ===========================================
  subroutine centerdif(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in )    :: WR, WL
    real(rt), intent(out)    :: dW

    dW = (WR+WL)/2.0d0

  end subroutine centerdif

  !================================== second order MC ==============================
  subroutine secondMC(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real
    use amrex_constants_module, only : ZERO, HALF, ONE, TWO

    implicit none

    real(rt), intent(in  )    :: WR, WL
    real(rt), intent(out )    :: dW

    real(rt)  :: dlim

    if (WR * WL .ge. ZERO) then
       dlim = TWO * min(abs(WR), abs(WL))
    else
       dlim = ZERO
    endif

    dW = min(HALF * abs(WR + WL), dlim ) * sign(ONE, WR + WL)


  end subroutine secondMC



  !================================================================
  subroutine slope(dW, WR, WL, flat)
    use amrex_fort_module, only : rt => amrex_real

    implicit none
    real(rt), intent(in )    :: flat
    real(rt), intent(in )    :: WR, WL
    real(rt), intent(out)    :: dW

    if  (mhd_plm_slope == 0)  then
       dW = 0.0
    elseif (mhd_plm_slope == 1) then
       call vanleer(dW,WR,WL)
    elseif (mhd_plm_slope == 2) then
       call centerdif(dW,WR,WL)
    elseif (mhd_plm_slope == 3) then
       call secondMC(dW,WR,WL)
    endif

    if (use_flattening == 1) then
        dW = flat * dW
    endif

  end subroutine slope

  !=========================================== Evals =========================================================

  subroutine evals(lam, Q, dir)

    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only: nspec
    implicit none

    real(rt), intent(in)  :: Q(NQ)
    real(rt), intent(out) :: lam(NEIGN) !7 waves
    integer, intent(in)   :: dir !Choose direction, 1 for x, 2 for y, 3 for z

    !The characteristic speeds of the system
    real(rt)  :: cfx, cfy, cfz, cax, cay, caz, csx, csy, csz, ca, as
    type(eos_t) :: eos_state

    !Speeeeeeeedssssss
    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % xn  = Q(QFS:QFS+nspec-1)
    eos_state % T   = Q(QTEMP)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)

    !Alfven
    ca  = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cax = (Q(QMAGX)**2)/Q(QRHO)
    cay = (Q(QMAGY)**2)/Q(QRHO)
    caz = (Q(QMAGZ)**2)/Q(QRHO)

    !Sloooooooooow
    csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
    csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
    csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))

    !Fassssst
    cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
    cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
    cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))

    if(dir.eq.1) then

       !Ax eigenvalues
       lam(1) = Q(QU) - sqrt(cfx)
       lam(2) = Q(QU) - sqrt(cax)
       lam(3) = Q(QU) - sqrt(csx)
       lam(4) = Q(QU)
       lam(5) = Q(QU) + sqrt(csx)
       lam(6) = Q(QU) + sqrt(cax)
       lam(7) = Q(QU) + sqrt(cfx)

    elseif(dir.eq.2) then

       !Ay eigenvalues
       lam(1) = Q(QV) - sqrt(cfy)
       lam(2) = Q(QV) - sqrt(cay)
       lam(3) = Q(QV) - sqrt(csy)
       lam(4) = Q(QV)
       lam(5) = Q(QV) + sqrt(csy)
       lam(6) = Q(QV) + sqrt(cay)
       lam(7) = Q(QV) + sqrt(cfy)

    else

       !Az eigenvalues
       lam(1) = Q(QW) - sqrt(cfz)
       lam(2) = Q(QW) - sqrt(caz)
       lam(3) = Q(QW) - sqrt(csz)
       lam(4) = Q(QW)
       lam(5) = Q(QW) + sqrt(csz)
       lam(6) = Q(QW) + sqrt(caz)
       lam(7) = Q(QW) + sqrt(cfz)
    endif

  end subroutine evals

  !====================================== Left Eigenvectors ===============================================

  !x direction
  subroutine lvecx(leig, Q)
    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ax
    real(rt), intent(in)  ::Q(NQ)
    real(rt), intent(out) ::leig(7,7)

    !The characteristic speeds of the system
    real(rt) :: cfx, cax, csx, ca, as, S, N
    real(rt) :: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz
    type (eos_t) :: eos_state

    !Speeeeeeeedssssss
    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % T   = Q(QTEMP)
    eos_state % xn  = Q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)
    !Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cax = (Q(QMAGX)**2)/Q(QRHO)

    !Sloooooooooow
    csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))

    !Fassssst
    cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))

    !useful constants
    alf = sqrt((as - csx)/(cfx - csx))
    als = sqrt((cfx - as)/(cfx - csx))
    if(cfx - as .lt. 0.d0) als = 0.d0
    if(as - csx .lt. 0.d0) alf = 0.d0

    if(abs(Q(QMAGY)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
       bety = 1.d0/sqrt(2.d0)
       betz = bety
    else
       bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
       betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
    endif

    cff = sqrt(cfx)*alf
    css = sqrt(csx)*als
    S = sign(1.0d0, Q(QMAGX))
    Qf = sqrt(cfx)*alf*S
    Qs = sqrt(csx)*als*S
    N = 0.5d0/as
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))


    leig(1,:) = (/0.d0, -N*Cff, N*Qs*bety, N*Qs*betz, N*alf/Q(QRHO), N*AAs*bety/Q(QRHO), N*AAs*betz/Q(QRHO)/) !u - cf
    leig(2,:) = (/0.d0,  0.d0, -0.5d0*betz, 0.5d0*bety, 0.d0, -0.5d0*S*betz/(sqrt(Q(QRHO))), 0.5d0*bety*S/(sqrt(Q(QRHO)))/) !u - cAx
    leig(3,:) = (/0.d0, -N*Css, -N*Qf*bety, -N*Qf*betz, N*als/Q(QRHO), -N*AAf*bety/Q(QRHO), -N*AAf*betz/Q(QRHO)/) !u - cs
    leig(4,:) = (/1.d0,  0.d0,  0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/) !u
    leig(5,:) = (/0.d0,  N*Css, N*Qf*bety, N*Qf*betz, N*als/Q(QRHO), -N*AAf*bety/Q(QRHO), -N*AAf*betz/Q(QRHO)/) !u + cs
    leig(6,:) = (/0.d0,  0.d0, 0.5d0*betz, -0.5d0*bety, 0.d0, -0.5d0*betz*S/(sqrt(Q(QRHO))), 0.5d0*bety*S/(sqrt(Q(QRHO)))/) !u + cAx
    leig(7,:) = (/0.d0, N*Cff, -N*Qs*bety, -N*Qs*betz, N*alf/Q(QRHO), N*AAs*bety/Q(QRHO), N*AAs*betz/Q(QRHO)/) !u + cf


  end subroutine lvecx

  !y direction
  subroutine lvecy(leig, Q)
    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ay
    real(rt), intent(in) ::Q(NQ)
    real(rt), intent(out) ::leig(7,7)

    !The characteristic speeds of the system
    real(rt) :: cfy, cay, csy, ca, as, S, N
    real(rt) :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

    type (eos_t) :: eos_state

    !Speeeeeeeedssssss
    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % T   = Q(QTEMP)
    eos_state % xn  = Q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)

    !Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cay = (Q(QMAGY)**2)/Q(QRHO)

    !Sloooooooooow
    csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))

    !Fassssst
    cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))

    !useful constants
    alf = sqrt((as - csy)/(cfy - csy))
    als = sqrt((cfy - as)/(cfy - csy))
    if(as - csy .lt. 0.d0) alf = 0.d0
    if(cfy - as .lt. 0.d0) als = 0.d0

    if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
       betx = 1.d0/sqrt(2.d0)
       betz = betx
    else
       betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
       betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
    endif

    cff = sqrt(cfy)*alf
    css = sqrt(csy)*als
    S = sign(1.0d0, Q(QMAGY))
    Qf = sqrt(cfy)*alf*S
    Qs = sqrt(csy)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))
    N = 0.5d0/as

    !Need to double check the rows
    leig(1,:) = (/0.d0, N*Qs*betx, -N*Cff , N*Qs*betz, N*alf/Q(QRHO), N*AAs*betx/Q(QRHO), N*AAs*betz/Q(QRHO)/) ! v - cf
    leig(2,:) = (/0.d0, -0.5d0*betz, 0.d0, 0.5d0*betx, 0.d0, -0.5d0*betz*S/(sqrt(Q(QRHO))), 0.5d0*betx*S/(sqrt(Q(QRHO)))/) ! v - cAy
    leig(3,:) = (/0.d0, -N*Qf*betx, -N*Css, -N*Qf*betz, N*als/Q(QRHO), -N*AAf*betx/Q(QRHO), -N*AAf*betz/Q(QRHO)/) ! v - cs
    leig(4,:) = (/1.d0,  0.d0, 0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/) ! v
    leig(5,:) = (/0.d0, N*Qf*betx, N*Css, N*Qf*betz, N*als/Q(QRHO), -N*AAf*betx/Q(QRHO), -N*AAf*betz/Q(QRHO)/) ! v + cs
    leig(6,:) = (/0.d0, 0.5d0*betz, 0.d0,  -0.5d0*betx, 0.d0, -0.5d0*betz*S/(sqrt(Q(QRHO))), 0.5d0*betx*S/(sqrt(Q(QRHO)))/) ! v + cAy
    leig(7,:) = (/0.d0, -N*Qs*betx, N*Cff, -N*Qs*betz, N*alf/Q(QRHO), N*AAs*betx/Q(QRHO), N*AAs*betz/Q(QRHO)/) ! v + cf


  end subroutine lvecy

  !z direction
  subroutine lvecz(leig, Q)
    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    !returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Az
    real(rt), intent(in)  ::Q(NQ)
    real(rt), intent(out) ::leig(7,7)

    !The characteristic speeds of the system
    real(rt)  :: cfz, caz, csz, ca, as, S, N
    real(rt)  :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

    type (eos_t) :: eos_state

    !Speeeeeeeedssssss

    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % T   = Q(QTEMP)
    eos_state % xn  = Q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)

    !Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    caz = (Q(QMAGZ)**2)/Q(QRHO)

    !Sloooooooooow
    csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))

    !Fassssst
    cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))

    !useful constants
    alf = sqrt((as - csz)/(cfz - csz))
    als = sqrt((cfz - as)/(cfz - csz))

    if(cfz - as .lt. 0.d0) als = 0.d0

    if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGY)).le. 1.d-14) then
       betx = 1.d0/sqrt(2.d0)
       bety = betx
    else
       betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
       bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
    endif

    cff = sqrt(cfz)*alf
    css = sqrt(csz)*als
    S = sign(1.0d0, Q(QMAGZ))
    Qf = sqrt(cfz)*alf*S
    Qs = sqrt(csz)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))
    N = 0.5d0/as

    !Need to double check the order
    leig(1,:) = (/0.d0, N*Qs*betx, N*Qs*bety, -N*Cff, N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO), N*AAs*bety/Q(QRHO)/) !w - cf
    leig(2,:) = (/0.d0, -0.5d0*bety, 0.5d0*betx, 0.d0, 0.d0, -0.5d0*S*bety/(sqrt(Q(QRHO))) , 0.5d0*betx*S/(sqrt(Q(QRHO)))/) !w - cAz
    leig(3,:) = (/0.d0, -N*Qf*betx, -N*Qf*bety, -N*Css, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO), -N*AAf*bety/Q(QRHO)/) !w - cs
    leig(4,:) = (/1.d0,  0.d0 ,  0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/) !w
    leig(5,:) = (/0.d0, N*Qf*betx, N*Qf*bety, N*Css, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO), -N*AAf*bety/Q(QRHO)/) !w + cs
    leig(6,:) = (/0.d0, 0.5d0*bety , -0.5d0*betx, 0.0d0, 0.d0 , -0.5d0*bety*S/(sqrt(Q(QRHO))) , 0.5d0*betx*S/(sqrt(Q(QRHO)))  /) !w + cAz
    leig(7,:) = (/0.d0, -N*Qs*betx, -N*Qs*bety, N*Cff , N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO) , N*AAs*bety/Q(QRHO) /) !w + cf
  end subroutine lvecz

  !====================================== Right Eigenvectors ===============================================
  !x direction
  subroutine rvecx(reig, Q)
    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    !returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ax
    real(rt), intent(in)   ::Q(NQ)
    real(rt), intent(out)  ::reig(7,7)

    !The characteristic speeds of the system
    real(rt)  :: cfx, cax, csx, ca, as, S
    real(rt)  :: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

    type(eos_t) :: eos_state

    !Speeeeeeeedssssss
    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % T   = Q(QTEMP)
    eos_state % xn  = Q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)

    !Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cax = (Q(QMAGX)**2)/Q(QRHO)

    !Sloooooooooow
    csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))

    !Fassssst
    cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))

    !useful constants
    alf = sqrt((as - csx)/(cfx - csx))
    als = sqrt((cfx - as)/(cfx - csx))

    if(cfx - as .lt. 0.d0) als = 0.d0

    if(as - csx .lt. 0.d0) alf = 0.d0

    if(abs(Q(QMAGY)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
       bety = 1.d0/sqrt(2.d0)
       betz = bety
    else
       bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
       betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
    endif

    cff = sqrt(cfx)*alf
    css = sqrt(csx)*als
    S = sign(1.0d0, Q(QMAGX))
    Qf = sqrt(cfx)*alf*S
    Qs = sqrt(csx)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))

    !   u - cf       u - Cax      u - cs   u    u + cs   u + Cax   u + cf
    reig(IEIGN_RHO,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(IEIGN_U,:) = (/-cff , 0.d0, -css, 0.d0, css, 0.d0, cff/)
    reig(IEIGN_V,:) = (/Qs*bety, -betz, -Qf*bety, 0.d0, Qf*bety, betz, -Qs*bety/)
    reig(IEIGN_W,:) = (/Qs*betz, bety, -Qf*betz, 0.d0, Qf*betz, -bety, -Qs*betz/)
    reig(IEIGN_P,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(IEIGN_BT,:) = (/AAs*bety, -betz*S*sqrt(Q(QRHO)), -AAf*bety, 0.d0  , -AAf*bety, -betz*S*sqrt(Q(QRHO)), AAs*bety/)
    reig(IEIGN_BTT,:) = (/AAs*betz, bety*S*sqrt(Q(QRHO)), -AAf*betz, 0.d0, -AAf*betz, bety*S*sqrt(Q(QRHO)), AAs*betz/)


  end subroutine rvecx

  !y direction
  subroutine rvecy(reig, Q)
    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    !returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ay
    real(rt), intent(in)  ::Q(NQ)
    real(rt), intent(out) ::reig(7,7)

    !The characteristic speeds of the system
    real(rt) :: cfy, cay, csy, ca, as, S
    real(rt) :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

    type (eos_t) :: eos_state

    !Speeeeeeeedssssss
    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % T   = Q(QTEMP)
    eos_state % xn  = Q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)
    !Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    cay = (Q(QMAGY)**2)/Q(QRHO)
    !Sloooooooooow
    csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
    !Fassssst
    cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
    !useful constants
    alf = sqrt((as - csy)/(cfy - csy))
    if(as - csy .lt. 0.d0) alf = 0.d0
    als = sqrt((cfy - as)/(cfy - csy))
    if(cfy - as .lt. 0.d0) als = 0.d0
    if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGZ)).le. 1.d-14) then
       betx = 1.d0/sqrt(2.d0)
       betz = betx
    else
       betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
       betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
    endif
    cff = sqrt(cfy)*alf
    css = sqrt(csy)*als
    S = sign(1.0d0, Q(QMAGY))
    Qf = sqrt(cfy)*alf*S
    Qs = sqrt(csy)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))

    !   v - cf   v - Cay   v - cs   v   v + cs   v + Cay   v + cf
    reig(IEIGN_RHO,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(IEIGN_V,:) = (/-cff, 0.d0, -css, 0.d0  , css, 0.d0, cff/)
    reig(IEIGN_U,:) = (/Qs*betx, -betz, -Qf*betx, 0.d0  , Qf*betx, betz, -Qs*betx/)
    reig(IEIGN_W,:) = (/Qs*betz, betx, -Qf*betz, 0.d0  , Qf*betz, -betx , -Qs*betz/)
    reig(IEIGN_P,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(IEIGN_BT,:) = (/AAs*betx, -betz*S*sqrt(Q(QRHO)), -AAf*betx, 0.d0  , -AAf*betx, -betz*S*sqrt(Q(QRHO)) , AAs*betx/)
    reig(IEIGN_BTT,:) = (/AAs*betz, betx*S*sqrt(Q(QRHO)), -AAf*betz, 0.d0, -AAf*betz, betx*S*sqrt(Q(QRHO)), AAs*betz/)


  end subroutine rvecy

  !z direction
  subroutine rvecz(reig, Q)
    use amrex_fort_module, only : rt => amrex_real
    use eos_module, only : eos
    use eos_type_module, only: eos_t, eos_input_rp
    use network, only : nspec

    implicit none

    !returnes reig, where the cols are the right eigenvectors of the characteristic matrix Az
    real(rt), intent(in)  ::Q(NQ)
    real(rt), intent(out) ::reig(7,7)

    !The characteristic speeds of the system
    real(rt) :: cfz, caz, csz, ca, as, S
    real(rt) :: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

    type(eos_t) :: eos_state

    !Speeeeeeeedssssss
    eos_state % rho = Q(QRHO)
    eos_state % p   = Q(QPRES)
    eos_state % T   = Q(QTEMP)
    eos_state % xn  = Q(QFS:QFS+nspec-1)

    call eos(eos_input_rp, eos_state)

    as = eos_state % gam1 * Q(QPRES)/Q(QRHO)

    !Alfven
    ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
    caz = (Q(QMAGZ)**2)/Q(QRHO)

    !Sloooooooooow
    csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))

    !Fassssst
    cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))

    !useful constants
    alf = sqrt((as - csz)/(cfz - csz))
    als = sqrt((cfz - as)/(cfz - csz))

    if(cfz - as .lt. 0.d0) als = 0.d0

    if(abs(Q(QMAGX)).le. 1.d-14 .and.abs(Q(QMAGY)).le. 1.d-14) then
       betx = 1.d0/sqrt(2.d0)
       bety = betx
    else
       betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
       bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
    endif

    cff = sqrt(cfz)*alf
    css = sqrt(csz)*als
    S = sign(1.0d0, Q(QMAGZ))
    Qf = sqrt(cfz)*alf*S
    Qs = sqrt(csz)*als*S
    AAf = sqrt(as)*alf*sqrt(Q(QRHO))
    AAs = sqrt(as)*als*sqrt(Q(QRHO))

    !   w - cf    w - Caz     w - cs     w    w + cs    w + Caz     w + cf
    reig(IEIGN_RHO,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(IEIGN_W,:) = (/-cff , 0.d0, -css, 0.d0, css, 0.d0 , cff/)
    reig(IEIGN_U,:) = (/Qs*betx, -bety, -Qf*betx, 0.d0, Qf*betx, bety, -Qs*betx/)
    reig(IEIGN_V,:) = (/Qs*bety, betx, -Qf*bety, 0.d0, Qf*bety, -betx , -Qs*bety/)
    reig(IEIGN_P,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(IEIGN_BT,:) = (/AAs*betx, -bety*S*sqrt(Q(QRHO)), -AAf*betx, 0.d0, -AAf*betx, -bety*S*sqrt(Q(QRHO)), AAs*betx/)
    reig(IEIGN_BTT,:) = (/AAs*bety, betx*S*sqrt(Q(QRHO)), -AAf*bety, 0.d0, -AAf*bety, betx*S*sqrt(Q(QRHO)), AAs*bety/)

  end subroutine rvecz

end module mhd_plm_module
