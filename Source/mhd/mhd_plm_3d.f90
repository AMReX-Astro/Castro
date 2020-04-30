module mhd_plm_module

  !Module that gives a piecewise linear interpolation for the primitive variables
  !They are projected onto the characteristic variables for tracing.

  use amrex_fort_module, only : rt => amrex_real
  use meth_params_module
  implicit none

  private centerdif,vanleer, lvecx, lvecy, lvecz, rvecx, rvecy, rvecz, evals, slope

  public plm

contains

  !
  ! characteristics based on u
  !
  !===========================================================================
  ! This is called from within threaded loops in advance_mhd_tile so *no* OMP here ...
  !===========================================================================

  subroutine plm(lo, hi, s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,&
                 flatn, & 
                 bx, bxlo, bxhi, &
                 by, bylo, byhi, &
                 bz, bzlo, bzhi, &
                 Ip,Im, ilo1,ilo2,ilo3,ihi1,ihi2,ihi3, &
                 srcQ, srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3, &
                 dx,dy,dz,dt)


    use network, only: nspec
    use eos_module
    use eos_type_module, only: eos_t, eos_input_rp

    implicit none

    integer , intent(in   ) ::  s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, lo(3), hi(3)
    integer , intent(in   ) ::  ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
    integer , intent(in   ) ::  bxlo(3), bxhi(3)
    integer , intent(in   ) ::  bylo(3), byhi(3)
    integer , intent(in   ) ::  bzlo(3), bzhi(3)
    integer , intent(in   ) ::  srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3

    real(rt), intent(in   ) ::  s(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,NQ) !Primitive Vars
    real(rt), intent(in   ) ::  bx(bxlo(1):bxhi(1), bxlo(2):bxhi(2), bxlo(3):bxhi(3))!Face Centered Magnetic Fields
    real(rt), intent(in   ) ::  by(bylo(1):byhi(1), bylo(2):byhi(2), bylo(3):byhi(3))
    real(rt), intent(in   ) ::  bz(bzlo(1):bzhi(1), bzlo(2):bzhi(2), bzlo(3):bzhi(3))
    real(rt), intent(in   ) ::  srcQ(srcq_l1:srcq_h1, srcq_l2:srcq_h2, srcq_l3:srcq_h3,NQSRC)
    real(rt), intent(in   ) ::  flatn(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3)

    real(rt), intent(out) :: Ip(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,NQ,3)
    real(rt), intent(out) :: Im(ilo1:ihi1,ilo2:ihi2,ilo3:ihi3,NQ,3)

    real(rt), intent(in   ) :: dx,dy,dz,dt
 
    real(rt) :: dQL(7), dQR(7), dW, dL, dR, leig(7,7), reig(7,7), lam(7), summ_p(7), summ_m(7)
    real(rt) :: smhd(7)
    real(rt) :: dt_over_a
    integer  :: ii,ibx,iby,ibz, i , j, k, n, idir

    type(eos_t) :: eos_state

    ! Ip and Im are the interface states in each dimension (the last index '3' is the
    ! direction

    ibx = 6
    iby = 7
    ibz = 8

    dt_over_a = dt
    Ip = 0.d0
    Im = 0.d0

       


       !=========================== PLM =========================================
    do k = s_l3+1,s_h3-1
       do j = s_l2+1, s_h2-1
          do i = s_l1+1, s_h1-1

             !=========================== X Direction ========================
             summ_p = 0.d0
             summ_m = 0.d0
             smhd = 0.d0
             dQL = 0.d0
             dQR = 0.d0
             dW = 0.d0
             reig = 0.d0
             leig = 0.d0
             lam = 0.d0
             !Skip Bx
             dQL(1) = s(i,j,k,QRHO) - s(i-1,j,k,QRHO)
             dQL(2) = s(i,j,k,QU) - s(i-1,j,k,QU)
             dQL(3) = s(i,j,k,QV) - s(i-1,j,k,QV)
             dQL(4) = s(i,j,k,QW) - s(i-1,j,k,QW)
             dQL(5) = s(i,j,k,QPRES)-s(i-1,j,k,QPRES)
             dQL(6) = s(i,j,k,QMAGY)-s(i-1,j,k,QMAGY)
             dQL(7) = s(i,j,k,QMAGZ)-s(i-1,j,k,QMAGZ)
             
             dQR(1) = s(i+1,j,k,QRHO) - s(i,j,k,QRHO)
             dQR(2) = s(i+1,j,k,QU) - s(i,j,k,QU)
             dQR(3) = s(i+1,j,k,QV) - s(i,j,k,QV)
             dQR(4) = s(i+1,j,k,QW) - s(i,j,k,QW)
             dQR(5) = s(i+1,j,k,QPRES)-s(i,j,k,QPRES)
             dQR(6) = s(i+1,j,k,QMAGY)-s(i,j,k,QMAGY)
             dQR(7) = s(i+1,j,k,QMAGZ)-s(i,j,k,QMAGZ)
       
            
             call evals(lam, s(i,j,k,:), 1) !!X dir eigenvalues
             call lvecx(leig,s(i,j,k,:))    !!left eigenvectors
             call rvecx(reig,s(i,j,k,:))    !!right eigenvectors
             !MHD Source Terms -- from the Miniati paper, Eq. 32 and 33
             smhd(2) = s(i,j,k,QMAGX)/s(i,j,k,QRHO)
             smhd(3) = s(i,j,k,QMAGY)/s(i,j,k,QRHO)
             smhd(4) = s(i,j,k,QMAGZ)/s(i,j,k,QRHO)
             smhd(5) = s(i,j,k,QMAGX)*s(i,j,k,QU) + s(i,j,k,QMAGY)*s(i,j,k,QV) + s(i,j,k,QMAGZ)*s(i,j,k,QW)
             smhd(6) = s(i,j,k,QV)
             smhd(7) = s(i,j,k,QW)
             smhd       = smhd*(bx(i+1,j,k) - bx(i,j,k))/dx !cross-talk of normal magnetic field direction
            !Interpolate
             !Plus
             !!Using HLLD so sum over all eigenvalues -- see the discussion after Eq. 31
             do ii = 1,7
                dL = dot_product(leig(ii,:),dQL)
                dR = dot_product(leig(ii,:),dQR)
                call slope(dW,dL,dR, flatn(i,j,k))
                summ_p(:) = summ_p(:) + (1 - dt_over_a/dx*lam(ii))*dW*reig(:,ii)
                summ_m(:) = summ_m(:) + (- 1 - dt_over_a/dx*lam(ii))*dW*reig(:,ii)
             enddo
             Ip(i,j,k,QRHO,1) = s(i,j,k,QRHO) + 0.5d0*summ_p(1) + 0.5d0*dt_over_a*smhd(1)
             Ip(i,j,k,QU,1) = s(i,j,k,QU) + 0.5d0*summ_p(2) + 0.5d0*dt_over_a*smhd(2)
             Ip(i,j,k,QV,1) = s(i,j,k,QV) + 0.5d0*summ_p(3) + 0.5d0*dt_over_a*smhd(3)
             Ip(i,j,k,QW,1) = s(i,j,k,QW) + 0.5d0*summ_p(4) + 0.5d0*dt_over_a*smhd(4)
             Ip(i,j,k,QPRES,1) = s(i,j,k,QPRES) + 0.5d0*summ_p(5) + 0.5d0*dt_over_a*smhd(5)

             Ip(i,j,k,QMAGX,1) = bx(i+1,j,k) !! Bx stuff
             Ip(i,j,k,QMAGY:QMAGZ,1)  = s(i,j,k,QMAGY:QMAGZ) + 0.5d0*summ_p(6:7) + 0.5d0*dt_over_a*smhd(6:7)
            
             !species
             do ii = QFS, QFS+nspec-1  
               dL = s(i,j,k,ii) - s(i-1,j,k,ii)
               dR = s(i+1,j,k,ii) - s(i,j,k,ii)
               call slope(dW,dL,dR, flatn(i,j,k))
               Ip(i,j,k,ii,1) = s(i,j,k,ii) + 0.5d0*(1-dt_over_a/dx*s(i,j,k,QU))*dW
               Im(i,j,k,ii,1) = s(i,j,k,ii) + 0.5d0*(- 1 - dt_over_a/dx*s(i,j,k,QU))*dW
             enddo 

             eos_state % rho = Ip(i,j,k,QRHO,1)
             eos_state % p   = Ip(i,j,k,QPRES,1)
             eos_state % T   = s(i,j,k,QTEMP) !some initial guess?
             eos_state % xn  = Ip(i,j,k,QFS:QFS+nspec-1,1)             

             call eos(eos_input_rp, eos_state)
             Ip(i,j,k,QREINT,1) = eos_state % e * eos_state % rho

             !Minus
             Im(i,j,k,QRHO,1) = s(i,j,k,QRHO) +0.5d0*summ_m(1) + 0.5d0*dt_over_a*smhd(1)
             Im(i,j,k,QU,1) = s(i,j,k,QU) +0.5d0*summ_m(2) + 0.5d0*dt_over_a*smhd(2)
             Im(i,j,k,QV,1) = s(i,j,k,QV) +0.5d0*summ_m(3) + 0.5d0*dt_over_a*smhd(3)
             Im(i,j,k,QW,1) = s(i,j,k,QW) +0.5d0*summ_m(4) + 0.5d0*dt_over_a*smhd(4)
             Im(i,j,k,QPRES,1) = s(i,j,k,QPRES) +0.5d0*summ_m(5) + 0.5d0*dt_over_a*smhd(5)

             Im(i,j,k,QMAGX,1)  = bx(i,j,k) !! Bx stuff
             Im(i,j,k,QMAGY:QMAGZ,1)  = s(i,j,k,QMAGY:QMAGZ) +0.5d0*summ_m(6:7) + 0.5d0*dt_over_a*smhd(6:7)
             
             
             eos_state % rho = Im(i,j,k,QRHO,1)
             eos_state % p   = Im(i,j,k,QPRES,1)
             eos_state % xn  = Im(i,j,k,QFS:QFS+nspec-1,1)

             call eos(eos_input_rp, eos_state)
             Im(i,j,k,QREINT,1) = eos_state % e * eos_state % rho



             !======================== Y Direction ===========================

             summ_p = 0.d0
             summ_m = 0.d0
             smhd = 0.d0
             dQL = 0.d0
             dQR = 0.d0
             dW = 0.d0
             reig = 0.d0
             leig = 0.d0
             lam = 0.d0
             !Skip By
             dQL(1) = s(i,j,k,QRHO) - s(i,j-1,k,QRHO)
             dQL(2) = s(i,j,k,QU) - s(i,j-1,k,QU)
             dQL(3) = s(i,j,k,QV) - s(i,j-1,k,QV)
             dQL(4) = s(i,j,k,QW) - s(i,j-1,k,QW)
             dQL(5) = s(i,j,k,QPRES)-s(i,j-1,k,QPRES)
             dQL(6) = s(i,j,k,QMAGX)-s(i,j-1,k,QMAGX)
             dQL(7) = s(i,j,k,QMAGZ)-s(i,j-1,k,QMAGZ)
             
             dQR(1) = s(i,j+1,k,QRHO) - s(i,j,k,QRHO)
             dQR(2) = s(i,j+1,k,QU) - s(i,j,k,QU)
             dQR(3) = s(i,j+1,k,QV) - s(i,j,k,QV)
             dQR(4) = s(i,j+1,k,QW) - s(i,j,k,QW)
             dQR(5) = s(i,j+1,k,QPRES)-s(i,j,k,QPRES)
             dQR(6) = s(i,j+1,k,QMAGX)-s(i,j,k,QMAGX)
             dQR(7) = s(i,j+1,k,QMAGZ)-s(i,j,k,QMAGZ)


             call evals(lam, s(i,j,k,:), 2) !!Y dir eigenvalues
             call lvecy(leig,s(i,j,k,:))    !!left eigenvectors
             call rvecy(reig,s(i,j,k,:))    !!right eigenvectors
             !!Using HLLD so sum over all eigenvalues
             do ii = 1,7
                dL = dot_product(leig(ii,:),dQL)
                dR = dot_product(leig(ii,:),dQR)
                call slope(dW,dL,dR,flatn(i,j,k))
                summ_p(:) = summ_p(:) + (1 - dt_over_a/dy*lam(ii))*dW*reig(:,ii)
                summ_m(:) = summ_m(:) + (- 1 - dt_over_a/dy*lam(ii))*dW*reig(:,ii)
             enddo
             !MHD Source Terms
             smhd(2) = s(i,j,k,QMAGX)/s(i,j,k,QRHO)
             smhd(3) = s(i,j,k,QMAGY)/s(i,j,k,QRHO)
             smhd(4) = s(i,j,k,QMAGZ)/s(i,j,k,QRHO)
             smhd(5) = s(i,j,k,QMAGX)*s(i,j,k,QU) + s(i,j,k,QMAGY)*s(i,j,k,QV) + s(i,j,k,QMAGZ)*s(i,j,k,QW)
             smhd(6) = s(i,j,k,QU)
             smhd(7) = s(i,j,k,QW)
             smhd    = smhd*(by(i,j+1,k) - by(i,j,k))/dy !cross-talk of normal magnetic field direction



             !Interpolate
             Ip(i,j,k,QRHO,2) = s(i,j,k,QRHO) +0.5d0*summ_p(1) + 0.5d0*dt_over_a*smhd(1) !!GAS
             Ip(i,j,k,QU,2)   = s(i,j,k,QU) +0.5d0*summ_p(2) + 0.5d0*dt_over_a*smhd(2)
             Ip(i,j,k,QV,2) = s(i,j,k,QV) +0.5d0*summ_p(3) + 0.5d0*dt_over_a*smhd(3)
             Ip(i,j,k,QW,2) = s(i,j,k,QW) +0.5d0*summ_p(4) + 0.5d0*dt_over_a*smhd(4)
             Ip(i,j,k,QPRES,2) = s(i,j,k,QPRES) +0.5d0*summ_p(5) + 0.5d0*dt_over_a*smhd(5)

             Ip(i,j,k,QMAGX,2)  = s(i,j,k,QMAGX) + 0.5d0*summ_p(6) + 0.5d0*dt_over_a*smhd(6)
             Ip(i,j,k,QMAGY,2)  = by(i,j+1,k) !! By stuff
             Ip(i,j,k,QMAGZ,2)  = s(i,j,k,QMAGZ) + 0.5d0*summ_p(7) + 0.5d0*dt_over_a*smhd(7)
             
             !species
             do ii = QFS, QFS+nspec-1  
               dL = s(i,j,k,ii) - s(i,j-1,k,ii)
               dR = s(i,j+1,k,ii) - s(i,j,k,ii)
               call slope(dW,dL,dR, flatn(i,j,k))
               Ip(i,j,k,ii,2) = s(i,j,k,ii) + 0.5d0*(1-dt_over_a/dy*s(i,j,k,QV))*dW
               Im(i,j,k,ii,2) = s(i,j,k,ii) + 0.5d0*(-1 - dt_over_a/dy*s(i,j,k,QV))*dW
             enddo

             
             eos_state % rho = Ip(i,j,k,QRHO,2)
             eos_state % p   = Ip(i,j,k,QPRES,2)
             eos_state % xn  = Ip(i,j,k,QFS:QFS+nspec-1,2)

             call eos(eos_input_rp, eos_state)
             Ip(i,j,k,QREINT,2) = eos_state % e * eos_state % rho


             Im(i,j,k,QRHO,2) = s(i,j,k,QRHO) + 0.5d0*summ_m(1) + 0.5d0*dt_over_a*smhd(1) !!GAS
             Im(i,j,k,QU,2)   = s(i,j,k,QU) + 0.5d0*summ_m(2) + 0.5d0*dt_over_a*smhd(2)
             Im(i,j,k,QV,2)   = s(i,j,k,QV) + 0.5d0*summ_m(3) + 0.5d0*dt_over_a*smhd(3)
             Im(i,j,k,QW,2) = s(i,j,k,QW) + 0.5d0*summ_m(4) + 0.5d0*dt_over_a*smhd(4)
             Im(i,j,k,QPRES,2) = s(i,j,k,QPRES) + 0.5d0*summ_m(5) + 0.5d0*dt_over_a*smhd(5)

             Im(i,j,k,QMAGX,2) = s(i,j,k,QMAGX) + 0.5d0*summ_m(6) + 0.5d0*dt_over_a*smhd(6)
             Im(i,j,k,QMAGY,2) = by(i,j,k) !! By stuff
             Im(i,j,k,QMAGZ,2) = s(i,j,k,QMAGZ) + 0.5d0*summ_m(7) + 0.5d0*dt_over_a*smhd(7)


             
             eos_state % rho = Im(i,j,k,QRHO,2)
             eos_state % p   = Im(i,j,k,QPRES,2)
             eos_state % xn  = Im(i,j,k,QFS:QFS+nspec-1,2)

             call eos(eos_input_rp, eos_state)
             Im(i,j,k,QREINT,2) = eos_state % e * eos_state % rho



             !======================= Z Direction ============================
             summ_p = 0.d0
             summ_m = 0.d0
             smhd = 0.d0
             dQL = 0.d0
             dQR = 0.d0
             dW = 0.d0
             reig = 0.d0
             leig = 0.d0
             lam = 0.d0
             !Skip Bz
             dQL(1) = s(i,j,k,QRHO) - s(i,j,k-1,QRHO)
             dQL(2) = s(i,j,k,QU) - s(i,j,k-1,QU)
             dQL(3) = s(i,j,k,QV) - s(i,j,k-1,QV)
             dQL(4) = s(i,j,k,QW) - s(i,j,k-1,QW)
             dQL(5) = s(i,j,k,QPRES)-s(i,j,k-1,QPRES)
             dQL(6) = s(i,j,k,QMAGX)-s(i,j,k-1,QMAGX)
             dQL(7) = s(i,j,k,QMAGY)-s(i,j,k-1,QMAGY)
             
             dQR(1) = s(i,j,k+1,QRHO) - s(i,j,k,QRHO)
             dQR(2) = s(i,j,k+1,QU) - s(i,j,k,QU)
             dQR(3) = s(i,j,k+1,QV) - s(i,j,k,QV)
             dQR(4) = s(i,j,k+1,QW) - s(i,j,k,QW)
             dQR(5) = s(i,j,k+1,QPRES)-s(i,j,k,QPRES)
             dQR(6) = s(i,j,k+1,QMAGX)-s(i,j,k,QMAGX)
             dQR(7) = s(i,j,k+1,QMAGY)-s(i,j,k,QMAGY)
                       
             call evals(lam, s(i,j,k,:), 3) !!Z dir eigenvalues
             call lvecz(leig,s(i,j,k,:))    !!left eigenvectors
             call rvecz(reig,s(i,j,k,:))    !!right eigenvectors

             !!Characteristic Tracing
             do ii = 1,7
                dL = dot_product(leig(ii,:),dQL)
                dR = dot_product(leig(ii,:),dQR)
                call slope(dW,dL,dR, flatn(i,j,k))
                summ_p(:) = summ_p(:) + (1 - dt_over_a/dz*lam(ii))*dW*reig(:,ii)
                summ_m(:) = summ_m(:) + (- 1 - dt_over_a/dz*lam(ii))*dW*reig(:,ii)
             enddo
             !MHD Source Terms
             smhd(2) = s(i,j,k,QMAGX)/s(i,j,k,QRHO)
             smhd(3) = s(i,j,k,QMAGY)/s(i,j,k,QRHO)
             smhd(4) = s(i,j,k,QMAGZ)/s(i,j,k,QRHO)
             smhd(5) = s(i,j,k,QMAGX)*s(i,j,k,QU) + s(i,j,k,QMAGY)*s(i,j,k,QV) + s(i,j,k,QMAGZ)*s(i,j,k,QW)
             smhd(6) = s(i,j,k,QU)
             smhd(7) = s(i,j,k,QV)
             smhd  = smhd*(bz(i,j,k+1) - bz(i,j,k))/dz !cross-talk of normal magnetic field direction
             !Interpolate
             Ip(i,j,k,QRHO,3) = s(i,j,k,QRHO) + 0.5d0*summ_p(1) + 0.5d0*dt_over_a*smhd(1) !!GAS
             Ip(i,j,k,QU,3)   = s(i,j,k,QU) + 0.5d0*summ_p(2) + 0.5d0*dt_over_a*smhd(2)
             Ip(i,j,k,QV,3)   = s(i,j,k,QV) + 0.5d0*summ_p(3) + 0.5d0*dt_over_a*smhd(3)
             Ip(i,j,k,QW,3)   = s(i,j,k,QW) + 0.5d0*summ_p(4) + 0.5d0*dt_over_a*smhd(4)
             Ip(i,j,k,QPRES,3) = s(i,j,k,QPRES) + 0.5d0*summ_p(5) + 0.5d0*dt_over_a*smhd(5)

             Ip(i,j,k,QMAGX:QMAGY,3)    = s(i,j,k,QMAGX:QMAGY) + 0.5d0*summ_p(6:7) + 0.5d0*dt_over_a*smhd(6:7)
             Ip(i,j,k,QMAGZ,3)          = bz(i,j,k+1) !! Bz stuff
             !species
             do ii = QFS, QFS+nspec-1  
               dL = s(i,j,k,ii) - s(i,j,k-1,ii)
               dR = s(i,j,k+1,ii) - s(i,j,k,ii)
               call slope(dW,dL,dR,flatn(i,j,k))
               Ip(i,j,k,ii,3) = s(i,j,k,ii) + 0.5d0*(1 - dt_over_a/dz*s(i,j,k,QW))*dW
               Im(i,j,k,ii,3) = s(i,j,k,ii) + 0.5d0*(-1 - dt_over_a/dz*s(i,j,k,QW))*dW
             enddo
     
             eos_state % rho = Ip(i,j,k,QRHO,3)
             eos_state % p   = Ip(i,j,k,QPRES,3)
             eos_state % xn  = Ip(i,j,k,QFS:QFS+nspec-1,3)

             call eos(eos_input_rp, eos_state)
             Ip(i,j,k,QREINT,3) = eos_state % e * eos_state % rho


             Im(i,j,k,QRHO,3) = s(i,j,k,QRHO) + 0.5d0*summ_m(1) + 0.5d0*dt_over_a*smhd(1) !!GAS
             Im(i,j,k,QU,3)   = s(i,j,k,QU) + 0.5d0*summ_m(2) + 0.5d0*dt_over_a*smhd(2)
             Im(i,j,k,QV,3)   = s(i,j,k,QV) + 0.5d0*summ_m(3) + 0.5d0*dt_over_a*smhd(3)
             Im(i,j,k,QW,3)   = s(i,j,k,QW) + 0.5d0*summ_m(4) + 0.5d0*dt_over_a*smhd(4)
             Im(i,j,k,QPRES,3) = s(i,j,k,QPRES) + 0.5d0*summ_m(5) + 0.5d0*dt_over_a*smhd(5)

             Im(i,j,k,QMAGX:QMAGY,3) = s(i,j,k,QMAGX:QMAGY) + 0.5d0*summ_m(6:7) + 0.5d0*dt_over_a*smhd(6:7)
             Im(i,j,k,QMAGZ,3) = bz(i,j,k) !! Bz stuff
             
            
             eos_state % rho = Im(i,j,k,QRHO,3)
             eos_state % p   = Im(i,j,k,QPRES,3)
             eos_state % xn  = Im(i,j,k,QFS:QFS+nspec-1, 3)

             call eos(eos_input_rp, eos_state)
             Im(i,j,k,QREINT,3) = eos_state % e * eos_state % rho



          enddo
       enddo
    enddo


    !Need to add source terms, heating cooling, gravity, etc.

     
  end subroutine plm

  !======================================== Minmod TVD slope limiter =========================================
  subroutine minmod(dW, WR, WL)

    use amrex_fort_module, only : rt => amrex_real

    implicit none

    real(rt), intent(in) :: WR, WL
    real(rt), intent(out) :: dW
    dW = 0.d0

    if(abs(WR).lt.abs(WL).and.WR*WL.gt. 0.d0) then
       dW = WR
    elseif(abs(WL).lt.abs(WR).and.WR*WL.gt.0.d0) then
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
       dlim = TWO * min(abs(WR),abs(WL))
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
    real(rt), intent(out) :: lam(7) !7 waves
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
    reig(1,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(2,:) = (/-cff , 0.d0, -css, 0.d0, css, 0.d0, cff/)
    reig(3,:) = (/Qs*bety, -betz, -Qf*bety, 0.d0, Qf*bety, betz, -Qs*bety/)
    reig(4,:) = (/Qs*betz, bety, -Qf*betz, 0.d0, Qf*betz, -bety, -Qs*betz/)
    reig(5,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(6,:) = (/AAs*bety, -betz*S*sqrt(Q(QRHO)), -AAf*bety, 0.d0  , -AAf*bety, -betz*S*sqrt(Q(QRHO)), AAs*bety/)
    reig(7,:) = (/AAs*betz, bety*S*sqrt(Q(QRHO)), -AAf*betz, 0.d0, -AAf*betz, bety*S*sqrt(Q(QRHO)), AAs*betz/)


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
    reig(1,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(3,:) = (/-cff, 0.d0, -css, 0.d0  , css, 0.d0, cff/)
    reig(2,:) = (/Qs*betx, -betz, -Qf*betx, 0.d0  , Qf*betx, betz, -Qs*betx/)
    reig(4,:) = (/Qs*betz, betx, -Qf*betz, 0.d0  , Qf*betz, -betx , -Qs*betz/)
    reig(5,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(6,:) = (/AAs*betx, -betz*S*sqrt(Q(QRHO)), -AAf*betx, 0.d0  , -AAf*betx, -betz*S*sqrt(Q(QRHO)) , AAs*betx/)
    reig(7,:) = (/AAs*betz, betx*S*sqrt(Q(QRHO)), -AAf*betz, 0.d0, -AAf*betz, betx*S*sqrt(Q(QRHO)), AAs*betz/)


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
    reig(1,:) = (/Q(QRHO)*alf, 0.d0, Q(QRHO)*als, 1.d0, Q(QRHO)*als, 0.d0, Q(QRHO)*alf/)
    reig(4,:) = (/-cff , 0.d0, -css, 0.d0, css, 0.d0 , cff/)
    reig(2,:) = (/Qs*betx, -bety, -Qf*betx, 0.d0, Qf*betx, bety, -Qs*betx/)
    reig(3,:) = (/Qs*bety, betx, -Qf*bety, 0.d0, Qf*bety, -betx , -Qs*bety/)
    reig(5,:) = (/Q(QRHO)*as*alf, 0.d0, Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*als, 0.d0, Q(QRHO)*as*alf/)
    reig(6,:) = (/AAs*betx, -bety*S*sqrt(Q(QRHO)), -AAf*betx, 0.d0, -AAf*betx, -bety*S*sqrt(Q(QRHO)), AAs*betx/)
    reig(7,:) = (/AAs*bety, betx*S*sqrt(Q(QRHO)), -AAf*bety, 0.d0, -AAf*bety, betx*S*sqrt(Q(QRHO)), AAs*bety/)


  end subroutine rvecz
end module mhd_plm_module
