module ppm_module


  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use amrex_constants_module, only: ZERO, HALF, ONE, TWO, SIXTH, &
                                    TWO3RD, THREE, SIX, SEVEN12TH, TWELFTH
  use prob_params_module, only : dg
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public ppm_reconstruct, ppm_int_profile, ppm_reconstruct_with_eos

contains

  subroutine ppm_reconstruct(s, s_lo, s_hi, ncomp, n, &
                             flatn, f_lo, f_hi, &
                             sxm, sxp, &
#if AMREX_SPACEDIM >= 2
                             sym, syp, &
#endif
#if AMREX_SPACEDIM == 3
                             szm, szp, &
#endif
                             sd_lo, sd_hi, &
                             lo, hi, dx)

    ! perform the ppm reconstruction on component n in the array s and
    ! store the limits of the parabola in s[xyz][mp] indexed with ic.

    ! Here, s[xyz][mp] is for only a single component
    ! s has ncomp components

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real

    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: ncomp, n
    integer, intent(in) ::  sd_lo(3),  sd_hi(3)
    integer, intent(in) ::  f_lo(3),  f_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt)        , intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3), ncomp)
    real(rt)        , intent(in) :: flatn( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3))
    real(rt)        , intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt)        , intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt)        , intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt)        , intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
    real(rt)        , intent(in) :: dx(3)

    ! local
    integer i, j, k

    real(rt)         dsl, dsr, dsc

    ! s_{\ib,+}, s_{\ib,-}
    real(rt)         :: sm, sp

    ! \delta s_{\ib}^{vL}
    real(rt)        , pointer :: dsvl(:,:,:)
    real(rt)         :: dsvlm, dsvl0, dsvlp

    ! s_{i+\half}^{H.O.}
    real(rt)        , pointer :: sedge(:,:,:)

#ifndef AMREX_USE_CUDA
    if (ppm_type .ne. 1) &
         call amrex_error("Should have ppm_type = 1 in ppm_type1")

    if (s_lo(1) .gt. lo(1)-3 .or. s_hi(1) .lt. hi(1)+3) then
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

#if AMREX_SPACEDIM >= 2
    if (s_lo(2) .gt. lo(2)-3 .or. s_hi(2) .lt. hi(2)+3) then
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if
#endif
#if AMREX_SPACEDIM == 3
    if (s_lo(3) .gt. lo(3)-3 .or. s_hi(3) .lt. hi(3)+3) then
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if
#endif
#endif
    ! cell-centered indexing w/extra ghost cell
    call bl_allocate(dsvl, lo(:)-2*dg(:), hi(:)+2*dg(:))

    ! edge-centered indexing
    call bl_allocate(sedge, lo(:)-dg(:), hi(:)+2*dg(:))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at x-edges

    ! compute van Leer slopes in x-direction

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i= lo(1)-2, hi(1)+2

             dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
             dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))

             if (dsl*dsr .gt. ZERO) then
                dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvl(i,j,k) = ZERO
             end if

          end do
       end do
    end do

    ! interpolate s to x-edges
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          !dir$ ivdep
          do i = lo(1)-1, hi(1)+2
             sedge(i,j,k) = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) &
                  - SIXTH*(dsvl(i,j,k)-dsvl(i-1,j,k))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k,n),s(i-1,j,k,n)))
             sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k,n),s(i-1,j,k,n)))
          end do
       end do
    end do

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sm = sedge(i  ,j,k)
             sp = sedge(i+1,j,k)

             ! flatten the parabola BEFORE doing the other
             ! monotonization -- this is the method that Flash does
             sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
             sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

             ! modify using quadratic limiters -- note this version of the limiting comes
             ! from Colella and Sekora (2008), not the original PPM paper.
             if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                sp = s(i,j,k,n)
                sm = s(i,j,k,n)

             else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                !else if (-(sp-sm)**2/SIX > &
                !     (sp - sm)*(s(i,j,k3d) - HALF*(sm + sp))) then
                sp = THREE*s(i,j,k,n) - TWO*sm

             else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                !else if ((sp-sm)*(s(i,j,k3d) - HALF*(sm + sp)) > &
                !     (sp - sm)**2/SIX) then
                sm = THREE*s(i,j,k,n) - TWO*sp
             end if

             sxp(i,j,k) = sp
             sxm(i,j,k) = sm

          end do
       end do
    end do

#if (AMREX_SPACEDIM >= 2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at y-edges

    ! compute van Leer slopes in y-direction
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-2, hi(2)+2
          do i=lo(1)-1, hi(1)+1

             dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
             dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))

             if (dsl*dsr .gt. ZERO) then
                dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvl(i,j,k) = ZERO
             end if

          end do
       end do
    end do

    ! interpolate s to y-edges
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-1, hi(2)+2

          !dir$ ivdep
          do i = lo(1)-1, hi(1)+1
             sedge(i,j,k) = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) &
                  - SIXTH*(dsvl(i,j,k)-dsvl(i,j-1,k))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k,n),s(i,j-1,k,n)))
             sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k,n),s(i,j-1,k,n)))
          end do
       end do
    end do

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sm = sedge(i,j  ,k)
             sp = sedge(i,j+1,k)

             ! flatten the parabola BEFORE doing the other
             ! monotonization -- this is the method that Flash does
             sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
             sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

             ! modify using quadratic limiters
             if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                sp = s(i,j,k,n)
                sm = s(i,j,k,n)

             else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                !else if (-(sp-sm)**2/SIX > &
                !     (sp - sm)*(s(i,j,k3d) - HALF*(sm + sp))) then
                sp = THREE*s(i,j,k,n) - TWO*sm

             else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                !else if ((sp-sm)*(s(i,j,k3d) - HALF*(sm + sp)) > &
                !     (sp - sm)**2/SIX) then
                sm = THREE*s(i,j,k,n) - TWO*sp
             end if

             syp(i,j,k) = sp
             sym(i,j,k) = sm

          end do
       end do
    end do
#endif

#if (AMREX_SPACEDIM == 3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at z-edges

    ! compute van Leer slopes in z-direction

    do k = lo(3)-2, hi(3)+2
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
             dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))

             if (dsl*dsr .gt. ZERO) then
                dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvl(i,j,k) = ZERO
             end if

          end do
       end do
    end do

    ! interpolate s to z-edges
    do k = lo(3)-1, hi(3)+2
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             ! interpolate to lo face
             sedge(i,j,k) = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) &
                  - SIXTH*(dsvl(i,j,k)-dsvl(i,j,k-1))
             ! make sure sedge lies in between adjacent cell-centered values
             sedge(i,j,k) = max(sedge(i,j,k),min(s(i,j,k,n),s(i,j,k-1,n)))
             sedge(i,j,k) = min(sedge(i,j,k),max(s(i,j,k,n),s(i,j,k-1,n)))
          end do
       end do
    end do

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             sm = sedge(i,j,k)
             sp = sedge(i,j,k+1)

             ! flatten the parabola BEFORE doing the other
             ! monotonization -- this is the method that Flash does
             sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
             sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

             ! modify using quadratic limiters
             if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then
                sp = s(i,j,k,n)
                sm = s(i,j,k,n)

             else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then
                !else if (-(sp-sm)**2/SIX > &
                !     (sp - sm)*(s(i,j,k3d) - HALF*(sm + sp))) then
                sp = THREE*s(i,j,k,n) - TWO*sm

             else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then
                !else if ((sp-sm)*(s(i,j,k3d) - HALF*(sm + sp)) > &
                !     (sp - sm)**2/SIX) then
                sm = THREE*s(i,j,k,n) - TWO*sp
             end if

             szp(i,j,k) = sp
             szm(i,j,k) = sm

          end do
       end do
    end do
#endif

    call bl_deallocate(dsvl)
    call bl_deallocate(sedge)

  end subroutine ppm_reconstruct


  subroutine ppm_int_profile(s, s_lo, s_hi, ncomp, n, &
                             q, qd_lo, qd_hi, &
                             qaux, qa_lo, qa_hi, &
                             sxm, sxp, &
#if AMREX_SPACEDIM >= 2
                             sym, syp, &
#endif
#if AMREX_SPACEDIM == 3
                             szm, szp, &
#endif
                             sd_lo, sd_hi, &
                             Ip, Im, I_lo, I_hi, icomp, ic, &
                             lo, hi, dx, dt)

    use meth_params_module, only : NQAUX, QC, NQ, QU, QV, QW

    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: ncomp, n
    integer, intent(in) :: icomp, ic
    integer, intent(in) :: qd_lo(3), qd_hi(3)
    integer, intent(in) :: qa_lo(3), qa_hi(3)
    integer, intent(in) :: sd_lo(3), sd_hi(3)
    integer, intent(in) ::  I_lo(3),  I_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3), ncomp)
    real(rt), intent(in) ::     q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3), NQ)
    real(rt), intent(in) ::  qaux(qa_lo(1):qa_hi(1),qa_lo(2):qa_hi(2),qa_lo(3):qa_hi(3), NQAUX)
    real(rt), intent(in) ::   sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(in) ::   sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(in) ::   sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(in) ::   syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(in) ::   szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(in) ::   szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
    real(rt), intent(inout) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3, icomp)
    real(rt), intent(inout) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3, icomp)

    real(rt), intent(in) :: dx(3), dt
    real(rt) :: speed

    ! local
    integer i,j,k

    real(rt)         dtdx, dtdy, dtdz
    real(rt)         sigma, s6
    real(rt)         :: sm, sp

    dtdx = dt/dx(1)
#if (AMREX_SPACEDIM >= 2)
    dtdy = dt/dx(2)
#endif
#if (AMREX_SPACEDIM == 3)
    dtdz = dt/dx(3)
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sp = sxp(i,j,k)
             sm = sxm(i,j,k)

             ! compute x-component of Ip and Im
             s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)

             ! Ip/m is the integral under the parabola for the extent
             ! that a wave can travel over a timestep
             !
             ! Ip integrates to the right edge of a cell
             ! Im integrates to the left edge of a cell

             ! u-c wave
             speed = q(i,j,k,QU)-qaux(i,j,k,QC)
             sigma = abs(speed)*dtdx

             ! if speed == ZERO, then either branch is the same
             if (speed <= ZERO) then
                Ip(i,j,k,1,1,ic) = sp
                Im(i,j,k,1,1,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,1,1,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,1,1,ic) = sm
             endif

             ! u wave
             speed = q(i,j,k,QU)
             sigma = abs(speed)*dtdx

             if (speed <= ZERO) then
                Ip(i,j,k,1,2,ic) = sp
                Im(i,j,k,1,2,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,1,2,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,1,2,ic) = sm
             endif

             ! u+c wave
             speed = q(i,j,k,QU)+qaux(i,j,k,QC)
             sigma = abs(speed)*dtdx

             if (speed <= ZERO) then
                Ip(i,j,k,1,3,ic) = sp
                Im(i,j,k,1,3,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,1,3,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,1,3,ic) = sm
             endif

          end do
       end do
    end do

#if (AMREX_SPACEDIM >= 2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             ! copy sedge into sp and sm
             sp = syp(i,j,k)
             sm = sym(i,j,k)

             ! compute y-component of Ip and Im
             s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)

             ! v-c wave
             speed = q(i,j,k,QV)-qaux(i,j,k,QC)
             sigma = abs(speed)*dtdy

             if (speed <= ZERO) then
                Ip(i,j,k,2,1,ic) = sp
                Im(i,j,k,2,1,ic) = sm + &
                  HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,2,1,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,2,1,ic) = sm
             endif

             ! v wave
             speed = q(i,j,k,QV)
             sigma = abs(speed)*dtdy

             if (speed <= ZERO) then
                Ip(i,j,k,2,2,ic) = sp
                Im(i,j,k,2,2,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,2,2,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,2,2,ic) = sm
             endif

             ! v+c wave
             speed = q(i,j,k,QV)+qaux(i,j,k,QC)
             sigma = abs(speed)*dtdy

             if (speed <= ZERO) then
                Ip(i,j,k,2,3,ic) = sp
                Im(i,j,k,2,3,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,2,3,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,2,3,ic) = sm
             endif

          end do
       end do
    end do
#endif

#if (AMREX_SPACEDIM == 3)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             sp = szp(i,j,k)
             sm = szm(i,j,k)

             ! compute z-component of Ip and Im
             s6 = SIX*s(i,j,k,n) - THREE*(sm+sp)

             ! w-c wave
             speed = q(i,j,k,QW)-qaux(i,j,k,QC)
             sigma = abs(speed)*dtdz

             if (speed <= ZERO) then
                Ip(i,j,k,3,1,ic) = sp
                Im(i,j,k,3,1,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,3,1,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,3,1,ic) = sm
             endif

             ! w wave
             speed = q(i,j,k,QW)
             sigma = abs(speed)*dtdz

             if (speed <= ZERO) then
                Ip(i,j,k,3,2,ic) = sp
                Im(i,j,k,3,2,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,3,2,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,3,2,ic) = sm
             endif

             ! w+c wave
             speed = q(i,j,k,QW)+qaux(i,j,k,QC)
             sigma = abs(speed)*dtdz

             if (speed <= ZERO) then
                Ip(i,j,k,3,3,ic) = sp
                Im(i,j,k,3,3,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             else
                Ip(i,j,k,3,3,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
                Im(i,j,k,3,3,ic) = sm
             endif

          end do
       end do
    end do
#endif

  end subroutine ppm_int_profile


  subroutine ca_ppm_reconstruct_cuda(lo, hi, &
                                     s, s_lo, s_hi, &
                                     flatn, f_lo, f_hi, &
                                     qm, qm_lo, qm_hi, &
                                     qp, qp_lo, qp_hi) bind(c,name='ca_ppm_reconstruct_cuda')

    use meth_params_module, only: NQ
#ifndef AMREX_USE_GPU
    use amrex_error_module, only: amrex_error
#endif

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: qm_lo(3), qm_hi(3)
    integer,  intent(in   ) :: qp_lo(3), qp_hi(3)

    real(rt), intent(in   ) :: s(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NQ)
    real(rt), intent(in   ) :: flatn(f_lo(1):f_hi(1), f_lo(2):f_hi(2), f_lo(3):f_hi(3))
    real(rt), intent(inout) :: qm(qm_lo(1):qm_hi(1),qm_lo(2):qm_hi(2),qm_lo(3):qm_hi(3),NQ,3)
    real(rt), intent(inout) :: qp(qp_lo(1):qp_hi(1),qp_lo(2):qp_hi(2),qp_lo(3):qp_hi(3),NQ,3)

    ! local
    integer :: i, j, k, n

    real(rt) :: dsl, dsr, dsc
    real(rt) :: dsvl_l, dsvl_r
    real(rt) :: sigma, s6

    ! s_{\ib,+}, s_{\ib,-}
    real(rt) :: sm, sp

    !$gpu

#ifndef AMREX_USE_GPU
    if (s_lo(1) .gt. lo(1)-3 .or. s_lo(2) .gt. lo(2)-3 .or. s_lo(3) .gt. lo(3)-3) then
         print *,'Low bounds of array: ',s_lo(1), s_lo(2),s_lo(3)
         print *,'Low bounds of  loop: ',lo(1),lo(2),lo(3)
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

    if (s_hi(1) .lt. hi(1)+3 .or. s_hi(2) .lt. hi(2)+3 .or. s_hi(3) .lt. hi(3)+3) then
         print *,'Hi  bounds of array: ',s_hi(1), s_hi(2), s_hi(3)
         print *,'Hi  bounds of  loop: ',hi(1),hi(2),hi(3)
         call amrex_error("Need more ghost cells on array in ppm_type1")
      end if
#endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, NQ
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Compute van Leer slopes

                dsl = TWO  * (s(i-1,j,k,n) - s(i-2,j,k,n))
                dsr = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i  ,j,k,n) - s(i-2,j,k,n))
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
                dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to x-edges

                sm = HALF*(s(i,j,k,n)+s(i-1,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sm = max(sm, min(s(i,j,k,n),s(i-1,j,k,n)))
                sm = min(sm, max(s(i,j,k,n),s(i-1,j,k,n)))

                ! Compute van Leer slopes

                dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
                dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsl = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))
                dsr = TWO  * (s(i+2,j,k,n) - s(i+1,j,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i+2,j,k,n) - s(i  ,j,k,n))
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to x-edges

                sp = HALF*(s(i+1,j,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sp = max(sp, min(s(i+1,j,k,n),s(i,j,k,n)))
                sp = min(sp, max(s(i+1,j,k,n),s(i,j,k,n)))

                ! Flatten the parabola
                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters -- note this version of the limiting comes
                ! from Colella and Sekora (2008), not the original PPM paper.
                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                qp(i  ,j,k,n,1) = sp
                qm(i+1,j,k,n,1) = sm

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! y-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, NQ
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Compute van Leer slopes

                dsl = TWO  * (s(i,j-1,k,n) - s(i,j-2,k,n))
                dsr = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j  ,k,n) - s(i,j-2,k,n))
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
                dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to y-edges

                sm = HALF*(s(i,j,k,n)+s(i,j-1,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sm = max(sm, min(s(i,j,k,n),s(i,j-1,k,n)))
                sm = min(sm, max(s(i,j,k,n),s(i,j-1,k,n)))

                ! Compute van Leer slopes

                dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
                dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsl = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))
                dsr = TWO  * (s(i,j+2,k,n) - s(i,j+1,k,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j+2,k,n) - s(i,j  ,k,n))
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to y-edges

                sp = HALF*(s(i,j+1,k,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sp = max(sp, min(s(i,j+1,k,n),s(i,j,k,n)))
                sp = min(sp, max(s(i,j+1,k,n),s(i,j,k,n)))

                ! Flatten the parabola

                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters

                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                qp(i,j  ,k,n,2) = sp
                qm(i,j+1,k,n,2) = sm

             end do
          end do
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! z-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n = 1, NQ
       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)

                ! Compute van Leer slopes

                dsl = TWO  * (s(i,j,k-1,n) - s(i,j,k-2,n))
                dsr = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j,k  ,n) - s(i,j,k-2,n))
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
                dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to z-edges

                sm = HALF*(s(i,j,k,n)+s(i,j,k-1,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sm = max(sm, min(s(i,j,k,n),s(i,j,k-1,n)))
                sm = min(sm, max(s(i,j,k,n),s(i,j,k-1,n)))

                ! Compute van Leer slopes

                dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
                dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
                   dsvl_l = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_l = ZERO
                end if

                dsl = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))
                dsr = TWO  * (s(i,j,k+2,n) - s(i,j,k+1,n))
                if (dsl*dsr .gt. ZERO) then
                   dsc = HALF * (s(i,j,k+2,n) - s(i,j,k  ,n))
                   dsvl_r = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
                else
                   dsvl_r = ZERO
                end if

                ! Interpolate s to z-edges

                sp = HALF*(s(i,j,k+1,n)+s(i,j,k,n)) - SIXTH*(dsvl_r - dsvl_l)

                ! Make sure sedge lies in between adjacent cell-centered values

                sp = max(sp, min(s(i,j,k+1,n),s(i,j,k,n)))
                sp = min(sp, max(s(i,j,k+1,n),s(i,j,k,n)))

                ! Flatten the parabola

                sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
                sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

                ! Modify using quadratic limiters

                if ((sp-s(i,j,k,n))*(s(i,j,k,n)-sm) .le. ZERO) then

                   sp = s(i,j,k,n)
                   sm = s(i,j,k,n)

                else if (abs(sp-s(i,j,k,n)) .ge. TWO*abs(sm-s(i,j,k,n))) then

                   sp = THREE*s(i,j,k,n) - TWO*sm

                else if (abs(sm-s(i,j,k,n)) .ge. TWO*abs(sp-s(i,j,k,n))) then

                   sm = THREE*s(i,j,k,n) - TWO*sp

                end if

                qp(i,j,k  ,n,3) = sp
                qm(i,j,k+1,n,3) = sm

             end do
          end do
       end do
    end do

  end subroutine ca_ppm_reconstruct_cuda


  subroutine ppm_reconstruct_with_eos(lo, hi, &
                                      Ip, Im, Ip_gc, Im_gc, I_lo, I_hi)

    use meth_params_module, only : NQ, QRHO, QTEMP, QPRES, QREINT, QFS, QFX
    use eos_type_module, only : eos_t, eos_input_rt
    use eos_module, only : eos
    use network, only : nspec, naux

    integer, intent(in) :: lo(3), hi(3)
    integer, intent(in) :: I_lo(3), I_hi(3)

    real(rt), intent(inout) :: Ip(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3, NQ)
    real(rt), intent(inout) :: Im(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3, NQ)
    real(rt), intent(inout) :: Ip_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3, 1)
    real(rt), intent(inout) :: Im_gc(I_lo(1):I_hi(1),I_lo(2):I_hi(2),I_lo(3):I_hi(3),1:AMREX_SPACEDIM,1:3, 1)

    integer :: iwave, idim, i, j, k

    type(eos_t) :: eos_state

    ! temperature-based PPM -- if desired, take the Ip(T)/Im(T)
    ! constructed above and use the EOS to overwrite Ip(p)/Im(p)
    ! get an edge-based gam1 here if we didn't get it from the EOS
    ! call above (for ppm_temp_fix = 1)

    do iwave = 1, 3
       do idim = 1, AMREX_SPACEDIM

          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                do i = lo(1), hi(1)

                   eos_state % rho = Ip(i,j,k,idim,iwave,QRHO)
                   eos_state % T   = Ip(i,j,k,idim,iwave,QTEMP)

                   eos_state % xn  = Ip(i,j,k,idim,iwave,QFS:QFS+nspec-1)
                   eos_state % aux = Ip(i,j,k,idim,iwave,QFX:QFX+naux-1)

                   call eos(eos_input_rt, eos_state)

                   Ip(i,j,k,idim,iwave,QPRES)  = eos_state % p
                   Ip(i,j,k,idim,iwave,QREINT) = Ip(i,j,k,idim,iwave,QRHO) * eos_state % e
                   Ip_gc(i,j,k,idim,iwave,1)   = eos_state % gam1


                   eos_state % rho = Im(i,j,k,idim,iwave,QRHO)
                   eos_state % T   = Im(i,j,k,idim,iwave,QTEMP)

                   eos_state % xn  = Im(i,j,k,idim,iwave,QFS:QFS+nspec-1)
                   eos_state % aux = Im(i,j,k,idim,iwave,QFX:QFX+naux-1)

                   call eos(eos_input_rt, eos_state)

                   Im(i,j,k,idim,iwave,QPRES)  = eos_state % p
                   Im(i,j,k,idim,iwave,QREINT) = Im(i,j,k,idim,iwave,QRHO) * eos_state % e
                   Im_gc(i,j,k,idim,iwave,1)   = eos_state % gam1
                end do
             end do
          end do

       end do
    end do

  end subroutine ppm_reconstruct_with_eos

end module ppm_module
