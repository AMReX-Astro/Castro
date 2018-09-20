module ppm_module


  ! this does the parabolic reconstruction on a variable and the (optional)
  ! integration under the characteristic domain of the parabola

  use amrex_constants_module, only: ZERO, HALF, ONE, TWO, SIXTH, &
                                    TWO3RD, THREE, SIX, SEVEN12TH, TWELFTH

  use prob_params_module, only : dg

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  private

  public ppm_reconstruct, ppm_int_profile

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
                             lo, hi, dx, &
                             force_type_in)

    ! perform the ppm reconstruction on component n in the array s and
    ! store the limits of the parabola in s[xyz][mp] indexed with ic.

    ! Here, s[xyz][mp] is for only a single component
    ! s has ncomp components

    use meth_params_module, only : ppm_type

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

    integer, intent(in), optional :: force_type_in

    integer :: ppm_type_to_use

    ppm_type_to_use = ppm_type
    if (present(force_type_in)) ppm_type_to_use = force_type_in

    if (ppm_type_to_use == 1) then

        call ppm_type1(s, s_lo, s_hi, ncomp, n, &
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

    else if (ppm_type_to_use == 2) then

        call ppm_type2(s, s_lo, s_hi, ncomp, n, &
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

    end if

  end subroutine ppm_reconstruct

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type1(s, s_lo, s_hi, ncomp, n, &
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

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : ppm_type

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: ncomp, n
    integer, intent(in) :: sd_lo(3), sd_hi(3)
    integer, intent(in) ::  f_lo(3),  f_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3), ncomp)
    real(rt), intent(in) :: flatn( f_lo(1): f_hi(1), f_lo(2): f_hi(2), f_lo(3): f_hi(3))
    real(rt), intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
    real(rt), intent(in) :: dx(3)

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

#if (AMREX_SPACEDIM >= 2)
    if (s_lo(2) .gt. lo(2)-3 .or. s_hi(2) .lt. hi(2)+3) then
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

             dsc = HALF * (s(i+1,j,k,n) - s(i-1,j,k,n))
             dsl = TWO  * (s(i  ,j,k,n) - s(i-1,j,k,n))
             dsr = TWO  * (s(i+1,j,k,n) - s(i  ,j,k,n))

             if (dsl*dsr .gt. ZERO) then
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
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-2, hi(2)+2
          do i=lo(1)-1, hi(1)+1

             dsc = HALF * (s(i,j+1,k,n) - s(i,j-1,k,n))
             dsl = TWO  * (s(i,j  ,k,n) - s(i,j-1,k,n))
             dsr = TWO  * (s(i,j+1,k,n) - s(i,j  ,k,n))

             if (dsl*dsr .gt. ZERO) then
                dsvl(i,j,k) = sign(ONE,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
             else
                dsvl(i,j,k) = ZERO
             end if

          end do
       end do
    end do

    ! interpolate s to y-edges
    do k = lo(3)-1, hi(3)+1
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

    do k = lo(3)-1, hi(3)+1
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

             dsc = HALF * (s(i,j,k+1,n) - s(i,j,k-1,n))
             dsl = TWO  * (s(i,j,k  ,n) - s(i,j,k-1,n))
             dsr = TWO  * (s(i,j,k+1,n) - s(i,j,k  ,n))

             if (dsl*dsr .gt. ZERO) then
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

  end subroutine ppm_type1

  ! :::
  ! ::: ----------------------------------------------------------------
  ! :::

  subroutine ppm_type2(s, s_lo, s_hi, ncomp, n, &
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

    use amrex_mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only : ppm_type

    use amrex_error_module
    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) ::  s_lo(3),  s_hi(3)
    integer, intent(in) :: ncomp, n
    integer, intent(in) :: sd_lo(3), sd_hi(3)
    integer, intent(in) ::  f_lo(3),  f_hi(3)
    integer, intent(in) :: lo(3), hi(3)

    real(rt), intent(in) ::     s( s_lo(1): s_hi(1), s_lo(2): s_hi(2), s_lo(3): s_hi(3), ncomp)
    real(rt), intent(in) :: flatn(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(inout) :: sxm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(inout) :: sxp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: sym( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(inout) :: syp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
#if AMREX_SPACEDIM == 3
    real(rt), intent(inout) :: szm( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
    real(rt), intent(inout) :: szp( sd_lo(1): sd_hi(1), sd_lo(2): sd_hi(2), sd_lo(3): sd_hi(3))
#endif
    real(rt), intent(in) :: dx(3)


    ! local
    integer i,j,k
    logical extremum, bigp, bigm

    real(rt)         D2, D2C, D2L, D2R, D2LIM, alphap, alpham
    real(rt)         sgn
    real(rt)         dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin
    real(rt)         dachkm, dachkp
    real(rt)         amax, delam, delap

    ! s_{\ib,+}, s_{\ib,-}
    real(rt)         :: sm, sp

    ! s_{i+\half}^{H.O.}
    real(rt)        , pointer :: sedge(:,:,:)

    ! constant used in Colella 2008
    real(rt)        , parameter :: C = 1.25e0_rt

    ! a constant used for testing extrema
    real(rt), parameter :: SMALL = 1.e-10_rt

#ifndef AMREX_USE_CUDA
    if (ppm_type .ne. 2) &
         call amrex_error("Should have ppm_type = 2 in ppm_type2")

    if (s_lo(1) .gt. lo(1)-3 .or. s_hi(1) .lt. hi(1)+3) then
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if

#if (AMREX_SPACEDIM >= 2)
    if (s_lo(2) .gt. lo(2)-3 .or. s_hi(2) .lt. hi(2)+3) then
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if
#endif

#if (AMREX_SPACEDIM == 3)
    if (s_lo(3) .gt. lo(3)-3 .or. s_hi(3) .lt. hi(3)+3) then
         call amrex_error("Need more ghost cells on array in ppm_type1")
    end if
#endif


#endif

    ! edge-centered indexing
    call bl_allocate(sedge, lo(:)-2*dg(:), hi(:)+3*dg(:))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! x-direction
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute s at x-edges

    ! interpolate s to x-edges
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-2, hi(1)+3
             sedge(i,j,k) = SEVEN12TH*(s(i-1,j,k,n)+s(i  ,j,k,n)) &
                  - TWELFTH*(s(i-2,j,k,n)+s(i+1,j,k,n))
             !
             ! limit sedge
             !
             if ((sedge(i,j,k)-s(i-1,j,k,n)) * (s(i,j,k,n)-sedge(i,j,k)) .lt. ZERO) then
                D2  = THREE*(s(i-1,j,k,n) - TWO*sedge(i,j,k) + s(i,j,k,n))
                D2L = s(i-2,j,k,n) - TWO*s(i-1,j,k,n) + s(i,j,k,n)
                D2R = s(i-1,j,k,n) - TWO*s(i,j,k,n) + s(i+1,j,k,n)
                sgn = sign(ONE, D2)
                D2LIM = sgn*max(min(C*sgn*D2L, C*sgn*D2R, sgn*D2), ZERO)
                sedge(i,j,k) = HALF*(s(i-1,j,k,n)+s(i,j,k,n)) - SIXTH*D2LIM
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do k = lo(3)-dg(3), hi(3)+dg(3)
       do j = lo(2)-dg(2), hi(2)+dg(2)
          do i = lo(1)-1, hi(1)+1

             alphap   = sedge(i+1,j,k) - s(i,j,k,n)
             alpham   = sedge(i  ,j,k) - s(i,j,k,n)
             bigp     = abs(alphap) .gt. TWO*abs(alpham)
             bigm     = abs(alpham) .gt. TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                !
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                !
                dafacem   = sedge(i,j,k) - sedge(i-1,j,k)
                dafacep   = sedge(i+2,j,k) - sedge(i+1,j,k)
                dabarm    = s(i,j,k,n) - s(i-1,j,k,n)
                dabarp    = s(i+1,j,k,n) - s(i,j,k,n)
                dafacemin = min(abs(dafacem),abs(dafacep))
                dabarmin  = min(abs(dabarm),abs(dabarp))

                if (dafacemin.ge.dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. ZERO)
             end if

             if (extremum) then
                D2     = SIX*(alpham + alphap)
                D2L    = s(i-2,j,k,n) - TWO*s(i-1,j,k,n) + s(i,j,k,n)
                D2R    = s(i,j,k,n) - TWO*s(i+1,j,k,n) + s(i+2,j,k,n)
                D2C    = s(i-1,j,k,n) - TWO*s(i,j,k,n) + s(i+1,j,k,n)
                sgn    = sign(ONE, D2)
                D2LIM  = max(min(sgn*D2, C*sgn*D2L, C*sgn*D2R, C*sgn*D2C), ZERO)
                alpham = alpham*D2LIM/max(abs(D2), SMALL)
                alphap = alphap*D2LIM/max(abs(D2), SMALL)
             else
                if (bigp) then
                   sgn   = sign(ONE, alpham)
                   amax  = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i-1,j,k,n) - s(i,j,k,n)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge. SMALL) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn   = sign(ONE, alphap)
                   amax  = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i+1,j,k,n) - s(i,j,k,n)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge. SMALL) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm = s(i,j,k,n) + alpham
             sp = s(i,j,k,n) + alphap

             ! flatten the parabola AFTER doing the monotonization
             sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
             sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

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

    ! interpolate s to y-edges
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-2, hi(2)+3
          do i = lo(1)-1, hi(1)+1
             sedge(i,j,k) = SEVEN12TH*(s(i,j-1,k,n)+s(i,j,k,n)) &
                  - TWELFTH*(s(i,j-2,k,n)+s(i,j+1,k,n))
             !
             ! limit sedge
             !
             if ((sedge(i,j,k)-s(i,j-1,k,n)) * (s(i,j,k,n)-sedge(i,j,k)) .lt. ZERO) then
                D2  = THREE*(s(i,j-1,k,n) - TWO*sedge(i,j,k) + s(i,j,k,n))
                D2L = s(i,j-2,k,n) - TWO*s(i,j-1,k,n) + s(i,j,k,n)
                D2R = s(i,j-1,k,n) - TWO*s(i,j,k,n) + s(i,j+1,k,n)
                sgn = sign(ONE, D2)
                D2LIM = sgn*max(min(C*sgn*D2L, C*sgn*D2R, sgn*D2), ZERO)
                sedge(i,j,k) = HALF*(s(i,j-1,k,n) + s(i,j,k,n)) - SIXTH*D2LIM
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             alphap   = sedge(i,j+1,k) - s(i,j,k,n)
             alpham   = sedge(i,j  ,k) - s(i,j,k,n)
             bigp     = abs(alphap) .gt. TWO*abs(alpham)
             bigm     = abs(alpham) .gt. TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                !
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                !
                dafacem   = sedge(i,j,k) - sedge(i,j-1,k)
                dafacep   = sedge(i,j+2,k) - sedge(i,j+1,k)
                dabarm    = s(i,j,k,n) - s(i,j-1,k,n)
                dabarp    = s(i,j+1,k,n) - s(i,j,k,n)
                dafacemin = min(abs(dafacem), abs(dafacep))
                dabarmin  = min(abs(dabarm), abs(dabarp))
                if (dafacemin .ge. dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. ZERO)
             end if

             if (extremum) then
                D2     = SIX*(alpham + alphap)
                D2L    = s(i,j-2,k,n) - TWO*s(i,j-1,k,n) + s(i,j,k,n)
                D2R    = s(i,j,k,n) - TWO*s(i,j+1,k,n) + s(i,j+2,k,n)
                D2C    = s(i,j-1,k,n) - TWO*s(i,j,k,n) + s(i,j+1,k,n)
                sgn    = sign(ONE, D2)
                D2LIM  = max(min(sgn*D2, C*sgn*D2L, C*sgn*D2R, C*sgn*D2C), ZERO)
                alpham = alpham*D2LIM/max(abs(D2), SMALL)
                alphap = alphap*D2LIM/max(abs(D2), SMALL)
             else
                if (bigp) then
                   sgn   = sign(ONE, alpham)
                   amax  = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j-1,k,n) - s(i,j,k,n)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge. SMALL) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn   = sign(ONE, alphap)
                   amax  = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i,j+1,k,n) - s(i,j,k,n)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge. SMALL) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm = s(i,j,k,n) + alpham
             sp = s(i,j,k,n) + alphap

             ! flatten the parabola AFTER doing the monotonization
             sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
             sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

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

    ! interpolate s to z-edges
    do k = lo(3)-2, hi(3)+3
       do j = lo(2)-1, hi(2)+1
          !dir$ ivdep
          do i = lo(1)-1, hi(1)+1
             sedge(i,j,k) = SEVEN12TH*(s(i,j,k-1,n)+s(i,j,k,n)) &
                  - TWELFTH*(s(i,j,k-2,n)+s(i,j,k+1,n))
             !
             ! limit sedge
             !
             if ((sedge(i,j,k)-s(i,j,k-1,n)) * (s(i,j,k,n)-sedge(i,j,k)) .lt. ZERO) then
                D2  = THREE*(s(i,j,k-1,n) - TWO*sedge(i,j,k) + s(i,j,k,n))
                D2L = s(i,j,k-2,n) - TWO*s(i,j,k-1,n) + s(i,j,k,n)
                D2R = s(i,j,k-1,n) - TWO*s(i,j,k,n) + s(i,j,k+1,n)
                sgn = sign(ONE, D2)
                D2LIM = sgn*max(min(C*sgn*D2L, C*sgn*D2R, sgn*D2), ZERO)
                sedge(i,j,k) = HALF*(s(i,j,k-1,n) + s(i,j,k,n)) - SIXTH*D2LIM
             end if
          end do
       end do
    end do
    !
    ! Use Colella 2008 limiters.
    !
    ! This is a new version of the algorithm to eliminate sensitivity to roundoff.
    !
    do k = lo(3)-1, hi(3)+1
       do j = lo(2)-1, hi(2)+1
          do i = lo(1)-1, hi(1)+1

             alphap   = sedge(i,j,k+1) - s(i,j,k,n)
             alpham   = sedge(i,j,k  ) - s(i,j,k,n)
             bigp     = abs(alphap) .gt. TWO*abs(alpham)
             bigm     = abs(alpham) .gt. TWO*abs(alphap)
             extremum = .false.

             if (alpham*alphap .ge. ZERO) then
                extremum = .true.
             else if (bigp .or. bigm) then
                !
                ! Possible extremum. We look at cell centered values and face
                ! centered values for a change in sign in the differences adjacent to
                ! the cell. We use the pair of differences whose minimum magnitude is the
                ! largest, and thus least susceptible to sensitivity to roundoff.
                !
                dafacem   = sedge(i,j,k) - sedge(i,j,k-1)
                dafacep   = sedge(i,j,k+2) - sedge(i,j,k+1)
                dabarm    = s(i,j,k,n) - s(i,j,k-1,n)
                dabarp    = s(i,j,k+1,n) - s(i,j,k,n)
                dafacemin = min(abs(dafacem), abs(dafacep))
                dabarmin  = min(abs(dabarm), abs(dabarp))
                if (dafacemin .ge. dabarmin) then
                   dachkm = dafacem
                   dachkp = dafacep
                else
                   dachkm = dabarm
                   dachkp = dabarp
                endif
                extremum = (dachkm*dachkp .le. ZERO)
             end if

             if (extremum) then
                D2     = SIX*(alpham + alphap)
                D2L    = s(i,j,k-2,n) - TWO*s(i,j,k-1,n) + s(i,j,k,n)
                D2R    = s(i,j,k,n) - TWO*s(i,j,k+1,n) + s(i,j,k+2,n)
                D2C    = s(i,j,k-1,n) - TWO*s(i,j,k,n) + s(i,j,k+1,n)
                sgn    = sign(ONE, D2)
                D2LIM  = max(min(sgn*D2, C*sgn*D2L, C*sgn*D2R, C*sgn*D2C), ZERO)
                alpham = alpham*D2LIM/max(abs(D2), SMALL)
                alphap = alphap*D2LIM/max(abs(D2), SMALL)
             else
                if (bigp) then
                   sgn   = sign(ONE, alpham)
                   amax  = -alphap**2 / (4*(alpham + alphap))
                   delam = s(i,j,k-1,n) - s(i,j,k,n)
                   if (sgn*amax .ge. sgn*delam) then
                      if (sgn*(delam - alpham).ge. SMALL) then
                         alphap = (-TWO*delam - TWO*sgn*sqrt(delam**2 - delam*alpham))
                      else
                         alphap = -TWO*alpham
                      endif
                   endif
                end if
                if (bigm) then
                   sgn   = sign(ONE, alphap)
                   amax  = -alpham**2 / (4*(alpham + alphap))
                   delap = s(i,j,k+1,n) - s(i,j,k,n)
                   if (sgn*amax .ge. sgn*delap) then
                      if (sgn*(delap - alphap).ge. SMALL) then
                         alpham = (-TWO*delap - TWO*sgn*sqrt(delap**2 - delap*alphap))
                      else
                         alpham = -TWO*alphap
                      endif
                   endif
                end if
             end if

             sm = s(i,j,k,n) + alpham
             sp = s(i,j,k,n) + alphap

             ! flatten the parabola AFTER doing the monotonization
             sm = flatn(i,j,k)*sm + (ONE-flatn(i,j,k))*s(i,j,k,n)
             sp = flatn(i,j,k)*sp + (ONE-flatn(i,j,k))*s(i,j,k,n)

             szp(i,j,k) = sp
             szm(i,j,k) = sm

          end do
       end do
    end do
#endif

    call bl_deallocate(sedge)

  end subroutine ppm_type2


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
             sigma = abs(q(i,j,k,QU)-qaux(i,j,k,QC))*dtdx

             if (q(i,j,k,QU)-qaux(i,j,k,QC) <= ZERO) then
                Ip(i,j,k,1,1,ic) = sp
             else
                Ip(i,j,k,1,1,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QU)-qaux(i,j,k,QC) >= ZERO) then
                Im(i,j,k,1,1,ic) = sm
             else
                Im(i,j,k,1,1,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! u wave
             sigma = abs(q(i,j,k,QU))*dtdx

             if (q(i,j,k,QU) <= ZERO) then
                Ip(i,j,k,1,2,ic) = sp
             else
                Ip(i,j,k,1,2,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QU) >= ZERO) then
                Im(i,j,k,1,2,ic) = sm
             else
                Im(i,j,k,1,2,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! u+c wave
             sigma = abs(q(i,j,k,QU)+qaux(i,j,k,QC))*dtdx

             if (q(i,j,k,QU)+qaux(i,j,k,QC) <= ZERO) then
                Ip(i,j,k,1,3,ic) = sp
             else
                Ip(i,j,k,1,3,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QU)+qaux(i,j,k,QC) >= ZERO) then
                Im(i,j,k,1,3,ic) = sm
             else
                Im(i,j,k,1,3,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
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
             sigma = abs(q(i,j,k,QV)-qaux(i,j,k,QC))*dtdy

             if (q(i,j,k,QV)-qaux(i,j,k,QC) <= ZERO) then
                Ip(i,j,k,2,1,ic) = sp
             else
                Ip(i,j,k,2,1,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QV)-qaux(i,j,k,QC) >= ZERO) then
                Im(i,j,k,2,1,ic) = sm
             else
                Im(i,j,k,2,1,ic) = sm + &
                  HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! v wave
             sigma = abs(q(i,j,k,QV))*dtdy

             if (q(i,j,k,QV) <= ZERO) then
                Ip(i,j,k,2,2,ic) = sp
             else
                Ip(i,j,k,2,2,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QV) >= ZERO) then
                Im(i,j,k,2,2,ic) = sm
             else
                Im(i,j,k,2,2,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! v+c wave
             sigma = abs(q(i,j,k,QV)+qaux(i,j,k,QC))*dtdy

             if (q(i,j,k,QV)+qaux(i,j,k,QC) <= ZERO) then
                Ip(i,j,k,2,3,ic) = sp
             else
                Ip(i,j,k,2,3,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QV)+qaux(i,j,k,QC) >= ZERO) then
                Im(i,j,k,2,3,ic) = sm
             else
                Im(i,j,k,2,3,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
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
             sigma = abs(q(i,j,k,QW)-qaux(i,j,k,QC))*dtdz

             if (q(i,j,k,QW)-qaux(i,j,k,QC) <= ZERO) then
                Ip(i,j,k,3,1,ic) = sp
             else
                Ip(i,j,k,3,1,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QW)-qaux(i,j,k,QC) >= ZERO) then
                Im(i,j,k,3,1,ic) = sm
             else
                Im(i,j,k,3,1,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! w wave
             sigma = abs(q(i,j,k,QW))*dtdz

             if (q(i,j,k,QW) <= ZERO) then
                Ip(i,j,k,3,2,ic) = sp
             else
                Ip(i,j,k,3,2,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QW) >= ZERO) then
                Im(i,j,k,3,2,ic) = sm
             else
                Im(i,j,k,3,2,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

             ! w+c wave
             sigma = abs(q(i,j,k,QW)+qaux(i,j,k,QC))*dtdz

             if (q(i,j,k,QW)+qaux(i,j,k,QC) <= ZERO) then
                Ip(i,j,k,3,3,ic) = sp
             else
                Ip(i,j,k,3,3,ic) = sp - &
                     HALF*sigma*(sp-sm-(ONE-TWO3RD*sigma)*s6)
             endif

             if (q(i,j,k,QW)+qaux(i,j,k,QC) >= ZERO) then
                Im(i,j,k,3,3,ic) = sm
             else
                Im(i,j,k,3,3,ic) = sm + &
                     HALF*sigma*(sp-sm+(ONE-TWO3RD*sigma)*s6)
             endif

          end do
       end do
    end do
#endif

  end subroutine ppm_int_profile

end module ppm_module
