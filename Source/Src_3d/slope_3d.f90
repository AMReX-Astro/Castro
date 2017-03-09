module slope_module
  
  use amrex_fort_module, only : rt => c_real
  implicit none

  private

  public uslope, pslope

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uslope(q,flatn,qd_lo,qd_hi, &
                        dqx,dqy,dqz,qpd_lo,qpd_hi, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,nv)

      use mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module
      use bl_constants_module

      use amrex_fort_module, only : rt => c_real
      implicit none

      integer          :: qd_lo(3), qd_hi(3)
      integer          :: qpd_lo(3),qpd_hi(3)
      integer          :: ilo1, ilo2, ihi1, ihi2, kc, k3d, nv

      real(rt)         :: q(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3),nv)
      real(rt)         :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
      real(rt)         :: dqx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      real(rt)         :: dqy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)
      real(rt)         :: dqz(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3),nv)

      integer i, j, k, n

      real(rt)         dlft, drgt, slop, dq1
      real(rt)         dm, dp, dc, ds, sl, dl, dfm, dfp

      integer ilo, ihi      
      
      real(rt)        , pointer::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

      ilo = MIN(ilo1,ilo2)
      ihi = MAX(ihi1,ihi2)

      call bl_allocate (dsgn, ilo-2,ihi+2,ilo-2,ihi+2)
      call bl_allocate (dlim, ilo-2,ihi+2,ilo-2,ihi+2)
      call bl_allocate (  df, ilo-2,ihi+2,ilo-2,ihi+2)
      call bl_allocate (dcen, ilo-2,ihi+2,ilo-2,ihi+2)

      if(plm_iorder.eq.1) then

         do n = 1, nv
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1
                  dqx(i,j,kc,n) = ZERO
                  dqy(i,j,kc,n) = ZERO
                  dqz(i,j,kc,n) = ZERO
               enddo
            enddo
         enddo

      else

         do n = 1, nv 

            ! Compute slopes in first coordinate direction
            do j = ilo2-1, ihi2+1

               ! First compute Fromm slopes
               do i = ilo1-2, ihi1+2
                  dlft = TWO*(q(i ,j,k3d,n) - q(i-1,j,k3d,n))
                  drgt = TWO*(q(i+1,j,k3d,n) - q(i ,j,k3d,n))
                  dcen(i,j) = FOURTH * (dlft+drgt)
                  dsgn(i,j) = sign(ONE, dcen(i,j))
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. ZERO) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = ZERO
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do i = ilo1-1, ihi1+1
                  dq1       = FOUR3RD*dcen(i,j) - SIXTH*(df(i+1,j) + df(i-1,j))
                  dqx(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo

            enddo

            ! Compute slopes in second coordinate direction
            do i = ilo1-1, ihi1+1
               ! First compute Fromm slopes for this column
               do j = ilo2-2, ihi2+2
                  dlft = TWO*(q(i,j ,k3d,n) - q(i,j-1,k3d,n))
                  drgt = TWO*(q(i,j+1,k3d,n) - q(i,j ,k3d,n))
                  dcen(i,j) = FOURTH * (dlft+drgt)
                  dsgn(i,j) = sign( ONE, dcen(i,j) )
                  slop = min( abs(dlft), abs(drgt) )
                  if (dlft*drgt .ge. ZERO) then
                     dlim(i,j) = slop
                  else
                     dlim(i,j) = ZERO
                  endif
                  df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
               enddo

               ! Now compute limited fourth order slopes
               do j = ilo2-1, ihi2+1
                  dq1 = FOUR3RD*dcen(i,j) - SIXTH*( df(i,j+1) + df(i,j-1) )
                  dqy(i,j,kc,n) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dq1))
               enddo
            enddo

            ! Compute slopes in third coordinate direction
            do j = ilo2-1, ihi2+1
               do i = ilo1-1, ihi1+1

                  ! Compute Fromm slope on slab below
                  k = k3d-1
                  dm = TWO*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = TWO*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = FOURTH*(dm+dp)
                  ds = sign( ONE, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. ZERO) then
                     dl = sl
                  else
                     dl = ZERO
                  endif
                  dfm = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on slab above
                  k = k3d+1
                  dm = TWO*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = TWO*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = FOURTH*(dm+dp)
                  ds = sign( ONE, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. ZERO) then
                     dl = sl
                  else
                     dl = ZERO
                  endif
                  dfp = ds*min(dl,abs(dc))

                  ! Compute Fromm slope on current slab
                  k = k3d
                  dm = TWO*(q(i,j,k ,n) - q(i,j,k-1,n))
                  dp = TWO*(q(i,j,k+1,n) - q(i,j,k, n))
                  dc = FOURTH*(dm+dp)
                  ds = sign( ONE, dc )
                  sl = min( abs(dm), abs(dp) )
                  if (dm*dp .ge. ZERO) then
                     dl = sl
                  else
                     dl = ZERO
                  endif

                  ! Now compute limited fourth order slopes
                  dq1 = FOUR3RD*dc - SIXTH*( dfp + dfm )
                  dqz(i,j,kc,n) = flatn(i,j,k3d)*ds*min(dl,abs(dq1))
               enddo
            enddo
         enddo

      endif

      call bl_deallocate (dsgn)
      call bl_deallocate (dlim)
      call bl_deallocate (  df)
      call bl_deallocate (dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pslope(p,rho,flatn,qd_lo,qd_hi, &
                        dpx,dpy,dpz,qpd_lo,qpd_hi, &
                        src,src_lo,src_hi, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,dx)
        
        use mempool_module, only : bl_allocate, bl_deallocate
        use meth_params_module
        use bl_constants_module

        use amrex_fort_module, only : rt => c_real
        implicit none

        integer          :: qd_lo(3), qd_hi(3)
        integer          :: qpd_lo(3),qpd_hi(3)
        integer          :: src_lo(3),src_hi(3)
        integer          :: ilo1, ilo2, ihi1, ihi2, kc, k3d

        real(rt)         :: p  (qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
        real(rt)         :: rho(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
        real(rt)         :: flatn(qd_lo(1):qd_hi(1),qd_lo(2):qd_hi(2),qd_lo(3):qd_hi(3))
        real(rt)         :: dpx(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3))
        real(rt)         :: dpy(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3))
        real(rt)         :: dpz(qpd_lo(1):qpd_hi(1),qpd_lo(2):qpd_hi(2),qpd_lo(3):qpd_hi(3))
        real(rt)         :: src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),QVAR)
        real(rt)         :: dx(3)

        integer i, j, k

        integer ilo,ihi        
        
        real(rt)         dlft, drgt, dp1
        real(rt)         dm, dp, dc, dl, dfm, dfp, ds

        !     Local arrays
        real(rt)        , pointer::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

        ilo = MIN(ilo1,ilo2)
        ihi = MAX(ihi1,ihi2)

        call bl_allocate (dsgn, ilo-2,ihi+2,ilo-2,ihi+2)
        call bl_allocate (dlim, ilo-2,ihi+2,ilo-2,ihi+2)
        call bl_allocate (  df, ilo-2,ihi+2,ilo-2,ihi+2)
        call bl_allocate (dcen, ilo-2,ihi+2,ilo-2,ihi+2)

        if(plm_iorder.eq.1) then

           do j = ilo2-1, ihi2+1
              do i = ilo1-1, ihi1+1
                 dpx(i,j,kc) = ZERO
                 dpy(i,j,kc) = ZERO
                 dpz(i,j,kc) = ZERO
              enddo
           enddo

        else
           ! Compute slopes in first coordinate direction
           do j = ilo2-1, ihi2+1

              ! First compute Fromm slopes
              do i = ilo1-2, ihi1+2

                 dlft = p(i  ,j,k3d) - p(i-1,j,k3d)
                 drgt = p(i+1,j,k3d) - p(i  ,j,k3d)

                 ! Subtract off (rho * acceleration) so as not to limit that part of the slope
                 dlft = dlft - FOURTH * &
                      (rho(i,j,k3d)+rho(i-1,j,k3d))*(src(i,j,k3d,QU)+src(i-1,j,k3d,QU))*dx(1)
                 drgt = drgt - FOURTH * &
                      (rho(i,j,k3d)+rho(i+1,j,k3d))*(src(i,j,k3d,QU)+src(i+1,j,k3d,QU))*dx(1)

                 dcen(i,j) = HALF*(dlft+drgt)
                 dsgn(i,j) = sign(ONE, dcen(i,j))
                 if (dlft*drgt .ge. ZERO) then
                    dlim(i,j) = TWO * min( abs(dlft), abs(drgt) )
                 else
                    dlim(i,j) = ZERO
                 endif
                 df(i,j) = dsgn(i,j)*min( dlim(i,j), abs(dcen(i,j)) )
              enddo

              ! Now limited fourth order slopes
              do i = ilo1-1, ihi1+1
                 dp1         = FOUR3RD*dcen(i,j) - SIXTH*(df(i+1,j) + df(i-1,j))
                 dpx(i,j,kc) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dp1))
                 dpx(i,j,kc) = dpx(i,j,kc) + rho(i,j,k3d)*src(i,j,k3d,QU)*dx(1)
              enddo
           enddo

           ! Compute slopes in second coordinate direction
           do i = ilo1-1, ihi1+1

              ! First compute Fromm slopes
              do j = ilo2-2, ihi2+2
                 dlft = p(i,j  ,k3d) - p(i,j-1,k3d)
                 drgt = p(i,j+1,k3d) - p(i,j  ,k3d)

                 ! Subtract off (rho * acceleration) so as not to limit that part of the slope
                 dlft = dlft - FOURTH * &
                      (rho(i,j,k3d)+rho(i,j-1,k3d))*(src(i,j,k3d,QV)+src(i,j-1,k3d,QV))*dx(2)
                 drgt = drgt - FOURTH * &
                      (rho(i,j,k3d)+rho(i,j+1,k3d))*(src(i,j,k3d,QV)+src(i,j+1,k3d,QV))*dx(2)

                 dcen(i,j) = HALF*(dlft+drgt)
                 dsgn(i,j) = sign( ONE, dcen(i,j) )
                 if (dlft*drgt .ge. ZERO) then
                    dlim(i,j) = TWO * min( abs(dlft), abs(drgt) )
                 else
                    dlim(i,j) = ZERO
                 endif
                 df(i,j) = dsgn(i,j)*min( dlim(i,j),abs(dcen(i,j)) )
              enddo

              ! Now limited fourth order slopes
              do j = ilo2-1, ihi2+1
                 dp1 = FOUR3RD*dcen(i,j) - SIXTH*( df(i,j+1) + df(i,j-1) )
                 dpy(i,j,kc) = flatn(i,j,k3d)*dsgn(i,j)*min(dlim(i,j),abs(dp1))
                 dpy(i,j,kc) = dpy(i,j,kc) + rho(i,j,k3d)*src(i,j,k3d,QV)*dx(2)
              enddo
           enddo

           ! Compute slopes in third coordinate direction
           do j = ilo2-1, ihi2+1
              do i = ilo1-1, ihi1+1

                 ! compute Fromm slopes on slab below
                 k = k3d-1
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - FOURTH * (rho(i,j,k)+rho(i,j,k-1))* &
                      (src(i,j,k,QW)+src(i,j,k-1,QW))*dx(3)
                 dp = dp - FOURTH * (rho(i,j,k)+rho(i,j,k+1))* &
                      (src(i,j,k,QW)+src(i,j,k+1,QW))*dx(3)
                 dc = HALF*(dm+dp)
                 ds = sign( ONE, dc )
                 if (dm*dp .ge. ZERO) then
                    dl = TWO * min( abs(dm), abs(dp) )
                 else
                    dl = ZERO
                 endif
                 dfm = ds*min(dl,abs(dc))

                 ! compute Fromm slopes on slab above
                 k = k3d+1
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - FOURTH * (rho(i,j,k)+rho(i,j,k-1))* &
                      (src(i,j,k,QW)+src(i,j,k-1,QW))*dx(3)
                 dp = dp - FOURTH * (rho(i,j,k)+rho(i,j,k+1))* &
                      (src(i,j,k,QW)+src(i,j,k+1,QW))*dx(3)
                 dc = HALF*(dm+dp)
                 ds = sign( ONE, dc )
                 if (dm*dp .ge. ZERO) then
                    dl = TWO * min( abs(dm), abs(dp) )
                 else
                    dl = ZERO
                 endif
                 dfp = ds*min(dl,abs(dc))

                 ! compute Fromm slopes on current slab
                 k = k3d
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - FOURTH * (rho(i,j,k)+rho(i,j,k-1))* &
                      (src(i,j,k,QW)+src(i,j,k-1,QW))*dx(3)
                 dp = dp - FOURTH * (rho(i,j,k)+rho(i,j,k+1))* &
                      (src(i,j,k,QW)+src(i,j,k+1,QW))*dx(3)
                 dc = HALF*(dm+dp)
                 ds = sign( ONE, dc )
                 if (dm*dp .ge. ZERO) then
                    dl = TWO * min( abs(dm), abs(dp) )
                 else
                    dl = ZERO
                 endif

                 ! now limited fourth order slopes
                 dp1 = FOUR3RD*dc - SIXTH*( dfp + dfm )
                 dpz(i,j,kc) = flatn(i,j,k3d)*ds*min(dl,abs(dp1))
                 dpz(i,j,kc) = dpz(i,j,kc) + rho(i,j,k3d)*src(i,j,k3d,QW)*dx(3)
              enddo
           enddo

        endif

        call bl_deallocate (dsgn)
        call bl_deallocate (dlim)
        call bl_deallocate (  df)
        call bl_deallocate (dcen)

      end subroutine pslope

end module slope_module
