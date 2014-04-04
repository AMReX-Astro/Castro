module slope_module
  
  implicit none

  private

  public uslope, pslope

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uslope(q,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dqx,dqy,dqz,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,nv)

      use meth_params_module
      use bl_constants_module

      implicit none

      integer ilo,ihi
      integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
      integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
      integer ilo1,ilo2,ihi1,ihi2,kc,k3d,nv

      double precision q(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,nv)
      double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
      double precision dqx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      double precision dqy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)
      double precision dqz(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,nv)

      integer i, j, k, n

      double precision dlft, drgt, slop, dq1
      double precision dm, dp, dc, ds, sl, dl, dfm, dfp

      double precision, allocatable::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

      ilo = MIN(ilo1,ilo2)
      ihi = MAX(ihi1,ihi2)

      allocate (dsgn(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dlim(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (  df(ilo-2:ihi+2,ilo-2:ihi+2))
      allocate (dcen(ilo-2:ihi+2,ilo-2:ihi+2))

      if(iorder.eq.1) then

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
            !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,slop,dq1)
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
            !$OMP END PARALLEL DO

            ! Compute slopes in second coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,slop,dq1)
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
            !$OMP END PARALLEL DO

            ! Compute slopes in third coordinate direction
            !$OMP PARALLEL DO PRIVATE(i,j,k,dm,dp,dc,ds,sl,dl,dfm,dfp,dq1)
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
            !$OMP END PARALLEL DO
         enddo

      endif

      deallocate(dsgn,dlim,df,dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pslope(p,rho,flatn,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                        dpx,dpy,dpz,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                        grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
                        ilo1,ilo2,ihi1,ihi2,kc,k3d,dx,dy,dz)
        
        use meth_params_module
        use bl_constants_module

        implicit none

        integer ilo,ihi
        integer qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
        integer qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
        integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
        integer ilo1,ilo2,ihi1,ihi2,kc,k3d

        double precision p  (qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
        double precision rho(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
        double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)
        double precision dpx(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3)
        double precision dpy(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3)
        double precision dpz(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3)
        double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,gv_l3:gv_h3,3)
        double precision dx,dy,dz

        integer i, j, k

        double precision dlft, drgt, dp1
        double precision dm, dp, dc, dl, dfm, dfp, ds

        !     Local arrays
        double precision, allocatable::dsgn(:,:),dlim(:,:),df(:,:),dcen(:,:)

        ilo = MIN(ilo1,ilo2)
        ihi = MAX(ihi1,ihi2)

        allocate (dsgn(ilo-2:ihi+2,ilo-2:ihi+2))
        allocate (dlim(ilo-2:ihi+2,ilo-2:ihi+2))
        allocate (  df(ilo-2:ihi+2,ilo-2:ihi+2))
        allocate (dcen(ilo-2:ihi+2,ilo-2:ihi+2))

        if(iorder.eq.1) then

           do j = ilo2-1, ihi2+1
              do i = ilo1-1, ihi1+1
                 dpx(i,j,kc) = ZERO
                 dpy(i,j,kc) = ZERO
                 dpz(i,j,kc) = ZERO
              enddo
           enddo

        else
           ! Compute slopes in first coordinate direction
           !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,dp1)
           do j = ilo2-1, ihi2+1

              ! First compute Fromm slopes
              do i = ilo1-2, ihi1+2

                 dlft = p(i  ,j,k3d) - p(i-1,j,k3d)
                 drgt = p(i+1,j,k3d) - p(i  ,j,k3d)

                 ! Subtract off (rho * grav) so as not to limit that part of the slope
                 dlft = dlft - FOURTH * &
                      (rho(i,j,k3d)+rho(i-1,j,k3d))*(grav(i,j,k3d,1)+grav(i-1,j,k3d,1))*dx
                 drgt = drgt - FOURTH * &
                      (rho(i,j,k3d)+rho(i+1,j,k3d))*(grav(i,j,k3d,1)+grav(i+1,j,k3d,1))*dx

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
                 dpx(i,j,kc) = dpx(i,j,kc) + rho(i,j,k3d)*grav(i,j,k3d,1)*dx
              enddo
           enddo
           !$OMP END PARALLEL DO

           ! Compute slopes in second coordinate direction
           !$OMP PARALLEL DO PRIVATE(i,j,dlft,drgt,dp1)
           do i = ilo1-1, ihi1+1

              ! First compute Fromm slopes
              do j = ilo2-2, ihi2+2
                 dlft = p(i,j  ,k3d) - p(i,j-1,k3d)
                 drgt = p(i,j+1,k3d) - p(i,j  ,k3d)

                 ! Subtract off (rho * grav) so as not to limit that part of the slope
                 dlft = dlft - FOURTH * &
                      (rho(i,j,k3d)+rho(i,j-1,k3d))*(grav(i,j,k3d,2)+grav(i,j-1,k3d,2))*dy
                 drgt = drgt - FOURTH * &
                      (rho(i,j,k3d)+rho(i,j+1,k3d))*(grav(i,j,k3d,2)+grav(i,j+1,k3d,2))*dy

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
                 dpy(i,j,kc) = dpy(i,j,kc) + rho(i,j,k3d)*grav(i,j,k3d,2)*dy
              enddo
           enddo
           !$OMP END PARALLEL DO

           ! Compute slopes in third coordinate direction
           !$OMP PARALLEL DO PRIVATE(i,j,k,dm,dp,dc,ds,dl,dfm,dfp,dp1)
           do j = ilo2-1, ihi2+1
              do i = ilo1-1, ihi1+1

                 ! compute Fromm slopes on slab below
                 k = k3d-1
                 dm = p(i,j,k  ) - p(i,j,k-1)
                 dp = p(i,j,k+1) - p(i,j,k  )
                 dm = dm - FOURTH * (rho(i,j,k)+rho(i,j,k-1))* &
                      (grav(i,j,k,3)+grav(i,j,k-1,3))*dz
                 dp = dp - FOURTH * (rho(i,j,k)+rho(i,j,k+1))* &
                      (grav(i,j,k,3)+grav(i,j,k+1,3))*dz
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
                      (grav(i,j,k,3)+grav(i,j,k-1,3))*dz
                 dp = dp - FOURTH * (rho(i,j,k)+rho(i,j,k+1))* &
                      (grav(i,j,k,3)+grav(i,j,k+1,3))*dz
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
                      (grav(i,j,k,3)+grav(i,j,k-1,3))*dz
                 dp = dp - FOURTH * (rho(i,j,k)+rho(i,j,k+1))* &
                      (grav(i,j,k,3)+grav(i,j,k+1,3))*dz
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
                 dpz(i,j,kc) = dpz(i,j,kc) + rho(i,j,k3d)*grav(i,j,k3d,3)*dz
              enddo
           enddo
           !$OMP END PARALLEL DO

        endif

        deallocate(dsgn,dlim,df,dcen)

      end subroutine pslope

end module slope_module
