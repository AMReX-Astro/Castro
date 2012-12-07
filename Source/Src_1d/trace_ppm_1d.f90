module trace_ppm_module

  implicit none

  private

  public trace_ppm

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine trace_ppm(q,dq,c,flatn,gamc,qd_l1,qd_h1, &
                           dloga,dloga_l1,dloga_h1, &
                           srcQ,src_l1,src_h1,&
                           grav,gv_l1,gv_h1, &
                           qxm,qxp,qpd_l1,qpd_h1, &
                           ilo,ihi,domlo,domhi,dx,dt)

      use network, only : nspec, naux
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QREINT, QPRES, QFA, QFS, QFX, & 
                                     nadv, small_dens, ppm_type, fix_mass_flux
      use prob_params_module, only : physbc_lo, physbc_hi, Outflow
      use ppm_module, only : ppm

      implicit none

      integer ilo,ihi
      integer domlo(1),domhi(1)
      integer    qd_l1,   qd_h1
      integer dloga_l1,dloga_h1
      integer   qpd_l1,  qpd_h1
      integer   src_l1,  src_h1
      integer    gv_l1,   gv_h1
      double precision dx, dt
      double precision     q( qd_l1: qd_h1,QVAR)
      double precision  srcQ(src_l1:src_h1,QVAR)
      double precision  gamc(qd_l1:qd_h1)
      double precision flatn(qd_l1:qd_h1)
      double precision     c(qd_l1:qd_h1)
      double precision dloga(dloga_l1:dloga_h1)

      double precision   dq( qpd_l1: qpd_h1,QVAR)
      double precision  qxm( qpd_l1: qpd_h1,QVAR)
      double precision  qxp( qpd_l1: qpd_h1,QVAR)
      double precision grav(  gv_l1:  gv_h1)

!     Local variables
      integer i
      integer n, iadv
      integer ns, ispec, iaux

      double precision hdt,dtdx
      double precision cc, csq, rho, u, p, rhoe
      double precision drho, du, dp, drhoe
      double precision drhop, dup, dpp, drhoep
      double precision drhom, dum, dpm, drhoem

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision apleft, amleft, azrleft, azeleft
      double precision acmprght, acmpleft
      double precision ascmprght, ascmpleft
      double precision sourcr,sourcp,source,courn,eta,dlogatmp

      logical :: fix_mass_flux_lo, fix_mass_flux_hi

      double precision, allocatable :: Ip(:,:,:)
      double precision, allocatable :: Im(:,:,:)

      fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) &
           .and. (ilo .eq. domlo(1))
      fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) &
           .and. (ihi .eq. domhi(1))

      if (ppm_type .eq. 0) then
        print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
        call bl_error("Error:: ppm_1d.f90 :: trace_ppm")
      end if

       hdt = 0.5d0 * dt
       dtdx = dt/dx

         allocate(Ip(ilo-1:ihi+1,3,QVAR))
         allocate(Im(ilo-1:ihi+1,3,QVAR))

         ! Compute Ip and Im
         do n=1,QVAR
            call ppm(q(:,n),qd_l1,qd_h1,q(:,QU),c,Ip(:,:,n),Im(:,:,n),ilo,ihi,dx,dt)
         end do

         ! Trace to left and right edges using upwind PPM
         do i = ilo-1, ihi+1

             cc = c(i)
             csq = cc**2
             rho = q(i,QRHO)
             u = q(i,QU)
             p = q(i,QPRES)
             rhoe = q(i,QREINT)
             enth = ( (rhoe+p)/rho )/csq

             ! plus state on face i
             drhom  = flatn(i)*(rho  - Im(i,1,QRHO))
             dum    = flatn(i)*(u    - Im(i,1,QU))
             dpm    = flatn(i)*(p    - Im(i,1,QPRES))
             drhoem = flatn(i)*(rhoe - Im(i,1,QREINT))
             
             drho  = flatn(i)*(rho  - Im(i,2,QRHO))
             du    = flatn(i)*(u    - Im(i,2,QU))
             dp    = flatn(i)*(p    - Im(i,2,QPRES))
             drhoe = flatn(i)*(rhoe - Im(i,2,QREINT))
             
             drhop  = flatn(i)*(rho  - Im(i,3,QRHO))
             dup    = flatn(i)*(u    - Im(i,3,QU))
             dpp    = flatn(i)*(p    - Im(i,3,QPRES))
             drhoep = flatn(i)*(rhoe - Im(i,3,QREINT))

             alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth

             if (u-cc .gt. 0.d0) then
                amright = 0.d0
             else if (u-cc .lt. 0.d0) then
                amright = -alpham
             else
                amright = -0.5d0*alpham
             endif
             if (u+cc .gt. 0.d0) then
                apright = 0.d0
             else if (u+cc .lt. 0.d0) then
                apright = -alphap
             else
                apright = -0.5d0*alphap
             endif
             if (u .gt. 0.d0) then
                azrright = 0.d0
                azeright = 0.d0
             else if (u .lt. 0.d0) then
                azrright = -alpha0r
                azeright = -alpha0e
             else
                azrright = -0.5d0*alpha0r
                azeright = -0.5d0*alpha0e        
             endif

             if (i .ge. ilo) then
                qxp(i,QRHO) = rho + apright + amright + azrright
                qxp(i,QRHO) = max(small_dens,qxp(i,QRHO))
                qxp(i,QU) = u + (apright - amright)*cc/rho
                qxp(i,QPRES) = p + (apright + amright)*csq
                qxp(i,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

                ! add non-gravitational source term
                qxp(i  ,QRHO  ) = qxp(i,QRHO  ) + hdt*srcQ(i,QRHO)
                qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
                qxp(i  ,QU    ) = qxp(i,QU    ) + hdt*srcQ(i,QU)
                qxp(i  ,QREINT) = qxp(i,QREINT) + hdt*srcQ(i,QREINT)
                qxp(i  ,QPRES ) = qxp(i,QPRES ) + hdt*srcQ(i,QPRES)
                
                ! add gravitational source term
                qxp(i  ,QU) = qxp(i,QU) + hdt*grav(i)
             end if

             ! minus state on face i+1
             drhom  = flatn(i)*(rho  - Ip(i,1,QRHO))
             dum    = flatn(i)*(u    - Ip(i,1,QU))
             dpm    = flatn(i)*(p    - Ip(i,1,QPRES))
             drhoem = flatn(i)*(rhoe - Ip(i,1,QREINT))

             drho  = flatn(i)*(rho  - Ip(i,2,QRHO))
             du    = flatn(i)*(u    - Ip(i,2,QU))
             dp    = flatn(i)*(p    - Ip(i,2,QPRES))
             drhoe = flatn(i)*(rhoe - Ip(i,2,QREINT))

             drhop  = flatn(i)*(rho  - Ip(i,3,QRHO))
             dup    = flatn(i)*(u    - Ip(i,3,QU))
             dpp    = flatn(i)*(p    - Ip(i,3,QPRES))
             drhoep = flatn(i)*(rhoe - Ip(i,3,QREINT))

             alpham = 0.5d0*(dpm/(rho*cc) - dum)*rho/cc
             alphap = 0.5d0*(dpp/(rho*cc) + dup)*rho/cc
             alpha0r = drho - dp/csq
             alpha0e = drhoe - dp*enth

             if (u-cc .gt. 0.d0) then
                amleft = -alpham
             else if (u-cc .lt. 0.d0) then
                amleft = 0.d0
             else
                amleft = -0.5d0*alpham
             endif
             if (u+cc .gt. 0.d0) then
                apleft = -alphap
             else if (u+cc .lt. 0.d0) then
                apleft = 0.d0
             else
                apleft = -0.5d0*alphap
             endif
             if (u .gt. 0.d0) then
                azrleft = -alpha0r
                azeleft = -alpha0e
             else if (u .lt. 0.d0) then
                azrleft = 0.d0
                azeleft = 0.d0
             else
                azrleft = -0.5d0*alpha0r
                azeleft = -0.5d0*alpha0e
             endif

             if (i .le. ihi) then
                qxm(i+1,QRHO) = rho + apleft + amleft + azrleft
                qxm(i+1,QRHO) = max(qxm(i+1,QRHO),small_dens)
                qxm(i+1,QU) = u + (apleft - amleft)*cc/rho
                qxm(i+1,QPRES) = p + (apleft + amleft)*csq
                qxm(i+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

                ! add non-gravitational source term
                qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + hdt*srcQ(i,QRHO)
                qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
                qxm(i+1,QU    ) = qxm(i+1,QU    ) + hdt*srcQ(i,QU)
                qxm(i+1,QREINT) = qxm(i+1,QREINT) + hdt*srcQ(i,QREINT)
                qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + hdt*srcQ(i,QPRES)

                ! add gravitational source term
                qxm(i+1,QU) = qxm(i+1,QU) + hdt*grav(i)
             end if

             if(dloga(i).ne.0)then
                courn = dtdx*(cc+abs(u))
                eta = (1.d0-courn)/(cc*dt*abs(dloga(i)))
                dlogatmp = min(eta,1.d0)*dloga(i)
                sourcr = -0.5d0*dt*rho*dlogatmp*u
                sourcp = sourcr*csq
                source = sourcp*enth
                if (i .le. ihi) then
                   qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
                   qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
                   qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
                   qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
                end if
                if (i .ge. ilo) then
                   qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
                   qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
                   qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
                   qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
                end if
             endif

         end do

         ! Enforce constant mass flux rate if specified
         if (fix_mass_flux_lo) then
             qxm(ilo,QRHO  ) = q(domlo(1)-1,QRHO)
             qxm(ilo,QU    ) = q(domlo(1)-1,QU  )
             qxm(ilo,QPRES ) = q(domlo(1)-1,QPRES)
             qxm(ilo,QREINT) = q(domlo(1)-1,QREINT)
         end if
   
         ! Enforce constant mass flux rate if specified
         if (fix_mass_flux_hi) then
             qxp(ihi+1,QRHO  ) = q(domhi(1)+1,QRHO)
             qxp(ihi+1,QU    ) = q(domhi(1)+1,QU  )
             qxp(ihi+1,QPRES ) = q(domhi(1)+1,QPRES)
             qxp(ihi+1,QREINT) = q(domhi(1)+1,QREINT)
         end if

         ! Now do the passively advected quantities
         do iadv = 1, nadv
            n = QFA + iadv - 1

            ! plus state on face i
            do i = ilo, ihi+1
               u = q(i,QU)
               if (u .gt. 0.d0) then
                  qxp(i,n) = q(i,n)
               else if (u .lt. 0.d0) then
                  qxp(i,n) = q(i,n) + flatn(i)*(Im(i,2,n) - q(i,n))
               else
                  qxp(i,n) = q(i,n) + 0.5d0*flatn(i)*(Im(i,2,n) - q(i,n))
               endif
            enddo
            if (fix_mass_flux_hi) qxp(ihi+1,n) = q(ihi+1,n)

            ! minus state on face i+1
            do i = ilo-1, ihi
               u = q(i,QU)
               if (u .gt. 0.d0) then
                  qxm(i+1,n) = q(i,n) + flatn(i)*(Ip(i,2,n) - q(i,n))
               else if (u .lt. 0.d0) then
                  qxm(i+1,n) = q(i,n)
               else
                  qxm(i+1,n) = q(i,n) + 0.5d0*flatn(i)*(Ip(i,2,n) - q(i,n))
               endif
            enddo
            if (fix_mass_flux_lo) qxm(ilo,n) = q(ilo-1,n)

         enddo

         do ispec = 1, nspec
            ns = QFS + ispec - 1

            ! plus state on face i
            do i = ilo, ihi+1
               u = q(i,QU)
               if (u .gt. 0.d0) then
                  qxp(i,ns) = q(i,ns)
               else if (u .lt. 0.d0) then
                  qxp(i,ns) = q(i,ns) + flatn(i)*(Im(i,2,ns) - q(i,ns))
               else
                  qxp(i,ns) = q(i,ns) + 0.5d0*flatn(i)*(Im(i,2,ns) - q(i,ns))
               endif
            enddo
            if (fix_mass_flux_hi) qxp(ihi+1,ns) = q(ihi+1,ns)

            ! minus state on face i+1
            do i = ilo-1, ihi
               u = q(i,QU)
               if (u .gt. 0.d0) then
                  qxm(i+1,ns) = q(i,ns) + flatn(i)*(Ip(i,2,ns) - q(i,ns))
               else if (u .lt. 0.d0) then
                  qxm(i+1,ns) = q(i,ns)
               else
                  qxm(i+1,ns) = q(i,ns) + 0.5d0*flatn(i)*(Ip(i,2,ns) - q(i,ns))
               endif
            enddo
            if (fix_mass_flux_lo) qxm(ilo,ns) = q(ilo-1,ns)

         enddo

         do iaux = 1, naux
            ns = QFX + iaux - 1

            ! plus state on face i
            do i = ilo, ihi+1
               u = q(i,QU)
               if (u .gt. 0.d0) then
                  qxp(i,ns) = q(i,ns)
               else if (u .lt. 0.d0) then
                  qxp(i,ns) = q(i,ns) + flatn(i)*(Im(i,2,ns) - q(i,ns))
               else
                  qxp(i,ns) = q(i,ns) + 0.5d0*flatn(i)*(Im(i,2,ns) - q(i,ns))
               endif
            enddo
            if (fix_mass_flux_hi) qxp(ihi+1,ns) = q(ihi+1,ns)

            ! minus state on face i+1
            do i = ilo-1, ihi
               u = q(i,QU)
               if (u .gt. 0.d0) then
                  qxm(i+1,ns) = q(i,ns) + flatn(i)*(Ip(i,2,ns) - q(i,ns))
               else if (u .lt. 0.d0) then
                  qxm(i+1,ns) = q(i,ns)
               else
                  qxm(i+1,ns) = q(i,ns) + 0.5d0*flatn(i)*(Ip(i,2,ns) - q(i,ns))
               endif
            enddo
            if (fix_mass_flux_lo) qxm(ilo,ns) = q(ilo-1,ns)

         enddo

         deallocate(Ip,Im)

     end subroutine trace_ppm

end module trace_ppm_module
