module trace_module

  implicit none

  private

  public trace

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine trace(q,dq,c,flatn,qd_l1,qd_h1, &
                       dloga,dloga_l1,dloga_h1, &
                       srcQ,src_l1,src_h1,&
                       grav,gv_l1,gv_h1, &
                       qxm,qxp,qpd_l1,qpd_h1, &
                       ilo,ihi,domlo,domhi,dx,dt)

      use network, only : nspec, naux
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QREINT, QPRES, QFA, QFS, QFX, & 
                                     nadv, small_dens, ppm_type, fix_mass_flux, use_pslope
      use prob_params_module, only : physbc_lo, physbc_hi, Outflow
      use slope_module, only : uslope, pslope
      implicit none

      integer domlo(1),domhi(1)
      integer ilo,ihi
      integer    qd_l1,   qd_h1
      integer dloga_l1,dloga_h1
      integer   qpd_l1,  qpd_h1
      integer   src_l1,  src_h1
      integer    gv_l1,   gv_h1
      double precision dx, dt
      double precision     q( qd_l1: qd_h1,QVAR)
      double precision  srcQ(src_l1:src_h1,QVAR)
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

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision apleft, amleft, azrleft, azeleft
      double precision acmprght, acmpleft
      double precision ascmprght, ascmpleft
      double precision sourcr,sourcp,source,courn,eta,dlogatmp

      logical :: fix_mass_flux_lo, fix_mass_flux_hi

      fix_mass_flux_lo = &
           (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
      fix_mass_flux_hi = &
           (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in trace with ppm_type != 0'
        call bl_error("Error:: Castro_1d.f90 :: trace")
      end if

      hdt = 0.5d0 * dt
      dtdx = dt/dx

      ! Compute slopes
      if (iorder .eq. 1) then

         dq(ilo-1:ihi+1,1:QVAR) = 0.d0

      else

         call uslope(q,flatn,qd_l1,qd_h1, &
              dq,qpd_l1,qpd_h1, &
              ilo,ihi,QVAR)

         if (use_pslope .eq. 1) &
            call pslope(q(:,QPRES),q(:,QRHO), &
                 flatn      , qd_l1, qd_h1, &
                 dq(:,QPRES),qpd_l1,qpd_h1, &
                 grav       , gv_l1, gv_h1, &
                 ilo,ihi,dx)

      endif

      ! Compute left and right traced states
      do i = ilo-1, ihi+1 

         cc = c(i)
         csq = cc**2
         rho = q(i,QRHO)
         u = q(i,QU)
         p = q(i,QPRES)
         rhoe = q(i,QREINT)
         enth = ( (rhoe+p)/rho )/csq

         drho  = dq(i,QRHO)
         du    = dq(i,QU)
         dp    = dq(i,QPRES) 
         drhoe = dq(i,QREINT)

         alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
         alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
         alpha0r = drho - dp/csq
         alpha0e = drhoe - dp*enth

         if (u-cc .gt. 0.d0) then
            spminus = -1.d0
         else
            spminus = (u-cc)*dtdx
         endif
         if (u+cc .gt. 0.d0) then
            spplus = -1.d0
         else
            spplus = (u+cc)*dtdx
         endif
         if (u .gt. 0.d0) then
            spzero = -1.d0
         else
            spzero = u*dtdx
         endif

         apright = 0.5d0*(-1.d0 - spplus )*alphap
         amright = 0.5d0*(-1.d0 - spminus)*alpham
         azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
         azeright = 0.5d0*(-1.d0 - spzero )*alpha0e

         if (i .ge. ilo) then
            qxp(i,QRHO  ) = rho + apright + amright + azrright
            qxp(i,QRHO  ) = max(small_dens,qxp(i,QRHO))
            qxp(i,QU    ) = u + (apright - amright)*cc/rho
            qxp(i,QPRES ) = p + (apright + amright)*csq
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

         if (u-cc .ge. 0.d0) then
            spminus = (u-cc)*dtdx
         else
            spminus = 1.d0
         endif
         if (u+cc .ge. 0.d0) then
            spplus = (u+cc)*dtdx
         else
            spplus = 1.d0
         endif
         if (u .ge. 0.d0) then
            spzero = u*dtdx
         else
            spzero = 1.d0
         endif

         apleft = 0.5d0*(1.d0 - spplus )*alphap
         amleft = 0.5d0*(1.d0 - spminus)*alpham
         azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
         azeleft = 0.5d0*(1.d0 - spzero )*alpha0e

         if (i .le. ihi) then
            qxm(i+1,QRHO  ) = rho + apleft + amleft + azrleft
            qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
            qxm(i+1,QU    ) = u + (apleft - amleft)*cc/rho
            qxm(i+1,QPRES ) = p + (apleft + amleft)*csq
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
      enddo

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

      do iadv = 1, nadv
         n = QFA + iadv - 1

         ! Right state
         do i = ilo,ihi+1
            u = q(i,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,n)
            qxp(i,n) = q(i,n) + acmprght
         enddo
         if (fix_mass_flux_hi) qxp(ihi+1,n) = q(ihi+1,n)

         ! Left state
         do i = ilo-1,ihi
            u = q(i,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,n)
            qxm(i+1,n) = q(i,n) + acmpleft
         enddo
         if (fix_mass_flux_lo) qxm(ilo,n) = q(ilo-1,n)
      enddo

      do ispec = 1, nspec
         ns = QFS + ispec - 1

         ! Right state
         do i = ilo,ihi+1
            u = q(i,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,ns)
            qxp(i,ns) = q(i,ns) + ascmprght
            if (fix_mass_flux_hi .and. i.eq.domhi(1)+1) & 
               qxp(i,ns) = q(i,ns)
         enddo
         if (fix_mass_flux_hi) qxp(ihi+1,ns) = q(ihi+1,ns)

         ! Left state
         do i = ilo-1,ihi
            u = q(i,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,ns)
            qxm(i+1,ns) = q(i,ns) + ascmpleft
         enddo
         if (fix_mass_flux_lo) qxm(ilo,ns) = q(ilo-1,ns)
      enddo

      do iaux = 1, naux
         ns = QFX + iaux - 1

         ! Right state
         do i = ilo,ihi+1
            u = q(i,QU)
            if (u .gt. 0.d0) then
               spzero = -1.d0
            else
               spzero = u*dtdx
            endif
            ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,ns)
            qxp(i,ns) = q(i,ns) + ascmprght
         enddo
         if (fix_mass_flux_hi) qxp(ihi+1,ns) = q(ihi+1,ns)

         ! Left state
         do i = ilo-1,ihi
            u = q(i,QU)
            if (u .ge. 0.d0) then
               spzero = u*dtdx
            else
               spzero = 1.d0
            endif
            ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,ns)
            qxm(i+1,ns) = q(i,ns) + ascmpleft
         enddo
         if (fix_mass_flux_lo) qxm(ilo,ns) = q(ilo-1,ns)
      end do

     end subroutine trace

end module trace_module
