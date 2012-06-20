
     ! characteristics based on u
     subroutine ppm(s,qd_l1,qd_h1,u,cspd,Ip,Im,ilo,ihi,dx,dt)
       
       use meth_params_module, only : ppm_type

       implicit none
       
       integer          qd_l1,qd_h1
       integer          ilo,ihi
       double precision s(qd_l1:qd_h1)
       double precision u(qd_l1:qd_h1)
       double precision cspd(qd_l1:qd_h1)
       double precision Ip(ilo-1:ihi+1,1:3)
       double precision Im(ilo-1:ihi+1,1:3)
       double precision dx,dt

       ! local
       integer i
       logical extremum, bigp, bigm

       double precision dsl, dsr, dsc, D2, D2C, D2L, D2R, D2LIM, C, alphap, alpham
       double precision sgn, sigma, s6, w0cc, amax, delam, delap
       double precision dafacem, dafacep, dabarm, dabarp, dafacemin, dabarmin, dachkm, dachkp

       ! s_{\ib,+}, s_{\ib,-}
       double precision, allocatable :: sp(:)
       double precision, allocatable :: sm(:)

       ! \delta s_{\ib}^{vL}
       double precision, allocatable :: dsvl(:)

       ! s_{i+\half}^{H.O.}
       double precision, allocatable :: sedge(:)

       ! cell-centered indexing
       allocate(sp(ilo-1:ihi+1))
       allocate(sm(ilo-1:ihi+1))

       ! constant used in Colella 2008
       C = 1.25d0

       ! cell-centered indexing w/extra x-ghost cell
       allocate(dsvl(ilo-2:ihi+2))

       ! edge-centered indexing for x-faces
       if (ppm_type .eq. 1) then
          allocate(sedge(ilo-1:ihi+2))
       else
          allocate(sedge(ilo-2:ihi+3))
       end if

    ! compute s at x-edges
    if (ppm_type .eq. 1) then

       ! compute van Leer slopes in x-direction
       dsvl = 0.d0
       do i=ilo-2,ihi+2
          dsc = 0.5d0 * (s(i+1) - s(i-1))
          dsl = 2.d0  * (s(i  ) - s(i-1))
          dsr = 2.d0  * (s(i+1) - s(i  ))
          if (dsl*dsr .gt. 0.d0) dsvl(i) = sign(1.d0,dsc)*min(abs(dsc),abs(dsl),abs(dsr))
       end do
       
       ! interpolate s to x-edges
       do i=ilo-1,ihi+2
          sedge(i) = 0.5d0*(s(i)+s(i-1)) - (1.d0/6.d0)*(dsvl(i)-dsvl(i-1))
          ! make sure sedge lies in between adjacent cell-centered values
          sedge(i) = max(sedge(i),min(s(i),s(i-1)))
          sedge(i) = min(sedge(i),max(s(i),s(i-1)))
       end do

       ! copy sedge into sp and sm
       do i=ilo-1,ihi+1
          sp(i) = sedge(i+1)
          sm(i) = sedge(i  )
       end do

       ! modify using quadratic limiters
       do i=ilo-1,ihi+1
          if ((sp(i)-s(i))*(s(i)-sm(i)) .le. 0.d0) then
             sp(i) = s(i)
             sm(i) = s(i)
          else if (abs(sp(i)-s(i)) .ge. 2.d0*abs(sm(i)-s(i))) then
             sp(i) = 3.d0*s(i) - 2.d0*sm(i)
          else if (abs(sm(i)-s(i)) .ge. 2.d0*abs(sp(i)-s(i))) then
             sm(i) = 3.d0*s(i) - 2.d0*sp(i)
          end if
       end do

    else if (ppm_type .eq. 2) then
       
       ! interpolate s to x-edges
       do i=ilo-2,ihi+3
          sedge(i) = (7.d0/12.d0)*(s(i-1)+s(i)) - (1.d0/12.d0)*(s(i-2)+s(i+1))
          ! limit sedge
          if ((sedge(i)-s(i-1))*(s(i)-sedge(i)) .lt. 0.d0) then
             D2  = 3.d0*(s(i-1)-2.d0*sedge(i)+s(i))
             D2L = s(i-2)-2.d0*s(i-1)+s(i)
             D2R = s(i-1)-2.d0*s(i)+s(i+1)
             sgn = sign(1.d0,D2)
             D2LIM = sgn*max(min(C*sgn*D2L,C*sgn*D2R,sgn*D2),0.d0)
             sedge(i) = 0.5d0*(s(i-1)+s(i)) - (1.d0/6.d0)*D2LIM
          end if
       end do

       ! use Colella 2008 limiters
       ! This is a new version of the algorithm 
       ! to eliminate sensitivity to roundoff.
       do i=ilo-1,ihi+1

          alphap = sedge(i+1)-s(i)
          alpham = sedge(i  )-s(i)
          bigp = abs(alphap).gt.2.d0*abs(alpham)
          bigm = abs(alpham).gt.2.d0*abs(alphap)
          extremum = .false.

          if (alpham*alphap .ge. 0.d0) then
             extremum = .true.
          else if (bigp .or. bigm) then
             ! Possible extremum. We look at cell centered values and face
             ! centered values for a change in sign in the differences adjacent to
             ! the cell. We use the pair of differences whose minimum magnitude is the
             ! largest, and thus least susceptible to sensitivity to roundoff.
             dafacem = sedge(i) - sedge(i-1)
             dafacep = sedge(i+2) - sedge(i+1)
             dabarm = s(i) - s(i-1)
             dabarp = s(i+1) - s(i)
             dafacemin = min(abs(dafacem),abs(dafacep))
             dabarmin= min(abs(dabarm),abs(dabarp))
             if (dafacemin.ge.dabarmin) then
                dachkm = dafacem
                dachkp = dafacep
             else
                dachkm = dabarm
                dachkp = dabarp
             endif
             extremum = (dachkm*dachkp .le. 0.d0)
          end if

          if (extremum) then
             D2  = 6.d0*(alpham + alphap)
             D2L = s(i-2)-2.d0*s(i-1)+s(i)
             D2R = s(i)-2.d0*s(i+1)+s(i+2)
             D2C = s(i-1)-2.d0*s(i)+s(i+1)
             sgn = sign(1.d0,D2)
             D2LIM = max(min(sgn*D2,C*sgn*D2L,C*sgn*D2R,C*sgn*D2C),0.d0)
             alpham = alpham*D2LIM/max(abs(D2),1.d-10)
             alphap = alphap*D2LIM/max(abs(D2),1.d-10)
          else
             if (bigp) then
                sgn = sign(1.d0,alpham)
                amax = -alphap**2 / (4*(alpham + alphap))
                delam = s(i-1) - s(i)
                if (sgn*amax .ge. sgn*delam) then
                   if (sgn*(delam - alpham).ge.1.d-10) then
                      alphap = (-2.d0*delam - 2.d0*sgn*sqrt(delam**2 - delam*alpham))
                   else 
                      alphap = -2.d0*alpham
                   endif
                endif
             end if
             if (bigm) then
                sgn = sign(1.d0,alphap)
                amax = -alpham**2 / (4*(alpham + alphap))
                delap = s(i+1) - s(i)
               if (sgn*amax .ge. sgn*delap) then
                  if (sgn*(delap - alphap).ge.1.d-10) then
                     alpham = (-2.d0*delap - 2.d0*sgn*sqrt(delap**2 - delap*alphap))
                  else
                     alpham = -2.d0*alphap
                  endif
               endif
             end if
          end if

          sm(i) = s(i) + alpham
          sp(i) = s(i) + alphap

       end do

    end if

        ! compute x-component of Ip and Im
       do i=ilo-1,ihi+1
          s6 = 6.0d0*s(i) - 3.0d0*(sm(i)+sp(i))
          sigma = abs(u(i)-cspd(i))*dt/dx
          Ip(i,1) = sp(i) - &
               (sigma/2.0d0)*(sp(i)-sm(i)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,1) = sm(i) + &
               (sigma/2.0d0)*(sp(i)-sm(i)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i))*dt/dx
          Ip(i,2) = sp(i) - &
               (sigma/2.0d0)*(sp(i)-sm(i)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,2) = sm(i) + &
               (sigma/2.0d0)*(sp(i)-sm(i)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          sigma = abs(u(i)+cspd(i))*dt/dx
          Ip(i,3) = sp(i) - &
               (sigma/2.0d0)*(sp(i)-sm(i)-(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
          Im(i,3) = sm(i) + &
               (sigma/2.0d0)*(sp(i)-sm(i)+(1.0d0-(2.0d0/3.0d0)*sigma)*s6)
       end do
       
       deallocate(sedge,dsvl,sp,sm)

     end subroutine ppm

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
