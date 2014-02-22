module CastroCG_module

  implicit none

  private

  public ctoprimcg, umeth1dcg

contains

!
! ::: ---------------------------------------------------------------
! ::: :: UMETH1D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  sound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux      <=  (modify) flux in X direction on X edges
! ::: ----------------------------------------------------------------

      subroutine umeth1dcg(q,c,gamc,csml,flatn,qd_l1,qd_h1, &
                         srcQ,src_l1,src_h1, &
                         grav, gv_l1, gv_h1, &
                         ilo,ihi,dx,dt, &
                         flux ,   fd_l1,   fd_h1, &
                         pgdnv,  pgd_l1,  pgd_h1, &
                         dloga,dloga_l1,dloga_h1)

      use meth_params_module, only : QVAR, NVAR

      implicit none
      integer dloga_l1,dloga_h1
      integer qd_l1,qd_h1
      integer src_l1,src_h1
      integer fd_l1,fd_h1
      integer pgd_l1,pgd_h1
      integer gv_l1,gv_h1
      integer ilo,ihi
      double precision dx, dt
      double precision     q(   qd_l1:qd_h1,QVAR)
      double precision  gamc(   qd_l1:qd_h1)
      double precision flatn(   qd_l1:qd_h1)
      double precision  csml(   qd_l1:qd_h1)
      double precision     c(   qd_l1:qd_h1)
      double precision dloga(dloga_l1:dloga_h1)
      double precision pgdnv(  pgd_l1:pgd_h1)
      double precision  flux(fd_l1   :fd_h1,NVAR)
      double precision  srcQ(src_l1  :src_h1,NVAR)
      double precision  grav(gv_l1   :gv_h1)

!     Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
      allocate ( dq(ilo-1:ihi+1,QVAR))
      allocate ( qm(ilo-1:ihi+1,QVAR+3))
      allocate ( qp(ilo-1:ihi+1,QVAR+3))

!     Trace to edges w/o transverse flux correction terms
      call tracecg(q,dq,c,flatn,qd_l1,qd_h1, &
                 dloga,dloga_l1,dloga_h1, &
                 srcQ,src_l1,src_h1, &
                 grav,gv_l1,gv_h1, &
                 qm,qp,ilo-1,ihi+1, &
                 ilo,ihi,dx,dt)

!     Solve Riemann problem, compute xflux from improved predicted states 
      call cmpflxcg(qm, qp, ilo-1,ihi+1, &
                  flux ,  fd_l1, fd_h1, &
                  pgdnv, pgd_l1,pgd_h1, &
                  gamc, csml,c,qd_l1,qd_h1,ilo,ihi)

      deallocate (dq,qm,qp)

      end subroutine umeth1dcg

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ctoprimcg(lo,hi,uin,uin_l1,uin_h1, &
                         q,c,gamc,csml,flatn,q_l1,q_h1,&
                         src,srcQ,src_l1,src_h1, &
                         courno,dx,dt,ngp,ngf,iflaten)

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
                                     QVAR, QRHO, QU, QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
                                     QGAMC,QGAME, nadv, small_temp, allow_negative_energy
      use advection_module, only : uflaten
      use bl_constants_module

      implicit none

      double precision, parameter:: small = 1.d-8

!     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
!     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
!     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
!     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
!     routine that computes flatn).  

      integer lo(1), hi(1)
      integer uin_l1,uin_h1
      integer q_l1,q_h1
      integer src_l1,src_h1
      double precision   uin(uin_l1:uin_h1,NVAR)
      double precision     q(  q_l1:  q_h1,QVAR)
      double precision     c(  q_l1:  q_h1)
      double precision  gamc(  q_l1:  q_h1)
      double precision  csml(  q_l1:  q_h1)
      double precision flatn(  q_l1:  q_h1)
      double precision   src(src_l1:src_h1,NVAR)
      double precision  srcQ(src_l1:src_h1,QVAR)
      double precision dx, dt, courno
      integer iflaten

      integer i
      double precision eken, courx, courmx
      integer ngp, ngf, loq(1), hiq(1)
      integer iadv, n, nq
      integer ispec, nqs
      integer iaux

      double precision, allocatable :: dpdrho(:), dpde(:)

      type (eos_t) :: eos_state

      loq(1) = lo(1)-ngp
      hiq(1) = hi(1)+ngp

      allocate(dpdrho(q_l1:q_h1))
      allocate(dpde  (q_l1:q_h1))

!     Make q (all but p), except put e in slot for rho.e, fix after eos call
      do i = loq(1),hiq(1)
         q(i,QRHO) = uin(i,URHO)
         q(i,QU) = uin(i,UMX)/uin(i,URHO)
         eken = HALF*q(i,QU)**2
!        q(i,QREINT) = uin(i,UEDEN)/q(i,QRHO) - eken
         q(i,QREINT) = uin(i,UEINT)/q(i,QRHO)
         q(i,QTEMP  ) = uin(i,UTEMP)
      enddo

!     Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nqs = QFA + iadv - 1
         q(loq(1):hiq(1),nqs) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
      enddo
      
!     Load species, c, into q, assuming they arrived in uin as rho.c
      do ispec = 1, nspec
         n = UFS + ispec - 1
         nqs = QFS + ispec - 1
         q(loq(1):hiq(1),nqs) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
      enddo

!     Load auxiliary variables which are needed in the EOS
      do iaux = 1, naux
         n = UFX + iaux - 1
         nqs = QFX + iaux - 1
         q(loq(1):hiq(1),nqs) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
      enddo

!     Get gamc, p, c, csml using q state
      do i = loq(1), hiq(1)
         
         eos_state % T   = q(i,QTEMP)
         eos_state % rho = q(i,QRHO)
         eos_state % xn  = q(i,QFS:QFS+nspec-1)

         ! If necessary, reset the energy using small_temp
         if (allow_negative_energy .eq. 0 .and. q(i,QREINT) .le. ZERO) then
            q(i,QTEMP) = small_temp
            eos_state % T = q(i,QTEMP)
            
            call eos(eos_input_rt, eos_state)
            q(i,QREINT) = eos_state % e

            if (q(i,QREINT) .lt. ZERO) then
               print *,'NEW E NEGATIVE IN CTOPRIM ',i,q(i,QREINT)
               call bl_error("Error:: CastroCG_1d.f90 :: ctoprimcg")
            end if
         end if

         eos_state % e = q(i,QREINT)

         call eos(eos_input_re, eos_state)

         q(i,QTEMP)  = eos_state % T
         q(i,QREINT) = eos_state % e
         q(i,QPRES)  = eos_state % p

         gamc(i)     = eos_state % gam1
         c(i)        = eos_state % cs
         dpdrho(i)   = eos_state % dpdr_e
         dpde(i)     = eos_state % dpde

         q(i,QGAMC)  = eos_state % gam1
         q(i,QGAME)  = ONE + q(i,QPRES) / (q(i,QRHO)*Q(i,QREINT))

         csml(i) = max(small, small * c(i))
      end do

      ! compute srcQ terms
      do i = lo(1)-1, hi(1)+1
         srcQ(i,QRHO   ) = src(i,URHO)
         srcQ(i,QU     ) = src(i,UMX) / q(i,QRHO)
         srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QRHO) * src(i,UMX)
         srcQ(i,QPRES  ) =    dpde(i)   * (src(i,UEDEN) - q(i,QU)*src(i,UMX)) / q(i,QRHO) &
                           +dpdrho(i) *  src(i,URHO)

         srcQ(i,QFS:QFS+nspec-1) = ZERO
         srcQ(i,QFX:QFX+naux -1) = ZERO

      end do

!     Make this "rho e" instead of "e"
      do i = loq(1),hiq(1)
         q(i,QREINT ) = q(i,QREINT )*q(i,QRHO)
      enddo

!     Compute Courant Number on grid interior, assume input courno reasonable
      courmx = courno
      do i = lo(1),hi(1)
         courx = ( c(i)+abs(q(i,QU)) ) * dt/dx
         courmx = max( courmx, courx)
      enddo
      courno = courmx

!     Compute flattening coef for slope calculations
      if (iflaten.eq.1) then
         loq(1)=lo(1)-ngf
         hiq(1)=hi(1)+ngf
         call uflaten(loq,hiq, &
                      q(q_l1,QPRES), &
                      q(q_l1,QU), &
                      flatn,q_l1,q_h1)


      endif
      end subroutine ctoprimcg

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
      subroutine tracecg(q,dq,c,flatn,qd_l1,qd_h1, &
                       dloga,dloga_l1,dloga_h1, &
                       srcQ,src_l1,src_h1,&
                       grav,gv_l1,gv_h1, &
                       qxm,qxp,qpd_l1,qpd_h1, &
                       ilo,ihi,dx,dt)

      use network, only : nspec, naux
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QREINT, QPRES, QFA, QFS, QFX, & 
                                     QGAMC,QGAME,  &
                                     nadv, small_dens
      use slope_module, only : uslope
      implicit none

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
      double precision  qxm( qpd_l1: qpd_h1,QVAR+3)
      double precision  qxp( qpd_l1: qpd_h1,QVAR+3)
      double precision grav(  gv_l1:  gv_h1)

!     Local variables
      integer i
      double precision hdt,dtdx
      double precision cc, csq, rho, u, p, rhoe
      double precision drho, du, dp, drhoe
      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision apleft, amleft, azrleft, azeleft

      integer n, iadv
      integer ns, ispec, iaux
      double precision sourcr,sourcp,source,courn,eta,dlogatmp
      double precision clag,t1,t2,t3,dgdlp,drght,drght0,dd,dd0,rtilde
      double precision rtr0r,dleft,dleft0,gam,dgam,ptr0r,rtr0l,ptr0l
      double precision gmtlde, betapr,beta0r, betaml, beta0l
      double precision gmin,gmax

       hdt = 0.5d0 * dt
      dtdx = dt/dx
 
!     Compute slopes
      if (iorder .eq. 1) then
         do n=1,QVAR
            do i = ilo-1, ihi+1 
               dq(i,n) = 0.d0
            enddo
         enddo
      else
         call uslope(q,flatn,qd_l1,qd_h1, &
                     dq,qpd_l1,qpd_h1, &
                     ilo,ihi,QVAR)

!        call pslope(q(:,QPRES),q(:,QRHO), &
!                    flatn      , qd_l1, qd_h1, &
!                    dq(:,QPRES),qpd_l1,qpd_h1, &
!                    grav       , gv_l1, gv_h1, &
!                    ilo,ihi,dx)
      endif
      
!     Compute left and right traced states
      do i = ilo-1, ihi+1 
            cc = c(i)
            csq = cc**2
            rho = q(i,QRHO)
            u = q(i,QU)
            p = q(i,QPRES)
            gam = q(i,QGAME)
            rhoe = q(i,QREINT)

            drho  = dq(i,QRHO)
            du    = dq(i,QU)
            dp    = dq(i,QPRES) 
            dgam  = dq(i,QGAME)
            drhoe = dq(i,QREINT)

            clag = rho*cc
            t1 = 1.d0/clag
            t2 = 0.5*dtdx/rho
            t3 = t1**2

            dgdlp = 2.d0*(1.d0-q(i,QGAME)/q(i,QGAMC))*(q(i,QGAME)-1.d0)

            qxp(i,QVAR+1) = dgdlp
            if(i.le.ihi)then
               qxm(i+1,QVAR+1) = dgdlp
            endif

            drght = -.5d0 - 0.5d0* min(0.d0,u - cc)*dtdx
            drght0 =  - .5*(1.+u*dtdx)
            dd = drght
            dd0 = drght0
            rtilde = rho + drght*drho
            rtr0r = rho + drght0*drho
            ptr0r = p + drght0*dp
            gmtlde = q(i,QGAME) + drght*dgam

            betapr =-(dp*t1+du)*t2
            beta0r = (dp*t3 -drho/(rtilde*rtr0r)) *0.5d0*dtdx*cc

            if (u+cc .ge. 0.0 ) then
                betapr = 0.
            endif
            if (u .ge. 0.0 ) then
                beta0r = 0.
            endif

            qxp(i,QRHO) = 1./(1./rtilde - beta0r - betapr)
            qxp(i,QPRES) = p + drght*dp +   clag*clag*betapr
            qxp(i,QU) = u + drght*du + clag*betapr
            qxp(i,QGAME) = gam + drght0*dgam + dgdlp*(qxp(i,QPRES) - ptr0r)/(qxp(i,QPRES) + ptr0r)
            if (u .ge. 0.0 ) then
              qxp(i,QGAME) = gmtlde
            endif

            enth = ( (rhoe+p)/rho )/csq
            alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
            alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
            alpha0r = drho - dp/csq
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

            qxp(i,QREINT) = rhoe + (apright + amright)*enth*csq + azeright

!           *****************************************************************************

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

            dleft = .5d0 - max(0.d0,u + cc)*0.5d0*dtdx
            dleft0 = .5d0*(1.d0-u*dtdx)
            dd0 = dleft0
            dd = dleft
            rtilde = rho + dleft*drho
            gmtlde = gam + dleft*dgam
            ptr0l = p + dleft0*dp
            rtr0l = rho + dleft0*drho
            betaml =  (dp*t1-du)*t2
            beta0l =-(dp*t3 - drho/(rtilde*rtr0l))*0.5d0*dtdx*cc

            if (u-cc .lt. 0.0 ) then
              betaml = 0.
            endif
            if (u .lt. 0.0 ) then
              beta0l = 0.
            endif

            if (i .le. ihi) then
               qxm(i+1,QRHO) = 1./(1./rtilde - beta0l  - betaml)
               qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
               qxm(i+1,QPRES) = p + dleft*dp + clag*clag*betaml
               qxm(i+1,QU) = u + dleft*du - clag*betaml
               qxm(i+1,QGAME) = gam + dleft0*dgam+dgdlp*(qxm(i+1,QPRES) - ptr0l)/(qxm(i+1,QPRES) + ptr0l)
               if (u .lt. 0.0 ) then
                 qxm(i+1,QGAME) = gmtlde
               endif

              qxm(i+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft

              ! add non-gravitational source term
              qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + hdt*srcQ(i,QRHO)
              qxm(i+1,QU    ) = qxm(i+1,QU    ) + hdt*srcQ(i,QU)
              qxm(i+1,QREINT) = qxm(i+1,QREINT) + hdt*srcQ(i,QREINT)
              qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + hdt*srcQ(i,QPRES)

              ! add gravitational source term
              qxm(i+1,QU) = qxm(i+1,QU) + hdt*grav(i)

            end if

            ! add non-gravitational source term
            qxp(i  ,QRHO  ) = qxp(i,QRHO  ) + hdt*srcQ(i,QRHO)
            qxp(i  ,QU    ) = qxp(i,QU    ) + hdt*srcQ(i,QU)
            qxp(i  ,QREINT) = qxp(i,QREINT) + hdt*srcQ(i,QREINT)
            qxp(i  ,QPRES ) = qxp(i,QPRES ) + hdt*srcQ(i,QPRES)

            ! add gravitational source term
            qxp(i  ,QU) = qxp(i,QU) + hdt*grav(i)

            if(dloga(i).ne.0)then
               courn = dtdx*(cc+abs(u))
               eta = (1.d0-courn)/(cc*dt*abs(dloga(i)))
               dlogatmp = min(eta,1.d0)*dloga(i)
               sourcr = -0.5d0*dt*rho*dlogatmp*u
               sourcp = sourcr*csq
               source = sourcp*enth
               if (i .le. ihi) then
                 qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
                 qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
                 qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
                 qxm(i+1,QGAME) =  qxm(i+1,QGAME) +  dgdlp*sourcp/(qxm(i+1,QPRES)-.5*sourcp)

               end if
               qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
               qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
               qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
               qxp(i,QGAME) =  qxp(i,QGAME) +  dgdlp*sourcp/(qxp(i,QPRES)-.5*sourcp)
            endif

            gmin = min(q(i-2,QGAME), q(i-1,QGAME), q(i,QGAME), q(i+1,QGAME))
            gmax = max(q(i-2,QGAME), q(i-1,QGAME), q(i,QGAME), q(i+1,QGAME))
            qxp(i,QGAME) = min(gmax,max(gmin,qxp(i,QGAME)))
            qxp(i,QVAR+2) = gmin
            qxp(i,QVAR+3) = gmax
            if (i.le.ihi)then
              gmin = min(q(i+2,QGAME), q(i-1,QGAME), q(i,QGAME), q(i+1,QGAME))
              gmax = max(q(i+2,QGAME), q(i-1,QGAME), q(i,QGAME), q(i+1,QGAME))
              qxm(i+1,QGAME) = min(gmax,max(gmin,qxm(i+1,QGAME)))
              qxm(i+1,QVAR+2) = gmin
              qxm(i+1,QVAR+3) = gmax
            endif
      enddo

      do iadv = 1, nadv
         n = QFA + iadv - 1

            ! Right state
            do i = ilo,ihi+1
               u = q(i,QU)
               cc = c(i)
               drght = -.5d0 - 0.5d0* min(0.d0,u - cc)*dtdx
               drght0 =  - .5*(1.+u*dtdx)
               qxp(i,n) = q(i,n) + drght0*dq(i,n)
               if (drght0.lt.-0.5d0) then
                 qxp(i,n) = q(i,n) + drght*dq(i,n)
               endif
            enddo
               
            ! Left state
            do i = ilo-1,ihi
               u = q(i,QU)
               cc = c(i)
               dleft = .5d0 - max(0.d0,u + cc)*0.5d0*dtdx
               dleft0 = .5d0*(1.d0-u*dtdx)
               qxm(i+1,n) = q(i,n) + dleft0*dq(i,n)
               if(dleft0.gt.0.5)then
                 qxm(i+1,n) = q(i,n) + dleft*dq(i,n)
               endif
            enddo
      enddo

      do ispec = 1, nspec
         ns = QFS + ispec - 1

           ! Right state
            do i = ilo,ihi+1
               u = q(i,QU)
               cc = c(i)
               drght = -.5d0 - 0.5d0* min(0.d0,u - cc)*dtdx
               drght0 =  - .5*(1.+u*dtdx)
               qxp(i,ns) = q(i,ns) + drght0*dq(i,ns)
               if (drght0.lt.-0.5d0) then
                 qxp(i,ns) = q(i,ns) + drght*dq(i,ns)
               endif
            enddo
               
            ! Left state
            do i = ilo-1,ihi
               u = q(i,QU)
               cc = c(i)
               dleft = .5d0 - max(0.d0,u + cc)*0.5d0*dtdx
               dleft0 = .5d0*(1.d0-u*dtdx)
               qxm(i+1,ns) = q(i,ns) + dleft0*dq(i,ns)
               if(dleft0.gt.0.5)then
                 qxm(i+1,ns) = q(i,ns) + dleft*dq(i,ns)
               endif
            enddo
      enddo

      do iaux = 1, naux
         ns = QFX + iaux - 1

           ! Right state
            do i = ilo,ihi+1
               u = q(i,QU)
               cc = c(i)
               drght = -.5d0 - 0.5d0* min(0.d0,u - cc)*dtdx
               drght0 =  - .5*(1.+u*dtdx)
               qxp(i,ns) = q(i,ns) + drght0*dq(i,ns)
               if (drght0.lt.-0.5d0) then
                 qxp(i,ns) = q(i,ns) + drght*dq(i,ns)
               endif
            enddo
               
            ! Left state
            do i = ilo-1,ihi
               u = q(i,QU)
               cc = c(i)
               dleft = .5d0 - max(0.d0,u + cc)*0.5d0*dtdx
               dleft0 = .5d0*(1.d0-u*dtdx)
               qxm(i+1,ns) = q(i,ns) + dleft0*dq(i,ns)
               if(dleft0.gt.0.5)then
                 qxm(i+1,ns) = q(i,ns) + dleft*dq(i,ns)
               endif
            enddo
      enddo


      end subroutine tracecg

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine cmpflxcg(qm,qp,qpd_l1,qpd_h1, &
                          flx,flx_l1,flx_h1, &
                          pgdnv,pg_l1,pg_h1, &
                          gamc,csml,c,qd_l1,qd_h1,ilo,ihi)

      use meth_params_module, only : QVAR, NVAR

      implicit none
      integer ilo,ihi
      integer qpd_l1,qpd_h1
      integer flx_l1, flx_h1
      integer  pg_l1, pg_h1
      integer  qd_l1,  qd_h1
      double precision    qm(qpd_l1:qpd_h1, QVAR+3)
      double precision    qp(qpd_l1:qpd_h1, QVAR+3)
      double precision   flx(flx_l1:flx_h1, NVAR)
      double precision pgdnv( pg_l1: pg_h1)
      double precision  gamc( qd_l1: qd_h1)
      double precision     c( qd_l1: qd_h1)
      double precision  csml( qd_l1: qd_h1)

!     Local variables
      integer i
      double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:), ugdnv(:)
      
      allocate ( smallc(ilo:ihi+1) )
      allocate ( cavg(ilo:ihi+1) )
      allocate ( gamcp(ilo:ihi+1) )
      allocate ( gamcm(ilo:ihi+1) )
      allocate ( ugdnv(ilo:ihi+1) )

      do i = ilo, ihi+1 
          smallc(i) = max( csml(i), csml(i-1) )
          cavg(i) = 0.5d0*( c(i) + c(i-1) )
          gamcm(i) = gamc(i-1)
          gamcp(i) = gamc(i)
      enddo

!     Solve Riemann problem (gdnv state passed back, but only (u,p) saved)
!     call riemannus(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
!                    gamcm, gamcp, flx, flx_l1, flx_h1, ugdnv,  &
!                    pgdnv, pg_l1, pg_h1, ilo, ihi )

      call riemanncg(qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
                     gamcm, gamcp, flx, flx_l1, flx_h1, ugdnv,  &
                     pgdnv, pg_l1, pg_h1, ilo, ihi )

      deallocate (smallc,cavg,gamcp,gamcm,ugdnv)

      end subroutine cmpflxcg

! ::: ------------------------------------------------------------------
! ::: ------------------------------------------------------------------
! ::: 

      subroutine riemanncg(ql,qr,qpd_l1,qpd_h1,smallc,cav, &
                           gamcl,gamcr,uflx,uflx_l1,uflx_h1,ugdnv, &
                           pgdnv,pg_l1,pg_h1,ilo,ihi)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QPRES, QREINT, QFA, QFS, QFX, &
                                     QGAMC,QGAME, &
                                     URHO, UMX, UEDEN, UFA, UFS, UFX, nadv, small_dens, small_pres

      implicit none

      double precision, parameter:: small = 1.d-8

      integer ilo,ihi
      integer  qpd_l1,  qpd_h1
      integer   pg_l1,   pg_h1
      integer uflx_l1, uflx_h1
      double precision ql(qpd_l1:qpd_h1, QVAR+3)
      double precision qr(qpd_l1:qpd_h1, QVAR+3)
      double precision   cav(ilo:ihi+1), smallc(ilo:ihi+1)
      double precision gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
      double precision  uflx(uflx_l1:uflx_h1, NVAR)
      double precision ugdnv(ilo:ihi+1)
      double precision pgdnv( pg_l1: pg_h1)

      double precision rgdnv, regdnv, ustar
      double precision rl, ul, pl, rel
      double precision rr, ur, pr, rer
      double precision wl, wr, scr
      double precision rstar, cstar, estar, pstar
      double precision ro, uo, po, reo, co, gamco, entho
      double precision sgnm, spin, spout, ushock, frac

      double precision wsmall, csmall
      double precision weakwv
        
      double precision taur,taul,tauo,clsql,clsqr,gambarl,gambarr
      double precision pstnm1,gamstar,wlsq,wrsq
      double precision ustarp, zp,dpditer, zm, gmin,gmax
      double precision ustarm, ustnm1,ustnp1, denom, gameo,gdot
      double precision wo,wosq, dpjmp,clsq,ushok,gamgdnv
      double precision gdotl,gdotr
      
      integer iadv, n, nq
      integer ispec, iaux, ns, nqs
      integer k
      integer idebug
      integer itno, iter

!     Solve Riemann Problem
!     NOTE: The calling routine will order velocity unknowns so that
!     for the purposes of this routine, the normal component is always
!     loaded in the QU slot.


      data itno/3/
      data weakwv/1.d-3/

      do k = ilo, ihi+1
         idebug = 1
         rl  = max( ql(k,QRHO), small_dens)
         ul  = ql(k,QU)
         pl  = max( ql(k,QPRES), small_pres)
         rel = ql(k,QREINT)
         rr  = max( qr(k,QRHO), small_dens)
         ur  = qr(k,QU)
         pr  = max (qr(k,QPRES), small_pres)
         rer = qr(k,QREINT)
         taur = 1.d0/rr
         taul = 1.d0/rl
         clsql = gamcl(k)*pl*rl
         clsqr = gamcr(k)*pr*rr
         gambarl = ql(k,QGAME)
         gambarr = qr(k,QGAME)
         gmin = min(ql(k,QVAR+2),qr(k,QVAR+2))
         gmax = max(ql(k,QVAR+3),qr(k,QVAR+3))
         gdotl = ql(k,QVAR+1)
         gdotr = qr(k,QVAR+1)
         
         csmall = smallc(k)
         wsmall = small_dens*csmall
         wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
         wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))
         pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
         pstar = max(pstar,small_pres)
         pstnm1 = pstar
         call wsqge(pl,taul,gambarl,gdotl,  &
                gamstar,pstar,wlsq,clsql,gmin,gmax)

         call wsqge(pr,taur,gambarr,gdotr,  &
                gamstar,pstar,wrsq,clsqr,gmin,gmax)

          wl = sqrt(wlsq)
          wr = sqrt(wrsq)
          ustarp = ul - (pstar-pl)/wl
          ustarm = ur + (pstar-pr)/wr
          pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
          pstar = max(pstar,small_pres)

          do iter =1,itno

            call wsqge(pl,taul,gambarl,gdotl,  &
                gamstar,pstar,wlsq,clsql,gmin,gmax)
            call wsqge(pr,taur,gambarr,gdotr,  &
                gamstar,pstar,wrsq,clsqr,gmin,gmax)

            wl = 1.d0 / sqrt(wlsq)
            wr = 1.d0 / sqrt(wrsq)
            ustnm1 = ustarm
            ustnp1 = ustarp
            ustarm = ur-(pr-pstar)*wr
            ustarp = ul+(pl-pstar)*wl
 
            dpditer=abs(pstnm1-pstar)
            zp=abs(ustarp-ustnp1)
            if(zp-weakwv*cav(k).le.0.d0)then
              zp = dpditer*wl
            endif
               
            zm=abs(ustarm-ustnm1)
            if(zm-weakwv*cav(k).le.0.d0)then
              zm = dpditer*wr
            endif
 
            denom=dpditer/max(zp+zm,small*(cav(k)))
            pstnm1 = pstar
            pstar=pstar-denom*(ustarm-ustarp)
            pstar=max(pstar,small_pres)
 
          enddo


           ustar = 0.5d0* ( ustarp + ustarm )

         if (ustar .ge. 0.d0) then
            ro = rl
            uo = ul
            po = pl
            tauo = taul
            reo = rel
            gamco = gamcl(k)
            gameo = gambarl
            gdot = ql(k,QVAR+1)
         else
            ro = rr
            uo = ur
            po = pr
            tauo = taur
            reo = rer
            gamco = gamcr(k)
            gameo = gambarr
            gdot = qr(k,QVAR+1)
         endif

         ro = max(small_dens,ro)
         co = sqrt(abs(gamco*po/ro))
         co = max(csmall,co)
         clsq = (co*ro)**2
         call wsqge(po,tauo,gameo,gdot,   &
                 gamstar,pstar,wosq,clsq,gmin,gmax)

         sgnm = sign(1.d0,ustar)

         wo = sqrt(wosq)
         dpjmp = pstar - po
         rstar=max(1.-ro*dpjmp/wosq,(gameo-1.)/(gameo+1.))
         rstar=ro/rstar
         rstar = max(small_dens,rstar)
         ushok = -sgnm*ustar + wo/rstar
         entho = (reo/ro + po/ro)/co**2
         estar = reo + (pstar - po)*entho


         cstar = sqrt(abs(gamco*pstar/rstar))
         cstar = max(cstar,csmall)

         spout = co - sgnm*uo
         spin = cstar - sgnm*ustar

         ushock = 0.5d0*(spin + spout)
         if (pstar-po .ge. 0.d0) then
            spin = ushock
            spout = ushock
         endif
         if (spout-spin .eq. 0.d0) then
            scr = small*cav(k)
         else
            scr = spout-spin
         endif
         frac = (1.d0 + (spout + spin)/scr)*0.5d0
         frac = max(0.d0,min(1.d0,frac))

         rgdnv = frac*rstar + (1.d0 - frac)*ro
         ugdnv(k) = frac*ustar + (1.d0 - frac)*uo
         pgdnv(k) = frac*pstar + (1.d0 - frac)*po
         gamgdnv =  frac*gamstar + (1.d0-frac)*gameo
         regdnv = frac*estar + (1.d0 - frac)*reo
         if (spout .lt. 0.d0) then
            rgdnv = ro
            ugdnv(k) = uo
            pgdnv(k) = po
            regdnv = reo
            gamgdnv = gameo
         endif
         if (spin .ge. 0.d0) then
            rgdnv = rstar
            ugdnv(k) = ustar
            pgdnv(k) = pstar
            regdnv = estar
            gamgdnv = gamstar
         endif

!     Compute fluxes, order as conserved state (not q)
         uflx(k,URHO) = rgdnv*ugdnv(k)
         uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 
!        uflx(k,UEDEN) = ugdnv(k)*pgdnv(k) +    &
!            uflx(k,URHO)*(0.5d0*ugdnv(k)**2 + pgdnv(k)/((gamgdnv-1.d0)))
         uflx(k,UEDEN) = ugdnv(k)*pgdnv(k) +    &
             0.5d0*uflx(k,URHO)*ugdnv(k)**2 + ugdnv(k)*regdnv

         do iadv = 1, nadv
            n = UFA + iadv - 1
            nq = QFA + iadv - 1
            if (ustar .ge. 0.d0) then
               uflx(k,n) = uflx(k,URHO)*ql(k,nq)
            else
               uflx(k,n) = uflx(k,URHO)*qr(k,nq)
            endif
         enddo

         do ispec = 1, nspec
            ns = UFS + ispec - 1
            nqs = QFS + ispec - 1
            if (ustar .ge. 0.d0) then
               uflx(k,ns) = uflx(k,URHO)*ql(k,nqs)
            else
               uflx(k,ns) = uflx(k,URHO)*qr(k,nqs)
            endif
         enddo

         do iaux = 1, naux
            ns = UFX + iaux - 1
            nqs = QFX + iaux - 1
            if (ustar .ge. 0.d0) then
               uflx(k,ns) = uflx(k,URHO)*ql(k,nqs)
            else
               uflx(k,ns) = uflx(k,URHO)*qr(k,nqs)
            endif
         enddo

      enddo
      end subroutine riemanncg

      subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)
 
      double precision p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax
      double precision smlp1,small,divide,temp
 
      data smlp1,small/.001,1.e-07/
 
      gstar=(pstar-p)*gdot/(pstar+p) + gam
      gstar=max(gmin,min(gmax,gstar))
         wsq = (0.5d0*(gstar-1.0d0)*(pstar+p)+pstar)
      temp = ((gstar-gam)/(gam-1.))
 
      if (pstar-p.eq.0.0d0) then
         divide=small
      else
         divide=pstar-p
      endif
 
      temp=temp/divide
      wsq = wsq/(v - temp*p*v)
      if (abs(pstar/p-1.d0)-smlp1 .lt. 0.0d0 ) then
        wsq = csq
      endif
      wsq=max(wsq,(.5d0*(gam-1.d0)/gam)*csq)
 
 
        return
      end subroutine wsqge

end module CastroCG_module
