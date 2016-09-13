      subroutine ca_workspace(nThermo, nGrowHyp, nMaxSpec, numadv,
     &     numspec, Density, Xmom, Eden, Temp, FirstAdv, FirstSpec,
     &     gas, poly, stell)
      implicit none
      include 'castro2d.fi'
      integer nThermo, nGrowHyp, nMaxSpec, numadv, numspec
      integer Density, Xmom, Eden, Temp, FirstAdv, FirstSpec
      integer gas, poly, stell
      
c     Set info going from C++ to Fortran
      nadv = numadv
      nspec = numspec
      QVAR = QTHERM + nadv + nspec
      QFA = QTHERM + 1
      QFS = QTHERM + nadv

      NVAR = NTHERM + nadv + nspec
      URHO = Density + 1
      UMX = Xmom + 1
      UMY = Xmom + 2
      UEDEN = Eden + 1
      UTEMP = Temp + 1
      UFA = FirstAdv + 1
      UFS = FirstSpec + 1

      SRHO = 1
      SMX = 2
      SMY = 3
      SEDEN = 4
      SPRES = 5
      SFA = 6
      SFS = SFA + nadv - 1
      SVAR = SPRES + nadv + nspec

c     Make sure that these point to something valid, even if useless
      if (numadv.lt.1)  then
         QFA=QRHO
         UFA=URHO
         SFA=SRHO
      endif
      if (numspec.lt.1) then
         QFS=QRHO
         UFS=URHO
         SFS=SRHO
      endif

      gasType = gas
      polytropic = poly
      stellar = stell
      
c     Set info going from Fortran to C++
      nThermo = NTHERM
      nGrowHyp = NHYP
      nMaxSpec = MAXSPEC
      end

      subroutine ca_umdrv(lo,hi,
     &     uin,uin_l1,uin_l2,uin_h1,uin_h2,
     &     uout,uout_l1,uout_l2,uout_h1,uout_h2,
     &     src,src_l1,src_l2,src_h1,src_h2,delta,dt,
     &     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2,
     &     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2,
     &     courno)
      implicit none
      include 'castro2d.fi'
      integer lo(2),hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
      integer src_l1,src_l2,src_h1,src_h2
      double precision   uin(  uin_l1:uin_h1,    uin_l2:uin_h2,  NVAR)
      double precision  uout( uout_l1:uout_h1,  uout_l2:uout_h2, NVAR)
      double precision   src(  src_l1:src_h1,    src_l2:src_h2,  SVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision delta(2),dt,courno

c     Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:)
      double precision, allocatable:: gamc(:,:)
      double precision, allocatable:: flatn(:,:)
      double precision, allocatable:: c(:,:)
      double precision, allocatable:: csml(:,:)
      double precision, allocatable:: div(:,:)

      double precision dx,dy
      integer nx,ny,ngq,ngf,iflaten

      allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,QVAR))
      allocate(  gamc(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate( flatn(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate(     c(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate(  csml(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate(   div(lo(1):hi(1)+1,lo(2):hi(2)+1))


c     NOTE: src ordered to correspond to QVAR variable ordering, where
c     the source term for pressure is in slot 5

      dx = delta(1)
      dy = delta(2)

      nx = hi(1) - lo(1) + 1
      ny = hi(2) - lo(2) + 1


      ngq = NHYP
      ngf = 1
      iflaten = 1

c     Translate to primitive variables, compute sound speeds
c     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
c       and set to correspond to coordinates of (lo:hi)
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2,
     &     q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_h1,uin_h2,
     &     courno,dx,dy,dt,ngq,ngf,iflaten)

c     Compute hyperbolic fluxes using unsplit Godunov
      call umeth2d(q,c,gamc,csml,flatn,src,nx,ny,dx,dy,dt,
     &     flux1,flux2)

c ::: compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi,q,uin_l1,uin_l2,uin_h1,uin_h2,
     &     dx,dy,div,lo(1),lo(2),hi(1)+1,hi(2)+1,
     &     QU,QV)

c ::: artificial viscosity, conservative update
      call consup(uin,uin_l1,uin_l2,uin_h1,uin_h2,
     &     uout,uout_l1,uout_l2,uout_h1,uout_h2,
     &     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2,
     &     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2,
     &     div,lo,hi,dx,dy,dt)

      end


c ::: ---------------------------------------------------------------
c ::: :: UMETH2D     Compute hyperbolic fluxes using unsplit second
c ::: ::               order Godunov integrator for polytropic gas dynamics.
c ::: :: 
c ::: :: inputs/outputs
c ::: :: q           => (const)  input state, primitives
c ::: :: c           => (const)  sound speed
c ::: :: gamc        => (const)  cound speed gamma
c ::: :: csml        => (const)  local small c val
c ::: :: flatn       => (const)  flattening parameter
c ::: :: src         => (const)  source
c ::: :: nx          => (const)  number of cells in X direction
c ::: :: ny          => (const)  number of cells in Y direction
c ::: :: dx          => (const)  grid spacing in X direction
c ::: :: dy          => (const)  grid spacing in Y direction
c ::: :: dt          => (const)  time stepsize
c ::: :: flux1      <=  (modify) flux in X direction on X edges
c ::: :: flux2      <=  (modify) flux in Y direction on Y edges
c ::: ----------------------------------------------------------------
      subroutine umeth2d(q,c,gamc,csml,flatn,src,nx,ny,dx,dy,dt,
     &     flux1,flux2)
      implicit none
      include 'castro2d.fi'

      integer nx, ny
      double precision dx, dy, dt
      double precision     q(1-NHYP:nx+NHYP,1-NHYP:ny+NHYP,QVAR)
      double precision  gamc(1-NHYP:nx+NHYP,1-NHYP:ny+NHYP)
      double precision flatn(1-NHYP:nx+NHYP,1-NHYP:ny+NHYP)
      double precision  csml(1-NHYP:nx+NHYP,1-NHYP:ny+NHYP)
      double precision     c(1-NHYP:nx+NHYP,1-NHYP:ny+NHYP)
      double precision flux1(nx+1,ny,NVAR)
      double precision flux2(nx,ny+1,NVAR)
      double precision src(0:nx+1,0:ny+1,SVAR)

c ::: left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:,:),  qm(:,:,:),   qp(:,:,:)
      double precision, allocatable::qxm(:,:,:),qym(:,:,:)
      double precision, allocatable::qxp(:,:,:),qyp(:,:,:)

c     Work arrays to hold 3 planes of riemann state and conservative fluxes
      double precision, allocatable::   fx(:,:,:),  fy(:,:,:)
      double precision, allocatable:: xtmp(:,:,:),ytmp(:,:,:)
      double precision, allocatable::  upx(:,:,:), upy(:,:,:)

c ::: local scalar variables
      integer ig1, ig2, jg1, jg2
      integer is1, is2, js1, js2
      double precision dtdx
      double precision hdtdx, hdt

      integer j

      allocate ( dq(0:nx+1,0:ny+1,QVAR))
      allocate ( qm(0:nx+1,0:ny+1,QVAR))
      allocate ( qp(0:nx+1,0:ny+1,QVAR))
      allocate (qxm(0:nx+1,0:ny+1,QVAR))
      allocate (qxp(0:nx+1,0:ny+1,QVAR))
      allocate (qym(0:nx+1,0:ny+1,QVAR))
      allocate (qyp(0:nx+1,0:ny+1,QVAR))

      allocate (  fx(1:nx+1,0:ny+1,NVAR))
      allocate (  fy(0:nx+1,1:ny+1,NVAR))
      allocate (xtmp(1:nx+1,0:ny+1,2))
      allocate (ytmp(0:nx+1,1:ny+1,2))
      allocate ( upx(1:nx+1,1:ny,2))
      allocate ( upy(1:nx,1:ny+1,2))

      ig1 = 1-NHYP
      ig2 = nx+NHYP
      jg1 = 1-NHYP
      jg2 = ny+NHYP

      is1 = 0
      is2 = nx+1
      js1 = 0
      js2 = ny+1

c ::: local constants
      dtdx = dt/dx
      hdtdx = 0.5d0*dtdx
      hdt = 0.5d0*dt

c     NOTE: Geometry terms need to be punched through

c     Trace to edges w/o transverse flux correction terms
      call trace(q,dq,c,flatn,gamc,nx,ny,qxm,qxp,qym,qyp,dx,dy,dt)

c     Solve Riemann problem, compute x fluxes from initial traces
      call cmpflx(qxm, 1, nx+2, 0, ny+1, qxp, 0, nx+1, 0, ny+1,
     &     1, nx+1, 0, ny+1, QU, QV, fx, 1, nx+1, 0, ny+1, xtmp,
     &     1, nx+1, 0, ny+1, q, gamc, csml, c, nx, ny)

c     Solve Riemann problem, compute y fluxes from initial traces
      call cmpflx(qym, 0, nx+1, 1, ny+2, qyp, 0, nx+1, 0, ny+1,
     &     0, nx+1, 1, ny+1, QV, QU, fy, 0, nx+1, 1, ny+1, ytmp,
     &     0, nx+1, 1, ny+1, q, gamc, csml, c, nx, ny)

c     Modify minus x predicted state with dG/dy and source term
      call transy(qxm, 0, nx+1, 0, ny+1, qm, 0, nx+1, 0, ny+1,
     &     fy, ytmp, 0, nx+1, 1, ny+1, gamc, ig1, ig2, jg1, jg2, 
     &     hdt, src, is1, is2, js1, js2, hdtdx, 0, nx, 1, ny)

c     Modify plus x predicted state with dG/dy and source term
      call transy(qxp, 0, nx+1, 0, ny+1, qp, 0, nx+1, 0, ny+1,
     &     fy, ytmp, 0, nx+1, 1, ny+1, gamc, ig1, ig2, jg1, jg2, 
     &     hdt, src, is1, is2, js1, js2, hdtdx, 1, nx+1, 1, ny)

c     Solve Riemann problem, compute xflux from improved predicted states 
      call cmpflx(qm, 1, nx+2, 0, ny+1, qp, 0, nx+1, 0, ny+1,
     &     1, nx+1, 1, ny, QU, QV, flux1, 1, nx+1, 1, ny, upx,
     &     1, nx+1, 1, ny, q, gamc, csml, c, nx, ny)

c     Modify minus y predicted state with dF/dx and source term
      call transx(qym, 0, nx+1, 0, ny+1, qm, 0, nx+1, 0, ny+1,
     &     fx, xtmp, 1, nx+1, 0, ny+1, gamc, ig1, ig2, jg1, jg2, 
     &     hdt, src, is1, is2, js1, js2, hdtdx, 1, nx, 0, ny)

c     Modify plus y predicted state with dF/dx and source term
      call transx(qyp, 0, nx+1, 0, ny+1, qp, 0, nx+1, 0, ny+1,
     &     fx, xtmp, 1, nx+1, 0, ny+1, gamc, ig1, ig2, jg1, jg2, 
     &     hdt, src, is1, is2, js1, js2, hdtdx, 1, nx, 1, ny+1)

c     Solve Riemann problem, compute yflux from improved predicted states 
      call cmpflx(qm, 0, nx+1, 1, ny+2, qp, 0, nx+1, 0, ny+1,
     &     1, nx, 1, ny+1, QV, QU, flux2, 1, nx, 1, ny+1, upy,
     &     1, nx, 1, ny+1, q, gamc, csml, c, nx, ny)

      end

c ::: 
c ::: ------------------------------------------------------------------
c ::: 
      subroutine ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2,
     &     q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2,
     &     courno,dx,dy,dt,ngp,ngf,iflaten)

c     Will give primitive variables (u,p,rho,gam) on lo-ngp:hi+ngp, and flatn
c     on lo-ngf:hi+ngf if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn
c     are given by DIMS(q).  This declared region is assumed to encompass
c     lo-ngp:hi+ngp.  Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is
c     used by the routine that computes flatn).  May reduce courno (=umax*dt/dx)
c     if necessary to ensure CFL stability.
      implicit none
      include 'castro2d.fi'
      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer q_l1,q_l2,q_h1,q_h2
      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision q(q_l1:q_h1,q_l2:q_h2,QVAR)
      double precision c(q_l1:q_h1,q_l2:q_h2)
      double precision gamc(q_l1:q_h1,q_l2:q_h2)
      double precision csml(q_l1:q_h1,q_l2:q_h2)
      double precision flatn(q_l1:q_h1,q_l2:q_h2)
      double precision dx, dy, dt, courno
      integer iflaten

      integer i, j
      double precision eken, courmx, courmy
      integer ngp, ngf, cnt, loq(2), hiq(2)
      integer iadv, n, nq
      integer ispec, ns, nqs

      do i=1,2
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo
      

c     Make q (all but p), compute gamma after eos call
      do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)
            q(i,j,QRHO) = uin(i,j,URHO)
            q(i,j,QU) = uin(i,j,UMX)/uin(i,j,URHO)
            q(i,j,QV) = uin(i,j,UMY)/uin(i,j,URHO)
            eken = 0.5d0*(q(i,j,QU)**2+q(i,j,QV)**2)
            q(i,j,QREINT) = uin(i,j,UEDEN)/q(i,j,QRHO) - eken
            q(i,j,QTEMP) = uin(i,j,UTEMP)
         enddo
      enddo

c     Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
            enddo
         enddo
      enddo
      
c     Load chemical species, c, into q, assuming they arrived in uin as rho.c
      do ispec = 1, nspec
         ns = UFS + ispec - 1
         nqs = QFS + ispec - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nqs) = uin(i,j,ns)/q(i,j,QRHO)
            enddo
         enddo
      enddo

c :::  Get gamc, p, c, csml using q state (where rho.e comp = e)
c     NOTE: Currently does not pass in advected fields
      call gp_eos(
     &     gamc,               q_l1,q_l2,q_h1,q_h2,
     &     q(q_l1,q_l2,QPRES), q_l1,q_l2,q_h1,q_h2,
     &     c,                  q_l1,q_l2,q_h1,q_h2,
     &     csml,               q_l1,q_l2,q_h1,q_h2,
     &     q(q_l1,q_l2,QRHO),  q_l1,q_l2,q_h1,q_h2,
     &     q(q_l1,q_l2,QREINT),q_l1,q_l2,q_h1,q_h2,
     &     q(q_l1,q_l2,QFS),   q_l1,q_l2,q_h1,q_h2,
     &     q(q_l1,q_l2,QTEMP), q_l1,q_l2,q_h1,q_h2,
     &     loq, hiq, nspec, cnt)

c     Set gam=(p/rho.e) + 1  over all q space
      do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)
            q(i,j,QGAM) =
     &           q(i,j,QPRES)/(q(i,j,QREINT)*q(i,j,QRHO)) + 1.d0
         enddo
      enddo

c     Compute Courant Number on grid interior, assume input courno reasonable
      courmx = courno*dx/dt
      courmy = courno*dy/dt
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            courmx = max( courmx, c(i,j)+abs(q(i,j,QU)) )
            courmy = max( courmy, c(i,j)+abs(q(i,j,QV)) )
         enddo
      enddo
      courmx = courmx*dt/dx
      courmy = courmy*dt/dy
      courno = max( courmx, courmy )

c ::: compute flattening coef for slope calculations
      if(iflaten.eq.1)then
         do n=1,2
            loq(n)=lo(n)-ngf
            hiq(n)=hi(n)+ngf
         enddo
         call uflaten(loq,hiq,
     &                q(q_l1,q_l2,QPRES),
     &                q(q_l1,q_l2,QU),
     &                q(q_l1,q_l2,QV),
     &                q(q_l1,q_l2,QRHO),
     &                c,flatn,q_l1,q_l2,q_h1,q_h2)

      endif
      end

c ::: 
c ::: ------------------------------------------------------------------
c ::: 
      subroutine trace(q,dq,c,flatn,gamc,nx,ny,qxm,qxp,qym,qyp,dx,dy,dt)
      implicit none
      include 'castro2d.fi'
      integer nx, ny
      double precision dx, dy, dt
      double precision     q(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP, QVAR)
      double precision  gamc(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)
      double precision flatn(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)
      double precision     c(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)

      double precision  dq(0:nx+1, 0:ny+1, QVAR)
      double precision qxm(0:nx+1, 0:ny+1, QVAR)
      double precision qxp(0:nx+1, 0:ny+1, QVAR)
      double precision qym(0:nx+1, 0:ny+1, QVAR)
      double precision qyp(0:nx+1, 0:ny+1, QVAR)

c :::  declare local variables
      integer i, j
      double precision dtdx, dtdy, dth

      integer n, iadv
      double precision clag,t1,t2,t3,dgdlp,ured,dleft,dleft0,rtilde
      double precision gmtlde,ptr0l,rtr0l,betaml,beta0l,drght,drght0
      double precision rtr0r,ptr0r,betapr,beta0r
      integer ns, ispec
      double precision ascmprght, ascmpleft, ascmpbot, ascmptop

      dth = 0.5d0 * dt

c ::: compute slopes
      if (iorder .eq. 1) then
         do n=1,QVAR
            do j = 0, ny+1 
               do i = 0, nx+1 
                  dq(i,j,n) = 0.d0
               enddo
            enddo
         enddo
      else
         call uslope(q,dq,flatn,nx,ny,QVAR,1)
      endif

      do j = 0, ny+1 
         do i = 0, nx+1 

            clag = c(i,j) * q(i,j,QRHO)
            t1 = 1.d0/clag
            t2 = dth/(q(i,j,QRHO)*dx)
            t3 = 1.d0/(c(i,j)*q(i,j,QRHO))**2
            dgdlp = 2.d0*(1.d0-q(i,j,QGAM)
     &           / gamc(i,j))*(q(i,j,QGAM)-1.d0)

            ured = q(i,j,QU)
            dleft  = .5d0 - max(0.d0,ured+c(i,j))*dth/dx
            dleft0 = .5d0*(1.d0-ured*dt/dx)
            rtilde = q(i,j,QRHO) + dleft*dq(i,j,QRHO)
            gmtlde = q(i,j,QGAM) + dleft*dq(i,j,QGAM)
            ptr0l  = q(i,j,QPRES) + dleft0*dq(i,j,QPRES)
            rtr0l  = q(i,j,QRHO) + dleft0*dq(i,j,QRHO)

            
            if (ured-c(i,j) .lt. 0.d0 ) then
               betaml = 0.d0
            else
               betaml = (dq(i,j,QPRES)*t1-dq(i,j,QU))*t2
            endif

            if (ured .lt. 0.d0 ) then
               beta0l = 0.d0
            else
               beta0l =-(dq(i,j,QPRES)*t3 - dq(i,j,QRHO)
     &              /(rtilde*rtr0l))*dth*c(i,j)/dx
            endif

            qxm(i,j,QRHO) = 1.d0/(1.d0/rtilde - beta0l  - betaml)
            qxm(i,j,QPRES) = q(i,j,QPRES) + dleft*dq(i,j,QPRES) +
     1           clag*clag*betaml
            qxm(i,j,QU) = q(i,j,QU) +  dleft*dq(i,j,QU) - clag*betaml

            if (ured .lt. 0.d0 ) then
               qxm(i,j,QV) = q(i,j,QV) + dleft*dq(i,j,QV)
               qxm(i,j,QGAM)= gmtlde
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qxm(i,j,n) = q(i,j,n) + dleft*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qxm(i,j,ns) = q(i,j,ns) + dleft*dq(i,j,ns)
               enddo
            else
               qxm(i,j,QV) = q(i,j,QV) + dleft0*dq(i,j,QV)
               qxm(i,j,QGAM) = q(i,j,QGAM) + dleft0*dq(i,j,QGAM)
     1              +dgdlp*(qxm(i,j,QPRES)-ptr0l)/(qxm(i,j,QPRES)+ptr0l)
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qxm(i,j,n) = q(i,j,n) + dleft0*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qxm(i,j,ns) = q(i,j,ns) + dleft0*dq(i,j,ns)
               enddo
            endif

            drght = -.5d0 - min(0.d0,ured - c(i,j))*dth/dx
            drght0 = -.5d0*(1.d0+ured*dt/dx)
            rtilde = q(i,j,QRHO) + drght*dq(i,j,QRHO)
            rtr0r = q(i,j,QRHO) + drght0*dq(i,j,QRHO)
            ptr0r = q(i,j,QPRES) + drght0*dq(i,j,QPRES)
            gmtlde = q(i,j,QGAM) + drght*dq(i,j,QGAM)

            if (ured+c(i,j) .ge. 0.0d0 ) then
               betapr = 0.d0
            else
               betapr =-(dq(i,j,QPRES)*t1+dq(i,j,QU))*t2
            endif

            if (ured .ge. 0.d0 ) then
               beta0r = 0.d0
            else
               beta0r = (dq(i,j,QPRES)*t3 -dq(i,j,QRHO)
     &              /(rtilde*rtr0r))*dth*c(i,j)/dx
            endif

            qxp(i,j,QRHO) = 1.d0/(1.d0/rtilde - beta0r - betapr)
            qxp(i,j,QPRES) = q(i,j,QPRES) + drght*dq(i,j,QPRES) +
     1           clag*clag*betapr
            qxp(i,j,QU) = q(i,j,QU) + drght*dq(i,j,QU) + clag*betapr

            if (ured .lt. 0.0 ) then
               qxp(i,j,QV) = q(i,j,QV)+drght0*dq(i,j,QV)
               qxp(i,j,QGAM) = q(i,j,QGAM) + drght0*dq(i,j,QGAM)
     1           + dgdlp*(qxp(i,j,QPRES)-ptr0r)/(qxp(i,j,QPRES)+ptr0r)
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qxp(i,j,n) = q(i,j,n) + drght0*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qxp(i,j,ns) = q(i,j,ns) + drght0*dq(i,j,ns)
               enddo
            else
               qxp(i,j,QV) = q(i,j,QV)+drght*dq(i,j,QV)
               qxp(i,j,QGAM) = gmtlde
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qxp(i,j,n) = q(i,j,n) + drght*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qxp(i,j,ns) = q(i,j,ns) + drght*dq(i,j,ns)
               enddo
            endif

         enddo
      enddo

c ::: compute slopes
      if (iorder .ne. 1) then
         call uslope(q,dq,flatn,nx,ny,QVAR,2)
      endif

      do j = 0, ny+1 
         do i = 0, nx+1 

            clag = c(i,j) * q(i,j,QRHO)
            t1 = 1.d0/clag
            t2 = dth/(q(i,j,QRHO)*dy)
            t3 = 1.d0/(c(i,j)*q(i,j,QRHO))**2
            dgdlp = 2.d0*(1.d0-q(i,j,QGAM)
     &           / gamc(i,j))*(q(i,j,QGAM)-1.d0)

            ured = q(i,j,QV)
            dleft  = .5d0 - max(0.d0,ured+c(i,j))*dth/dy
            dleft0 = .5d0*(1.d0-ured*dt/dy)
            rtilde = q(i,j,QRHO) + dleft*dq(i,j,QRHO)
            gmtlde = q(i,j,QGAM) + dleft*dq(i,j,QGAM)
            ptr0l  = q(i,j,QPRES) + dleft0*dq(i,j,QPRES)
            rtr0l  = q(i,j,QRHO) + dleft0*dq(i,j,QRHO)

            
            if (ured-c(i,j) .lt. 0.d0 ) then
               betaml = 0.d0
            else
               betaml = (dq(i,j,QPRES)*t1-dq(i,j,QV))*t2
            endif

            if (ured .lt. 0.d0 ) then
               beta0l = 0.d0
            else
               beta0l =-(dq(i,j,QPRES)*t3 - dq(i,j,QRHO)
     &              /(rtilde*rtr0l))*dth*c(i,j)/dy
            endif

            qym(i,j,QRHO) = 1.d0/(1.d0/rtilde - beta0l  - betaml)
            qym(i,j,QPRES) = q(i,j,QPRES) + dleft*dq(i,j,QPRES) +
     1           clag*clag*betaml
            qym(i,j,QV) = q(i,j,QV) +  dleft*dq(i,j,QV) - clag*betaml

            if (ured .lt. 0.d0 ) then
               qym(i,j,QU) = q(i,j,QU) + dleft*dq(i,j,QU)
               qym(i,j,QGAM)= gmtlde
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qym(i,j,n) = q(i,j,n) + dleft*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qym(i,j,ns) = q(i,j,ns) + dleft*dq(i,j,ns)
               enddo
            else
               qym(i,j,QU) = q(i,j,QU) + dleft0*dq(i,j,QU)
               qym(i,j,QGAM) = q(i,j,QGAM) + dleft0*dq(i,j,QGAM)
     1              +dgdlp*(qym(i,j,QPRES)-ptr0l)/(qym(i,j,QPRES)+ptr0l)
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qym(i,j,n) = q(i,j,n) + dleft0*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qym(i,j,ns) = q(i,j,ns) + dleft0*dq(i,j,ns)
               enddo
            endif

            drght = -.5d0 - min(0.d0,ured - c(i,j))*dth/dy
            drght0 = -.5d0*(1.d0+ured*dt/dy)
            rtilde = q(i,j,QRHO) + drght*dq(i,j,QRHO)
            rtr0r = q(i,j,QRHO) + drght0*dq(i,j,QRHO)
            ptr0r = q(i,j,QPRES) + drght0*dq(i,j,QPRES)
            gmtlde = q(i,j,QGAM) + drght*dq(i,j,QGAM)

            if (ured+c(i,j) .ge. 0.0d0 ) then
               betapr = 0.d0
            else
               betapr =-(dq(i,j,QPRES)*t1+dq(i,j,QV))*t2
            endif

            if (ured .ge. 0.d0 ) then
               beta0r = 0.d0
            else
               beta0r = (dq(i,j,QPRES)*t3 -dq(i,j,QRHO)
     &              /(rtilde*rtr0r))*dth*c(i,j)/dy
            endif

            qyp(i,j,QRHO) = 1.d0/(1.d0/rtilde - beta0r - betapr)
            qyp(i,j,QPRES) = q(i,j,QPRES) + drght*dq(i,j,QPRES) +
     1           clag*clag*betapr
            qyp(i,j,QV) = q(i,j,QV) + drght*dq(i,j,QV) + clag*betapr

            if (ured .lt. 0.0 ) then
               qyp(i,j,QU) = q(i,j,QU)+drght0*dq(i,j,QU)
               qyp(i,j,QGAM) = q(i,j,QGAM) + drght0*dq(i,j,QGAM)
     1           + dgdlp*(qyp(i,j,QPRES)-ptr0r)/(qyp(i,j,QPRES)+ptr0r)
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qyp(i,j,n) = q(i,j,n) + drght0*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qyp(i,j,ns) = q(i,j,ns) + drght0*dq(i,j,ns)
               enddo
            else
               qyp(i,j,QU) = q(i,j,QU)+drght*dq(i,j,QU)
               qyp(i,j,QGAM) = gmtlde
               do iadv = 1, nadv
                  n = QFA + iadv - 1
                  qyp(i,j,n) = q(i,j,n) + drght*dq(i,j,n)
               enddo
               do ispec = 1, nspec
                  ns = QFS + ispec - 1
                  qyp(i,j,ns) = q(i,j,ns) + drght*dq(i,j,ns)
               enddo
            endif

         enddo
      enddo

      end

c ::: 
c ::: ------------------------------------------------------------------
c ::: 
      subroutine consup(uin,uin_l1,uin_l2,uin_h1,uin_h2,
     &     uout,uout_l1,uout_l2,uout_h1,uout_h2,
     &     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2,
     &     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2,
     &     div,lo,hi,dx,dy,dt)
      implicit none
      include 'castro2d.fi'
      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
      double precision dx, dy, dt

      integer i, j, n
      double precision div1, dvol

c ::: add diffusive flux where flow is compressing
c ::: and convert to extensive quantities
c     Skip flux mod and update for Temperature (u comp = UTEMP)
c     ...in fact, hard code the flux to zero so that AMR is happy/consistent
      do n = 1, NVAR
         if (n.eq.UTEMP) then
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  flux1(i,j,n) = 0.d0
               enddo
            enddo

            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  flux2(i,j,n) = 0.d0
               enddo
            enddo
         else
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  div1 = .5d0*(div(i,j) + div(i,j+1))
                  div1 = difmag*min(0.d0,div1)
                  flux1(i,j,n) = flux1(i,j,n)
     &                 + dx*div1*(uin(i,j,n) - uin(i-1,j,n))
                  flux1(i,j,n) = flux1(i,j,n) * dy*dt
               enddo
            enddo

            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  div1 = .5d0*(div(i,j) + div(i+1,j))
                  div1 = difmag*min(0.d0,div1)
                  flux2(i,j,n) = flux2(i,j,n)
     &                 + dy*div1*(uin(i,j,n) - uin(i,j-1,n))
                  flux2(i,j,n) = flux2(i,j,n) * dx*dt
               enddo
            enddo
         endif
      enddo

c ::: conservative update  ..  if not polytropic, what is update for auxiliary variables, such as Temp?
      dvol = 1.d0/(dx*dy)
      if (gasType.eq.polytropic) then
         do n = 1, NVAR
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)
                  uout(i,j,n) = uin(i,j,n)
     &                 + dvol*( flux1(i,j,n) - flux1(i+1,j,n)
     &                 +        flux2(i,j,n) - flux2(i,j+1,n) )
               enddo
            enddo
         enddo
      else
         call bl_abort('consup for stellar?')
      endif

      end

      subroutine cmpflx(qm,im1,im2,jm1,jm2,
     &                  qp,ip1,ip2,jp1,jp2,
     &                  ilo,ihi,jlo,jhi,ln,lt,
     &                  flx,if1,if2,jf1,jf2,
     &                  tmp,it1,it2,jt1,jt2,
     &                  q,gamc,csml,c,nx,ny)
      implicit none
      include 'castro2d.fi'
      integer ln, lt, nx, ny
      integer im1,im2,jm1,jm2
      integer ip1,ip2,jp1,jp2
      integer if1,if2,jf1,jf2
      integer it1,it2,jt1,jt2
      integer ilo,ihi,jlo,jhi
      double precision qm(im1:im2, jm1:jm2, QVAR)
      double precision qp(ip1:ip2, jp1:jp2, QVAR)
      double precision flx(if1:if2, jf1:jf2, NVAR)
      double precision tmp(it1:it2, jt1:jt2, 2)
      double precision    q(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP, QVAR)
      double precision gamc(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)
      double precision    c(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)
      double precision csml(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)

c     Local variables
      integer nedge, is, ie, i, j, ii, jj, n, nto, nfr
      double precision, allocatable :: sm(:,:),sp(:,:),gmin(:),gmax(:)
      double precision, allocatable :: sflx(:,:), ugdnv(:), pgdnv(:)
      double precision, allocatable :: gamm(:), gamp(:)
      double precision, allocatable :: clagm(:), clagp(:), smallc(:)
      double precision, allocatable :: gamcm(:), gamcp(:)
      double precision, allocatable :: cm(:), cp(:)
      double precision, allocatable :: dgdlpm(:), dgdlpp(:)

      if (ln .eq. QU) then
         is = ilo
         ie = ihi
      else
         is = jlo
         ie = jhi
      endif

      nedge = ie - is + 1

      allocate ( sm(nedge,QVAR) )
      allocate ( sp(nedge,QVAR) )
      allocate ( gmin(nedge) )
      allocate ( gmax(nedge) )
      allocate ( clagm(nedge) )
      allocate ( clagp(nedge) )
      allocate ( smallc(nedge) )
      allocate ( gamcm(nedge) )
      allocate ( gamcp(nedge) )
      allocate ( cm(nedge) )
      allocate ( cp(nedge) )
      allocate ( dgdlpm(nedge) )
      allocate ( dgdlpp(nedge) )
      allocate ( sflx(nedge,NVAR) )
      allocate ( ugdnv(nedge) )
      allocate ( pgdnv(nedge) )


      if (ln .eq. QU) then
         do j = jlo, jhi
            do n = 1, QVAR 
               nfr = n
               if (n .eq. QU) then
                  nfr = ln
               endif
               if (n .eq. QV) then
                  nfr = lt
               endif
               do i = ilo, ihi 
                  ii = i-ilo+1
                  sm(ii,n) = qm(i,j,nfr)
                  sp(ii,n) = qp(i,j,nfr)
               enddo
            enddo

            do i = ilo, ihi 
               ii = i-ilo+1
               gmin(ii) = min(q(i-2,j,QGAM),q(i-1,j,QGAM)
     &              ,         q(i  ,j,QGAM),q(i+1,j,QGAM))
               gmax(ii) = max(q(i-2,j,QGAM),q(i-1,j,QGAM)
     &              ,         q(i  ,j,QGAM),q(i+1,j,QGAM))
               gamcm(ii) = gamc(i-1,j)
               gamcp(ii) = gamc(i  ,j)
               clagm(ii) = c(i-1,j)*q(i-1,j,QRHO)
               clagp(ii) = c(i  ,j)*q(i  ,j,QRHO)
               dgdlpm(ii) = (q(i-1,j,QGAM)-1.d0)
     &              * 2.d0*(1.d0-q(i-1,j,QGAM)/gamc(i-1,j))
               dgdlpp(ii) = (q(i  ,j,QGAM)-1.d0)
     &              * 2.d0*(1.d0-q(i  ,j,QGAM)/gamc(i  ,j))
               smallc(ii) = max( csml(i,j), csml(i-1,j) )
               cm(ii) = c(i-1,j)
               cp(ii) = c(i  ,j)
            enddo
            
c     solve Riemann problem (gdnv state passed back, but only (u,p) saved)
            call plmdeRP(sm, sp, nedge, gmin, gmax,
     &           gamcm, gamcp, clagm, clagp, cm, cp, smallc,
     &           dgdlpm, dgdlpp, sflx, ugdnv, pgdnv )
            
            do n = 1, NVAR 
               nto = n
               if (n .eq. QU) then
                  nto = ln
               endif
               if (n .eq. QV) then
                  nto = lt
               endif
               do i = ilo, ihi 
                  ii = i-ilo+1
                  flx(i,j,nto) = sflx(ii,n)
               enddo
            enddo
            
            do i = ilo, ihi 
               ii = i-ilo+1
               tmp(i,j,1) = pgdnv(ii)
               tmp(i,j,2) = ugdnv(ii)
            enddo
         enddo
      else
         do i = ilo, ihi
            do n = 1, QVAR 
               nfr = n
               if (n .eq. QU) then
                  nfr = ln
               endif
               if (n .eq. QV) then
                  nfr = lt
               endif
               do j = jlo, jhi 
                  ii = j-jlo+1
                  sm(ii,n) = qm(i,j,nfr)
                  sp(ii,n) = qp(i,j,nfr)
               enddo
            enddo

            do j = jlo, jhi 
               jj = j-jlo+1
               gmin(jj) = min(q(i,j-2,QGAM),q(i,j-1,QGAM)
     &              ,         q(i  ,j,QGAM),q(i,j+1,QGAM))
               gmax(jj) = max(q(i,j-2,QGAM),q(i,j-1,QGAM)
     &              ,         q(i  ,j,QGAM),q(i,j+1,QGAM))
               gamcm(jj) = gamc(i,j-1)
               gamcp(jj) = gamc(i  ,j)
               clagm(jj) = c(i,j-1)*q(i,j-1,QRHO)
               clagp(jj) = c(i  ,j)*q(i  ,j,QRHO)
               dgdlpm(jj) = (q(i,j-1,QGAM)-1.d0)
     &              * 2.d0*(1.d0-q(i,j-1,QGAM)/gamc(i,j-1))
               dgdlpp(jj) = (q(i  ,j,QGAM)-1.d0)
     &              * 2.d0*(1.d0-q(i  ,j,QGAM)/gamc(i  ,j))
               smallc(jj) = max( csml(i,j), csml(i,j-1) )
               cm(jj) = c(i,j-1)
               cp(jj) = c(i,j)
            enddo
            
c     solve Riemann problem (gdnv state passed back, but only (u,p) saved)
            call plmdeRP(sm, sp, nedge, gmin, gmax,
     &           gamcm, gamcp, clagm, clagp, cm, cp, smallc,
     &           dgdlpm, dgdlpp, sflx, ugdnv, pgdnv )
            
            do n = 1, NVAR 
               nto = n
               if (n .eq. QU) then
                  nto = ln
               endif
               if (n .eq. QV) then
                  nto = lt
               endif
               do j = jlo, jhi 
                  ii = j-jlo+1
                  flx(i,j,nto) = sflx(ii,n)
               enddo
            enddo
            
            do j = jlo, jhi 
               ii = j-jlo+1
               tmp(i,j,1) = pgdnv(ii)
               tmp(i,j,2) = ugdnv(ii)
            enddo
         enddo
      endif

      end

      subroutine wsqge(p,v,gam,gdot,gstar,pstar,wsq,csq,gmin,gmax)
      implicit none
      double precision pstar,p,v,gam, gdot, wsq, gstar
      double precision csq,gmax,gmin,smlp1,small,temp,divide
      data smlp1,small/.001,1.e-07/

      gstar=(pstar-p)*gdot/(pstar+p) + gam
      gstar=max(gmin,min(gmax,gstar))
      wsq = (0.5d0*(gstar-1.d0)*(pstar+p)+pstar)
      temp = ((gstar-gam)/(gam-1.))
      if (pstar-p.eq.0.d0) then
         divide=small
      else
         divide=pstar-p
      endif
      
      temp=temp/divide
      wsq = wsq/(v - temp*p*v)
      if (abs(pstar/p-1.d0)-smlp1 .lt. 0.d0 ) then
         wsq = csq
      endif
      wsq=max(wsq,(.5d0*(gam-1.d0)/gam)*csq)
      end


      subroutine plmdeRP(qm, qp, nedge, gmin, gmax,
     &         gamcm, gamcp, clagm, clagp, cm, cp, smallc,
     &         dgdlpm, dgdlpp, uflx, ugdnv, pgdnv )

      implicit none
      include 'castro2d.fi'
      integer offset, nedge
      double precision qm(nedge, QVAR)
      double precision qp(nedge, QVAR)
      double precision gmin(*), gmax(*)
      double precision gamcm(*), gamcp(*)
      double precision clagm(*), clagp(*), cm(*), cp(*), smallc(*)
      double precision dgdlpm(*), dgdlpp(*)
      double precision uflx(nedge, NVAR)
      double precision ugdnv(nedge), pgdnv(nedge)
      double precision adgdnv(MAXADV),spgdnv(MAXSPEC)

      double precision pbarl,pbarr,rbarl,rbarr,clsql,clsqr,vbarl,vbarr
      double precision ubarl,ubarr,utbarl,utbarr,pstar,pstnm1,cstar
      double precision ustarm,ustarp,ustar,wr,wl,sflag,wlsq,wrsq
      double precision dgdlpl,dgdlpr,gamcl,gamcr
      double precision gamstar,wsi,wso,rstar,ushok,wtr,wtl,ctr,ctl,dpjmp
      double precision cffan,uutr,ptr,rhotr,vtr,gamtr,gambarr,gambarl
      double precision gamctr,ustnp1,ustnm1,dpditer,zp,zm,denom,t3
      double precision utgdnv,rgdnv,gamgdnv,wtrsq

      integer iadv, n, nq
      integer ispec, ns, nqs
      integer k, iter, itno
      double precision ff,a1,b1,c1,d1,e1,f1,weakwv
      data itno /3/
      data weakwv/1.e-03/

      ff(a1,b1,c1,d1,e1,f1)=(a1*f1 + c1*d1 - c1*f1*(e1-b1))/(c1+f1)

c ::: Solve Riemann Problem
c     NOTE: The calling routine will order velocity unknowns so that
c     for the purposes of this routine, the normal component is always
c     loaded in the QU slot.
      do k = 1, nedge 
         pbarl = max(qm(k,QPRES),smallp)
         pbarr = max(qp(k,QPRES),smallp)
         rbarl = max(qm(k,QRHO),smallr)
         rbarr = max(qp(k,QRHO),smallr)
         clsql = gamcm(k)*pbarl*rbarl
         clsqr = gamcp(k)*pbarr*rbarr
         vbarl = 1.d0/rbarl
         vbarr = 1.d0/rbarr
         ubarl = qm(k,QU)
         ubarr = qp(k,QU)
         gambarl = min(gmax(k),max(gmin(k),qm(k,QGAM)))
         gambarr = min(gmax(k),max(gmin(k),qp(k,QGAM)))
         utbarl = qm(k,QV)
         utbarr = qp(k,QV)
         dgdlpl = dgdlpm(k)
         dgdlpr = dgdlpp(k)
         gamcl = gamcm(k)
         gamcr = gamcp(k)

         wl = clagm(k)
         wr = clagp(k)

         pstar =(wl*pbarr+wr*pbarl-wr*wl*(ubarr-ubarl))/(wl+wr)
         pstar=max(pstar,smallp)
         pstnm1 = pstar
         
         call wsqge(pbarl,vbarl,gambarl,dgdlpm(k),
     1        gamstar,pstar,wlsq,clsql,gmin(k),gmax(k))
         call wsqge(pbarr,vbarr,gambarr,dgdlpp(k),
     1        gamstar,pstar,wrsq,clsqr,gmin(k),gmax(k))

         wl = sqrt(wlsq)
         wr = sqrt(wrsq)
         ustarp = ubarl - (pstar-pbarl)/wl
         ustarm = ubarr + (pstar-pbarr)/wr
         pstar = ff(pbarl,ubarl,wl,pbarr,ubarr,wr)
         pstar = max(pstar,smallp)


         do iter = 1,itno
            
            call wsqge(pbarl,vbarl,gambarl,dgdlpm(k),
     1           gamstar,pstar,wlsq,clsql,gmin(k),gmax(k))
            call wsqge(pbarr,vbarr,gambarr,dgdlpp(k),
     1           gamstar,pstar,wrsq,clsqr,gmin(k),gmax(k))
            
            wl = 1.d0/sqrt(wlsq)
            wr = 1.d0/sqrt(wrsq)
            
            ustnm1=ustarm
            ustnp1=ustarp
            
            ustarm = ubarr - (pbarr-pstar)*wr
            ustarp = ubarl + (pbarl-pstar)*wl
            
            dpditer=abs(pstnm1-pstar)
            zp=abs(ustarp-ustnp1)
            if (zp-weakwv*cm(k) .lt. 0.d0 ) then
               zp = dpditer*wl
            endif
            zm=abs(ustarm-ustnm1)
            if (zm-weakwv*cp(k) .lt. 0.d0 ) then
               zm = dpditer*wr
            endif
            
            denom=dpditer/max(zp+zm,smallc(k))
            pstnm1 = pstar
            pstar = pstar - denom*(ustarm-ustarp)
            pstar = max(pstar,smallp)

            if (abs(dpditer)/pstar .lt. 1.e-8) goto 100
   
         enddo

 100     ustar = .5d0*(ustarp + ustarm)
         sflag = -sign(1.d0,ustar)

         if (sflag .ge. 0.0 ) then
            uutr = ubarr
            ptr = pbarr
            rhotr = rbarr
            vtr = vbarr
            gamtr = gambarr
            gamctr = gamcr
            t3 = dgdlpr
            utgdnv = utbarr
            do iadv = 1, nadv
               n = QFA + iadv - 1
               adgdnv(iadv) = qp(k,n)
            enddo
            do ispec = 1, nspec
               ns = QFS + ispec - 1
               spgdnv(ispec) = qp(k,ns)
            enddo
         else
            uutr = ubarl
            ptr = pbarl
            rhotr = rbarl
            vtr = vbarl
            gamtr = gambarl
            gamctr = gamcl
            t3 = dgdlpl
            utgdnv = utbarl
            do iadv = 1, nadv
               n = QFA + iadv - 1
               adgdnv(iadv) = qm(k,n)
            enddo
            do ispec = 1, nspec
               ns = QFS + ispec - 1
               spgdnv(ispec) = qm(k,ns)
            enddo
         endif

         ctr = sqrt(gamctr*ptr/rhotr)
         clsql = (ctr*rhotr)**2
         call wsqge(ptr,vtr,gamtr,t3,gamstar,pstar,wtrsq,
     &        clsql,gmin(k),gmax(k))

         wtr = sqrt(wtrsq)
         dpjmp = pstar - ptr
         rstar = max(1.d0-rhotr*dpjmp/wtrsq,
     1               (gamtr-1.d0)/(gamtr+1.d0))
         rstar = rhotr/rstar
         ushok = sflag*ustar + wtr/rstar
         cstar = sqrt(gamctr*pstar/rstar)
         if (dpjmp .ge. 0.d0 ) then
            wsi = ushok
         else
            wsi=sflag*ustar+cstar
         endif
         if (dpjmp .ge. 0.d0 ) then
            wso = ushok
         else
            wso=sflag*uutr+ctr
         endif
         
         cffan=(wso+wsi)/max(wso-wsi,wso+wsi,small)
         cffan=.5d0*(cffan+1.d0)
         cffan=max(0.d0,min(1.d0,cffan))
         
         rgdnv = rhotr + cffan*(rstar-rhotr)
         ugdnv(k)=uutr + cffan*(ustar-uutr)
         pgdnv(k)=ptr  + cffan*(pstar-ptr)
         gamgdnv = gamtr + cffan*(gamstar - gamtr)

         if (wso .lt. 0.d0 ) then
            pgdnv(k) = ptr
            ugdnv(k) = uutr
            rgdnv = rhotr
            gamgdnv = gamtr
         endif

         if (wsi .ge. 0.d0 ) then
            pgdnv(k) = pstar
            ugdnv(k) = ustar
            rgdnv = rstar
            gamgdnv = gamstar
         endif

c     Compute fluxes, order as conserved state (not q)
         uflx(k,URHO) = ugdnv(k)*rgdnv
         uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) + pgdnv(k)
         uflx(k,UMY) = uflx(k,URHO)*utgdnv
         uflx(k,UEDEN)=uflx(k,URHO)*(.5d0*(ugdnv(k)**2 + utgdnv**2)
     &        + pgdnv(k)/((gamgdnv-1.d0)*rgdnv))
     &        + ugdnv(k)*pgdnv(k)

         do iadv = 1, nadv
            n = UFA + iadv - 1
            nq = QFA + iadv - 1
            uflx(k,n) = uflx(k,URHO)*adgdnv(iadv)
         enddo

         do ispec = 1, nspec
            ns = UFS + ispec - 1
            nqs = QFS + ispec - 1
            uflx(k,ns) = uflx(k,URHO)*spgdnv(ispec)
         enddo

      enddo
      end

c ::: 
c ::: ------------------------------------------------------------------
c ::: 
      subroutine transx(q,iq1,iq2,jq1,jq2, qo,io1,io2,jo1,jo2, fx,xtmp,
     & if1,if2,jf1,jf2, gamc,ig1,ig2,jg1,jg2, 
     & hdt, src,is1,is2,js1,js2, cdtdx,ilo,ihi,jlo,
     & jhi)
      implicit none
      include 'castro2d.fi'
      integer iq1,iq2,jq1,jq2
      integer io1,io2,jo1,jo2
      integer if1,if2,jf1,jf2
      integer ig1,ig2,jg1,jg2
      integer is1,is2,js1,js2
      integer ilo,ihi,jlo,jhi
      double precision cdtdx
      double precision q(iq1:iq2, jq1:jq2, QVAR)
      double precision qo(io1:io2, jo1:jo2, QVAR)
      double precision gamc(ig1:ig2, jg1:jg2)
      double precision fx(if1:if2, jf1:jf2, NVAR)
      double precision xtmp(if1:if2, jf1:jf2, 2)
      double precision  hdt
      double precision  src(is1:is2, js1:js2, SVAR)

      integer i, j
      double precision rr, ru, rv, re, eken, rhoeken
      double precision rrnew, runew, rvnew, renew
      double precision pgp, pgm, ugp, ugm, dup, pav, du, pnew
      integer n, nq, iadv
      double precision compu, compn
      integer ns, nqs, ispec, nss
      double precision comps, compsn

      do iadv = 1, nadv
          n = UFA + iadv - 1
          nq = QFA + iadv - 1
          ns = SFA + iadv - 1
          do j = jlo, jhi 
              do i = ilo, ihi 
                  rr = q(i,j,QRHO)
                  compu = rr*q(i,j,nq)
                  rrnew = rr - cdtdx*(fx(i+1,j,URHO) - fx(i,j,URHO))
     &                 + hdt*src(i,j,SRHO)
                  compn = compu - cdtdx*(fx(i+1,j,n)-fx(i,j,n))
     &                 + hdt*src(i,j,ns)
                  rrnew = max(rrnew,smallr)
                  qo(i,j,nq) = compn/rrnew
              enddo
          enddo
      enddo

      do ispec = 1, nspec
          ns = UFS + ispec - 1
          nqs = QFS + ispec - 1
          nss = SFS + ispec - 1
          do j = jlo, jhi 
              do i = ilo, ihi 
                  rr = q(i,j,QRHO)
                  comps = rr*q(i,j,nqs)
                  rrnew = rr - cdtdx*(fx(i+1,j,URHO) - fx(i,j,URHO))
     &                 + hdt*src(i,j,SRHO)
                  compsn = comps - cdtdx*(fx(i+1,j,ns)-fx(i,j,ns))
     &                 + hdt*src(i,j,nss)
                  rrnew = max(rrnew,smallr)
                  qo(i,j,nqs) = compsn/rrnew
              enddo
          enddo
      enddo

      do j = jlo, jhi 
          do i = ilo, ihi 
c ::: . convert to conservation form
              rr = q(i,j,QRHO)
              ru = rr*q(i,j,QU)
              rv = rr*q(i,j,QV)
              eken = 0.5d0*rr*(q(i,j,QU)**2 + q(i,j,QV)**2)
              re = eken + q(i,j,QPRES)/(q(i,j,QGAM) - 1.d0)
c ::: . add transverse predictor
              rrnew = rr - cdtdx*(fx(i+1,j,URHO) - fx(i,j,URHO))
     &             + hdt*src(i,j,SRHO)
              runew = ru - cdtdx*(fx(i+1,j,UMX)  - fx(i,j,UMX))
     &             + hdt*src(i,j,SMX)
              rvnew = rv - cdtdx*(fx(i+1,j,UMY)  - fx(i,j,UMY))
     &             + hdt*src(i,j,SMY)
              renew = re - cdtdx*(fx(i+1,j,UEDEN)- fx(i,j,UEDEN))
     &             + hdt*src(i,j,SEDEN)

              pgp = xtmp(i+1,j,1)
              pgm = xtmp(i  ,j,1)
              ugp = xtmp(i+1,j,2)
              ugm = xtmp(i  ,j,2)
              dup = pgp*ugp - pgm*ugm
              pav = 0.5d0*(pgp+pgm)
              du = ugp-ugm
              pnew = q(i,j,QPRES) -
     &              cdtdx*(dup + pav*du*(gamc(i,j)-1.d0))
c ::: . convert back to non-conservation form
              qo(i,j,QRHO) = max( smallr, rrnew )
              qo(i,j,QU) = runew/qo(i,j,QRHO)
              qo(i,j,QV) = rvnew/qo(i,j,QRHO)
              rhoeken = 0.5d0*(runew**2+rvnew**2)/qo(i,j,QRHO)
              qo(i,j,QGAM)= pnew / (renew - rhoeken) + 1.d0
              qo(i,j,QPRES) = pnew
          enddo
      enddo
      end

c ::: 
c ::: ------------------------------------------------------------------
c ::: 
      subroutine transy(q ,iq1,iq2,jq1,jq2,
     &                  qo,io1,io2,jo1,jo2,
     &                  fy,ytmp,if1,if2,jf1,jf2,
     &                  gamc,ig1,ig2,jg1,jg2,
     &                  hdt, src, is1, is2, js1, js2, cdtdy,
     &                  ilo,ihi,jlo,jhi)
      implicit none
      include 'castro2d.fi'
      integer iq1,iq2,jq1,jq2
      integer io1,io2,jo1,jo2
      integer if1,if2,jf1,jf2
      integer ig1,ig2,jg1,jg2
      integer is1,is2,js1,js2
      integer ilo,ihi,jlo,jhi
      double precision cdtdy
      double precision q(iq1:iq2, jq1:jq2, QVAR)
      double precision qo(io1:io2, jo1:jo2, QVAR)
      double precision gamc(ig1:ig2, jg1:jg2)
      double precision fy(if1:if2, jf1:jf2, NVAR)
      double precision ytmp(if1:if2, jf1:jf2, 2)
      double precision  hdt
      double precision  src(is1:is2, js1:js2, SVAR)

      integer i, j
      double precision rr, ru, rv, re, eken, rhoeken
      double precision rrnew, runew, rvnew, renew
      double precision pgp, pgm, ugp, ugm, dup, pav, du, pnew

      integer n, nq, iadv
      double precision compu, compn
      integer ns, nqs, ispec, nss
      double precision comps, compsn

      do iadv = 1, nadv
          n = UFA + iadv - 1
          nq = QFA + iadv - 1
          ns = SFA + iadv - 1
          do j = jlo, jhi 
              do i = ilo, ihi 
                  rr = q(i,j,QRHO)
                  compu = rr*q(i,j,nq)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
     &                 + hdt*src(i,j,SRHO)
                  compn = compu - cdtdy*(fy(i,j+1,n)-fy(i,j,n))
     &                 + hdt*src(i,j,ns)
                  rrnew = max(rrnew,smallr)
                  qo(i,j,nq) = compn/rrnew
              enddo
          enddo
      enddo

      do ispec = 1, nspec 
          ns = UFS + ispec - 1
          nqs = QFS + ispec - 1
          nss = SFS + ispec - 1
          do j = jlo, jhi 
              do i = ilo, ihi 
                  rr = q(i,j,QRHO)
                  comps = rr*q(i,j,nqs)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO))
     &                 + hdt*src(i,j,SRHO)
                  compsn = comps - cdtdy*(fy(i,j+1,ns)-fy(i,j,ns))
     &                 + hdt*src(i,j,nss)
                  rrnew = max(rrnew,smallr)
                  qo(i,j,nqs) = compsn/rrnew
              enddo
          enddo
      enddo

      do j = jlo, jhi 
          do i = ilo, ihi 
c ::: . convert to conservation form
              rr = q(i,j,QRHO)
              ru = rr*q(i,j,QU)
              rv = rr*q(i,j,QV)
              eken = 0.5d0*rr*(q(i,j,QU)**2 + q(i,j,QV)**2)
              re = eken + q(i,j,QPRES)/(q(i,j,QGAM) - 1.d0)

c ::: . add transverse predictor
              rrnew = rr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO))
     &             + hdt*src(i,j,SRHO)
              runew = ru - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX))
     &             + hdt*src(i,j,SMX)
              rvnew = rv - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY))
     &             + hdt*src(i,j,SMY)
              renew = re - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN))
     &             + hdt*src(i,j,SEDEN)

              pgp = ytmp(i,j+1,1)
              pgm = ytmp(i,j  ,1)
              ugp = ytmp(i,j+1,2)
              ugm = ytmp(i,j  ,2)
              dup = pgp*ugp - pgm*ugm
              pav = 0.5d0*(pgp+pgm)
              du = ugp-ugm
              pnew = q(i,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-1.d0))

c ::: . convert back to non-conservation form
              qo(i,j,QRHO) = max( smallr, rrnew )
              qo(i,j,QU) = runew/qo(i,j,QRHO)
              qo(i,j,QV) = rvnew/qo(i,j,QRHO)
              rhoeken = 0.5d0*(runew**2+rvnew**2)/qo(i,j,QRHO)
              qo(i,j,QGAM)= pnew / (renew - rhoeken) + 1.d0
              qo(i,j,QPRES) = pnew
          enddo
      enddo
      end

c ::: 
c ::: ------------------------------------------------------------------
c ::: 
      subroutine uslope(q,dq,flatn,nx,ny,nv,idir)
      implicit none
      include 'castro.fi'
      integer nx, ny, idir, nv
      double precision q(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP, nv)
      double precision dq(0:nx+1, 0:ny+1, nv)
      double precision flatn(1-NHYP:nx+NHYP, 1-NHYP:ny+NHYP)

c     Automatic local arrays
      double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

      integer i, j, n, nmax
      double precision dlft, drgt, slop, dq1
      double precision four3rd, sixth

      four3rd = 4.d0/3.d0
      sixth = 1.d0/6.d0

      nmax = MAX(nx,ny)
      allocate (dsgn(1-NHYP:nmax+NHYP))
      allocate (dlim(1-NHYP:nmax+NHYP))
      allocate (  df(1-NHYP:nmax+NHYP))
      allocate (dcen(1-NHYP:nmax+NHYP))

      do n = 1, nv 
          if (idir .eq. 1) then

c :::  slopes in first coordinate direction
              do j = 0, ny+1 

c ::: ..::::: first compute Fromm slopes
                  do i = -1, nx+2 
                      dlft = 2.d0*(q(i  ,j,n) - q(i-1,j,n))
                      drgt = 2.d0*(q(i+1,j,n) - q(i  ,j,n))
                      dcen(i) = .25d0 * (dlft+drgt)
                      dsgn(i) = sign(1.d0, dcen(i))
                      slop = min( abs(dlft), abs(drgt) )
c                      dlim(i) = cvmgp( slop, 0.d0, dlft*drgt )
                      if (dlft*drgt .ge. 0.0) then
                         dlim(i) = slop
                      else
                         dlim(i) = 0.d0
                      endif
                      df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
                  enddo

c ::: ..::::: now limited fourth order slopes
                  do i = 0, nx+1 
                      dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                      dq(i,j,n) = flatn(i,j)*
     &                            dsgn(i)*min(dlim(i),abs(dq1))
                  enddo
              enddo

          else


c :::  compute slopes in second coordinate direction
              do i = 0, nx+1 
c ::: ..::::: first compute Fromm slopes for this column
                  do j = -1, ny+2 
                      dlft = 2.d0*(q(i,j  ,n) - q(i,j-1,n))
                      drgt = 2.d0*(q(i,j+1,n) - q(i,j  ,n))
                      dcen(j) = .25d0 * (dlft+drgt)
                      dsgn(j) = sign( 1.d0, dcen(j) )
                      slop = min( abs(dlft), abs(drgt) )
c                      dlim(j) = cvmgp( slop, 0.d0, dlft*drgt )
                      if (dlft*drgt .ge. 0.0) then
                         dlim(j) = slop
                      else
                         dlim(j) = 0.d0
                      endif
                      df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
                  enddo

c ::: ..::::: now compute limited fourth order slopes
                  do j = 0, ny+1 
                      dq1 = four3rd*dcen(j) -
     &                      sixth*( df(j+1) + df(j-1) )
                      dq(i,j,n) = flatn(i,j)*
     &                            dsgn(j)*min(dlim(j),abs(dq1))
                  enddo
              enddo

          endif
      enddo
      end

c ::: 
c ::: ------------------------------------------------------------------
c -
c ::: 
      subroutine uflaten(lo,hi,p,u,v,rho,c,flatn,
     &     q_l1,q_l2,q_h1,q_h2)
      implicit none
      include 'castro.fi'
      integer lo(2),hi(2)
      integer q_l1,q_l2,q_h1,q_h2
      double precision p(q_l1:q_h1,q_l2:q_h2)
      double precision u(q_l1:q_h1,q_l2:q_h2)
      double precision v(q_l1:q_h1,q_l2:q_h2)
      double precision rho(q_l1:q_h1,q_l2:q_h2)
      double precision c(q_l1:q_h1,q_l2:q_h2)
      double precision flatn(q_l1:q_h1,q_l2:q_h2)

c ::: Local arrays
      double precision, allocatable :: dp(:), z(:), chi(:)

      integer i, j, idx, ishft
      double precision shktst, zcut1, zcut2, dzcut
      double precision denom, zeta, tst, tmp, ftmp
      integer nx,ny,nmax

c :::  knobs for detection of strong shock
      data shktst /0.33d0/
      data zcut1 /0.75d0/
      data zcut2 /0.85d0/

      nx = hi(1)-lo(1)+3
      ny = hi(2)-lo(2)+3
      nmax = max(nx,ny)
      allocate(dp(0:nmax-1),z(0:nmax-1),chi(0:nmax-1))

      dzcut = 1.d0/(zcut2-zcut1)

      if (iorder .eq. 3) then
         do j = lo(2),hi(2) 
            do i = lo(1),hi(1) 
               flatn(i,j) = 1.d0
            enddo
          enddo
          return
      endif

c :::   x-direction flattening coef
      do j = lo(2),hi(2) 
         do i = lo(1)-1,hi(1)+1
            idx = i-lo(1)+1
            dp(idx) = p(i+1,j) - p(i-1,j)
            denom = max(smallp,abs(p(i+2,j)-p(i-2,j)))
            zeta = abs(dp(idx))/denom
            z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
c           tst = cvmgp( 1.d0, 0.d0, u(i-1,j) - u(i+1,j) )
            if (u(i-1,j)-u(i+1,j) .ge. 0.0) then
               tst = 1.d0
            else
               tst = 0.d0
            endif
            tmp = min(p(i+1,j),p(i-1,j))
c           chi(idx) = cvmgt(tst,0.d0,(abs(dp(idx))/tmp).gt.shktst)
            if ((abs(dp(idx))/tmp).gt.shktst) then
               chi(idx) = tst
            else
               chi(idx) = 0.d0
            endif
         enddo
         do i = lo(1),hi(1)
            idx = i-lo(1)+1
            if(dp(idx).gt.0.d0)then
               ishft = 1
            else
               ishft = -1
            endif
            flatn(i,j) = 1.d0 -
     &           max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
         enddo
      enddo
      
c :::   y-direction flattening coef
      do i = lo(1),hi(1)
         do j = lo(2)-1,hi(2)+1
            idx = j-lo(2)+1
            dp(idx) = p(i,j+1) - p(i,j-1)
            denom = max(smallp,abs(p(i,j+2)-p(i,j-2)))
            zeta = abs(dp(idx))/denom
            z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
c           tst = cvmgp( 1.d0, 0.d0, v(i,j-1) - v(i,j+1) )
            if (v(i,j-1)-v(i,j+1) .ge. 0.0) then
               tst = 1.d0
            else
               tst = 0.d0
            endif
            tmp = min(p(i,j+1),p(i,j-1))
c           chi(idx) = cvmgt(tst,0.d0,(abs(dp(idx))/tmp).gt.shktst)
            if ((abs(dp(idx))/tmp).gt.shktst) then
               chi(idx) = tst
            else
               chi(idx) = 0.d0
            endif
         enddo
         do j = lo(2),hi(2)
            idx = j-lo(2)+1
            if(dp(idx).gt.0.d0)then
               ishft = 1
            else
               ishft = -1
            endif
            ftmp = 1.d0 -
     &           max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
            flatn(i,j) = min( flatn(i,j), ftmp )
         enddo
      enddo

      end

      subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,dx,dy,
     &     div,div_l1,div_l2,div_h1,div_h2,QU,QV)
      implicit none
      integer lo(2),hi(2)
      integer q_l1,q_l2,q_h1,q_h2
      integer div_l1,div_l2,div_h1,div_h2
      double precision q(q_l1:q_h1,q_l2:q_h2,*)
      double precision div(div_l1:div_h1,div_l2:div_h2)
      double precision dx, dy
      integer QU,QV

      integer i, j
      double precision ux, vy

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            ux=.5d0*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/dx
            vy=.5d0*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/dy
            div(i,j) = ux + vy
         enddo
      enddo
      end

      subroutine ca_gsrc(lo,hi,
     &     gnew,gnew_l1,gnew_l2,gnew_h1,gnew_h2,
     &     state,state_l1,state_l2,state_h1,state_h2,
     &     src,src_l1,src_l2,src_h1,src_h2,dt)
      implicit none
      include 'castro2d.fi'
      integer lo(2),hi(2)
      integer gnew_l1,gnew_l2,gnew_h1,gnew_h2
      integer state_l1,state_l2,state_h1,state_h2
      integer src_l1,src_l2,src_h1,src_h2
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,2)
      double precision  state(state_l1:state_h1,state_l2:state_h2,SVAR)
      double precision  src(src_l1:src_h1,src_l2:src_h2,SVAR)
      double precision dt

      integer i,j,n
      double precision rho, Up, Vp, gx, gy
      double precision SrU, SrV, SrE
      double precision pres, presnew

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            
            rho = state(i,j,URHO)
            Up = state(i,j,UMX) / rho
            Vp = state(i,j,UMY) / rho
            
            gx = gnew(i,j,1)
            gy = gnew(i,j,2)
            
c     Form predicted momentum sources (constant from old-time data)
            SrU = rho * gx
            SrV = rho * gy

c     Form predicted Eden source (computed directly from corrected momentum sources)
            SrE = SrU*(Up + SrU*dt/(2*rho))
     &           +SrV*(Vp + SrV*dt/(2*rho))

            src(i,j,1) = SrU
            src(i,j,2) = SrV
            src(i,j,3) = SrE

         enddo
      enddo
      end

      subroutine ca_corrgsrc(lo,hi,
     &     sold,sold_l1,sold_l2,sold_h1,sold_h2,
     &     gnew,gnew_l1,gnew_l2,gnew_h1,gnew_h2,
     &     uold,uold_l1,uold_l2,uold_h1,uold_h2,
     &     unew,unew_l1,unew_l2,unew_h1,unew_h2,dt)
      implicit none
      include 'castro2d.fi'
      integer lo(2),hi(2)
      integer sold_l1,sold_l2,sold_h1,sold_h2
      integer gnew_l1,gnew_l2,gnew_h1,gnew_h2
      integer uold_l1,uold_l2,uold_h1,uold_h2
      integer unew_l1,unew_l2,unew_h1,unew_h2
      double precision   sold(sold_l1:sold_h1,sold_l2:sold_h2,2+1)
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,2)
      double precision  uold(uold_l1:uold_h1,uold_l2:uold_h2,SVAR)
      double precision  unew(unew_l1:unew_h1,unew_l2:unew_h2,SVAR)
      double precision dt

      integer i,j,n
      double precision SrU_new, SrV_new, SrU_old, SrV_old, SrE_old
      double precision rhon, Upn, Vpn
      double precision gx, gy
      double precision SrUcorr,SrVcorr,SrEcorr

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            
            rhon = unew(i,j,URHO)
            Upn = unew(i,j,UMX) / rhon
            Vpn = unew(i,j,UMY) / rhon
            
            gx = gnew(i,j,1)
            gy = gnew(i,j,2)
            
            SrU_old = sold(i,j,1)
            SrV_old = sold(i,j,2)
            SrE_old = sold(i,j,3)
            
            SrU_new = rhon * gx
            SrV_new = rhon * gy
            
            SrUcorr = 0.5d0*(SrU_new - SrU_old)
            SrVcorr = 0.5d0*(SrV_new - SrV_old)
            SrEcorr = SrUcorr*(Upn + SrUcorr*dt/(2*rhon))
     &               +SrVcorr*(Vpn + SrVcorr*dt/(2*rhon))

c     Correct state by removing old source, adding new
            unew(i,j,UMX)   = unew(i,j,UMX)   + SrUcorr*dt
            unew(i,j,UMY)   = unew(i,j,UMY)   + SrVcorr*dt
            unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr*dt
            
         enddo
      enddo
      end

      subroutine ca_estdt(u,u_l1,u_l2,u_h1,u_h2,lo,hi,dx,dt)
      implicit none
      include 'castro2d.fi'
      integer u_l1,u_l2,u_h1,u_h2
      integer lo(2), hi(2)
      double precision u(u_l1:u_h1,u_l2:u_h2,NVAR)
      double precision dx(2), dt

      double precision, allocatable ::    p(:,:)
      double precision, allocatable ::    e(:,:)
      double precision, allocatable :: gamc(:,:)
      double precision, allocatable ::    c(:,:)
      double precision, allocatable :: csml(:,:)
      double precision, allocatable ::    T(:,:)
      double precision, allocatable ::    Y(:,:,:)

      double precision rhoInv,ux,uy,dt1,dt2
      integer i,j,n,cnt

      allocate(   p(lo(1):hi(1),lo(2):hi(2)))
      allocate(   e(lo(1):hi(1),lo(2):hi(2)))
      allocate(gamc(lo(1):hi(1),lo(2):hi(2)))
      allocate(   c(lo(1):hi(1),lo(2):hi(2)))
      allocate(csml(lo(1):hi(1),lo(2):hi(2)))
      allocate(   T(lo(1):hi(1),lo(2):hi(2)))
      allocate(   Y(lo(1):hi(1),lo(2):hi(2),nspec))

c     Translate to primitive variables, compute sound speed (call eos), get dtmax
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = 1.d0/u(i,j,URHO)
            ux = u(i,j,UMX)*rhoInv
            uy = u(i,j,UMY)*rhoInv
            e(i,j)=u(i,j,UEDEN)*rhoInv-0.5d0*(ux**2+uy**2)
            T(i,j)=u(i,j,UTEMP)
            do n=1,nspec
               Y(i,j,n)=u(i,j,UFS+n-1)*rhoInv
            enddo
         enddo
      enddo

      call gp_eos(
     &     gamc, lo(1),lo(2),hi(1),hi(2),
     &     p,    lo(1),lo(2),hi(1),hi(2),
     &     c,    lo(1),lo(2),hi(1),hi(2),
     &     csml, lo(1),lo(2),hi(1),hi(2),
     &     u,    u_l1, u_l2, u_h1, u_h2,
     &     e,    lo(1),lo(2),hi(1),hi(2),
     &     Y,    lo(1),lo(2),hi(1),hi(2),
     &     T,    lo(1),lo(2),hi(1),hi(2),
     &     lo, hi, nspec, cnt)

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = 1.d0/u(i,j,URHO)
            dt1 = dx(1)/(c(i,j) + abs(u(i,j,UMX)*rhoInv))
            dt2 = dx(2)/(c(i,j) + abs(u(i,j,UMY)*rhoInv))
            dt = min(dt,dt1,dt2)
         enddo
      enddo
      end

      subroutine ca_derpres(p,p_l1,p_l2,p_h1,p_h2,np,
     &     u,u_l1,u_l2,u_h1,u_h2,nc,lo,hi,domlo,
     &     domhi,dx,xlo,time,dt,bc,level,grid_no)
      implicit none
      include 'castro2d.fi'
      integer p_l1,p_l2,p_h1,p_h2,np
      integer u_l1,u_l2,u_h1,u_h2,nc
      integer lo(2), hi(2), domlo(2), domhi(2)
      double precision p(p_l1:p_h1,p_l2:p_h2,np)
      double precision u(u_l1:u_h1,u_l2:u_h2,nc)
      double precision dx(2), xlo(2), time, dt
      integer bc(2,2,nc), level, grid_no

      double precision, allocatable ::    e(:,:)
      double precision, allocatable :: gamc(:,:)
      double precision, allocatable ::    c(:,:)
      double precision, allocatable :: csml(:,:)
      double precision, allocatable ::    T(:,:)
      double precision, allocatable ::    Y(:,:,:)

      double precision rhoInv,ux,uy
      integer i,j,n,cnt

      allocate(   e(lo(1):hi(1),lo(2):hi(2)))
      allocate(gamc(lo(1):hi(1),lo(2):hi(2)))
      allocate(   c(lo(1):hi(1),lo(2):hi(2)))
      allocate(csml(lo(1):hi(1),lo(2):hi(2)))
      allocate(   T(lo(1):hi(1),lo(2):hi(2)))
      allocate(   Y(lo(1):hi(1),lo(2):hi(2),nspec))

c     Translate to primitive variables, compute sound speed (call eos), get dtmax
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            rhoInv = 1.d0/u(i,j,URHO)
            ux = u(i,j,UMX)*rhoInv
            uy = u(i,j,UMY)*rhoInv
            e(i,j)=u(i,j,UEDEN)*rhoInv-0.5d0*(ux**2+uy**2)
            T(i,j)=u(i,j,UTEMP)
            do n=1,nspec
               Y(i,j,n)=u(i,j,UFS+n-1)*rhoInv
            enddo
         enddo
      enddo

      call gp_eos(
     &     gamc, lo(1),lo(2),hi(1),hi(2),
     &     p,    lo(1),lo(2),hi(1),hi(2),
     &     c,    lo(1),lo(2),hi(1),hi(2),
     &     csml, lo(1),lo(2),hi(1),hi(2),
     &     u,    u_l1, u_l2, u_h1, u_h2,
     &     e,    lo(1),lo(2),hi(1),hi(2),
     &     Y,    lo(1),lo(2),hi(1),hi(2),
     &     T,    lo(1),lo(2),hi(1),hi(2),
     &     lo, hi, nspec, cnt)

      end

c :: ----------------------------------------------------------
c :: Volume-weight average the fine grid data onto the coarse
c :: grid.  Overlap is given in coarse grid coordinates.
c ::
c :: INPUTS / OUTPUTS:
c ::  crse      <=  coarse grid data
c ::  clo,chi    => index limits of crse array interior
c ::  ngc        => number of ghost cells in coarse array
c ::  nvar	 => number of components in arrays
c ::  fine       => fine grid data
c ::  flo,fhi    => index limits of fine array interior
c ::  ngf        => number of ghost cells in fine array
c ::  rfine      => (ignore) used in 2-D RZ calc
c ::  lo,hi      => index limits of overlap (crse grid)
c ::  lrat       => refinement ratio
c ::
c :: NOTE:
c ::  Assumes all data cell centered
c :: ----------------------------------------------------------
c ::
      subroutine ca_avgdown (crse,c_l1,c_l2,c_h1,c_h2,nvar,
     &     cv,cv_l1,cv_l2,cv_h1,cv_h2,
     &     fine,f_l1,f_l2,f_h1,f_h2,
     &     fv,fv_l1,fv_l2,fv_h1,fv_h2,lo,hi,lrat)
      implicit none
      integer c_l1,c_l2,c_h1,c_h2
      integer cv_l1,cv_l2,cv_h1,cv_h2
      integer f_l1,f_l2,f_h1,f_h2
      integer fv_l1,fv_l2,fv_h1,fv_h2
      integer lo(2), hi(2)
      integer nvar, lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2,nvar)
      double precision cv(cv_l1:cv_h1,cv_l2:cv_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2,nvar)
      double precision fv(fv_l1:fv_h1,fv_l2:fv_h2,nvar)

      integer i, j, n, ic, jc, ioff, joff
      integer lenx, leny, mxlen
      integer lratx, lraty
      double precision   volfrac

      lratx = lrat(1)
      lraty = lrat(2)
      lenx = hi(1)-lo(1)+1
      leny = hi(2)-lo(2)+1
      mxlen = max(lenx,leny)
      volfrac = 1.d0/float(lrat(1)*lrat(2))

      if (lenx .eq. mxlen) then
         do n = 1, nvar
c
c         ::::: set coarse grid to zero on overlap
c
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo
c
c         ::::: sum fine data
c
            do joff = 0, lraty-1
               do jc = lo(2), hi(2)
                  j = jc*lraty + joff
                  do ioff = 0, lratx-1
                     do ic = lo(1), hi(1)
                        i = ic*lratx + ioff
                        crse(ic,jc,n) = crse(ic,jc,n) + fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo
c            
c         ::::: divide out by volume weight
c
            do jc = lo(2), hi(2)
               do ic = lo(1), hi(1)
                  crse(ic,jc,n) = volfrac*crse(ic,jc,n)
               enddo
            enddo
            
         enddo

      else

         do n = 1, nvar
c
c         ::::: set coarse grid to zero on overlap
c
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = 0.d0
               enddo
            enddo
c
c         ::::: sum fine data
c
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  do joff = 0, lraty-1
                     do jc = lo(2), hi(2)
                        j = jc*lraty + joff
                        crse(ic,jc,n) = crse(ic,jc,n) + fine(i,j,n)
                     enddo
                  enddo
               enddo
            enddo
c            
c         ::::: divide out by volume weight
c
            do ic = lo(1), hi(1)
               do jc = lo(2), hi(2)
                  crse(ic,jc,n) = volfrac*crse(ic,jc,n)
               enddo
            enddo
            
         enddo

      end if

      end

c :: ----------------------------------------------------------
c :: Average the fine grid phi onto the coarse
c :: grid.  Overlap is given in coarse grid coordinates.
c :: Note this differs from ca_avgdown in that there is no volume weighting.
c ::
c :: INPUTS / OUTPUTS:
c ::  crse      <=  coarse grid data
c ::  clo,chi    => index limits of crse array interior
c ::  fine       => fine grid data
c ::  flo,fhi    => index limits of fine array interior
c ::  rfine      => (ignore) used in 2-D RZ calc
c ::  lo,hi      => index limits of overlap (crse grid)
c ::  lrat       => refinement ratio
c ::
c :: NOTE:
c ::  Assumes all data cell centered
c :: ----------------------------------------------------------
c ::
      subroutine ca_avgdown_phi (crse,c_l1,c_l2,c_h1,c_h2,
     &                           fine,f_l1,f_l2,f_h1,f_h2,
     &                           lo,hi,lrat)
      implicit none
      integer c_l1,c_l2,c_h1,c_h2
      integer f_l1,f_l2,f_h1,f_h2
      integer lo(2), hi(2)
      integer lrat(2)
      double precision crse(c_l1:c_h1,c_l2:c_h2)
      double precision fine(f_l1:f_h1,f_l2:f_h2)

      integer i, j, n, ic, jc, ioff, joff
      integer lratx, lraty
      double precision volfrac

      lratx = lrat(1)
      lraty = lrat(2)
      volfrac = 1.d0/float(lrat(1)*lrat(2))
c
c     ::::: set coarse grid to zero on overlap
c
      do jc = lo(2), hi(2)
         do ic = lo(1), hi(1)
            crse(ic,jc) = 0.d0
         enddo
      enddo
c
c         ::::: sum fine data
c
      do joff = 0, lraty-1
         do jc = lo(2), hi(2)
            j = jc*lraty + joff
            do ioff = 0, lratx-1
               do ic = lo(1), hi(1)
                  i = ic*lratx + ioff
                  crse(ic,jc) = crse(ic,jc) + fine(i,j)
               enddo
            enddo
         enddo
      enddo

      do ic = lo(1), hi(1)
         do jc = lo(2), hi(2)
            crse(ic,jc) = volfrac*crse(ic,jc)
         enddo
      enddo

      end


c :: ----------------------------------------------------------
c :: summass
c ::             MASS = sum{ vol(i,j)*rho(i,j) }
c ::
c :: INPUTS / OUTPUTS:
c ::  rho        => density field
c ::  rlo,rhi    => index limits of rho aray
c ::  lo,hi      => index limits of grid interior
c ::  delta      => cell size
c ::  mass      <=  total mass
c ::  r          => radius at cell center
c ::  irlo,hi    => index limits of r array
c ::  rz_flag    => == 1 if R_Z coords
c ::  tmp        => temp column array
c :: ----------------------------------------------------------
c ::
       subroutine ca_summass(rho,r_l1,r_l2,r_h1,r_h2,lo,hi,
     &     delta,mass,r,irlo,irhi,rz_flag,tmp,tlo,thi)
       implicit none
       integer irlo, irhi, rz_flag
       integer r_l1,r_l2,r_h1,r_h2
       integer lo(2), hi(2)
       double precision mass, delta(2)
       double precision rho(r_l1:r_h1,r_l2:r_h2)
       double precision r(irlo:irhi)
       integer tlo, thi
       double precision  tmp(tlo:thi)

       integer i, j
       double precision vol

       vol = delta(1)*delta(2)
       do j = lo(2),hi(2)
          tmp(j) = 0.d0
       enddo

       do i = lo(1), hi(1)
          do j = lo(2), hi(2)
             tmp(j) = tmp(j) + vol*rho(i,j)
          enddo
       enddo

       mass = 0.d0
       do j = lo(2), hi(2)
          mass = mass + tmp(j)
       enddo

       end

c ::
c :: ----------------------------------------------------------
c ::

       subroutine ca_grad_phi(
     &     gp,gp_l1,gp_l2,gp_h1,gp_h2,
     &     po,po_l1,po_l2,po_h1,po_h2,
     &     lo,hi,dx,ng)
       implicit none
       integer gp_l1,gp_l2,gp_h1,gp_h2
       integer po_l1,po_l2,po_h1,po_h2
       double precision gp(gp_l1:gp_h1,gp_l2:gp_h2,2)
       double precision po(po_l1:po_h1,po_l2:po_h2)
       integer lo(2), hi(2)
       double precision dx(2)
       integer ng

       integer i, j, n
       double precision dxInv(2)

       dxInv(1) = .5d0 / dx(1)
       dxInv(2) = .5d0 / dx(2)

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             gp(i,j,1) = (po(i+1,j)-po(i-1,j)) * dxInv(1)
             gp(i,j,2) = (po(i,j+1)-po(i,j-1)) * dxInv(2)
          enddo
       enddo

c     Do constant extrap for ng cells around the grid (skip corners)
       do j = lo(2), hi(2)
          do n = 1,ng
             gp(lo(1)-n,j,1) = gp(lo(1),j,1)
             gp(hi(1)+n,j,1) = gp(hi(1),j,1)
             gp(lo(1)-n,j,2) = gp(lo(1),j,2)
             gp(hi(1)+n,j,2) = gp(hi(1),j,2)
          enddo
       enddo
       do i = lo(1), hi(1)
          do n = 1,ng
             gp(i,lo(2)-n,1) = gp(i,lo(2),1)
             gp(i,hi(2)+n,1) = gp(i,hi(2),1)
             gp(i,lo(2)-n,2) = gp(i,lo(2),2)
             gp(i,hi(2)+n,2) = gp(i,hi(2),2)
          enddo
       enddo

       end

       subroutine ca_grad_phi_on_edges(
     &     xflux,xflux_l1,xflux_l2,xflux_h1,xflux_h2,
     &     yflux,yflux_l1,yflux_l2,yflux_h1,yflux_h2,
     &     po,po_l1,po_l2,po_h1,po_h2,
     &     lo,hi,dx)
       implicit none
       integer xflux_l1,xflux_l2,xflux_h1,xflux_h2
       integer yflux_l1,yflux_l2,yflux_h1,yflux_h2
       integer po_l1,po_l2,po_h1,po_h2
       double precision xflux(xflux_l1:xflux_h1,xflux_l2:xflux_h2)
       double precision yflux(yflux_l1:yflux_h1,yflux_l2:yflux_h2)
       double precision po(po_l1:po_h1,po_l2:po_h2)
       integer lo(2), hi(2)
       double precision dx(2)

       integer i, j
       double precision dxInv(2)

       dxInv(1) = 1.d0 / dx(1)
       dxInv(2) = 1.d0 / dx(2)

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1
             xflux(i,j) = (po(i,j)-po(i-1,j)) * dxInv(1) * dx(2)
          enddo
       enddo

       do j = lo(2), hi(2)+1
          do i = lo(1), hi(1)
             yflux(i,j) = (po(i,j)-po(i,j-1)) * dxInv(2) * dx(1)
          enddo
       enddo

       end

c-----------------------------------------------------------------------

      subroutine ca_edge_interp(flo, fhi, nc, ratio, dir,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer flo(0:2-1), fhi(0:2-1), nc, ratio(0:2-1), dir
      integer fine_l0, fine_l1, fine_h0, fine_h1
      double precision fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,ii,jj,n,P,M,clo(0:2-1),chi(0:2-1)
      double precision val, df

c     Do linear in dir, pc transverse to dir, leave alone the fine values
c     lining up with coarse edges--assume these have been set to hold the 
c     values you want to interpolate to the rest.
      if (dir.eq.0) then
         do n=1,nc
            do j=flo(1),fhi(1),ratio(1)
               do i=flo(0),fhi(0)-ratio(dir),ratio(0)
                  df = fine(i+ratio(dir),j,n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(j,flo(1)),MIN(j+ratio(1)-1,fhi(1))
                        fine(i+M,P,n) = val
                     enddo
                  enddo                     
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=flo(1),fhi(1)-ratio(dir),ratio(1)
               do i=flo(0),fhi(0)
                  df = fine(i,j+ratio(dir),n)-fine(i,j,n)
                  do M=1,ratio(dir)-1
                     val = fine(i,j,n) + df*dble(M)/dble(ratio(dir))
                     do P=MAX(i,flo(0)),MIN(i+ratio(0)-1,fhi(0))
                        fine(P,j+M,n) = val
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      end

c-----------------------------------------------------------------------

      subroutine ca_pc_edge_interp(lo, hi, nc, ratio, dir,
     &     crse, crse_l0, crse_l1, crse_h0, crse_h1,
     &     fine, fine_l0, fine_l1, fine_h0, fine_h1)
      implicit none
      integer lo(2),hi(2), nc, ratio(0:2-1), dir
      integer crse_l0, crse_l1, crse_h0, crse_h1
      integer fine_l0, fine_l1, fine_h0, fine_h1
      double precision crse(crse_l0:crse_h0,crse_l1:crse_h1,nc)
      double precision fine(fine_l0:fine_h0,fine_l1:fine_h1,nc)
      integer i,j,ii,jj,n,L
      double precision val, dc

c     For edge-based data, fill fine values with piecewise-constant interp of coarse data.
c     Operate only on faces that overlap--ie, only fill the fine faces that make up each
c     coarse face, leave the in-between faces alone.
      if (dir.eq.0) then
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(1)-1
                     fine(ii,jj+L,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      else
         do n=1,nc
            do j=lo(2),hi(2)
               jj = ratio(1)*j
               do i=lo(1),hi(1)
                  ii = ratio(0)*i
                  do L=0,ratio(0)-1
                     fine(ii+L,jj,n) = crse(i,j,n)
                  enddo
               enddo
            enddo
         enddo
      endif

      end

c-----------------------------------------------------------------------

      subroutine ca_avg_ec_to_cc(lo, hi,
     &     cc, ccl1, ccl2, cch1, cch2,
     &     ecx, ecxl1, ecxl2, ecxh1, ecxh2,
     &     ecy, ecyl1, ecyl2, ecyh1, ecyh2)
      implicit none
      integer lo(2),hi(2)
      integer ccl1, ccl2, cch1, cch2
      integer ecxl1, ecxl2, ecxh1, ecxh2
      integer ecyl1, ecyl2, ecyh1, ecyh2
      double precision cc(ccl1:cch1,ccl2:cch2,2)
      double precision ecx(ecxl1:ecxh1,ecxl2:ecxh2)
      double precision ecy(ecyl1:ecyh1,ecyl2:ecyh2)
      integer i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            cc(i,j,1) = 0.5d0 * ( ecx(i+1,j) + ecx(i,j) )
         enddo
      enddo

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            cc(i,j,2) = 0.5d0 * ( ecy(i,j+1) + ecy(i,j) )
         enddo
      enddo

      end
c-----------------------------------------------------------------------

      subroutine ca_test_residual(lo, hi,
     &     rhs, rhl1, rhl2, rhh1, rhh2,
     &     ecx, ecxl1, ecxl2, ecxh1, ecxh2,
     &     ecy, ecyl1, ecyl2, ecyh1, ecyh2,
     &     dx)
      implicit none
      integer lo(2),hi(2)
      integer rhl1, rhl2, rhh1, rhh2
      integer ecxl1, ecxl2, ecxh1, ecxh2
      integer ecyl1, ecyl2, ecyh1, ecyh2
      double precision rhs(rhl1:rhh1,rhl2:rhh2)
      double precision ecx(ecxl1:ecxh1,ecxl2:ecxh2)
      double precision ecy(ecyl1:ecyh1,ecyl2:ecyh2)
      double precision dx(2)
      integer i,j

      double precision lapphi

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            lapphi = (ecx(i+1,j)-ecx(i,j)) / dx(1) +
     $               (ecy(i,j+1)-ecy(i,j)) / dx(2)
            rhs(i,j) = rhs(i,j) - lapphi
         enddo
      enddo

      end

c-----------------------------------------------------------------------
      subroutine ca_average_ec (
     &     fx, fxl1, fxl2, fxh1, fxh2,
     &     fy, fyl1, fyl2, fyh1, fyh2,
     &     cx, cxl1, cxl2, cxh1, cxh2,
     &     cy, cyl1, cyl2, cyh1, cyh2,
     $     lo, hi, rr)
c
      integer lo(2),hi(2)
      integer fxl1, fxl2, fxh1, fxh2
      integer fyl1, fyl2, fyh1, fyh2
      integer cxl1, cxl2, cxh1, cxh2
      integer cyl1, cyl2, cyh1, cyh2
      double precision fx(fxl1:fxh1,fxl2:fxh2)
      double precision fy(fyl1:fyh1,fyl2:fyh2)
      double precision cx(cxl1:cxh1,cxl2:cxh2)
      double precision cy(cyl1:cyh1,cyl2:cyh2)
      integer rr(2)
c
      double precision facx,facy
      integer i,j,n

      facx = dble(rr(1))
      facy = dble(rr(2))

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)+1
            cx(i,j) = 0.d0
            do n = 0,facy-1
              cx(i,j) = cx(i,j) + fx(facx*i,facy*j+n)
            end do
            cx(i,j) = cx(i,j) / facy
         end do
      end do

      do i = lo(1), hi(1)
         do j = lo(2), hi(2)+1
            cy(i,j) = 0.d0
            do n = 0,facx-1
              cy(i,j) = cy(i,j) + fy(facx*i+n,facy*j)
            end do
            cy(i,j) = cy(i,j) / facx
         end do
      end do
c
      end

      subroutine ca_dervel(vel,vel_l1,vel_l2,vel_h1,vel_h2,nv,
     &                     dat,dat_l1,dat_l2,dat_h1,dat_h2,nc,lo,hi,domlo,
     &                     domhi,delta,xlo,time,dt,bc,level,grid_no)
c
c     This routine will derive kinetic energy from density
c     and the velocity field.
c
      integer          lo(2), hi(2)
      integer          vel_l1,vel_l2,vel_h1,vel_h2,nv
      integer          dat_l1,dat_l2,dat_h1,dat_h2,nc
      integer          domlo(2), domhi(2)
      integer          bc(2,2,nc)
      double precision delta(2), xlo(2), time, dt
      double precision vel(vel_l1:vel_h1,vel_l2:vel_h2,nv)
      double precision dat(dat_l1:dat_h1,dat_l2:dat_h2,nc)
      integer    level, grid_no
 
      integer    i,j
 
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            vel(i,j,1) = dat(i,j,2) / dat(i,j,1)
         end do
      end do
 
      end

