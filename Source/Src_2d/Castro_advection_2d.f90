
! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_umdrv(is_finest_level,time,lo,hi,domlo,domhi, &
                          uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                          uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                          ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                          ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                          src,src_l1,src_l2,src_h1,src_h2, &
                          grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                          delta,dt, &
                          flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                          flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                          area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                          area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                          dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                          vol,vol_l1,vol_l2,vol_h1,vol_h2,courno,verbose)

      use meth_params_module, only : URHO, QVAR, NVAR, NHYP, &
                                     do_sponge, normalize_species, allow_negative_energy

      implicit none

      integer is_finest_level
      integer lo(2),hi(2),verbose
      integer domlo(2),domhi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
      integer area1_l1,area1_l2,area1_h1,area1_h2
      integer area2_l1,area2_l2,area2_h1,area2_h2
      integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
      integer vol_l1,vol_l2,vol_h1,vol_h2
      integer src_l1,src_l2,src_h1,src_h2
      integer gv_l1,gv_l2,gv_h1,gv_h2

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision src(src_l1:src_h1,src_l2:src_h2,NVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
      double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
      double precision delta(2),dt,time,courno

!     Automatic arrays for workspace
      double precision, allocatable:: q(:,:,:)
      double precision, allocatable:: gamc(:,:)
      double precision, allocatable:: flatn(:,:)
      double precision, allocatable:: c(:,:)
      double precision, allocatable:: csml(:,:)
      double precision, allocatable:: div(:,:)
      double precision, allocatable:: pgdx(:,:)
      double precision, allocatable:: pgdy(:,:)
      double precision, allocatable:: srcQ(:,:,:)
      double precision, allocatable:: pdivu(:,:)

      integer ngq,ngf,iflaten
!     integer i_c,j_c

      double precision dx,dy

      allocate(     q(uin_l1:uin_h1,uin_l2:uin_h2,QVAR))
      allocate(  gamc(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate( flatn(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate(     c(uin_l1:uin_h1,uin_l2:uin_h2))
      allocate(  csml(uin_l1:uin_h1,uin_l2:uin_h2))

      allocate(  srcQ(src_l1:src_h1,src_l2:src_h2,QVAR))

      allocate(   div(lo(1)  :hi(1)+1,lo(2)  :hi(2)+1))
      allocate( pdivu(lo(1)  :hi(1)  ,lo(2)  :hi(2)))
      allocate(  pgdx(lo(1)  :hi(1)+1,lo(2)-1:hi(2)+1))
      allocate(  pgdy(lo(1)-1:hi(1)+1,lo(2)  :hi(2)+1))

      dx = delta(1)
      dy = delta(2)

      ngq = NHYP
      ngf = 1
      iflaten = 1

!     Translate to primitive variables, compute sound speeds
!     Note that (q,c,gamc,csml,flatn) are all dimensioned the same
!       and set to correspond to coordinates of (lo:hi)
      call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                   q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_h1,uin_h2, &
                   src,srcQ,src_l1,src_l2,src_h1,src_h2, &
                   courno,dx,dy,dt,ngq,ngf,iflaten)

!     Compute hyperbolic fluxes using unsplit Godunov
      call umeth2d(q,c,gamc,csml,flatn,uin_l1,uin_l2,uin_h1,uin_h2, &
                   srcQ, src_l1, src_l2, src_h1, src_h2,  &
                   grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                   lo(1),lo(2),hi(1),hi(2),dx,dy,dt, &
                   flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                   flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                   pgdx, lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
                   pgdy, lo(1)-1, lo(2), hi(1)+1, hi(2)+1, &
                   ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                   ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                   area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                   area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                   pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                   dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2)

!     Compute divergence of velocity field (on surroundingNodes(lo,hi))
      call divu(lo,hi,q,uin_l1,uin_l2,uin_h1,uin_h2, &
                delta,div,lo(1),lo(2),hi(1)+1,hi(2)+1)

!     Conservative update
      call consup(uin,    uin_l1,  uin_l2,  uin_h1,  uin_h2, &
                  uout,  uout_l1, uout_l2, uout_h1, uout_h2, &
                  pgdx,  lo(1), lo(2)-1, hi(1)+1, hi(2)+1, &
                  pgdy,lo(1)-1,   lo(2), hi(1)+1, hi(2)+1, &
                  src,    src_l1,  src_l2,  src_h1,  src_h2, &
                  grav,    gv_l1,   gv_l2,   gv_h1,   gv_h2, &
                  flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                  flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                  area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                  area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                  vol,    vol_l1,  vol_l2,  vol_h1,  vol_h2, &
                  div,pdivu,lo,hi,dx,dy,dt)

      ! Enforce the density >= small_dens.
      call enforce_minimum_density( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                                   uout,uout_l1,uout_l2,uout_h1,uout_h2,&
                                   lo,hi,verbose)

      ! Enforce the species >= 0
      call ca_enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      ! Normalize the species 
      if (normalize_species .eq. 1) &
         call normalize_new_species(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi)

      if (do_sponge .eq. 1) &
           call sponge(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi, &
                       time,dt, &
                       dx,dy,domlo,domhi)

      deallocate(q,gamc,flatn,c,csml,div,pgdx,pgdy,srcQ,pdivu)

      end subroutine ca_umdrv

! ::: ---------------------------------------------------------------
! ::: :: UMETH2D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: :: 
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: c           => (const)  sound speed
! ::: :: gamc        => (const)  cound speed gamma
! ::: :: csml        => (const)  local small c val
! ::: :: flatn       => (const)  flattening parameter
! ::: :: src         => (const)  source
! ::: :: nx          => (const)  number of cells in X direction
! ::: :: ny          => (const)  number of cells in Y direction
! ::: :: dx          => (const)  grid spacing in X direction
! ::: :: dy          => (const)  grid spacing in Y direction
! ::: :: dt          => (const)  time stepsize
! ::: :: flux1      <=  (modify) flux in X direction on X edges
! ::: :: flux2      <=  (modify) flux in Y direction on Y edges
! ::: ----------------------------------------------------------------

      subroutine umeth2d(q, c, gamc, csml, flatn, qd_l1, qd_l2, qd_h1, qd_h2,&
                         srcQ, src_l1, src_l2, src_h1, src_h2, &
                         grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                         ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
                         flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                         flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                         pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                         pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                         ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
                         ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
                         area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                         area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                         pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                         dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, ppm_type

      implicit none

      integer qd_l1, qd_l2, qd_h1, qd_h2
      integer dloga_l1, dloga_l2, dloga_h1, dloga_h2
      integer src_l1, src_l2, src_h1, src_h2
      integer gv_l1, gv_l2, gv_h1, gv_h2
      integer fd1_l1, fd1_l2, fd1_h1, fd1_h2
      integer fd2_l1, fd2_l2, fd2_h1, fd2_h2
      integer pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2
      integer pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2
      integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
      integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
      integer area1_l1, area1_l2, area1_h1, area1_h2
      integer area2_l1, area2_l2, area2_h1, area2_h2
      integer vol_l1, vol_l2, vol_h1, vol_h2
      integer ilo1, ilo2, ihi1, ihi2

      double precision dx, dy, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision  srcQ(src_l1:src_h1,src_l2:src_h2)
      double precision  grav( gv_l1: gv_h1, gv_l2: gv_h2)
      double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
      double precision pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
      double precision pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
      double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
      double precision pdivu(ilo1:ihi1,ilo2:ihi2)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)

!     Left and right state arrays (edge centered, cell centered)
      double precision, allocatable:: dq(:,:,:),  qm(:,:,:),   qp(:,:,:)
      double precision, allocatable::qxm(:,:,:),qym(:,:,:)
      double precision, allocatable::qxp(:,:,:),qyp(:,:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
      double precision, allocatable::   fx(:,:,:),  fy(:,:,:)
      double precision, allocatable::   pgdxtmp(:,:) ,  ugdxtmp(:,:)

!     Local scalar variables
      double precision :: dtdx
      double precision :: hdtdx, hdt, hdtdy
      integer          :: i,j

      allocate ( pgdxtmp(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2))
      allocate ( ugdxtmp(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2))
      allocate ( dq(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QVAR) )
      allocate ( fx(ilo1:ihi1+1,ilo2-1:ihi2+1,NVAR))
      allocate ( fy(ilo1-1:ihi1+1,ilo2:ihi2+1,NVAR))

!     Local constants
      dtdx = dt/dx
      hdtdx = 0.5d0*dtdx
      hdtdy = 0.5d0*dt/dy
      hdt = 0.5d0*dt

!     NOTE: Geometry terms need to be punched through

!     Trace to edges w/o transverse flux correction terms
      if (ppm_type .eq. 0) then
         call trace(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                    dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                    dq,qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                    grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                    ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
      else
         call trace_ppm(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                        qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
      end if

      call cmpflx(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                  pgdxtmp, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                  ugdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                  gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                  1, ilo1, ihi1, ilo2-1, ihi2+1)

      call cmpflx(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                  pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                  gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                  2, ilo1-1, ihi1+1, ilo2, ihi2)

      call transy(qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
                  pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                  gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                  srcQ, src_l1, src_l2, src_h1, src_h2, &
                  grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                  hdt, hdtdy, &
                  ilo1-1, ihi1+1, ilo2, ihi2)

      call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                  pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                  ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                  gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                  1, ilo1, ihi1, ilo2, ihi2)
      
      call transx(qym, qm,qyp,qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
                  pgdxtmp, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                  ugdxtmp, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                  gamc, qd_l1, qd_l2, qd_h1, qd_h2, &
                  srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
                  grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                  hdt, hdtdx, &
                  area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                  vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                  ilo1, ihi1, ilo2-1, ihi2+1)

      call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
                  flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                  pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                  ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                  gamc, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
                  2, ilo1, ihi1, ilo2, ihi2)
      

      do j = ilo2,ihi2
      do i = ilo1,ihi1
         pdivu(i,j) = 0.5d0 * &
              ((pgdx(i+1,j)+pgdx(i,j))*(ugdx(i+1,j)*area1(i+1,j)-ugdx(i,j)*area1(i,j)) &
              +(pgdy(i,j+1)+pgdy(i,j))*(ugdy(i,j+1)*area2(i,j+1)-ugdy(i,j)*area2(i,j)) ) / vol(i,j)
      end do
      end do

      deallocate(dq,qm,qp,qxm,qxp,qym,qyp)
      deallocate(fx,fy)
      deallocate(pgdxtmp,ugdxtmp)

      end subroutine umeth2d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ctoprim(lo,hi, &
                         uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                         q,c,gamc,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
                         src,srcQ,src_l1,src_l2,src_h1,src_h2, &
                         courno,dx,dy,dt,ngp,ngf,iflaten)

!     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
!     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
!     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
!     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
!     routine that computes flatn).  

      use network, only : nspec, naux
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
                                     QVAR, QRHO, QU, QV, QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
                                     nadv, allow_negative_energy, small_temp

      implicit none

      double precision, parameter:: small = 1.d-8

      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer q_l1,q_l2,q_h1,q_h2
      integer src_l1,src_l2,src_h1,src_h2
      integer iflaten

      double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
      double precision :: c(q_l1:q_h1,q_l2:q_h2)
      double precision :: gamc(q_l1:q_h1,q_l2:q_h2)
      double precision :: csml(q_l1:q_h1,q_l2:q_h2)
      double precision :: flatn(q_l1:q_h1,q_l2:q_h2)
      double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
      double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
      double precision :: dx, dy, dt, courno

      double precision, allocatable :: dpdrho(:,:)
      double precision, allocatable :: dpde(:,:)
      double precision, allocatable :: dpdX_er(:,:,:)

      integer          :: i, j
      integer          :: pt_index(2)
      integer          :: ngp, ngf, loq(2), hiq(2)
      integer          :: iadv, ispec, iaux, n, nq
      double precision :: courx, coury, courmx, courmy

      allocate(     dpdrho(q_l1:q_h1,q_l2:q_h2))
      allocate(     dpde(q_l1:q_h1,q_l2:q_h2))
      allocate(  dpdX_er(q_l1:q_h1,q_l2:q_h2,nspec))

      do i=1,2
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo

!     Make q (all but p), except put e in slot for rho.e, fix after eos call
!     The temperature is used as an initial guess for the eos call and will be overwritten
      do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)

            if (uin(i,j,URHO) .le. 0.d0) then
               print *,'   '
               print *,'>>> Error: Castro_2d::ctoprim ',i,j
               print *,'>>> ... negative density ',uin(i,j,URHO)
               print *,'    '
               call bl_error("Error:: Castro_2d.f90 :: ctoprim")
            end if

            q(i,j,QRHO) = uin(i,j,URHO)
            q(i,j,QU) = uin(i,j,UMX)/uin(i,j,URHO)
            q(i,j,QV) = uin(i,j,UMY)/uin(i,j,URHO)
            q(i,j,QREINT ) = uin(i,j,UEINT)/q(i,j,QRHO)
            q(i,j,QTEMP  ) = uin(i,j,UTEMP)
         enddo
      enddo

!    Load advected quatities, c, into q, assuming they arrived in uin as rho.c
     do iadv = 1, nadv
         n  = UFA + iadv - 1
         nq = QFA + iadv - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
            enddo
         enddo
      enddo
      
!     Load chemical species, c, into q, assuming they arrived in uin as rho.c
      do ispec = 1, nspec
         n  = UFS + ispec - 1
         nq = QFS + ispec - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
            enddo
         enddo
      enddo
      
!     Load auxiliary variables which are needed in the EOS
      do iaux = 1, naux
         n  = UFX + iaux - 1
         nq = QFX + iaux - 1
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
            enddo
         enddo
      enddo

!     Get gamc, p, T, c, csml using q state 
      do j = loq(2), hiq(2)
        do i = loq(1), hiq(1)

           ! If necessary, reset the energy using small_temp
           if ((allow_negative_energy .eq. 0) .and. (q(i,j,QREINT) .lt. 0)) then
              q(i,j,QTEMP) = small_temp
              call eos_given_RTX(q(i,j,QREINT),q(i,j,QPRES),q(i,j,QRHO), &
                                 q(i,j,QTEMP),q(i,j,QFS:))
              if (q(i,j,QREINT) .lt. 0.d0) then
                 print *,'   '
                 print *,'>>> Error: Castro_2d::ctoprim ',i,j
                 print *,'>>> ... new e from eos_given_RTX call is negative ',q(i,j,QREINT)
                 print *,'    '
                 call bl_error("Error:: Castro_2d.f90 :: ctoprim")
              end if
           end if

           pt_index(1) = i
           pt_index(2) = j
           call eos_given_ReX(gamc(i,j), q(i,j,QPRES), c(i,j), q(i,j,QTEMP), &
                              dpdrho(i,j), dpde(i,j), &
                              q(i,j,QRHO), q(i,j,QREINT), q(i,j,QFS:), &
                              pt_index=pt_index)!, &
!                              dpdX_er=dpdX_er(i,j,:))
           csml(i,j) = max(small, small * c(i,j))
        end do
      end do

!     Make this "rho e" instead of "e"
      do j = loq(2),hiq(2)
         do i = loq(1),hiq(1)
            q(i,j,QREINT) = q(i,j,QREINT)*q(i,j,QRHO)
         enddo
      enddo

!     Compute sources in terms of Q
      do j = lo(2)-1, hi(2)+1
        do i = lo(1)-1, hi(1)+1

           srcQ(i,j,QRHO  ) = src(i,j,URHO)
           srcQ(i,j,QU    ) = (src(i,j,UMX) - q(i,j,QU) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
           srcQ(i,j,QV    ) = (src(i,j,UMY) - q(i,j,QV) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
           srcQ(i,j,QREINT) = src(i,j,UEDEN) - q(i,j,QU) *src(i,j,UMX)   &
                                             - q(i,j,QV) *src(i,j,UMY) + &
                      0.5d0 * (q(i,j,QU)**2 + q(i,j,QV)**2) * srcQ(i,j,QRHO)
           srcQ(i,j,QPRES ) = dpde(i,j) * &
                (srcQ(i,j,QREINT) - q(i,j,QREINT)*srcQ(i,j,QRHO)/q(i,j,QRHO))/q(i,j,QRHO) + &
                dpdrho(i,j) * srcQ(i,j,QRHO)! + &
!                sum(dpdX_er(i,j,:)*(src(i,j,UFS:UFS+nspec-1) - &
!                    q(i,j,QFS:QFS+nspec-1)*srcQ(i,j,QRHO))) / q(i,j,QRHO)

           do ispec = 1,nspec
              srcQ(i,j,QFS+ispec-1) = ( src(i,j,UFS+ispec-1) - q(i,j,QFS+ispec-1) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
           enddo

           do iaux = 1,naux
              srcQ(i,j,QFX+iaux-1) = ( src(i,j,UFX+iaux-1) - q(i,j,QFX+iaux-1) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
           enddo

           do iadv = 1,nadv
              srcQ(i,j,QFA+iadv-1) = ( src(i,j,UFA+iadv-1) - q(i,j,QFA+iadv-1) * srcQ(i,j,QRHO) ) / q(i,j,QRHO)
           enddo

        end do
      end do

!     Compute running max of Courant number over grids
      courmx = courno
      courmy = courno
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            courx =  ( c(i,j)+abs(q(i,j,QU)) ) * dt/dx
            coury =  ( c(i,j)+abs(q(i,j,QV)) ) * dt/dy
            courmx = max( courmx, courx )
            courmy = max( courmy, coury )

            if (courx .gt. 1.d0) then
               print *,'   '
               call bl_warning("Warning:: Castro_2d.f90 :: CFL violation in ctoprim")
               print *,'>>> ... (u+c) * dt / dx > 1 ', courx
               print *,'>>> ... at cell (i,j)     : ',i,j
               print *,'>>> ... u, c                ',q(i,j,QU), c(i,j)
               print *,'>>> ... density             ',q(i,j,QRHO)
            end if

            if (coury .gt. 1.d0) then
               print *,'   '
               call bl_warning("Warning:: Castro_2d.f90 :: CFL violation in ctoprim")
               print *,'>>> ... (v+c) * dt / dx > 1 ', coury
               print *,'>>> ... at cell (i,j)     : ',i,j
               print *,'>>> ... v, c                ',q(i,j,QV), c(i,j)
               print *,'>>> ... density             ',q(i,j,QRHO)
            end if

         enddo
      enddo
      courno = max( courmx, courmy )

!     Compute flattening coef for slope calculations
      if(iflaten.eq.1)then
         do n=1,2
            loq(n)=lo(n)-ngf
            hiq(n)=hi(n)+ngf
         enddo
         call uflaten(loq,hiq, &
                      q(q_l1,q_l2,QPRES), &
                      q(q_l1,q_l2,QU), &
                      q(q_l1,q_l2,QV), &
                      flatn,q_l1,q_l2,q_h1,q_h2)
      else
         flatn = 1.d0
      endif

      deallocate(dpdrho,dpde)

      end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine trace(q,c,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                       dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                       dq,qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                       grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                       ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

      use network, only : nspec, naux
      use meth_params_module, only : iorder, QVAR, QRHO, QU, QV, &
                                     QREINT, QPRES, QFA, QFS, QFX, &
                                     nadv, small_dens, ppm_type

      implicit none

      integer ilo1,ilo2,ihi1,ihi2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
      integer gv_l1,gv_l2,gv_h1,gv_h2

      double precision dx, dy, dt
      double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
      double precision  dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)

      ! Local variables
      integer i, j
      integer n, iadv, ispec, iaux

      double precision dtdx, dtdy
      double precision cc, csq, rho, u, v, p, rhoe
      double precision drho, du, dv, dp, drhoe

      double precision enth, alpham, alphap, alpha0r, alpha0e
      double precision alpha0u, alpha0v
      double precision spminus, spplus, spzero
      double precision apright, amright, azrright, azeright
      double precision azu1rght, azv1rght
      double precision apleft, amleft, azrleft, azeleft
      double precision azu1left, azv1left
      double precision acmprght, acmpleft, acmpbot, acmptop
      double precision ascmprght, ascmpleft, ascmpbot, ascmptop
      double precision sourcr,sourcp,source,courn,eta,dlogatmp

      if (ppm_type .ne. 0) then
        print *,'Oops -- shouldnt be in trace with ppm_type != 0'
        call bl_error("Error:: Castro_2d.f90 :: trace")
      end if

      dtdx = dt/dx
      dtdy = dt/dy

         ! Compute slopes
         if (iorder .eq. 1) then

            dq(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:QVAR) = 0.d0

         else
            
            call uslope(q, &
                        flatn, qd_l1, qd_l2, qd_h1, qd_h2, &
                        dq   ,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        ilo1,ilo2,ihi1,ihi2,QVAR,1)

            call pslope(q(:,:,QPRES),q(:,:,QRHO),  &
                        flatn        , qd_l1, qd_l2, qd_h1, qd_h2, &
                        dq(:,:,QPRES),qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        grav         , gv_l1, gv_l2, gv_h1, gv_h2, &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,1)

         endif
      
         ! Compute left and right traced states
         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1
               
               cc = c(i,j)
               csq = cc**2
               rho = q(i,j,QRHO)
               u = q(i,j,QU)
               v = q(i,j,QV)
               p = q(i,j,QPRES)
               rhoe = q(i,j,QREINT)
               enth = ( (rhoe+p)/rho )/csq
               
               drho = dq(i,j,QRHO)
               du = dq(i,j,QU)
               dv = dq(i,j,QV)
               dp = dq(i,j,QPRES)
               drhoe = dq(i,j,QREINT)
               
               alpham = 0.5d0*(dp/(rho*cc) - du)*rho/cc
               alphap = 0.5d0*(dp/(rho*cc) + du)*rho/cc
               alpha0r = drho - dp/csq
               alpha0e = drhoe - dp*enth
               alpha0v = dv
               
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
               azv1rght = 0.5d0*(-1.d0 - spzero )*alpha0v

               if (i .ge. ilo1) then
                  qxp(i,j,QRHO) = rho + apright + amright + azrright
                  qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
                  qxp(i,j,QU) = u + (apright - amright)*cc/rho
                  qxp(i,j,QV) = v + azv1rght
                  qxp(i,j,QPRES) = p + (apright + amright)*csq
                  qxp(i,j,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
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
               azv1left = 0.5d0*(1.d0 - spzero )*alpha0v

               if (i .le. ihi1) then
                  qxm(i+1,j,QRHO) = rho + apleft + amleft + azrleft
                  qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                  qxm(i+1,j,QU) = u + (apleft - amleft)*cc/rho
                  qxm(i+1,j,QV) = v + azv1left
                  qxm(i+1,j,QPRES) = p + (apleft + amleft)*csq
                  qxm(i+1,j,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
               end if

               if(dloga(i,j).ne.0)then
                  courn = dtdx*(cc+abs(u))
                  eta = (1.d0-courn)/(cc*dt*abs(dloga(i,j)))
                  dlogatmp = min(eta,1.d0)*dloga(i,j)
                  sourcr = -0.5d0*dt*rho*dlogatmp*u
                  sourcp = sourcr*csq
                  source = sourcp*enth
                  if (i .le. ihi1) then
                     qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
                     qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
                     qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
                     qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
                  end if
                  if (i .ge. ilo1) then
                     qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
                     qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
                     qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
                     qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
                  end if
               endif

            enddo
         enddo

         do iadv = 1, nadv
            n = QFA + iadv - 1
            do j = ilo2-1, ihi2+1

               ! Right state
               do i = ilo1, ihi1+1
                  u = q(i,j,QU)
                  if (u .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = u*dtdx
                  endif
                  acmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qxp(i,j,n) = q(i,j,n) + acmprght
               enddo

               ! Left state
               do i = ilo1-1, ihi1
                  u = q(i,j,QU)
                  if (u .ge. 0.d0) then
                     spzero = u*dtdx
                  else
                     spzero = 1.d0
                  endif
                  acmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qxm(i+1,j,n) = q(i,j,n) + acmpleft
               enddo

            enddo
         enddo

         do ispec = 1, nspec
            n = QFS + ispec - 1
            do j = ilo2-1, ihi2+1

               ! Right state
               do i = ilo1, ihi1+1
                  u = q(i,j,QU)
                  if (u .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = u*dtdx
                  endif
                  ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qxp(i,j,n) = q(i,j,n) + ascmprght
               enddo

               ! Left state
               do i = ilo1-1, ihi1
                  u = q(i,j,QU)
                  if (u .ge. 0.d0) then
                     spzero = u*dtdx
                  else
                     spzero = 1.d0
                  endif
                  ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qxm(i+1,j,n) = q(i,j,n) + ascmpleft
               enddo

            enddo
         enddo

         do iaux = 1, naux
            n = QFX + iaux - 1
            do j = ilo2-1, ihi2+1

               ! Right state
               do i = ilo1, ihi1+1
                  u = q(i,j,QU)
                  if (u .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = u*dtdx
                  endif
                  ascmprght = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qxp(i,j,n) = q(i,j,n) + ascmprght
               enddo

               ! Left state
               do i = ilo1-1, ihi1
                  u = q(i,j,QU)
                  if (u .ge. 0.d0) then
                     spzero = u*dtdx
                  else
                     spzero = 1.d0
                  endif
                  ascmpleft = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qxm(i+1,j,n) = q(i,j,n) + ascmpleft
               enddo

            enddo
         enddo

         ! Compute slopes
         if (iorder .ne. 1) then

            call uslope(q,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        ilo1,ilo2,ihi1,ihi2,QVAR,2)

            call pslope(q(:,:,QPRES),q(:,:,QRHO),  &
                        flatn        , qd_l1, qd_l2, qd_h1, qd_h2, &
                        dq(:,:,QPRES),qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        grav         , gv_l1, gv_l2, gv_h1, gv_h2, &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,2)

         endif

         do j = ilo2-1, ihi2+1
            do i = ilo1-1, ihi1+1

               cc = c(i,j)
               csq = cc**2
               rho = q(i,j,QRHO)
               u = q(i,j,QU)
               v = q(i,j,QV)
               p = q(i,j,QPRES)
               rhoe = q(i,j,QREINT)
               enth = ( (rhoe+p)/rho )/csq

               drho = dq(i,j,QRHO)
               du = dq(i,j,QU)
               dv = dq(i,j,QV)
               dp = dq(i,j,QPRES)
               drhoe = dq(i,j,QREINT)

               alpham = 0.5d0*(dp/(rho*cc) - dv)*rho/cc
               alphap = 0.5d0*(dp/(rho*cc) + dv)*rho/cc
               alpha0r = drho - dp/csq
               alpha0e = drhoe - dp*enth
               alpha0u = du

               if (v-cc .gt. 0.d0) then
                  spminus = -1.d0
               else
                  spminus = (v-cc)*dtdy
               endif
               if (v+cc .gt. 0.d0) then
                  spplus = -1.d0
               else
                  spplus = (v+cc)*dtdy
               endif
               if (v .gt. 0.d0) then
                  spzero = -1.d0
               else
                  spzero = v*dtdy
               endif

               apright = 0.5d0*(-1.d0 - spplus )*alphap
               amright = 0.5d0*(-1.d0 - spminus)*alpham
               azrright = 0.5d0*(-1.d0 - spzero )*alpha0r
               azeright = 0.5d0*(-1.d0 - spzero )*alpha0e
               azu1rght = 0.5d0*(-1.d0 - spzero )*alpha0u

               if (j .ge. ilo2) then
                  qyp(i,j,QRHO) = rho + apright + amright + azrright
                  qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
                  qyp(i,j,QV) = v + (apright - amright)*cc/rho
                  qyp(i,j,QU) = u + azu1rght
                  qyp(i,j,QPRES) = p + (apright + amright)*csq
                  qyp(i,j,QREINT) = rhoe + (apright + amright)*enth*csq + azeright
               end if

               if (v-cc .ge. 0.d0) then
                  spminus = (v-cc)*dtdy
               else
                  spminus = 1.d0
               endif
               if (v+cc .ge. 0.d0) then
                  spplus = (v+cc)*dtdy
               else
                  spplus = 1.d0
               endif
               if (v .ge. 0.d0) then
                  spzero = v*dtdy
               else
                  spzero = 1.d0
               endif

               apleft = 0.5d0*(1.d0 - spplus )*alphap
               amleft = 0.5d0*(1.d0 - spminus)*alpham
               azrleft = 0.5d0*(1.d0 - spzero )*alpha0r
               azeleft = 0.5d0*(1.d0 - spzero )*alpha0e
               azu1left = 0.5d0*(1.d0 - spzero )*alpha0u

               if (j .le. ihi2) then
                  qym(i,j+1,QRHO) = rho + apleft + amleft + azrleft
                  qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
                  qym(i,j+1,QV) = v + (apleft - amleft)*cc/rho
                  qym(i,j+1,QU) = u + azu1left
                  qym(i,j+1,QPRES) = p + (apleft + amleft)*csq
                  qym(i,j+1,QREINT) = rhoe + (apleft + amleft)*enth*csq + azeleft
               end if

            enddo
         enddo

         do iadv = 1, nadv
            n = QFA + iadv - 1
            do i = ilo1-1, ihi1+1

               ! Top state
               do j = ilo2, ihi2+1
                  v = q(i,j,QV)
                  if (v .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = v*dtdy
                  endif
                  acmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qyp(i,j,n) = q(i,j,n) + acmptop
               enddo

               ! Bottom state
               do j = ilo2-1, ihi2
                  v = q(i,j,QV)
                  if (v .ge. 0.d0) then
                     spzero = v*dtdy
                  else
                     spzero = 1.d0
                  endif
                  acmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qym(i,j+1,n) = q(i,j,n) + acmpbot
               enddo

            enddo
         enddo

         do ispec = 1, nspec
            n = QFS + ispec - 1
            do i = ilo1-1, ihi1+1

               ! Top state
               do j = ilo2, ihi2+1
                  v = q(i,j,QV)
                  if (v .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = v*dtdy
                  endif
                  ascmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qyp(i,j,n) = q(i,j,n) + ascmptop
               enddo

               ! Bottom state
               do j = ilo2-1, ihi2
                  v = q(i,j,QV)
                  if (v .ge. 0.d0) then
                     spzero = v*dtdy
                  else
                     spzero = 1.d0
                  endif
                  ascmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qym(i,j+1,n) = q(i,j,n) + ascmpbot
               enddo

            enddo
         enddo

         do iaux = 1, naux
            n = QFX + iaux - 1
            do i = ilo1-1, ihi1+1

               ! Top state
               do j = ilo2, ihi2+1
                  v = q(i,j,QV)
                  if (v .gt. 0.d0) then
                     spzero = -1.d0
                  else
                     spzero = v*dtdy
                  endif
                  ascmptop = 0.5d0*(-1.d0 - spzero )*dq(i,j,n)
                  qyp(i,j,n) = q(i,j,n) + ascmptop
               enddo

               ! Bottom state
               do j = ilo2-1, ihi2
                  v = q(i,j,QV)
                  if (v .ge. 0.d0) then
                     spzero = v*dtdy
                  else
                     spzero = 1.d0
                  endif
                  ascmpbot = 0.5d0*(1.d0 - spzero )*dq(i,j,n)
                  qym(i,j+1,n) = q(i,j,n) + ascmpbot
               enddo

            enddo
         enddo

      end subroutine trace

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine consup( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                        uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                        pgdx,pgdx_l1,pgdx_l2,pgdx_h1,pgdx_h2, &
                        pgdy,pgdy_l1,pgdy_l2,pgdy_h1,pgdy_h2, &
                        src , src_l1, src_l2, src_h1, src_h2, &
                        grav,  gv_l1,  gv_l2,  gv_h1,  gv_h2, &
                        flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                        flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                        area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                        area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                        vol,vol_l1,vol_l2,vol_h1,vol_h2, &
                        div,pdivu,lo,hi,dx,dy,dt)

      use eos_module
      use network, only : nspec, naux
      use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                     normalize_species

      implicit none

      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer uout_l1,uout_l2,uout_h1,uout_h2
      integer pgdx_l1,pgdx_l2,pgdx_h1,pgdx_h2
      integer pgdy_l1,pgdy_l2,pgdy_h1,pgdy_h2
      integer   src_l1,  src_l2,  src_h1,  src_h2
      integer    gv_l1,   gv_l2,   gv_h1,   gv_h2
      integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
      integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
      integer area1_l1,area1_l2,area1_h1,area1_h2
      integer area2_l1,area2_l2,area2_h1,area2_h2
      integer vol_l1,vol_l2,vol_h1,vol_h2

      double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
      double precision pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
      double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,NVAR)
      double precision  grav(   gv_l1:   gv_h1,   gv_l2:   gv_h2,2)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
      double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
      double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
      double precision pdivu(lo(1):hi(1),lo(2):hi(2))
      double precision dx, dy, dt

      integer i, j, n

      double precision div1
      double precision SrU, SrV
      double precision rho, Up, Vp, SrE

      ! Normalize the species fluxes
      if (normalize_species .eq. 1) &
         call normalize_species_fluxes( &
                  flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                  flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                  lo,hi)

      do n = 1, NVAR
         if ( n.eq.UTEMP) then
            flux1(:,:,n) = 0.d0
            flux2(:,:,n) = 0.d0
         else 
            do j = lo(2),hi(2)
               do i = lo(1),hi(1)+1
                  div1 = .5d0*(div(i,j) + div(i,j+1))
                  div1 = difmag*min(0.d0,div1)
                  flux1(i,j,n) = flux1(i,j,n) &
                       + dx*div1*(uin(i,j,n) - uin(i-1,j,n))
                  flux1(i,j,n) = area1(i,j)*flux1(i,j,n)
               enddo
            enddo

            do j = lo(2),hi(2)+1
               do i = lo(1),hi(1)
                  div1 = .5d0*(div(i,j) + div(i+1,j))
                  div1 = difmag*min(0.d0,div1)
                  flux2(i,j,n) = flux2(i,j,n) &
                       + dy*div1*(uin(i,j,n) - uin(i,j-1,n))
                  flux2(i,j,n) = area2(i,j)*flux2(i,j,n)
               enddo
            enddo
         endif
      enddo
         
      do n = 1, NVAR
         if (n .eq. UTEMP) then
            uout(lo(1):hi(1),lo(2):hi(2),n) = uin(lo(1):hi(1),lo(2):hi(2),n)
         else 
            do j = lo(2),hi(2)
            do i = lo(1),hi(1)
               uout(i,j,n) = uin(i,j,n) + dt * &
                      ( flux1(i,j,n) - flux1(i+1,j,n) &
                    +   flux2(i,j,n) - flux2(i,j+1,n) ) / vol(i,j) &
                    +   dt * src(i,j,n)
            enddo
            enddo
         end if
      enddo

      ! Add gradp term to momentum equation
      do j = lo(2),hi(2)
      do i = lo(1),hi(1)
!         uout(i,j,UMX) = uout(i,j,UMX)+ 0.5d0*(area1(i,j)+area1(i+1,j))* &
!            dt * ( pgdx(i,j)-pgdx(i+1,j) )/vol(i,j)
!         uout(i,j,UMY) = uout(i,j,UMY)+ 0.5d0*(area2(i,j)+area2(i,j+1))* &
!            dt * ( pgdy(i,j)-pgdy(i,j+1) )/vol(i,j)

          uout(i,j,UMX) = uout(i,j,UMX) - dt * (pgdx(i+1,j)-pgdx(i,j))/ dx
          uout(i,j,UMY) = uout(i,j,UMY) - dt * (pgdy(i,j+1)-pgdy(i,j))/ dy
      enddo
      enddo

      ! Add source term to (rho e)
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            uout(i,j,UEINT) = uout(i,j,UEINT)  - dt * pdivu(i,j)
         enddo
      enddo

      ! Add gravitational source terms to momentum and energy equations
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

            rho = uin(i,j,URHO)
            Up  = uin(i,j,UMX) / rho
            Vp  = uin(i,j,UMY) / rho

            SrU = rho * grav(i,j,1)
            SrV = rho * grav(i,j,2)

            ! This doesn't work (in 1-d)
            ! SrE = SrU*(Up + SrU*dt/(2*rho)) &
            !      +SrV*(Vp + SrV*dt/(2*rho))

            ! This does work (in 1-d)
            SrE = uin(i,j,UMX) * grav(i,j,1) + uin(i,j,UMY) * grav(i,j,2)

            uout(i,j,UMX)   = uout(i,j,UMX)   + dt * SrU
            uout(i,j,UMY)   = uout(i,j,UMY)   + dt * SrV
            uout(i,j,UEDEN) = uout(i,j,UEDEN) + dt * SrE

         enddo
      enddo

      do j = lo(2),hi(2)
      do i = lo(1),hi(1)+1
        flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR)
        flux1(i,j,   UMX) = flux1(i,j,UMX) + dt*area1(i,j)*pgdx(i,j)
      enddo
      enddo

      do j = lo(2),hi(2)+1 
      do i = lo(1),hi(1)
        flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR)
        flux2(i,j,UMY) = flux2(i,j,UMY) + dt*area2(i,j)*pgdy(i,j)
      enddo
      enddo

      end subroutine consup

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine cmpflx(qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        flx,flx_l1,flx_l2,flx_h1,flx_h2, &
                        pgd,pgd_l1,pgd_l2,pgd_h1,pgd_h2, &
                        ugd,ugd_l1,ugd_l2,ugd_h1,ugd_h2, &
                        gamc,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
                        idir,ilo,ihi,jlo,jhi)

      use meth_params_module, only : QVAR, NVAR

      implicit none

      integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
      integer flx_l1,flx_l2,flx_h1,flx_h2
      integer pgd_l1,pgd_l2,pgd_h1,pgd_h2
      integer ugd_l1,ugd_l2,ugd_h1,ugd_h2
      integer qd_l1,qd_l2,qd_h1,qd_h2
      integer idir,ilo,ihi,jlo,jhi

      double precision    qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision    qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision   flx(flx_l1:flx_h1,flx_l2:flx_h2,NVAR)
      double precision pgd(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
      double precision ugd(ugd_l1:ugd_h1,ugd_l2:ugd_h2)
      double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
      double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2)

!     Local variables
      integer i, j

      double precision, allocatable :: smallc(:,:), cavg(:,:)
      double precision, allocatable :: gamcm(:,:), gamcp(:,:)

      allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
      allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
      allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
      allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )

      if(idir.eq.1) then
         do j = jlo, jhi
            do i = ilo, ihi+1
               smallc(i,j) = max( csml(i,j), csml(i-1,j) )
               cavg(i,j) = 0.5d0*( c(i,j) + c(i-1,j) )
               gamcm(i,j) = gamc(i-1,j)
               gamcp(i,j) = gamc(i,j)
            enddo
         enddo
      else
         do j = jlo, jhi+1
            do i = ilo, ihi
               smallc(i,j) = max( csml(i,j), csml(i,j-1) )
               cavg(i,j) = 0.5d0*( c(i,j) + c(i,j-1) )
               gamcm(i,j) = gamc(i,j-1)
               gamcp(i,j) = gamc(i,j)
            enddo
         enddo
      endif

!     Solve Riemann problem (godunov state passed back, but only (u,p) saved
      call riemannus(qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                     gamcm, gamcp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
                     flx, flx_l1, flx_l2, flx_h1, flx_h2, &
                     pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                     ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                     idir, ilo, ihi, jlo, jhi)

      deallocate(smallc,cavg,gamcm,gamcp)

      end subroutine cmpflx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine riemannus(ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
                           gamcl, gamcr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
                           uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
                           pgdnv, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
                           ugdnv, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
                           idir, ilo1, ihi1, ilo2, ihi2)

      use network, only : nspec, naux
      use prob_params_module, only : physbc_lo,physbc_hi,Symmetry
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv, &
                                     small_dens, small_pres

      implicit none

      double precision, parameter:: small = 1.d-8

      integer qpd_l1, qpd_l2, qpd_h1, qpd_h2
      integer gd_l1, gd_l2, gd_h1, gd_h2
      integer uflx_l1, uflx_l2, uflx_h1, uflx_h2
      integer pgd_l1, pgd_l2, pgd_h1, pgd_h2
      integer ugd_l1, ugd_l2, ugd_h1, ugd_h2
      integer idir, ilo1, ihi1, ilo2, ihi2
      integer ilo,ihi,jlo,jhi

      double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QVAR)
      double precision gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision cav(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision smallc(gd_l1:gd_h1,gd_l2:gd_h2)
      double precision uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
      double precision pgdnv(pgd_l1:pgd_h1,pgd_l2:pgd_h2)
      double precision ugdnv(ugd_l1:ugd_h1,ugd_l2:ugd_h2)

      integer iadv, ispec, iaux, n, nq
      integer i, j

      double precision rgd, vgd, regd, ustar
      double precision rl, ul, vl, pl, rel
      double precision rr, ur, vr, pr, rer
      double precision wl, wr, rhoetot, scr
      double precision rstar, cstar, estar, pstar
      double precision ro, uo, po, reo, co, gamco, entho
      double precision sgnm, spin, spout, ushock, frac
      double precision wsmall, csmall,qavg

!************************************************************
!  set min/max based on normal direction
      if(idir.eq.1) then
         ilo = ilo1
         ihi = ihi1 + 1
         jlo = ilo2
         jhi = ihi2
      else
         ilo = ilo1
         ihi = ihi1
         jlo = ilo2
         jhi = ihi2+1
      endif

!     Solve Riemann Problem
!     NOTE: The calling routine will order velocity unknowns so that
!     for the purposes of this routine, the normal component is always
!     loaded in the QU slot.
      do j = jlo, jhi
         do i = ilo, ihi

            rl = ql(i,j,QRHO)

!  pick left velocities based on direction
            if(idir.eq.1) then
               ul = ql(i,j,QU)
               vl = ql(i,j,QV)
            else
               ul = ql(i,j,QV)
               vl = ql(i,j,QU)
            endif

            pl = ql(i,j,QPRES)
            rel = ql(i,j,QREINT)

            rr = qr(i,j,QRHO)

!  pick right velocities based on direction
            if(idir.eq.1) then
               ur = qr(i,j,QU)
               vr = qr(i,j,QV)
            else
               ur = qr(i,j,QV)
               vr = qr(i,j,QU)
            endif

            pr = qr(i,j,QPRES)
            rer = qr(i,j,QREINT)

            csmall = smallc(i,j)
            wsmall = small_dens*csmall
            wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
            wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

            pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
            pstar = max(pstar,small_pres)
            ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

            if (ustar .gt. 0.d0) then
               ro = rl
               uo = ul
               po = pl
               reo = rel
               gamco = gamcl(i,j)
            else if (ustar .lt. 0.d0) then
               ro = rr
               uo = ur
               po = pr
               reo = rer
               gamco = gamcr(i,j)
            else
               ro = 0.5d0*(rl+rr)
               uo = 0.5d0*(ul+ur)
               po = 0.5d0*(pl+pr)
               reo = 0.5d0*(rel+rer)
               gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))               
            endif
            ro = max(small_dens,ro)

            co = sqrt(abs(gamco*po/ro))
            co = max(csmall,co)
            entho = (reo/ro + po/ro)/co**2
            rstar = ro + (pstar - po)/co**2
            rstar = max(small_dens,rstar)
            estar = reo + (pstar - po)*entho
            cstar = sqrt(abs(gamco*pstar/rstar))
            cstar = max(cstar,csmall)

            sgnm = sign(1.d0,ustar)
            spout = co - sgnm*uo
            spin = cstar - sgnm*ustar
            ushock = 0.5d0*(spin + spout)
            if (pstar-po .ge. 0.d0) then
               spin = ushock
               spout = ushock
            endif
            if (spout-spin .eq. 0.d0) then
               scr = small*cav(i,j)
            else
               scr = spout-spin
            endif
            frac = (1.d0 + (spout + spin)/scr)*0.5d0
            frac = max(0.d0,min(1.d0,frac))

            if (ustar .gt. 0.d0) then
               vgd = vl
            else if (ustar .lt. 0.d0) then
               vgd = vr
            else
               vgd = 0.5d0*(vl+vr)
            endif
            rgd = frac*rstar + (1.d0 - frac)*ro

            ugdnv(i,j) = frac*ustar + (1.d0 - frac)*uo
            pgdnv(i,j) = frac*pstar + (1.d0 - frac)*po

            regd = frac*estar + (1.d0 - frac)*reo
            if (spout .lt. 0.d0) then
               rgd = ro
               ugdnv(i,j) = uo
               pgdnv(i,j) = po
               regd = reo
            endif
            if (spin .ge. 0.d0) then
               rgd = rstar
               ugdnv(i,j) = ustar
               pgdnv(i,j) = pstar
               regd = estar
            endif

            ! Enforce that fluxes through a symmetry plane are hard zero.
            if (i.eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) ugdnv(i,j) = 0.d0
            if (j.eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) ugdnv(i,j) = 0.d0

            ! Compute fluxes, order as conserved state (not q)
            uflx(i,j,URHO) = rgd*ugdnv(i,j)
            if(idir.eq.1) then
               uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j)
               uflx(i,j,UMY) = uflx(i,j,URHO)*vgd
            else
               uflx(i,j,UMX) = uflx(i,j,URHO)*vgd
               uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j)
            endif

            rhoetot = regd + 0.5d0*rgd*(ugdnv(i,j)**2 + vgd**2)
            uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv(i,j))
            uflx(i,j,UEINT) = ugdnv(i,j)*regd

            do iadv = 1, nadv
               n = UFA + iadv - 1
               nq = QFA + iadv - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  uflx(i,j,n) = uflx(i,j,URHO)*qavg
               endif
            enddo

            do ispec = 1, nspec
               n  = UFS + ispec - 1
               nq = QFS + ispec - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  uflx(i,j,n) = uflx(i,j,URHO)*qavg
               endif
            enddo

            do iaux = 1, naux
               n  = UFX + iaux - 1
               nq = QFX + iaux - 1
               if (ustar .gt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*ql(i,j,nq)
               else if (ustar .lt. 0.d0) then
                  uflx(i,j,n) = uflx(i,j,URHO)*qr(i,j,nq)
               else 
                  qavg = 0.5d0 * (ql(i,j,nq) + qr(i,j,nq))
                  uflx(i,j,n) = uflx(i,j,URHO)*qavg
               endif
            enddo

         enddo
      enddo
      end subroutine riemannus

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transx(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                        fx, fx_l1, fx_l2, fx_h1, fx_h2, &
                        pgdx, pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
                        ugdx, ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2, &
                        gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                        srcQ, src_l1, src_l2, src_h1, src_h2, &
                        grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                        hdt, cdtdx,  &
                        area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                        vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                        ilo, ihi, jlo, jhi)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, &
                                     nadv
      implicit none

      integer qd_l1, qd_l2, qd_h1, qd_h2
      integer gc_l1, gc_l2, gc_h1, gc_h2
      integer fx_l1, fx_l2, fx_h1, fx_h2
      integer pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2
      integer ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2
      integer src_l1, src_l2, src_h1, src_h2
      integer gv_l1, gv_l2, gv_h1, gv_h2
      integer area1_l1, area1_l2, area1_h1, area1_h2
      integer vol_l1, vol_l2, vol_h1, vol_h2
      integer ilo, ihi, jlo, jhi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision fx(fx_l1:fx_h1,fx_l2:fx_h2,NVAR)
      double precision pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
      double precision ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
      double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
      double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
      double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
      double precision hdt, cdtdx

      integer i, j
      integer n, nq, iadv
      integer ispec, iaux

      double precision rr, rrnew, compo, compn
      double precision rrr, rur, rvr, rer, ekinr, rhoekinr
      double precision rrnewr, runewr, rvnewr, renewr
      double precision rrl, rul, rvl, rel, ekinl, rhoekinl
      double precision rrnewl, runewl, rvnewl, renewl
      double precision pgp, pgm, ugp, ugm, dup, pav, du, pnewl,pnewr
      double precision rhotmp

      ! NOTE: it is better *not* to protect against small density in this routine

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i,j,  QRHO)
                  rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                    area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

                  compo = rr*qp(i,j  ,nq)
                  compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                       area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

                  qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i,j+1,QRHO)
                  rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                    area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

                  compo = rr*qm(i,j+1,nq)
                  compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                       area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

                  qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do ispec = 1, nspec
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i  ,j,QRHO)
                  rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                    area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

                  compo = rr*qp(i,j  ,nq)
                  compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                       area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

                  qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i,j+1,QRHO)
                  rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                    area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

                  compo  = rr*qm(i,j+1,nq)
                  compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                       area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

                  qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do iaux = 1, naux
          n  = UFX + iaux - 1
          nq = QFX + iaux - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i  ,j,QRHO)
                  rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                    area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

                  compo = rr*qp(i,j  ,nq)
                  compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                       area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

                  qpo(i,j  ,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i,j+1,QRHO)
                  rrnew = rr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                    area1(i  ,j)*fx(i  ,j,URHO))/vol(i,j) 

                  compo = rr*qm(i,j+1,nq)
                  compn = compo - hdt*(area1(i+1,j)*fx(i+1,j,n)- &
                                       area1(i  ,j)*fx(i  ,j,n))/vol(i,j) 

                  qmo(i,j+1,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do j = jlo, jhi 
         do i = ilo, ihi 

            pgp = pgdx(i+1,j)
            pgm = pgdx(i,j)
            ugp = ugdx(i+1,j)
            ugm = ugdx(i,j)

!           Convert to conservation form
            rrr = qp(i,j,QRHO)
            rur = rrr*qp(i,j,QU)
            rvr = rrr*qp(i,j,QV)
            ekinr = 0.5d0*rrr*(qp(i,j,QU)**2 + qp(i,j,QV)**2)
            rer = qp(i,j,QREINT) + ekinr

            rrl = qm(i,j+1,QRHO)
            rul = rrl*qm(i,j+1,QU)
            rvl = rrl*qm(i,j+1,QV)
            ekinl = 0.5d0*rrl*(qm(i,j+1,QU)**2 + qm(i,j+1,QV)**2)
            rel = qm(i,j+1,QREINT) + ekinl

!           Add transverse predictor
            rrnewr = rrr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                area1(i,j)*fx(i,j,URHO))/vol(i,j) 
            runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                 -0.5*hdt*(area1(i+1,j)+area1(i,j))*(pgp-pgm)/vol(i,j) 
            rvnewr = rvr - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                                area1(i,j)*fx(i,j,UMY))/vol(i,j) 
            renewr = rer - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                                area1(i,j)*fx(i,j,UEDEN))/vol(i,j) 

            rrnewl = rrl - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
                                area1(i,j)*fx(i,j,URHO))/vol(i,j) 
            runewl = rul - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
                                area1(i,j)*fx(i,j,UMX))/vol(i,j) &
                   -0.5*hdt*(area1(i+1,j)+area1(i,j))*(pgp-pgm)/vol(i,j) 
            rvnewl = rvl - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
                                area1(i,j)*fx(i,j,UMY))/vol(i,j) 
            renewl = rel - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
                                area1(i,j)*fx(i,j,UEDEN))/vol(i,j) 

            dup = pgp*ugp - pgm*ugm
            pav = 0.5d0*(pgp+pgm)
            du = ugp-ugm

            pnewr = qp(i,j,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j)-1.d0))

!           Convert back to non-conservation form
            rhotmp = rrnewr
            qpo(i,j,QRHO) = rhotmp        + hdt*srcQ(i,j,QRHO)
            qpo(i,j,QU  ) = runewr/rhotmp + hdt*srcQ(i,j,QU)  + hdt*grav(i,j,1)
            qpo(i,j,QV  ) = rvnewr/rhotmp + hdt*srcQ(i,j,QV)  + hdt*grav(i,j,2)
            rhoekinr = 0.5d0*(runewr**2+rvnewr**2)/rhotmp
            qpo(i,j,QREINT)= renewr - rhoekinr + hdt*srcQ(i,j,QREINT)
            qpo(i,j,QPRES) =  pnewr            + hdt*srcQ(i,j,QPRES)

            pnewl = qm(i,j+1,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j)-1.d0))

!           Convert back to non-conservation form
            rhotmp = rrnewl
            qmo(i,j+1,QRHO) = rhotmp         + hdt*srcQ(i,j,QRHO)
            qmo(i,j+1,QU  ) = runewl/rhotmp  + hdt*srcQ(i,j,QU)  + hdt*grav(i,j,1)
            qmo(i,j+1,QV  ) = rvnewl/rhotmp  + hdt*srcQ(i,j,QV)  + hdt*grav(i,j,2)
            rhoekinl = 0.5d0*(runewl**2+rvnewl**2)/rhotmp
            qmo(i,j+1,QREINT)= renewl - rhoekinl +hdt*srcQ(i,j,QREINT)
            qmo(i,j+1,QPRES) = pnewl +hdt*srcQ(i,j,QPRES)

        enddo
      enddo

      end subroutine transx

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine transy(qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
                        fy,fy_l1,fy_l2,fy_h1,fy_h2, &
                        pgdy, pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
                        ugdy, ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2, &
                        gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
                        srcQ, src_l1, src_l2, src_h1, src_h2, &
                        grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                        hdt, cdtdy, ilo, ihi, jlo, jhi)

      use network, only : nspec, naux
      use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
                                     URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, &
                                     nadv
      implicit none

      integer qd_l1, qd_l2, qd_h1, qd_h2
      integer gc_l1, gc_l2, gc_h1, gc_h2
      integer fy_l1, fy_l2, fy_h1, fy_h2
      integer pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2
      integer ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2
      integer src_l1, src_l2, src_h1, src_h2
      integer gv_l1, gv_l2, gv_h1, gv_h2
      integer ilo, ihi, jlo, jhi

      double precision qm(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qp(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QVAR)
      double precision fy(fy_l1:fy_h1,fy_l2:fy_h2,NVAR)
      double precision pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
      double precision ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
      double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
      double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
      double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
      double precision hdt, cdtdy

      integer i, j
      integer n, nq, iadv, ispec, iaux

      double precision rr,rrnew
      double precision pgp, pgm, ugp, ugm, dup, pav, du, pnewr,pnewl
      double precision rrr, rur, rvr, rer, ekinr, rhoekinr
      double precision rrnewr, runewr, rvnewr, renewr
      double precision rrl, rul, rvl, rel, ekinl, rhoekinl
      double precision rrnewl, runewl, rvnewl, renewl
      double precision rhotmp
      double precision compo, compn

      ! NOTE: it is better *not* to protect against small density in this routine

      do iadv = 1, nadv
          n  = UFA + iadv - 1
          nq = QFA + iadv - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i,j,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

                  compo = rr*qp(i,j,nq)
                  compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

                  qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i+1,j,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

                  compo = rr*qm(i+1,j,nq)
                  compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

                  qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do ispec = 1, nspec 
          n  = UFS + ispec - 1
          nq = QFS + ispec - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i,j,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

                  compo = rr*qp(i,j,nq)
                  compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

                  qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i+1,j,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

                  compo = rr*qm(i+1,j,nq)
                  compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

                  qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do iaux = 1, naux 
          n  = UFX + iaux - 1
          nq = QFX + iaux - 1
          do j = jlo, jhi 
              do i = ilo, ihi 

                  rr = qp(i,j,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

                  compo = rr*qp(i,j,nq)
                  compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

                  qpo(i,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

                  rr = qm(i+1,j,QRHO)
                  rrnew = rr - cdtdy*(fy(i,j+1,URHO)-fy(i,j,URHO)) 

                  compo = rr*qm(i+1,j,nq)
                  compn = compo - cdtdy*(fy(i,j+1,n)-fy(i,j,n)) 

                  qmo(i+1,j,nq) = compn/rrnew + hdt*srcQ(i,j,nq)

              enddo
          enddo
      enddo

      do j = jlo, jhi 
         do i = ilo, ihi 

            pgp = pgdy(i,j+1)
            pgm = pgdy(i,j)
            ugp = ugdy(i,j+1)
            ugm = ugdy(i,j)

!           Convert to conservation form
            rrr = qp(i,j,QRHO)
            rur = rrr*qp(i,j,QU)
            rvr = rrr*qp(i,j,QV)
            ekinr = 0.5d0*rrr*(qp(i,j,QU)**2 + qp(i,j,QV)**2)
            rer = qp(i,j,QREINT) + ekinr

            rrl = qm(i+1,j,QRHO)
            rul = rrl*qm(i+1,j,QU)
            rvl = rrl*qm(i+1,j,QV)
            ekinl = 0.5d0*rrl*(qm(i+1,j,QU)**2 + qm(i+1,j,QV)**2)
            rel = qm(i+1,j,QREINT) + ekinl

!           Add transverse predictor
            rrnewr = rrr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO)) 

            runewr = rur - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX)) 
            rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                 -cdtdy*(pgp-pgm) 
            renewr = rer - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) 

            rrnewl = rrl - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO)) 
            runewl = rul - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX)) 
            rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
                 -cdtdy*(pgp-pgm) 
            renewl = rel - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) 

            dup = pgp*ugp - pgm*ugm
            pav = 0.5d0*(pgp+pgm)
            du = ugp-ugm
            pnewr = qp(i  ,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-1.d0))
            pnewl = qm(i+1,j,QPRES)-cdtdy*(dup + pav*du*(gamc(i,j)-1.d0))

!           convert back to non-conservation form
            rhotmp =  rrnewr
            qpo(i,j,QRHO  ) = rhotmp           + hdt*srcQ(i,j,QRHO)
            qpo(i,j,QU    ) = runewr/rhotmp    + hdt*srcQ(i,j,QU) + hdt*grav(i,j,1)
            qpo(i,j,QV    ) = rvnewr/rhotmp    + hdt*srcQ(i,j,QV) + hdt*grav(i,j,2)
            rhoekinr = 0.5d0*(runewr**2+rvnewr**2)/rhotmp
            qpo(i,j,QREINT) = renewr - rhoekinr + hdt*srcQ(i,j,QREINT)
            qpo(i,j,QPRES ) =  pnewr            + hdt*srcQ(i,j,QPRES)

            rhotmp =  rrnewl
            qmo(i+1,j,QRHO  ) = rhotmp            + hdt*srcQ(i,j,QRHO)
            qmo(i+1,j,QU    ) = runewl/rhotmp     + hdt*srcQ(i,j,QU) + hdt*grav(i,j,1)
            qmo(i+1,j,QV    ) = rvnewl/rhotmp     + hdt*srcQ(i,j,QV) + hdt*grav(i,j,2)
            rhoekinl = 0.5d0*(runewl**2+rvnewl**2)/rhotmp
            qmo(i+1,j,QREINT) = renewl - rhoekinl + hdt*srcQ(i,j,QREINT)
            qmo(i+1,j,QPRES ) = pnewl             + hdt*srcQ(i,j,QPRES)

         enddo
      enddo

      end subroutine transy

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uslope(q,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dq,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        ilo1,ilo2,ihi1,ihi2,nv,idir)
        
        implicit none

        integer ilo,ihi
        integer qd_l1,qd_l2,qd_h1,qd_h2
        integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
        integer ilo1,ilo2,ihi1,ihi2,nv,idir

        double precision     q( qd_l1: qd_h1, qd_l2: qd_h2,nv)
        double precision flatn( qd_l1: qd_h1, qd_l2: qd_h2)
        double precision    dq(qpd_l1:qpd_h1,qpd_l2:qpd_h2,nv)

        ! local
        double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

        integer i, j, n
        double precision dlft, drgt, dq1
        double precision four3rd, sixth

        four3rd = 4.d0/3.d0
        sixth = 1.d0/6.d0

        ilo = MIN(ilo1,ilo2)
        ihi = MAX(ihi1,ihi2)

        allocate (dsgn(ilo-2:ihi+2))
        allocate (dlim(ilo-2:ihi+2))
        allocate (  df(ilo-2:ihi+2))
        allocate (dcen(ilo-2:ihi+2))

        do n = 1, nv 

           if (idir .eq. 1) then

              ! Slopes in first coordinate direction
              do j = ilo2-1, ihi2+1
                 ! First compute Fromm slopes
                 do i = ilo1-2, ihi1+2
                    dlft = q(i  ,j,n) - q(i-1,j,n)
                    drgt = q(i+1,j,n) - q(i  ,j,n)

                    dcen(i) = 0.5d0*(dlft+drgt)
                    dsgn(i) = sign(1.d0, dcen(i))

                    if (dlft*drgt .ge. 0.d0) then
                       dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
                    else
                       dlim(i) = 0.d0
                    endif
                    df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
                 enddo

                 ! now limited fourth order slopes
                 do i = ilo1-1, ihi1+1
                    dq1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                    dq(i,j,n) = flatn(i,j)*dsgn(i)*min(dlim(i),abs(dq1))
                 enddo
              enddo

           else

              ! Slopes in second coordinate direction
              do i = ilo1-1, ihi1+1
                 ! First compute Fromm slopes
                 do j = ilo2-2, ihi2+2
                    dlft = q(i,j  ,n) - q(i,j-1,n)
                    drgt = q(i,j+1,n) - q(i,j  ,n)

                    dcen(j) = 0.5d0*(dlft+drgt)
                    dsgn(j) = sign( 1.d0, dcen(j) )

                    if (dlft*drgt .ge. 0.d0) then
                       dlim(j) = 2.d0 * min( abs(dlft), abs(drgt) )
                    else
                       dlim(j) = 0.d0
                    endif
                    df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
                 enddo

                 ! now limited fourth order slopes
                 do j = ilo2-1, ihi2+1
                    dq1 = four3rd*dcen(j) - sixth*(df(j+1) + df(j-1))
                    dq(i,j,n) = flatn(i,j)*dsgn(j)*min(dlim(j),abs(dq1))
                 enddo
              enddo

           endif

        enddo

        deallocate(dsgn,dlim,df,dcen)

      end subroutine uslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine pslope(p,rho,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
                        grav,gv_l1,gv_l2,gv_h1,gv_h2, &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,idir)
        
        implicit none

        integer ilo,ihi
        integer qd_l1,qd_l2,qd_h1,qd_h2
        integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
        integer gv_l1,gv_l2,gv_h1,gv_h2
        integer ilo1,ilo2,ihi1,ihi2,idir

        double precision, intent(in   ) ::      p( qd_l1: qd_h1, qd_l2: qd_h2)
        double precision, intent(in   ) ::    rho( qd_l1: qd_h1, qd_l2: qd_h2)
        double precision, intent(in   ) ::  flatn( qd_l1: qd_h1, qd_l2: qd_h2)
        double precision, intent(  out) ::     dp(qpd_l1:qpd_h1,qpd_l2:qpd_h2)
        double precision, intent(in   ) ::   grav( gv_l1: gv_h1, gv_l2: gv_h2,2)
        double precision, intent(in   ) ::  dx,dy

        ! local
        double precision, allocatable::dsgn(:),dlim(:),df(:),dcen(:)

        integer i, j
        double precision dlft, drgt, dp1
        double precision four3rd, sixth

        four3rd = 4.d0/3.d0
        sixth = 1.d0/6.d0

        ilo = MIN(ilo1,ilo2)
        ihi = MAX(ihi1,ihi2)

        allocate (dsgn(ilo-2:ihi+2))
        allocate (dlim(ilo-2:ihi+2))
        allocate (  df(ilo-2:ihi+2))
        allocate (dcen(ilo-2:ihi+2))

        if (idir .eq. 1) then

           ! Slopes in first coordinate direction
           do j = ilo2-1, ihi2+1
              ! First compute Fromm slopes
              do i = ilo1-2, ihi1+2
                 dlft = p(i  ,j) - p(i-1,j)
                 drgt = p(i+1,j) - p(i  ,j)

                 ! Subtract off (rho * grav) so as not to limit that part of the slope
                 dlft = dlft - 0.25d0 * (rho(i,j)+rho(i-1,j))*(grav(i,j,1)+grav(i-1,j,1))*dx
                 drgt = drgt - 0.25d0 * (rho(i,j)+rho(i+1,j))*(grav(i,j,1)+grav(i+1,j,1))*dx

                 dcen(i) = 0.5d0*(dlft+drgt)
                 dsgn(i) = sign(1.d0, dcen(i))

                 if (dlft*drgt .ge. 0.d0) then
                    dlim(i) = 2.d0 * min( abs(dlft), abs(drgt) )
                 else
                    dlim(i) = 0.d0
                 endif
                 df(i) = dsgn(i)*min( dlim(i), abs(dcen(i)) )
              enddo

              ! now limited fourth order slopes
              do i = ilo1-1, ihi1+1
                 dp1 = four3rd*dcen(i) - sixth*(df(i+1) + df(i-1))
                 dp(i,j) = flatn(i,j)*dsgn(i)*min(dlim(i),abs(dp1))
                 dp(i,j) = dp(i,j) + rho(i,j)*grav(i,j,1)*dx
              enddo
           enddo

        else

           ! Slopes in second coordinate direction
           do i = ilo1-1, ihi1+1
              ! First compute Fromm slopes
              do j = ilo2-2, ihi2+2
                 dlft = p(i,j  ) - p(i,j-1)
                 drgt = p(i,j+1) - p(i,j  )

                 ! Subtract off (rho * grav) so as not to limit that part of the slope
                 dlft = dlft - 0.25d0 * (rho(i,j)+rho(i,j-1))*(grav(i,j,2)+grav(i,j-1,2))*dy
                 drgt = drgt - 0.25d0 * (rho(i,j)+rho(i,j+1))*(grav(i,j,2)+grav(i,j+1,2))*dy

                 dcen(j) = 0.5d0*(dlft+drgt)
                 dsgn(j) = sign( 1.d0, dcen(j) )

                 if (dlft*drgt .ge. 0.d0) then
                    dlim(j) = 2.d0 * min( abs(dlft), abs(drgt) )
                 else
                    dlim(j) = 0.d0
                 endif
                 df(j) = dsgn(j)*min( dlim(j),abs(dcen(j)) )
              enddo

              ! now limited fourth order slopes
              do j = ilo2-1, ihi2+1
                 dp1 = four3rd*dcen(j) - sixth*(df(j+1) + df(j-1))
                 dp(i,j) = flatn(i,j)*dsgn(j)*min(dlim(j),abs(dp1))
                 dp(i,j) = dp(i,j) + rho(i,j)*grav(i,j,2)*dy
              enddo
           enddo

        endif

        deallocate(dsgn,dlim,df,dcen)

      end subroutine pslope

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine uflaten(lo,hi,p,u,v,flatn, &
                         q_l1,q_l2,q_h1,q_h2)

      use meth_params_module, only : iorder, small_pres

      implicit none

      integer lo(2),hi(2)
      integer q_l1,q_l2,q_h1,q_h2
      double precision p(q_l1:q_h1,q_l2:q_h2)
      double precision u(q_l1:q_h1,q_l2:q_h2)
      double precision v(q_l1:q_h1,q_l2:q_h2)
      double precision flatn(q_l1:q_h1,q_l2:q_h2)

!     Local arrays
      double precision, allocatable :: dp(:), z(:), chi(:)

      integer i, j, idx, ishft
      double precision shktst, zcut1, zcut2, dzcut
      double precision denom, zeta, tst, tmp, ftmp
      integer nx,ny,nmax

!     Knobs for detection of strong shock
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
      
!     x-direction flattening coef
      do j = lo(2),hi(2) 
         do i = lo(1)-1,hi(1)+1
            idx = i-lo(1)+1
            dp(idx) = p(i+1,j) - p(i-1,j)
            denom = max(small_pres,abs(p(i+2,j)-p(i-2,j)))
            zeta = abs(dp(idx))/denom
            z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
            if (u(i-1,j)-u(i+1,j) .ge. 0.d0) then
               tst = 1.d0
            else
               tst = 0.d0
            endif
            tmp = min(p(i+1,j),p(i-1,j))
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
            flatn(i,j) = 1.d0 - &
                 max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
         enddo
      enddo

!     y-direction flattening coef
      do i = lo(1),hi(1)
         do j = lo(2)-1,hi(2)+1
            idx = j-lo(2)+1
            dp(idx) = p(i,j+1) - p(i,j-1)
            denom = max(small_pres,abs(p(i,j+2)-p(i,j-2)))
            zeta = abs(dp(idx))/denom
            z(idx) = min( 1.d0, max( 0.d0, dzcut*(zeta - zcut1) ) )
            if (v(i,j-1)-v(i,j+1) .ge. 0.d0) then
               tst = 1.d0
            else
               tst = 0.d0
            endif
            tmp = min(p(i,j+1),p(i,j-1))
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
            ftmp = 1.d0 - &
                 max(chi(idx-ishft)*z(idx-ishft),chi(idx)*z(idx))
            flatn(i,j) = min( flatn(i,j), ftmp )
         enddo
      enddo

      deallocate(dp,z,chi)

      end subroutine uflaten

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine divu(lo,hi,q,q_l1,q_l2,q_h1,q_h2,dx, &
                      div,div_l1,div_l2,div_h1,div_h2)

      use prob_params_module, only : coord_type
      use meth_params_module, only : QU, QV

      implicit none

      integer          :: lo(2),hi(2)
      integer          :: q_l1,q_l2,q_h1,q_h2
      integer          :: div_l1,div_l2,div_h1,div_h2
      double precision :: q(q_l1:q_h1,q_l2:q_h2,*)
      double precision :: div(div_l1:div_h1,div_l2:div_h2)
      double precision :: dx(2)

      integer          :: i, j
      double precision :: rl, rr, rc, ul, ur
      double precision :: vb, vt
      double precision :: ux,vy

      if (coord_type .eq. 0) then
         do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            ux = 0.5d0*(q(i,j,QU)-q(i-1,j,QU)+q(i,j-1,QU)-q(i-1,j-1,QU))/dx(1)
            vy = 0.5d0*(q(i,j,QV)-q(i,j-1,QV)+q(i-1,j,QV)-q(i-1,j-1,QV))/dx(2)
            div(i,j) = ux + vy
         enddo
         enddo
      else
         do i=lo(1),hi(1)+1

            if (i.eq.0) then

               div(i,lo(2):hi(2)+1) = 0.d0

            else 

               rl = (dble(i)-0.5d0) * dx(1)
               rr = (dble(i)+0.5d0) * dx(1)
               rc = (dble(i)      ) * dx(1)
   
               do j=lo(2),hi(2)+1
                  ! These are transverse averages in the y-direction
                  ul = 0.5d0 * (q(i-1,j,QU)+q(i-1,j-1,QU))
                  ur = 0.5d0 * (q(i  ,j,QU)+q(i  ,j-1,QU))
   
                  ! Take 1/r d/dr(r*u)
                  div(i,j) = (rr*ur - rl*ul) / dx(1) / rc
  
                  ! These are transverse averages in the x-direction
                  vb = 0.5d0 * (q(i,j-1,QV)+q(i-1,j-1,QV))
                  vt = 0.5d0 * (q(i,j  ,QV)+q(i-1,j  ,QV))
 
                  div(i,j) = div(i,j) + (vt - vb) / dx(2)
               enddo
   
            end if
         enddo
      end if

      end subroutine divu

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_corrgsrc(lo,hi, &
                             gold,gold_l1,gold_l2,gold_h1,gold_h2, &
                             gnew,gnew_l1,gnew_l2,gnew_h1,gnew_h2, &
                             uold,uold_l1,uold_l2,uold_h1,uold_h2, &
                             unew,unew_l1,unew_l2,unew_h1,unew_h2,dt)

      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN

      implicit none

      integer lo(2),hi(2)
      integer gold_l1,gold_l2,gold_h1,gold_h2
      integer gnew_l1,gnew_l2,gnew_h1,gnew_h2
      integer uold_l1,uold_l2,uold_h1,uold_h2
      integer unew_l1,unew_l2,unew_h1,unew_h2
      double precision   gold(gold_l1:gold_h1,gold_l2:gold_h2,2)
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,2)
      double precision  uold(uold_l1:uold_h1,uold_l2:uold_h2,NVAR)
      double precision  unew(unew_l1:unew_h1,unew_l2:unew_h2,NVAR)
      double precision dt

      integer i,j
      double precision SrU_old, SrV_old
      double precision SrU_new, SrV_new
      double precision SrUcorr, SrVcorr
      double precision SrEcorr

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

!           Define old source terms
            SrU_old = uold(i,j,URHO) * gold(i,j,1)
            SrV_old = uold(i,j,URHO) * gold(i,j,2)
            
!           Define new source terms
            SrU_new = unew(i,j,URHO) * gnew(i,j,1)
            SrV_new = unew(i,j,URHO) * gnew(i,j,2)
            
!           Define corrections to source terms
            SrUcorr = 0.5d0*(SrU_new - SrU_old)
            SrVcorr = 0.5d0*(SrV_new - SrV_old)

            ! This doesn't work (in 1-d)
            ! Upn = unew(i,j,UMX) / unew(i,j,URHO)
            ! Vpn = unew(i,j,UMY) / unew(i,j,URHO)
            ! SrEcorr = SrUcorr*(Upn + SrUcorr*dt/(2*unew(i,j,URHO))) &
            !          +SrVcorr*(Vpn + SrVcorr*dt/(2*unew(i,j,URHO)))

            ! This does work (in 1-d)
            SrEcorr = 0.5d0 * ( (unew(i,j,UMX) * gnew(i,j,1) - uold(i,j,UMX) * gold(i,j,1)) + &
                                (unew(i,j,UMY) * gnew(i,j,2) - uold(i,j,UMY) * gold(i,j,2)) )

!           Correct state with correction terms
            unew(i,j,UMX)   = unew(i,j,UMX)   + SrUcorr*dt
            unew(i,j,UMY)   = unew(i,j,UMY)   + SrVcorr*dt
            unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr*dt

         enddo
      enddo

      end subroutine ca_corrgsrc

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_syncgsrc(lo,hi, &
                             gphi,gphi_l1,gphi_l2,gphi_h1,gphi_h2, &
                             gdphi,gdphi_l1,gdphi_l2,gdphi_h1,gdphi_h2, &
                             state,state_l1,state_l2,state_h1,state_h2, &
                             dstate,dstate_l1,dstate_l2,dstate_h1,dstate_h2, &
                             sync_src,src_l1,src_l2,src_h1,src_h2, &
                             dt)

      use meth_params_module, only : NVAR, URHO, UMX, UMY

      implicit none

      integer lo(2),hi(2)
      integer gphi_l1,gphi_l2,gphi_h1,gphi_h2
      integer gdphi_l1,gdphi_l2,gdphi_h1,gdphi_h2
      integer state_l1,state_l2,state_h1,state_h2
      integer dstate_l1,dstate_l2,dstate_h1,dstate_h2
      integer src_l1,src_l2,src_h1,src_h2
      double precision     gphi(gphi_l1:gphi_h1,gphi_l2:gphi_h2,2)
      double precision    gdphi(gdphi_l1:gdphi_h1,gdphi_l2:gdphi_h2,2)
      double precision    state(state_l1:state_h1,state_l2:state_h2,NVAR)
      double precision   dstate(dstate_l1:dstate_h1,dstate_l2:dstate_h2,2+1)
      double precision sync_src(src_l1:src_h1,src_l2:src_h2,2+1)
      double precision dt

!     Note that dstate is drho and drhoU, state is the entire state, and src
!     is S_rhoU and S_rhoE

      integer i,j
      double precision rho_pre, rhoU_pre, rhoV_pre
      double precision gx, gy,dgx, dgy, SrU, SrV, SrE

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            
            rho_pre  = state(i,j,URHO) - dstate(i,j,1)
            rhoU_pre = state(i,j,UMX)  - dstate(i,j,2)
            rhoV_pre = state(i,j,UMY)  - dstate(i,j,3)
            
            gx  = gphi(i,j,1)
            gy  = gphi(i,j,2)

            dgx = gdphi(i,j,1)
            dgy = gdphi(i,j,2)

            SrU = dstate(i,j,1)*gx + rho_pre*dgx
            SrV = dstate(i,j,1)*gy + rho_pre*dgy

            SrE = SrU * (rhoU_pre + 0.5*SrU*dt)/rho_pre &
                 +SrV * (rhoV_pre + 0.5*SrV*dt)/rho_pre
            
            sync_src(i,j,1) = SrU
            sync_src(i,j,2) = SrV
            sync_src(i,j,3) = SrE

         enddo
      enddo
      end subroutine ca_syncgsrc
