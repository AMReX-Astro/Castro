module rad_advection_module

  implicit none

  private

  public umeth2d_rad, ctoprim_rad, consup_rad

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
subroutine ctoprim_rad(lo,hi, &
     uin,uin_l1,uin_l2,uin_h1,uin_h2, &
     Erin,Erin_l1,Erin_l2,Erin_h1,Erin_h2, &
     lam,lam_l1,lam_l2,lam_h1,lam_h2, &
     q,c,cg,gamc,gamcg,csml,flatn,q_l1,q_l2,q_h1,q_h2, &
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
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, comoving, &
       flatten_pp_threshold, first_order_hydro
  use rad_params_module, only : ngroups
  use flatten_module, only : uflaten
  use fluxlimiter_module, only : Edd_factor

  implicit none

  double precision, parameter:: small = 1.d-8

  integer lo(2), hi(2)
  integer uin_l1,uin_l2,uin_h1,uin_h2
  integer Erin_l1,Erin_l2,Erin_h1,Erin_h2
  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer q_l1,q_l2,q_h1,q_h2
  integer src_l1,src_l2,src_h1,src_h2
  integer iflaten

  double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  double precision :: Erin(Erin_l1:Erin_h1,Erin_l2:Erin_h2,0:ngroups-1)
  double precision :: lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision :: q(q_l1:q_h1,q_l2:q_h2,QRADVAR)
  double precision ::    c (q_l1:q_h1,q_l2:q_h2)
  double precision ::    cg(q_l1:q_h1,q_l2:q_h2)
  double precision :: gamc (q_l1:q_h1,q_l2:q_h2)
  double precision :: gamcg(q_l1:q_h1,q_l2:q_h2)
  double precision ::  csml(q_l1:q_h1,q_l2:q_h2)
  double precision :: flatn(q_l1:q_h1,q_l2:q_h2)
  double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
  double precision :: srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
  double precision :: dx, dy, dt, courno

  double precision, allocatable :: dpdrho(:,:)
  double precision, allocatable :: dpde(:,:)
  double precision, allocatable :: flatg(:,:)

  integer          :: i, j, g
  integer          :: ngp, ngf, loq(2), hiq(2)
  integer          :: iadv, ispec, iaux, n, nq
  double precision :: courx, coury, courmx, courmy

  double precision :: csrad2, prad, Eddf, gamr

  type(eos_t) :: eos_state

  allocate(dpdrho(q_l1:q_h1,q_l2:q_h2))
  allocate(dpde  (q_l1:q_h1,q_l2:q_h2))
  allocate(flatg (q_l1:q_h1,q_l2:q_h2))

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
        q(i,j,qrad:qradhi) = Erin(i,j,:)
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

        eos_state % rho = q(i,j,QRHO)
        eos_state % T   = q(i,j,QTEMP)
        eos_state % e   = q(i,j,QREINT)
        eos_state % xn  = q(i,j,QFS:QFS+nspec-1)
        eos_state % aux = q(i,j,QFX:QFX+naux-1)

        ! If necessary, reset the energy using small_temp
        if ((allow_negative_energy .eq. 0) .and. (q(i,j,QREINT) .lt. 0)) then
           q(i,j,QTEMP) = small_temp
           eos_state % T = q(i,j,QTEMP)
           call eos(eos_input_rt, eos_state)
           q(i,j,QPRES ) = eos_state % p
           q(i,j,QREINT) = eos_state % e
           if (q(i,j,QREINT) .lt. 0.d0) then
              print *,'   '
              print *,'>>> Error: ctoprim ',i,j
              print *,'>>> ... new e from eos call is negative ',q(i,j,QREINT)
              print *,'    '
              call bl_error("Error:: ctoprim")
           end if
        end if

        call eos(eos_input_re, eos_state)
        
        q(i,j,QTEMP) = eos_state % T
        q(i,j,QPRES) = eos_state % p
        dpdrho(i,j)  = eos_state % dpdr_e
        dpde(i,j)    = eos_state % dpde
        gamcg(i,j)   = eos_state % gam1
        cg(i,j)      = eos_state % cs

        csrad2 = 0.d0
        prad = 0.d0
        do g=0, ngroups-1
           if (comoving) then
              Eddf = Edd_factor(lam(i,j,g))
              gamr = (3.d0-Eddf)/2.d0
           else
              gamr = lam(i,j,g) + 1.d0
           end if
           prad = prad + lam(i,j,g)*q(i,j,qrad+g)
           csrad2 = csrad2 + gamr * (lam(i,j,g)*q(i,j,qrad+g)) / q(i,j,QRHO)
        end do

        q(i,j,qptot) = q(i,j,QPRES) + prad
        c(i,j) = cg(i,j)**2 + csrad2
        gamc(i,j) = c(i,j) * q(i,j,QRHO) / q(i,j,qptot)
        c(i,j) = sqrt(c(i,j))
        csml(i,j) = max(small, small * c(i,j))

!     Make this "rho e" instead of "e"
        q(i,j,QREINT) = q(i,j,QREINT)*q(i,j,QRHO)
        q(i,j,qreitot) = q(i,j,QREINT) + sum(q(i,j,qrad:qradhi)) 
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
        srcQ(i,j,QPRES ) =   dpde(i,j) * &
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
  if (first_order_hydro) then
     flatn = 0.d0
  else if(iflaten.eq.1)then
     do n=1,2
        loq(n)=lo(n)-ngf
        hiq(n)=hi(n)+ngf
     enddo
     call uflaten(loq,hiq, &
          q(:,:,qptot), &
          q(:,:,QU), &
          q(:,:,QV), &
          flatn,q_l1,q_l2,q_h1,q_h2)
     call uflaten(loq,hiq, &
          q(:,:,qpres), &
          q(:,:,QU), &
          q(:,:,QV), &
          flatg,q_l1,q_l2,q_h1,q_h2)
     flatn = flatn * flatg

     if (flatten_pp_threshold > 0.d0) then
        call ppflaten(loq,hiq, &
             flatn, q, q_l1,q_l2,q_h1,q_h2)
     end if
  else
     flatn = 1.d0
  endif

  deallocate(dpdrho,dpde, flatg)

end subroutine ctoprim_rad


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

subroutine umeth2d_rad(q, c,cg, gamc,gamcg, csml, flatn, qd_l1, qd_l2, qd_h1, qd_h2,&
     lam, lam_l1, lam_l2, lam_h1, lam_h2, &
     srcQ, src_l1, src_l2, src_h1, src_h2, &
     grav, gv_l1, gv_l2, gv_h1, gv_h2, &
     ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
     flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
     flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
     rflux1, rfd1_l1, rfd1_l2, rfd1_h1, rfd1_h2, &
     rflux2, rfd2_l1, rfd2_l2, rfd2_h1, rfd2_h2, &
     pgdx,pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2, &
     pgdy,pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2, &
     ergdx, ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2, &
     ergdy, ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2, &
     lmgdx, lmgdx_l1, lmgdx_l2, lmgdx_h1, lmgdx_h2, &
     lmgdy, lmgdy_l1, lmgdy_l2, lmgdy_h1, lmgdy_h2, &
     ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
     ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
     area1, area1_l1, area1_l2, area1_h1, area1_h2, &
     area2, area2_l1, area2_l2, area2_h1, area2_h2, &
     pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
     uy_xfc, ux_yfc, &
     dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2)

  use network, only : nspec, naux
  use meth_params_module, only : QVAR, NVAR, ppm_type 
  use radhydro_params_module, only : QRADVAR
  use rad_params_module, only : ngroups

  implicit none

  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2
  integer ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2
  integer lmgdx_l1, lmgdx_l2, lmgdx_h1, lmgdx_h2
  integer lmgdy_l1, lmgdy_l2, lmgdy_h1, lmgdy_h2
  integer rfd1_l1, rfd1_l2, rfd1_h1, rfd1_h2
  integer rfd2_l1, rfd2_l2, rfd2_h1, rfd2_h2
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
  double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision  gamc(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision gamcg(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision  csml(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision    cg(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision  srcQ(src_l1:src_h1,src_l2:src_h2)
  double precision  grav( gv_l1: gv_h1, gv_l2: gv_h2)
  double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  double precision  pgdx(pgdx_l1:pgdx_h1,pgdx_l2:pgdx_h2)
  double precision  ugdx(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
  double precision  pgdy(pgdy_l1:pgdy_h1,pgdy_l2:pgdy_h2)
  double precision  ugdy(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)
  double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
  double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
  double precision pdivu(ilo1:ihi1,ilo2:ihi2)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)

  double precision uy_xfc(ugdx_l1:ugdx_h1,ugdx_l2:ugdx_h2)
  double precision ux_yfc(ugdy_l1:ugdy_h1,ugdy_l2:ugdy_h2)

  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision ergdx(ergdx_l1:ergdx_h1,ergdx_l2:ergdx_h2,0:ngroups-1)
  double precision ergdy(ergdy_l1:ergdy_h1,ergdy_l2:ergdy_h2,0:ngroups-1)
  double precision lmgdx(lmgdx_l1:lmgdx_h1,lmgdx_l2:lmgdx_h2,0:ngroups-1)
  double precision lmgdy(lmgdy_l1:lmgdy_h1,lmgdy_l2:lmgdy_h2,0:ngroups-1)
  double precision rflux1(rfd1_l1:rfd1_h1,rfd1_l2:rfd1_h2,0:ngroups-1)
  double precision rflux2(rfd2_l1:rfd2_h1,rfd2_l2:rfd2_h2,0:ngroups-1)

!     Left and right state arrays (edge centered, cell centered)
  double precision, allocatable:: dq(:,:,:),  qm(:,:,:),   qp(:,:,:)
  double precision, allocatable::qxm(:,:,:),qym(:,:,:)
  double precision, allocatable::qxp(:,:,:),qyp(:,:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
  double precision, allocatable::   fx(:,:,:),  fy(:,:,:)
  double precision, allocatable::  rfx(:,:,:), rfy(:,:,:)
  double precision, allocatable::   pgdxtmp(:,:) ,  ugdxtmp(:,:),  ergdxtmp(:,:,:)

!     Local scalar variables
  double precision :: dtdx
  double precision :: hdtdx, hdt, hdtdy
  integer          :: i,j
  double precision :: pggdx, pggdy

  allocate ( pgdxtmp( pgdx_l1: pgdx_h1, pgdx_l2: pgdx_h2))
  allocate (ergdxtmp(ergdx_l1:ergdx_h1,ergdx_l2:ergdx_h2,0:ngroups-1))
  allocate ( ugdxtmp( ugdx_l1: ugdx_h1, ugdx_l2: ugdx_h2))
  allocate (  dq(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate (  qm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate (  qp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate ( qxm(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate ( qxp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate ( qym(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate ( qyp(ilo1-1:ihi1+2,ilo2-1:ihi2+2,QRADVAR) )
  allocate ( fx(ilo1  :ihi1+1,ilo2-1:ihi2+1,NVAR))
  allocate ( fy(ilo1-1:ihi1+1,ilo2  :ihi2+1,NVAR))
  allocate (rfx(ilo1  :ihi1+1,ilo2-1:ihi2+1,0:ngroups-1))
  allocate (rfy(ilo1-1:ihi1+1,ilo2  :ihi2+1,0:ngroups-1))

!     Local constants
  dtdx = dt/dx
  hdtdx = 0.5d0*dtdx
  hdtdy = 0.5d0*dt/dy
  hdt = 0.5d0*dt

!     NOTE: Geometry terms need to be punched through

! !     Trace to edges w/o transverse flux correction terms
  if (ppm_type .eq. 0) then
     call bl_error("ppm_type <=0 is not supported in umeth2d_rad")
  else
     call trace_ppm_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
          q,c,cg,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
          dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
          qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
          ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
  end if

  call cmpflx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
       qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
        fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
       rfx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
        pgdxtmp,  pgdx_l1,  pgdx_l2,  pgdx_h1,  pgdx_h2, &
       ergdxtmp, ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2, &
          lmgdx, lmgdx_l1, lmgdx_l2, lmgdx_h1, lmgdx_h2, &
        ugdxtmp,  ugdx_l1,  ugdx_l2,  ugdx_h1,  ugdx_h2, &
        uy_xfc, &
       gamc,gamcg, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
       1, ilo1, ihi1, ilo2-1, ihi2+1)

  call cmpflx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
       qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
        fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
       rfy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
        pgdy,  pgdy_l1,  pgdy_l2,  pgdy_h1,  pgdy_h2, &
       ergdy, ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2, &
       lmgdy, lmgdy_l1, lmgdy_l2, lmgdy_h1, lmgdy_h2, &
        ugdy,  ugdy_l1,  ugdy_l2,  ugdy_h1,  ugdy_h2, &
        ux_yfc, &
       gamc,gamcg, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
       2, ilo1-1, ihi1+1, ilo2, ihi2)

  call transy_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
       qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
        fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
       rfy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
        pgdy,  pgdy_l1 , pgdy_l2,  pgdy_h1,  pgdy_h2, &
       ergdy, ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2, &
        ugdy,  ugdy_l1,  ugdy_l2,  ugdy_h1,  ugdy_h2, &
       gamcg, qd_l1, qd_l2, qd_h1, qd_h2, &
       srcQ, src_l1, src_l2, src_h1, src_h2, &
       grav, gv_l1, gv_l2, gv_h1, gv_h2, &
       hdt, hdtdy, &
       ilo1-1, ihi1+1, ilo2, ihi2)

  call cmpflx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
       qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
        flux1,  fd1_l1,  fd1_l2,  fd1_h1,  fd1_h2, &
       rflux1, rfd1_l1, rfd1_l2, rfd1_h1, rfd1_h2, &
        pgdx,  pgdx_l1,  pgdx_l2,  pgdx_h1,  pgdx_h2, &
       ergdx, ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2, &
       lmgdx, lmgdx_l1, lmgdx_l2, lmgdx_h1, lmgdx_h2, &
        ugdx,  ugdx_l1,  ugdx_l2,  ugdx_h1,  ugdx_h2, &
        uy_xfc, &
       gamc,gamcg, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
       1, ilo1, ihi1, ilo2, ihi2)

  call transx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, & 
       qym, qm,qyp,qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
        fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
       rfx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
        pgdxtmp,  pgdx_l1,  pgdx_l2,  pgdx_h1,  pgdx_h2, &
       ergdxtmp, ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2, &
        ugdxtmp,  ugdx_l1,  ugdx_l2,  ugdx_h1,  ugdx_h2, &
       gamcg, qd_l1, qd_l2, qd_h1, qd_h2, &
       srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
       grav, gv_l1, gv_l2, gv_h1, gv_h2, &
       hdt, hdtdx, &
       area1, area1_l1, area1_l2, area1_h1, area1_h2, &
       vol, vol_l1, vol_l2, vol_h1, vol_h2, &
       ilo1, ihi1, ilo2-1, ihi2+1)
  
  call cmpflx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
       qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
        flux2,  fd2_l1,  fd2_l2,  fd2_h1,  fd2_h2, &
       rflux2, rfd2_l1, rfd2_l2, rfd2_h1, rfd2_h2, &
        pgdy,  pgdy_l1,  pgdy_l2,  pgdy_h1,  pgdy_h2, &
       ergdy, ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2, &
       lmgdy, lmgdy_l1, lmgdy_l2, lmgdy_h1, lmgdy_h2, &
        ugdy,  ugdy_l1,  ugdy_l2,  ugdy_h1,  ugdy_h2, &
        ux_yfc, &
       gamc,gamcg, csml, c, qd_l1, qd_l2, qd_h1, qd_h2, &
       2, ilo1, ihi1, ilo2, ihi2)

  do j = ilo2,ihi2
     do i = ilo1,ihi1
        pggdx = 0.5d0*(pgdx(i+1,j  )+pgdx(i,j)) 
        pggdy = 0.5d0*(pgdy(i  ,j+1)+pgdy(i,j)) 
        pdivu(i,j) = (pggdx*(ugdx(i+1,j  )*area1(i+1,j  )-ugdx(i,j)*area1(i,j)) &
             +        pggdy*(ugdy(i  ,j+1)*area2(i  ,j+1)-ugdy(i,j)*area2(i,j)) ) / vol(i,j)
     end do
  end do

  deallocate(dq,qm,qp,qxm,qxp,qym,qyp)
  deallocate(fx,fy)
  deallocate(rfx,rfy)
  deallocate(pgdxtmp,ugdxtmp,ergdxtmp)

end subroutine umeth2d_rad


! ::: 
! ::: ------------------------------------------------------------------
! ::: 
subroutine trace_ppm_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &   
     q,c,cg,flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
     dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
     qxm,qxp,qym,qyp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
     ilo1,ilo2,ihi1,ihi2,dx,dy,dt)

  use network, only : nspec, naux
  use meth_params_module, only : QVAR, QRHO, QU, QV, &
       QREINT, QPRES, QFA, QFS, QFX, &
       nadv, small_dens, ppm_type
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
  use rad_params_module, only : ngroups
  use ppm_module, only : ppm

  implicit none

  integer ilo1,ilo2,ihi1,ihi2
  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer qd_l1,qd_l2,qd_h1,qd_h2
  integer dloga_l1,dloga_l2,dloga_h1,dloga_h2
  integer qpd_l1,qpd_l2,qpd_h1,qpd_h2

  double precision dx, dy, dt
  double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision     c(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision    cg(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  double precision qxm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision qxp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision qym(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision qyp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)

  ! Local variables
  integer i, j, g
  integer n, iadv
  integer ns, ispec, iaux

  double precision dtdx, dtdy
  
  double precision, dimension(0:ngroups-1) :: er, der, alphar, sourcer, qrtmp,hr
  double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
  double precision cc, csq, rho, u, v, p, ptot, rhoe, enth, cgassq
  double precision dum, dvm, dptotm
  double precision drho, du, dv, drhoe, dptot
  double precision dup, dvp, dptotp
  
  double precision alpham, alphap, alpha0, alphae, alphau, alphav
  double precision sourcr,sourcp,source,courn,eta,dlogatmp

  double precision rhoe_g, h_g, alphae_g, drhoe_g

  double precision, allocatable :: Ip(:,:,:,:,:)
  double precision, allocatable :: Im(:,:,:,:,:)

  double precision :: er_foo

  dtdx = dt/dx
  dtdy = dt/dy

  if (ppm_type .eq. 0) then
     print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
     call bl_error("Error:: ppm_2d.f90 :: trace_ppm")
  end if
  
  allocate(Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QRADVAR))
  allocate(Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,2,3,QRADVAR))

  ! Compute Ip and Im
  do n=1,QRADVAR
     call ppm(q(:,:,n),qd_l1,qd_l2,qd_h1,qd_h2, &
          q(:,:,QU:), c, qd_l1,qd_l2,qd_h1,qd_h2,&
          flatn, &
          Ip(:,:,:,:,n),Im(:,:,:,:,n), &
          ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
  end do

  ! Trace to left and right edges using upwind PPM
  do j = ilo2-1, ihi2+1
     do i = ilo1-1, ihi1+1

        do g=0, ngroups-1
           lam0(g) = lam(i,j,g)
           lamp(g) = lam(i,j,g)
           lamm(g) = lam(i,j,g)
        end do

        cgassq = cg(i,j)**2
        cc = c(i,j)
        csq = cc**2

        rho = q(i,j,QRHO)
        u = q(i,j,QU)
        v = q(i,j,QV)
        p = q(i,j,QPRES)
        rhoe_g = q(i,j,QREINT)
        h_g = (p+rhoe_g) / rho
        er(:) = q(i,j,qrad:qradhi)
        hr(:) = (lam0+1.d0)*er/rho
        ptot = q(i,j,qptot)
        rhoe = q(i,j,qreitot)
        enth = ( (rhoe+ptot)/rho )/csq

        ! plus state on face i
        dum    = flatn(i,j)*(u    - Im(i,j,1,1,QU))
        dptotm = flatn(i,j)*(ptot - Im(i,j,1,1,qptot))

        drho    = flatn(i,j)*(rho    - Im(i,j,1,2,QRHO))
        dv      = flatn(i,j)*(v      - Im(i,j,1,2,QV))
        dptot   = flatn(i,j)*(ptot   - Im(i,j,1,2,qptot))
        drhoe   = flatn(i,j)*(rhoe   - Im(i,j,1,2,qreitot))
        drhoe_g = flatn(i,j)*(rhoe_g - Im(i,j,1,2,QREINT))
        der(:)  = flatn(i,j)*(er(:)  - Im(i,j,1,2,qrad:qradhi))

        dup    = flatn(i,j)*(u    - Im(i,j,1,3,QU))
        dptotp = flatn(i,j)*(ptot - Im(i,j,1,3,qptot))

        alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
        alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
        alpha0 = drho - dptot/csq
        alphae = drhoe - dptot*enth
        alphae_g = drhoe_g - dptot/csq*h_g
        alphav = dv
        alphar(:) = der(:) - dptot/csq*hr

        if (u-cc .gt. 0.d0) then
           alpham = 0.d0
        else if (u-cc .lt. 0.d0) then
           alpham = -alpham
        else
           alpham = -0.5d0*alpham
        endif
        if (u+cc .gt. 0.d0) then
           alphap = 0.d0
        else if (u+cc .lt. 0.d0) then
           alphap = -alphap
        else
           alphap = -0.5d0*alphap
        endif
        if (u .gt. 0.d0) then
           alpha0 = 0.d0
           alphav = 0.d0
           alphae = 0.d0 
           alphae_g = 0.d0 
           alphar(:) = 0.d0
        else if (u .lt. 0.d0) then
           alpha0 = -alpha0
           alphav = -alphav
           alphae = -alphae
           alphae_g = -alphae_g
           alphar(:) = -alphar(:)
        else
           alpha0 = -0.5d0*alpha0
           alphav = -0.5d0*alphav
           alphae = -0.5d0*alphae
           alphae_g = -0.5d0*alphae_g
           alphar(:) = -0.5d0*alphar(:)
        endif

        if (i .ge. ilo1) then
           qxp(i,j,QRHO) = rho + alphap + alpham + alpha0
           qxp(i,j,QRHO) = max(small_dens,qxp(i,j,QRHO))
           qxp(i,j,QU) = u + (alphap - alpham)*cc/rho
           qxp(i,j,QV) = v + alphav 
           qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
           qxp(i,j,qrad:qradhi) = qrtmp
           qxp(i,j,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
           qxp(i,j,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
           qxp(i,j,qptot) = ptot + (alphap + alpham)*csq 
           qxp(i,j,qreitot) = qxp(i,j,QREINT) + sum(qrtmp)

           do g=0, ngroups-1
              if (qxp(i,j,qrad+g) < 0.d0) then
                 er_foo = - qxp(i,j,qrad+g)
                 qxp(i,j,qrad+g) = 0.d0
                 qxp(i,j,qptot) = qxp(i,j,qptot) + lamp(g) * er_foo
                 qxp(i,j,qreitot) = qxp(i,j,qreitot) + er_foo
              end if
           end do

           if (qxp(i,j,QPRES) < 0.d0) then
              qxp(i,j,QPRES) = p
           end if
        end if

        ! minus state on face i+1
        dum    = flatn(i,j)*(u    - Ip(i,j,1,1,QU))
        dptotm = flatn(i,j)*(ptot - Ip(i,j,1,1,qptot))

        drho    = flatn(i,j)*(rho    - Ip(i,j,1,2,QRHO))
        dv      = flatn(i,j)*(v      - Ip(i,j,1,2,QV))
        dptot   = flatn(i,j)*(ptot   - Ip(i,j,1,2,qptot))
        drhoe   = flatn(i,j)*(rhoe   - Ip(i,j,1,2,qreitot))
        drhoe_g = flatn(i,j)*(rhoe_g - Ip(i,j,1,2,QREINT))
        der(:)  = flatn(i,j)*(er(:)  - Ip(i,j,1,2,qrad:qradhi))

        dup    = flatn(i,j)*(u    - Ip(i,j,1,3,QU))
        dptotp = flatn(i,j)*(ptot - Ip(i,j,1,3,qptot))

        alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
        alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
        alpha0 = drho - dptot/csq
        alphae = drhoe - dptot*enth
        alphae_g = drhoe_g - dptot/csq*h_g
        alphav = dv
        alphar(:) = der(:)- dptot/csq*hr

        if (u-cc .gt. 0.d0) then
           alpham = -alpham
        else if (u-cc .lt. 0.d0) then
           alpham = 0.d0
        else
           alpham = -0.5d0*alpham
        endif
        if (u+cc .gt. 0.d0) then
           alphap = -alphap
        else if (u+cc .lt. 0.d0) then
           alphap = 0.d0
        else
           alphap = -0.5d0*alphap
        endif
        if (u .gt. 0.d0) then
           alpha0 = -alpha0
           alphav = -alphav
           alphae = -alphae
           alphae_g = -alphae_g
           alphar(:) = -alphar(:)
        else if (u .lt. 0.d0) then
           alpha0 = 0.d0
           alphav = 0.d0
           alphae = 0.d0
           alphae_g = 0.d0
           alphar(:) = 0.d0
        else
           alpha0 = -0.5d0*alpha0
           alphav = -0.5d0*alphav
           alphae = -0.5d0*alphae
           alphae_g = -0.5d0*alphae_g
           alphar(:) = -0.5d0*alphar(:)
        endif

        if (i .le. ihi1) then
           qxm(i+1,j,QRHO) = rho + alphap + alpham + alpha0
           qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
           qxm(i+1,j,QU) = u + (alphap - alpham)*cc/rho
           qxm(i+1,j,QV) = v + alphav
           qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
           qxm(i+1,j,qrad:qradhi) = qrtmp
           qxm(i+1,j,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
           qxm(i+1,j,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
           qxm(i+1,j,qptot) = ptot + (alphap + alpham)*csq
           qxm(i+1,j,qreitot) = qxm(i+1,j,QREINT) + sum(qrtmp)

           do g=0, ngroups-1
              if (qxm(i+1,j,qrad+g) < 0.d0) then
                 er_foo = - qxm(i+1,j,qrad+g)
                 qxm(i+1,j,qrad+g) = 0.d0
                 qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + lamm(g) * er_foo
                 qxm(i+1,j,qreitot) = qxm(i+1,j,qreitot) + er_foo
              end if
           end do

           if (qxm(i+1,j,QPRES) < 0.d0) then
              qxm(i+1,j,QPRES) = p
           end if
        end if

        if(dloga(i,j).ne.0)then
           courn = dtdx*(cc+abs(u))
           eta = (1.d0-courn)/(cc*dt*abs(dloga(i,j)))
           dlogatmp = min(eta,1.d0)*dloga(i,j)
           sourcr = -0.5d0*dt*rho*dlogatmp*u
           sourcp = sourcr*cgassq
           source = sourcr*h_g
           sourcer(:) = -0.5d0*dt*dlogatmp*u*(lam0(:)+1.d0)*er(:)
           if (i .le. ihi1) then
              qxm(i+1,j,QRHO) = qxm(i+1,j,QRHO) + sourcr
              qxm(i+1,j,QRHO) = max(qxm(i+1,j,QRHO),small_dens)
              qxm(i+1,j,QPRES) = qxm(i+1,j,QPRES) + sourcp
              qxm(i+1,j,QREINT) = qxm(i+1,j,QREINT) + source
              qxm(i+1,j,qrad:qradhi) = qxm(i+1,j,qrad:qradhi) + sourcer(:)
!              qxm(i+1,j,qptot ) = sum(lamm(:)*qxm(i+1,j,qrad:qradhi)) + qxm(i+1,j,QPRES)
              qxm(i+1,j,qptot) = qxm(i+1,j,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
              qxm(i+1,j,qreitot) = sum(qxm(i+1,j,qrad:qradhi))  + qxm(i+1,j,QREINT)
           end if
           if (i .ge. ilo1) then
              qxp(i,j,QRHO) = qxp(i,j,QRHO) + sourcr
              qxp(i,j,QRHO) = max(qxp(i,j,QRHO),small_dens)
              qxp(i,j,QPRES) = qxp(i,j,QPRES) + sourcp
              qxp(i,j,QREINT) = qxp(i,j,QREINT) + source
              qxp(i,j,qrad:qradhi) = qxp(i,j,qrad:qradhi) + sourcer(:)
!              qxp(i,j,qptot ) = sum(lamp(:)*qxp(i,j,qrad:qradhi)) + qxp(i,j,QPRES)
              qxp(i,j,qptot) = qxp(i,j,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
              qxp(i,j,qreitot) = sum(qxp(i,j,qrad:qradhi))  + qxp(i,j,QREINT)
           end if
        endif
              
     end do
  end do

  ! Now do the passively advected quantities
  do iadv = 1, nadv
     n = QFA + iadv - 1
     do j = ilo2-1, ihi2+1

        ! plus state on face i
        do i = ilo1, ihi1+1
           u = q(i,j,QU)
           if (u .gt. 0.d0) then
              qxp(i,j,n) = q(i,j,n)
           else if (u .lt. 0.d0) then
              qxp(i,j,n) = q(i,j,n) + flatn(i,j)*(Im(i,j,1,2,n) - q(i,j,n))
           else
              qxp(i,j,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Im(i,j,1,2,n) - q(i,j,n))
           endif
        enddo
        
        ! minus state on face i+1
        do i = ilo1-1, ihi1
           u = q(i,j,QU)
           if (u .gt. 0.d0) then
              qxm(i+1,j,n) = q(i,j,n) + flatn(i,j)*(Ip(i,j,1,2,n) - q(i,j,n))
           else if (u .lt. 0.d0) then
              qxm(i+1,j,n) = q(i,j,n)
           else
              qxm(i+1,j,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,n) - q(i,j,n))
           endif
        enddo
        
     enddo
  enddo
  
  do ispec = 1, nspec
     ns = QFS + ispec - 1
     do j = ilo2-1, ihi2+1
        
        ! plus state on face i
        do i = ilo1, ihi1+1
           u = q(i,j,QU)
           if (u .gt. 0.d0) then
              qxp(i,j,ns) = q(i,j,ns)
           else if (u .lt. 0.d0) then
              qxp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
           else
              qxp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
           endif
        enddo
        
        ! minus state on face i+1
        do i = ilo1-1, ihi1
           u = q(i,j,QU)
           if (u .gt. 0.d0) then
              qxm(i+1,j,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
           else if (u .lt. 0.d0) then
              qxm(i+1,j,ns) = q(i,j,ns)
           else
              qxm(i+1,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
           endif
        enddo
        
     enddo
  enddo
  
  do iaux = 1, naux
     ns = QFX + iaux - 1
     do j = ilo2-1, ihi2+1
        
        ! plus state on face i
        do i = ilo1, ihi1+1
           u = q(i,j,QU)
           if (u .gt. 0.d0) then
              qxp(i,j,ns) = q(i,j,ns)
           else if (u .lt. 0.d0) then
              qxp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
           else
              qxp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,1,2,ns) - q(i,j,ns))
           endif
        enddo
        
        ! minus state on face i+1
        do i = ilo1-1, ihi1
           u = q(i,j,QU)
           if (u .gt. 0.d0) then
              qxm(i+1,j,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
           else if (u .lt. 0.d0) then
              qxm(i+1,j,ns) = q(i,j,ns)
           else
              qxm(i+1,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,1,2,ns) - q(i,j,ns))
           endif
        enddo
        
     enddo
  enddo

  ! Trace to bottom and top edges using upwind PPM
  do j = ilo2-1, ihi2+1
     do i = ilo1-1, ihi1+1

        do g=0, ngroups-1
           lam0(g) = lam(i,j,g)
           lamp(g) = lam(i,j,g)
           lamm(g) = lam(i,j,g)
        end do

        cgassq = cg(i,j)**2
        cc = c(i,j)
        csq = cc**2

        rho = q(i,j,QRHO)
        u = q(i,j,QU)
        v = q(i,j,QV)
        p = q(i,j,QPRES)
        rhoe_g = q(i,j,QREINT)
        h_g = (p+rhoe_g) / rho
        er(:) = q(i,j,qrad:qradhi)
        hr(:) = (lam0+1.d0)*er/rho
        ptot = q(i,j,qptot)
        rhoe = q(i,j,qreitot)
        enth = ( (rhoe+ptot)/rho )/csq

        ! plus state on face j
        dvm    = flatn(i,j)*(v    - Im(i,j,2,1,QV))
        dptotm = flatn(i,j)*(ptot - Im(i,j,2,1,qptot))

        drho    = flatn(i,j)*(rho    - Im(i,j,2,2,QRHO))
        du      = flatn(i,j)*(u      - Im(i,j,2,2,QU))
        dptot   = flatn(i,j)*(ptot   - Im(i,j,2,2,qptot))
        drhoe   = flatn(i,j)*(rhoe   - Im(i,j,2,2,qreitot))
        drhoe_g = flatn(i,j)*(rhoe_g - Im(i,j,2,2,QREINT))
        der(:)  = flatn(i,j)*(er(:)  - Im(i,j,2,2,qrad:qradhi))

        dvp    = flatn(i,j)*(v    - Im(i,j,2,3,QV))
        dptotp = flatn(i,j)*(ptot - Im(i,j,2,3,qptot))

        alpham = 0.5d0*(dptotm/(rho*cc) - dvm)*rho/cc
        alphap = 0.5d0*(dptotp/(rho*cc) + dvp)*rho/cc
        alpha0 = drho - dptot/csq
        alphae = drhoe - dptot*enth
        alphae_g = drhoe_g - dptot/csq*h_g
        alphau = du
        alphar(:) = der(:) - dptot/csq*hr

        if (v-cc .gt. 0.d0) then
           alpham = 0.d0
        else if (v-cc .lt. 0.d0) then 
           alpham = -alpham
        else
           alpham = -0.5d0*alpham
        endif
        if (v+cc .gt. 0.d0) then
           alphap = 0.d0
        else if (v+cc .lt. 0.d0) then
           alphap = -alphap
        else
           alphap = -0.5d0*alphap
        endif
        if (v .gt. 0.d0) then
           alpha0 = 0.d0
           alphau = 0.d0
           alphae = 0.d0
           alphae_g = 0.d0 
           alphar(:) = 0.d0
        else if (v .lt. 0.d0) then
           alpha0 = -alpha0
           alphau = -alphau
           alphae = -alphae
           alphae_g = -alphae_g
           alphar(:) = -alphar(:)
        else
           alpha0 = -0.5d0*alpha0
           alphau = -0.5d0*alphau
           alphae = -0.5d0*alphae
           alphae_g = -0.5d0*alphae_g
           alphar(:) = -0.5d0*alphar(:)
        endif

        if (j .ge. ilo2) then
           qyp(i,j,QRHO) = rho + alphap + alpham + alpha0
           qyp(i,j,QRHO) = max(small_dens, qyp(i,j,QRHO))
           qyp(i,j,QV) = v + (alphap - alpham)*cc/rho
           qyp(i,j,QU) = u + alphau
           qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
           qyp(i,j,qrad:qradhi) = qrtmp
           qyp(i,j,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
           qyp(i,j,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
           qyp(i,j,qptot) = ptot + (alphap + alpham)*csq 
           qyp(i,j,qreitot) = qyp(i,j,QREINT) + sum(qrtmp)

           do g=0, ngroups-1
              if (qyp(i,j,qrad+g) < 0.d0) then
                 er_foo = - qyp(i,j,qrad+g)
                 qyp(i,j,qrad+g) = 0.d0
                 qyp(i,j,qptot) = qyp(i,j,qptot) + lamp(g) * er_foo
                 qyp(i,j,qreitot) = qyp(i,j,qreitot) + er_foo
              end if
           end do

           if (qyp(i,j,QPRES) < 0.d0) then
              qyp(i,j,QPRES) = p
           end if
        end if

        ! minus state on face j+1
        dvm    = flatn(i,j)*(v    - Ip(i,j,2,1,QV))
        dptotm = flatn(i,j)*(ptot - Ip(i,j,2,1,qptot))

        drho    = flatn(i,j)*(rho    - Ip(i,j,2,2,QRHO))
        du      = flatn(i,j)*(u      - Ip(i,j,2,2,QU))
        dptot   = flatn(i,j)*(ptot   - Ip(i,j,2,2,qptot))
        drhoe   = flatn(i,j)*(rhoe   - Ip(i,j,2,2,qreitot))
        drhoe_g = flatn(i,j)*(rhoe_g - Ip(i,j,2,2,QREINT))
        der(:)  = flatn(i,j)*(er(:)  - Ip(i,j,2,2,qrad:qradhi))

        dvp    = flatn(i,j)*(v    - Ip(i,j,2,3,QV))
        dptotp = flatn(i,j)*(ptot - Ip(i,j,2,3,qptot))

        alpham = 0.5d0*(dptotm/(rho*cc) - dvm)*rho/cc
        alphap = 0.5d0*(dptotp/(rho*cc) + dvp)*rho/cc
        alpha0 = drho - dptot/csq 
        alphae = drhoe - dptot*enth
        alphae_g = drhoe_g - dptot/csq*h_g
        alphau = du
        alphar(:) = der(:)- dptot/csq*hr

        if (v-cc .gt. 0.d0) then
           alpham = -alpham
        else if (v-cc .lt. 0.d0) then
           alpham = 0.d0
        else
           alpham = -0.5d0*alpham
        endif
        if (v+cc .gt. 0.d0) then
           alphap = -alphap
        else if (v+cc .lt. 0.d0) then
           alphap = 0.d0
        else
           alphap = -0.5d0*alphap
        endif
        if (v .gt. 0.d0) then
           alpha0 = -alpha0
           alphau = -alphau
           alphae = -alphae
           alphae_g = -alphae_g
           alphar(:) = -alphar(:)
        else if (v .lt. 0.d0) then
           alpha0 = 0.d0
           alphau = 0.d0
           alphae = 0.d0
           alphae_g = 0.d0
           alphar(:) = 0.d0
        else
           alpha0 = -0.5d0*alpha0
           alphau = -0.5d0*alphau
           alphae = -0.5d0*alphae
           alphae_g = -0.5d0*alphae_g
           alphar(:) = -0.5d0*alphar(:)
        endif
        
        if (j .le. ihi2) then
           qym(i,j+1,QRHO) = rho + alphap + alpham + alpha0
           qym(i,j+1,QRHO) = max(small_dens, qym(i,j+1,QRHO))
           qym(i,j+1,QV) = v + (alphap - alpham)*cc/rho
           qym(i,j+1,QU) = u + alphau
           qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
           qym(i,j+1,qrad:qradhi) = qrtmp
           qym(i,j+1,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
           qym(i,j+1,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
           qym(i,j+1,qptot) = ptot + (alphap + alpham)*csq
           qym(i,j+1,qreitot) = qym(i,j+1,QREINT) + sum(qrtmp)

           do g=0, ngroups-1
              if (qym(i,j+1,qrad+g) < 0.d0) then
                 er_foo = - qym(i,j+1,qrad+g)
                 qym(i,j+1,qrad+g) = 0.d0
                 qym(i,j+1,qptot) = qym(i,j+1,qptot) + lamm(g) * er_foo
                 qym(i,j+1,qreitot) = qym(i,j+1,qreitot) + er_foo
              end if
           end do

           if (qym(i,j+1,QPRES) < 0.d0) then
              qym(i,j+1,QPRES) = p
           end if
        end if
        
     end do
  end do

  ! Now do the passively advected quantities
  do iadv = 1, nadv
     n = QFA + iadv - 1
     do i = ilo1-1, ihi1+1

        ! plus state on face j
        do j = ilo2, ihi2+1
           v = q(i,j,QV)
           if (v .gt. 0.d0) then
              qyp(i,j,n) = q(i,j,n)
           else if (v .lt. 0.d0) then
              qyp(i,j,n) = q(i,j,n) + flatn(i,j)*(Im(i,j,2,2,n) - q(i,j,n))
           else
              qyp(i,j,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Im(i,j,2,2,n) - q(i,j,n))
           endif
        enddo
        
        ! minus state on face j+1
        do j = ilo2-1, ihi2
           v = q(i,j,QV)
           if (v .gt. 0.d0) then
              qym(i,j+1,n) = q(i,j,n) + flatn(i,j)*(Ip(i,j,2,2,n) - q(i,j,n))
           else if (v .lt. 0.d0) then
              qym(i,j+1,n) = q(i,j,n)
           else
              qym(i,j+1,n) = q(i,j,n) + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,n) - q(i,j,n))
           endif
        enddo
        
     enddo
  enddo
  
  do ispec = 1, nspec
     ns = QFS + ispec - 1
     do i = ilo1-1, ihi1+1
        
        ! plus state on face j
        do j = ilo2, ihi2+1
           v = q(i,j,QV)
           if (v .gt. 0.d0) then
              qyp(i,j,ns) = q(i,j,ns)
           else if (v .lt. 0.d0) then
              qyp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
           else
              qyp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
           endif
        enddo
        
        ! minus state on face j+1
        do j = ilo2-1, ihi2
           v = q(i,j,QV)
           if (v .gt. 0.d0) then
              qym(i,j+1,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
           else if (v .lt. 0.d0) then
              qym(i,j+1,ns) = q(i,j,ns)
           else
              qym(i,j+1,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
           endif
        enddo
        
     enddo
  enddo
  
  do iaux = 1, naux
     ns = QFX + iaux - 1
     do i = ilo1-1, ihi1+1
        
        ! plus state on face j
        do j = ilo2, ihi2+1
           v = q(i,j,QV)
           if (v .gt. 0.d0) then
              qyp(i,j,ns) = q(i,j,ns)
           else if (v .lt. 0.d0) then
              qyp(i,j,ns) = q(i,j,ns) + flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
           else
              qyp(i,j,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Im(i,j,2,2,ns) - q(i,j,ns))
           endif
        enddo
        
        ! minus state on face j+1
        do j = ilo2-1, ihi2
           v = q(i,j,QV)
           if (v .gt. 0.d0) then
              qym(i,j+1,ns) = q(i,j,ns) + flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
           else if (v .lt. 0.d0) then
              qym(i,j+1,ns) = q(i,j,ns)
           else
              qym(i,j+1,ns) = q(i,j,ns) + 0.5d0*flatn(i,j)*(Ip(i,j,2,2,ns) - q(i,j,ns))
           endif
        enddo
        
     enddo
  enddo
  
  deallocate(Ip,Im)

end subroutine trace_ppm_rad

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine cmpflx_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2,&
     qm,qp,qpd_l1,qpd_l2,qpd_h1,qpd_h2, &
      flx, flx_l1, flx_l2, flx_h1, flx_h2, &
     rflx,rflx_l1,rflx_l2,rflx_h1,rflx_h2, &
      pgd, pgd_l1, pgd_l2, pgd_h1, pgd_h2, &
     ergd,ergd_l1,ergd_l2,ergd_h1,ergd_h2, &
     lmgd,lmgd_l1,lmgd_l2,lmgd_h1,lmgd_h2, &
      ugd, ugd_l1, ugd_l2, ugd_h1, ugd_h2, &
      vgd, &
     gamc,gamcg,csml,c,qd_l1,qd_l2,qd_h1,qd_h2, &
     idir,ilo,ihi,jlo,jhi)

  use meth_params_module, only : QVAR, NVAR
  use radhydro_params_module, only : QRADVAR 
  use rad_params_module, only : ngroups

  implicit none

  integer lam_l1, lam_l2, lam_h1, lam_h2
  integer qpd_l1,qpd_l2,qpd_h1,qpd_h2
  integer flx_l1,flx_l2,flx_h1,flx_h2
  integer rflx_l1,rflx_l2,rflx_h1,rflx_h2
  integer pgd_l1,pgd_l2,pgd_h1,pgd_h2
  integer ergd_l1,ergd_l2,ergd_h1,ergd_h2
  integer lmgd_l1,lmgd_l2,lmgd_h1,lmgd_h2
  integer ugd_l1,ugd_l2,ugd_h1,ugd_h2
  integer qd_l1,qd_l2,qd_h1,qd_h2
  integer idir,ilo,ihi,jlo,jhi

  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision    qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision    qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision   flx( flx_l1: flx_h1, flx_l2: flx_h2,NVAR)
  double precision  rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,0:ngroups-1)
  double precision  pgd( pgd_l1: pgd_h1, pgd_l2: pgd_h2)
  double precision ergd(ergd_l1:ergd_h1,ergd_l2:ergd_h2,0:ngroups-1)
  double precision lmgd(lmgd_l1:lmgd_h1,lmgd_l2:lmgd_h2,0:ngroups-1)
  double precision  ugd( ugd_l1: ugd_h1, ugd_l2: ugd_h2)
  double precision  vgd( ugd_l1: ugd_h1, ugd_l2: ugd_h2)
  double precision  gamc (qd_l1:qd_h1,qd_l2:qd_h2)
  double precision  gamcg(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision     c (qd_l1:qd_h1,qd_l2:qd_h2)
  double precision   csml(qd_l1:qd_h1,qd_l2:qd_h2)

  !     Local variables
  integer i, j
  
  double precision, allocatable :: smallc(:,:), cavg(:,:)
  double precision, allocatable :: gamcm(:,:), gamcp(:,:)
  double precision, allocatable :: gamcgm(:,:), gamcgp(:,:)
  
  allocate ( smallc(ilo-1:ihi+1,jlo-1:jhi+1) )
  allocate (   cavg(ilo-1:ihi+1,jlo-1:jhi+1) )
  allocate (  gamcm(ilo-1:ihi+1,jlo-1:jhi+1) )
  allocate (  gamcp(ilo-1:ihi+1,jlo-1:jhi+1) )
  allocate ( gamcgm(ilo-1:ihi+1,jlo-1:jhi+1) )
  allocate ( gamcgp(ilo-1:ihi+1,jlo-1:jhi+1) )

  if(idir.eq.1) then
     do j = jlo, jhi
        do i = ilo, ihi+1
           smallc(i,j) = max( csml(i,j), csml(i-1,j) )
           cavg(i,j) = 0.5d0*( c(i,j) + c(i-1,j) )
           gamcm(i,j) = gamc(i-1,j)
           gamcp(i,j) = gamc(i,j)
           gamcgm(i,j) = gamcg(i-1,j)
           gamcgp(i,j) = gamcg(i,j)
        enddo
     enddo
  else
     do j = jlo, jhi+1
        do i = ilo, ihi
           smallc(i,j) = max( csml(i,j), csml(i,j-1) )
           cavg(i,j) = 0.5d0*( c(i,j) + c(i,j-1) )
           gamcm(i,j) = gamc(i,j-1)
           gamcp(i,j) = gamc(i,j)
           gamcgm(i,j) = gamcg(i,j-1)
           gamcgp(i,j) = gamcg(i,j)
        enddo
     enddo
  endif

!     Solve Riemann problem (godunov state passed back, but only (u,p,Er) saved
  call riemannus_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2,&
       qm, qp, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
       gamcm, gamcp, gamcgm, gamcgp, cavg, smallc, ilo-1, jlo-1, ihi+1, jhi+1, &
        flx,  flx_l1,  flx_l2,  flx_h1,  flx_h2, &
       rflx, rflx_l1, rflx_l2, rflx_h1, rflx_h2, &
        pgd,  pgd_l1,  pgd_l2,  pgd_h1,  pgd_h2, &
       ergd, ergd_l1, ergd_l2, ergd_h1, ergd_h2, &
       lmgd, lmgd_l1, lmgd_l2, lmgd_h1, lmgd_h2, &
        ugd,  ugd_l1,  ugd_l2,  ugd_h1,  ugd_h2, &
        vgd, &
       idir, ilo, ihi, jlo, jhi)
  
  deallocate (smallc,cavg,gamcp,gamcm,gamcgp,gamcgm)
  
end subroutine cmpflx_rad


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine riemannus_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, & 
     ql, qr, qpd_l1, qpd_l2, qpd_h1, qpd_h2, &
     gamcl, gamcr, gamcgl, gamcgr, cav, smallc, gd_l1, gd_l2, gd_h1, gd_h2, &
     uflx, uflx_l1, uflx_l2, uflx_h1, uflx_h2, &
     rflx, rflx_l1, rflx_l2, rflx_h1, rflx_h2, &
      pgdnv,  pgd_l1,  pgd_l2,  pgd_h1,  pgd_h2, &
     ergdnv, ergd_l1, ergd_l2, ergd_h1, ergd_h2, &
     lmgdnv, lmgd_l1, lmgd_l2, lmgd_h1, lmgd_h2, &
      ugdnv,  ugd_l1,  ugd_l2,  ugd_h1,  ugd_h2, &
      vgdnv, &
     idir, ilo1, ihi1, ilo2, ihi2)

  use network, only : nspec, naux
  use prob_params_module, only : physbc_lo,physbc_hi,Symmetry
  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
       URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv, &
       small_dens, small_pres
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor

  implicit none

  double precision, parameter:: small = 1.d-8

  integer lam_l1, lam_l2, lam_h1, lam_h2
  integer qpd_l1, qpd_l2, qpd_h1, qpd_h2
  integer gd_l1, gd_l2, gd_h1, gd_h2
  integer uflx_l1, uflx_l2, uflx_h1, uflx_h2
  integer rflx_l1, rflx_l2, rflx_h1, rflx_h2
  integer pgd_l1, pgd_l2, pgd_h1, pgd_h2
  integer ergd_l1, ergd_l2, ergd_h1, ergd_h2
  integer lmgd_l1, lmgd_l2, lmgd_h1, lmgd_h2
  integer ugd_l1, ugd_l2, ugd_h1, ugd_h2
  integer idir, ilo1, ihi1, ilo2, ihi2
  integer ilo,ihi,jlo,jhi

  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision ql(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision qr(qpd_l1:qpd_h1,qpd_l2:qpd_h2,QRADVAR)
  double precision gamcl(gd_l1:gd_h1,gd_l2:gd_h2)
  double precision gamcr(gd_l1:gd_h1,gd_l2:gd_h2)
  double precision gamcgl(gd_l1:gd_h1,gd_l2:gd_h2)
  double precision gamcgr(gd_l1:gd_h1,gd_l2:gd_h2)
  double precision cav(gd_l1:gd_h1,gd_l2:gd_h2)
  double precision smallc(gd_l1:gd_h1,gd_l2:gd_h2)
  double precision uflx(uflx_l1:uflx_h1,uflx_l2:uflx_h2,NVAR)
  double precision rflx(rflx_l1:rflx_h1,rflx_l2:rflx_h2,0:ngroups-1)
  double precision  pgdnv( pgd_l1: pgd_h1, pgd_l2: pgd_h2)
  double precision ergdnv(ergd_l1:ergd_h1,ergd_l2:ergd_h2,0:ngroups-1)
  double precision lmgdnv(lmgd_l1:lmgd_h1,lmgd_l2:lmgd_h2,0:ngroups-1)
  double precision  ugdnv( ugd_l1: ugd_h1, ugd_l2: ugd_h2)
  double precision  vgdnv( ugd_l1: ugd_h1, ugd_l2: ugd_h2)

  double precision rgdnv, ustar
  double precision, dimension(0:ngroups-1) :: erl, err
  double precision rl, ul, vl, pl, rel, pl_g, rel_g, wl
  double precision rr, ur, vr, pr, rer, pr_g, rer_g, wr
  double precision rstar, cstar, pstar
  double precision ro, uo, po, reo, co, gamco, reo_g, po_g, co_g, gamco_g
  double precision sgnm, spin, spout, ushock, frac
  double precision rhoetot, scr
  
  double precision wsmall, csmall

  double precision :: regdnv_g, pgdnv_g, pgdnv_t
  double precision :: drho, estar_g, pstar_g
  double precision, dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r
  double precision :: eddf, f1
  integer :: i, j, g

  integer iadv, ispec, iaux, n, nq
  double precision :: qavg

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

        pl = ql(i,j,qptot)
        rel = ql(i,j,qreitot)
        erl(:) = ql(i,j,qrad:qradhi)
        pl_g = ql(i,j,QPRES)
        rel_g = ql(i,j,QREINT)

        rr = qr(i,j,QRHO)

!  pick right velocities based on direction
        if(idir.eq.1) then
           ur = qr(i,j,QU)
           vr = qr(i,j,QV)
        else
           ur = qr(i,j,QV)
           vr = qr(i,j,QU)
        endif

        pr  = qr(i,j,qptot)
        rer = qr(i,j,qreitot)
        err(:) = qr(i,j,qrad:qradhi)
        pr_g = qr(i,j,QPRES)
        rer_g = qr(i,j,QREINT) 

        csmall = smallc(i,j)
        wsmall = small_dens*csmall
        wl = max(wsmall,sqrt(abs(gamcl(i,j)*pl*rl)))
        wr = max(wsmall,sqrt(abs(gamcr(i,j)*pr*rr)))

        pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
        pstar = max(pstar,small_pres)
        ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

        if (ustar .gt. 0.d0) then
           if (idir.eq.1) then
              lambda(:) = lam(i-1,j,:)
           else
              lambda(:) = lam(i,j-1,:)
           end if
           ro = rl
           uo = ul
           po = pl
           po_g = pl_g
           po_r(:) = erl(:) * lambda(:)
           reo = rel
           reo_r(:) = erl(:)
           reo_g = rel_g
           gamco = gamcl(i,j)
           gamco_g = gamcgl(i,j)
        else if (ustar .lt. 0.d0) then
           lambda(:) = lam(i,j,:)
           ro = rr
           uo = ur
           po = pr
           po_g = pr_g
           po_r(:) = err(:) * lambda(:)
           reo = rer
           reo_r(:) = err(:)
           reo_g = rer_g
           gamco = gamcr(i,j)
           gamco_g = gamcgr(i,j)
        else
           if (idir.eq.1) then
              do g=0, ngroups-1
                 lambda(g) = 2.0d0*(lam(i-1,j,g)*lam(i,j,g))/(lam(i-1,j,g)+lam(i,j,g)+1.d-50)
              end do
           else
              do g=0, ngroups-1              
                 lambda(g) = 2.0d0*(lam(i,j-1,g)*lam(i,j,g))/(lam(i,j-1,g)+lam(i,j,g)+1.d-50)
              end do
           end if
           ro = 0.5d0*(rl+rr)
           uo = 0.5d0*(ul+ur)
           reo = 0.5d0*(rel+rer)
           reo_r(:) = 0.5d0*(erl(:)+err(:))
           reo_g = 0.5d0*(rel_g+rer_g)
           po = 0.5d0*(pl+pr)
           po_r(:) = lambda(:) * reo_r(:)
           po_g = 0.5*(pr_g+pl_g)
           gamco = 0.5d0*(gamcl(i,j)+gamcr(i,j))
           gamco_g = 0.5d0*(gamcgl(i,j)+gamcgr(i,j))
        endif
        ro = max(small_dens,ro)

        co = sqrt(abs(gamco*po/ro))
        co = max(csmall,co)
        drho = (pstar - po)/co**2
        rstar = ro + drho
        rstar = max(small_dens,rstar)
        estar_g = reo_g + drho*(reo_g + po_g)/ro
        co_g = sqrt(abs(gamco_g*po_g/ro))
        co_g = max(csmall,co_g)
        pstar_g = po_g + drho*co_g**2
        pstar_g = max(pstar_g,small_pres)
        estar_r(:) = reo_r(:) + drho*(reo_r(:) + po_r(:))/ro
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
           vgdnv(i,j) = vl
        else if (ustar .lt. 0.d0) then
           vgdnv(i,j) = vr
        else
           vgdnv(i,j) = 0.5d0*(vl+vr)
        endif

        rgdnv = frac*rstar + (1.d0 - frac)*ro
        ugdnv(i,j) = frac*ustar + (1.d0 - frac)*uo
        pgdnv_t = frac*pstar + (1.d0 - frac)*po
        pgdnv_g = frac*pstar_g + (1.d0 - frac)*po_g
        regdnv_g = frac*estar_g + (1.d0 - frac)*reo_g
        regdnv_r(:) = frac*estar_r(:) + (1.d0 - frac)*reo_r(:)

        if (spout .lt. 0.d0) then
           rgdnv = ro
           ugdnv(i,j) = uo
           pgdnv_t = po
           pgdnv_g = po_g
           regdnv_g = reo_g
           regdnv_r(:) = reo_r(:)
        endif
        if (spin .ge. 0.d0) then
           rgdnv = rstar
           ugdnv(i,j) = ustar
           pgdnv_t = pstar
           pgdnv_g = pstar_g
           regdnv_g = estar_g
           regdnv_r(:) = estar_r(:)
        endif
        
        ! Enforce that fluxes through a symmetry plane are hard zero.
        if (i.eq.0 .and. physbc_lo(1) .eq. Symmetry .and. idir .eq. 1) ugdnv(i,j) = 0.d0
        if (j.eq.0 .and. physbc_lo(2) .eq. Symmetry .and. idir .eq. 2) ugdnv(i,j) = 0.d0

        do g=0, ngroups-1
           ergdnv(i,j,g) = max(regdnv_r(g), 0.d0)
        end do

        pgdnv(i,j) = pgdnv_g

        lmgdnv(i,j,:) = lambda

        ! Compute fluxes, order as conserved state (not q)
        uflx(i,j,URHO) = rgdnv*ugdnv(i,j)
        if(idir.eq.1) then
           uflx(i,j,UMX) = uflx(i,j,URHO)*ugdnv(i,j)
           uflx(i,j,UMY) = uflx(i,j,URHO)*vgdnv(i,j)
        else
           uflx(i,j,UMX) = uflx(i,j,URHO)*vgdnv(i,j)
           uflx(i,j,UMY) = uflx(i,j,URHO)*ugdnv(i,j)
        endif

        rhoetot = regdnv_g + 0.5d0*rgdnv*(ugdnv(i,j)**2 + vgdnv(i,j)**2)
        uflx(i,j,UEDEN) = ugdnv(i,j)*(rhoetot + pgdnv_g)
        uflx(i,j,UEINT) = ugdnv(i,j)*regdnv_g

        if (fspace_type.eq.1) then
           do g=0,ngroups-1
              eddf = Edd_factor(lambda(g))
              f1 = 0.5d0*(1.d0-eddf)
              rflx(i,j,g) = (1.d0+f1) * ergdnv(i,j,g) * ugdnv(i,j)
           end do
        else ! type 2
           do g=0,ngroups-1
              rflx(i,j,g) = ergdnv(i,j,g) * ugdnv(i,j)
           end do
        end if

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

     end do
  end do

end subroutine riemannus_rad

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine transy_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &   
     qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
      fy, fy_l1, fy_l2, fy_h1, fy_h2, &
     rfy,rfy_l1,rfy_l2,rfy_h1,rfy_h2, &
      pgdy,  pgdy_l1,  pgdy_l2,  pgdy_h1,  pgdy_h2, &
     ergdy, ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2, &
      ugdy,  ugdy_l1,  ugdy_l2,  ugdy_h1,  ugdy_h2, &
     gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
     srcQ, src_l1, src_l2, src_h1, src_h2, &
     grav, gv_l1, gv_l2, gv_h1, gv_h2, &
     hdt, cdtdy, ilo, ihi, jlo, jhi)
  
  use network, only : nspec, naux
  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
       URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
       fspace_type, comoving
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor

  implicit none

  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer qd_l1, qd_l2, qd_h1, qd_h2
  integer gc_l1, gc_l2, gc_h1, gc_h2
  integer fy_l1, fy_l2, fy_h1, fy_h2
  integer rfy_l1, rfy_l2, rfy_h1, rfy_h2
  integer pgdy_l1, pgdy_l2, pgdy_h1, pgdy_h2
  integer ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2
  integer ugdy_l1, ugdy_l2, ugdy_h1, ugdy_h2
  integer src_l1, src_l2, src_h1, src_h2
  integer gv_l1, gv_l2, gv_h1, gv_h2
  integer ilo, ihi, jlo, jhi

  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision qm (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision qp (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision  fy( fy_l1: fy_h1, fy_l2: fy_h2,NVAR)
  double precision rfy(rfy_l1:rfy_h1,rfy_l2:rfy_h2,0:ngroups-1)
  double precision  pgdy( pgdy_l1: pgdy_h1, pgdy_l2: pgdy_h2)
  double precision ergdy(ergdy_l1:ergdy_h1,ergdy_l2:ergdy_h2,0:ngroups-1)
  double precision  ugdy( ugdy_l1: ugdy_h1, ugdy_l2: ugdy_h2)
  double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
  double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
  double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
  double precision hdt, cdtdy

  integer i, j, g
  integer n, nq, iadv, ispec, iaux
  
  double precision rr,rrnew
  double precision ugp, ugm, dup, pav, du, pnewr,pnewl
  double precision rrr, rur, rvr, rer, ekinr, rhoekinr
  double precision rrnewr, runewr, rvnewr, renewr
  double precision rrl, rul, rvl, rel, ekinl, rhoekinl
  double precision rrnewl, runewl, rvnewl, renewl
  double precision rhotmp
  double precision compo, compn

  double precision pggp, pggm, dre, dmom
  double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
       der, ernewr, ernewl
  double precision eddf, f1, ugc

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
        lambda(:) = lam(i,j,:)
        pggp = pgdy(i,j+1)
        pggm = pgdy(i,j)
        ugp = ugdy(i,j+1)
        ugm = ugdy(i,j)
        ugc = 0.5d0*(ugp+ugm)
        ergp(:) = ergdy(i,j+1,:)
        ergm(:) = ergdy(i,j,:)

!           Convert to conservation form
        rrr = qp(i,j,QRHO)
        rur = rrr*qp(i,j,QU)
        rvr = rrr*qp(i,j,QV)
        ekinr = 0.5d0*rrr*(qp(i,j,QU)**2 + qp(i,j,QV)**2)
        rer = qp(i,j,QREINT) + ekinr
        err(:) = qp(i,j,qrad:qradhi)

        rrl = qm(i+1,j,QRHO)
        rul = rrl*qm(i+1,j,QU)
        rvl = rrl*qm(i+1,j,QV)
        ekinl = 0.5d0*rrl*(qm(i+1,j,QU)**2 + qm(i+1,j,QV)**2)
        rel = qm(i+1,j,QREINT) + ekinl
        erl(:) = qm(i+1,j,qrad:qradhi)

!           Add transverse predictor   
        rrnewr = rrr - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO)) 

        runewr = rur - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX)) 
        lamge(:) = lambda(:) * (ergp(:)-ergm(:))
        dmom = -cdtdy*((pggp-pggm) + sum(lamge))
        rvnewr = rvr - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
             + dmom
        luge(:) = ugc * lamge(:)
        dre = -cdtdy*sum(luge)
        renewr = rer - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) &
             + dre

        if (fspace_type .eq. 1 .and. comoving) then
           do g=0, ngroups-1
              eddf = Edd_factor(lambda(g))
              f1 = 0.5d0*(1.d0-eddf)
              der(g) = cdtdy * ugc * f1 * (ergp(g) - ergm(g))
           end do
        else if (fspace_type .eq. 2) then
           do g=0, ngroups-1
              eddf = Edd_factor(lambda(g))
              f1 = 0.5d0*(1.d0-eddf)
              der(g) = cdtdy * f1 * 0.5d0*(ergp(g)+ergm(g)) * (ugm-ugp)
           end do
        else ! mixed frame
           der(:) = cdtdy * luge
        end if

        ernewr(:) = err(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
             + der(:)

        rrnewl = rrl - cdtdy*(fy(i,j+1,URHO) - fy(i,j,URHO)) 
        runewl = rul - cdtdy*(fy(i,j+1,UMX)  - fy(i,j,UMX)) 
        rvnewl = rvl - cdtdy*(fy(i,j+1,UMY)  - fy(i,j,UMY)) &
             + dmom 
        renewl = rel - cdtdy*(fy(i,j+1,UEDEN)- fy(i,j,UEDEN)) &
             + dre
        ernewl(:) = erl(:) - cdtdy*(rfy(i,j+1,:) - rfy(i,j,:)) &
             + der(:)

        dup = pggp*ugp - pggm*ugm
        pav = 0.5d0*(pggp+pggm)
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
        qpo(i,j,qrad:qradhi) = ernewr(:)
        qpo(i,j,qptot  ) = sum(lambda*ernewr) + qpo(i,j,QPRES)
        qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)
        
        rhotmp =  rrnewl
        qmo(i+1,j,QRHO  ) = rhotmp            + hdt*srcQ(i,j,QRHO)
        qmo(i+1,j,QU    ) = runewl/rhotmp     + hdt*srcQ(i,j,QU) + hdt*grav(i,j,1)
        qmo(i+1,j,QV    ) = rvnewl/rhotmp     + hdt*srcQ(i,j,QV) + hdt*grav(i,j,2)
        rhoekinl = 0.5d0*(runewl**2+rvnewl**2)/rhotmp
        qmo(i+1,j,QREINT) = renewl - rhoekinl + hdt*srcQ(i,j,QREINT)
        qmo(i+1,j,QPRES ) = pnewl             + hdt*srcQ(i,j,QPRES)
        qmo(i+1,j,qrad:qradhi) = ernewl(:)
        qmo(i+1,j,qptot  ) = sum(lambda*ernewl) + qmo(i+1,j,QPRES)
        qmo(i+1,j,qreitot) = sum(qmo(i+1,j,qrad:qradhi)) + qmo(i+1,j,QREINT)
        
     enddo
  enddo

end subroutine transy_rad

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine transx_rad(lam, lam_l1, lam_l2, lam_h1, lam_h2, &   
     qm, qmo, qp, qpo, qd_l1, qd_l2, qd_h1, qd_h2, &
      fx,  fx_l1,  fx_l2,  fx_h1,  fx_h2, &
     rfx, rfx_l1, rfx_l2, rfx_h1, rfx_h2, &
      pgdx,  pgdx_l1,  pgdx_l2,  pgdx_h1,  pgdx_h2, &
     ergdx, ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2, &
      ugdx,  ugdx_l1,  ugdx_l2,  ugdx_h1,  ugdx_h2, &
     gamc, gc_l1, gc_l2, gc_h1, gc_h2, &
     srcQ, src_l1, src_l2, src_h1, src_h2, &
     grav, gv_l1, gv_l2, gv_h1, gv_h2, &
     hdt, cdtdx,  &
     area1, area1_l1, area1_l2, area1_h1, area1_h2, &
     vol, vol_l1, vol_l2, vol_h1, vol_h2, &
     ilo, ihi, jlo, jhi)
  
  use network, only : nspec, naux
  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QV, QPRES, QREINT, QFA, QFS, QFX, &
       URHO, UMX, UMY, UEDEN, UEINT, UFA, UFS, UFX, nadv
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, &
       fspace_type, comoving
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor

  implicit none

  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer qd_l1, qd_l2, qd_h1, qd_h2
  integer gc_l1, gc_l2, gc_h1, gc_h2
  integer fx_l1, fx_l2, fx_h1, fx_h2
  integer rfx_l1, rfx_l2, rfx_h1, rfx_h2
  integer pgdx_l1, pgdx_l2, pgdx_h1, pgdx_h2
  integer ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2
  integer ugdx_l1, ugdx_l2, ugdx_h1, ugdx_h2
  integer src_l1, src_l2, src_h1, src_h2
  integer gv_l1, gv_l2, gv_h1, gv_h2
  integer area1_l1, area1_l2, area1_h1, area1_h2
  integer vol_l1, vol_l2, vol_h1, vol_h2
  integer ilo, ihi, jlo, jhi
  
  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision qm (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision qmo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision qp (qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision qpo(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision  fx( fx_l1: fx_h1, fx_l2: fx_h2,NVAR)
  double precision rfx(rfx_l1:rfx_h1,rfx_l2:rfx_h2,0:ngroups-1)
  double precision  pgdx( pgdx_l1: pgdx_h1, pgdx_l2: pgdx_h2)
  double precision ergdx(ergdx_l1:ergdx_h1,ergdx_l2:ergdx_h2,0:ngroups-1)
  double precision  ugdx( ugdx_l1: ugdx_h1, ugdx_l2: ugdx_h2)
  double precision gamc(gc_l1:gc_h1,gc_l2:gc_h2)
  double precision srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
  double precision grav(gv_l1:gv_h1,gv_l2:gv_h2,2)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
  double precision hdt, cdtdx
  
  integer i, j, g
  integer n, nq, iadv
  integer ispec, iaux
  
  double precision rr, rrnew, compo, compn
  double precision rrr, rur, rvr, rer, ekinr, rhoekinr
  double precision rrnewr, runewr, rvnewr, renewr
  double precision rrl, rul, rvl, rel, ekinl, rhoekinl
  double precision rrnewl, runewl, rvnewl, renewl
  double precision ugp, ugm, dup, pav, du, pnewl,pnewr
  double precision rhotmp
  
  double precision pggp, pggm, dre, dmom
  double precision, dimension(0:ngroups-1) :: lambda, ergp, ergm, err, erl, lamge, luge, &
       der, ernewr, ernewl
  double precision eddf, f1, ugc, divu
  
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
        lambda(:) = lam(i,j,:)
        pggp = pgdx(i+1,j)
        pggm = pgdx(i,j)
        ugp = ugdx(i+1,j)
        ugm = ugdx(i,j)
        ugc = 0.5d0*(ugp+ugm)
        ergp(:) = ergdx(i+1,j,:)
        ergm(:) = ergdx(i,j,:)

!           Convert to conservation form
        rrr = qp(i,j,QRHO)
        rur = rrr*qp(i,j,QU)
        rvr = rrr*qp(i,j,QV)
        ekinr = 0.5d0*rrr*(qp(i,j,QU)**2 + qp(i,j,QV)**2)
        rer = qp(i,j,QREINT) + ekinr
        err(:) = qp(i,j,qrad:qradhi)

        rrl = qm(i,j+1,QRHO)
        rul = rrl*qm(i,j+1,QU)
        rvl = rrl*qm(i,j+1,QV)
        ekinl = 0.5d0*rrl*(qm(i,j+1,QU)**2 + qm(i,j+1,QV)**2)
        rel = qm(i,j+1,QREINT) + ekinl
        erl(:) = qm(i,j+1,qrad:qradhi)

!           Add transverse predictor  
        rrnewr = rrr - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
             area1(i,j)*fx(i,j,URHO))/vol(i,j) 
        lamge(:) = lambda(:) * (ergp(:)-ergm(:))
        dmom = -0.5d0*hdt*(area1(i+1,j)+area1(i,j))*((pggp-pggm)+sum(lamge))/vol(i,j) 
        runewr = rur - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
             area1(i,j)*fx(i,j,UMX))/vol(i,j) &
             + dmom
        rvnewr = rvr - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
             area1(i,j)*fx(i,j,UMY))/vol(i,j) 
        luge(:) = ugc * lamge(:)
        dre = -cdtdx*sum(luge)
        renewr = rer - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
             area1(i,j)*fx(i,j,UEDEN))/vol(i,j) &
             + dre

        if (fspace_type .eq. 1 .and. comoving) then
           do g=0, ngroups-1
              eddf = Edd_factor(lambda(g))
              f1 = 0.5d0*(1.d0-eddf)
              der(g) = cdtdx * ugc * f1 * (ergp(g) - ergm(g))
           end do
        else if (fspace_type .eq. 2) then
           divu = (area1(i+1,j)*ugp-area1(i,j)*ugm)/vol(i,j)
           do g=0, ngroups-1
              eddf = Edd_factor(lambda(g))
              f1 = 0.5d0*(1.d0-eddf)
              der(g) = -hdt * f1 * 0.5d0*(ergp(g)+ergm(g)) * divu
           end do
        else ! mixed frame
           der(:) = cdtdx * luge(:)
        end if

        ernewr(:) = err(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
             area1(i,j)*rfx(i,j,:))/vol(i,j) &
             + der(:)

        rrnewl = rrl - hdt*(area1(i+1,j)*fx(i+1,j,URHO) -  &
             area1(i,j)*fx(i,j,URHO))/vol(i,j) 
        runewl = rul - hdt*(area1(i+1,j)*fx(i+1,j,UMX)  -  &
             area1(i,j)*fx(i,j,UMX))/vol(i,j) &
             + dmom 
        rvnewl = rvl - hdt*(area1(i+1,j)*fx(i+1,j,UMY)  -  &
             area1(i,j)*fx(i,j,UMY))/vol(i,j) 
        renewl = rel - hdt*(area1(i+1,j)*fx(i+1,j,UEDEN)-  &
             area1(i,j)*fx(i,j,UEDEN))/vol(i,j) &
             + dre
        ernewl(:) = erl(:) - hdt*(area1(i+1,j)*rfx(i+1,j,:)-  &
             area1(i,j)*rfx(i,j,:))/vol(i,j) &
             + der(:)

        dup = pggp*ugp - pggm*ugm
        pav = 0.5d0*(pggp+pggm)
        du = ugp-ugm

        pnewr = qp(i,j  ,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j)-1.d0))
        pnewl = qm(i,j+1,QPRES) - cdtdx*(dup + pav*du*(gamc(i,j)-1.d0))

        !           Convert back to non-conservation form
        rhotmp = rrnewr
        qpo(i,j,QRHO) = rhotmp        + hdt*srcQ(i,j,QRHO)
        qpo(i,j,QU  ) = runewr/rhotmp + hdt*srcQ(i,j,QU)  + hdt*grav(i,j,1)
        qpo(i,j,QV  ) = rvnewr/rhotmp + hdt*srcQ(i,j,QV)  + hdt*grav(i,j,2)
        rhoekinr = 0.5d0*(runewr**2+rvnewr**2)/rhotmp
        qpo(i,j,QREINT)= renewr - rhoekinr + hdt*srcQ(i,j,QREINT)
        qpo(i,j,QPRES) =  pnewr            + hdt*srcQ(i,j,QPRES)
        qpo(i,j,qrad:qradhi) = ernewr(:)
        qpo(i,j,qptot)   = sum(lambda*ernewr) + qpo(i,j,QPRES)
        qpo(i,j,qreitot) = sum(qpo(i,j,qrad:qradhi)) + qpo(i,j,QREINT)
                 
        !           Convert back to non-conservation form
        rhotmp = rrnewl
        qmo(i,j+1,QRHO) = rhotmp         + hdt*srcQ(i,j,QRHO)
        qmo(i,j+1,QU  ) = runewl/rhotmp  + hdt*srcQ(i,j,QU)  + hdt*grav(i,j,1)
        qmo(i,j+1,QV  ) = rvnewl/rhotmp  + hdt*srcQ(i,j,QV)  + hdt*grav(i,j,2)
        rhoekinl = 0.5d0*(runewl**2+rvnewl**2)/rhotmp
        qmo(i,j+1,QREINT)= renewl - rhoekinl +hdt*srcQ(i,j,QREINT)
        qmo(i,j+1,QPRES) = pnewl +hdt*srcQ(i,j,QPRES)
        qmo(i,j+1,qrad:qradhi) = ernewl(:)
        qmo(i,j+1,qptot)   = sum(lambda*ernewl) + qmo(i,j+1,QPRES)
        qmo(i,j+1,qreitot) = sum(qmo(i,j+1,qrad:qradhi)) + qmo(i,j+1,QREINT)
 
     enddo
  enddo

end subroutine transx_rad

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine consup_rad( uin, uin_l1, uin_l2, uin_h1, uin_h2, &
     uout,uout_l1,uout_l2,uout_h1,uout_h2, &
      Erin, Erin_l1, Erin_l2, Erin_h1, Erin_h2, &
     Erout,Erout_l1,Erout_l2,Erout_h1,Erout_h2, &
     pgdx,pgdx_l1,pgdx_l2,pgdx_h1,pgdx_h2, &
     pgdy,pgdy_l1,pgdy_l2,pgdy_h1,pgdy_h2, &
     ergdx,ergdx_l1,ergdx_l2,ergdx_h1,ergdx_h2, &
     ergdy,ergdy_l1,ergdy_l2,ergdy_h1,ergdy_h2, &
     lmgdx,lmgdx_l1,lmgdx_l2,lmgdx_h1,lmgdx_h2, &
     lmgdy,lmgdy_l1,lmgdy_l2,lmgdy_h1,lmgdy_h2, &
     ugdx,ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2, &
     ugdy,ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2, &
     src , src_l1, src_l2, src_h1, src_h2, &
     grav,  gv_l1,  gv_l2,  gv_h1,  gv_h2, &
     flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
     flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
     rflux1,rflux1_l1,rflux1_l2,rflux1_h1,rflux1_h2, &
     rflux2,rflux2_l1,rflux2_l2,rflux2_h1,rflux2_h2, &
     area1,area1_l1,area1_l2,area1_h1,area1_h2, &
     area2,area2_l1,area2_l2,area2_h1,area2_h2, &
     vol,vol_l1,vol_l2,vol_h1,vol_h2, &
     div,pdivu, uy_xfc, ux_yfc, &
     lo,hi,dx,dy,dt, nstep_fsp)

  use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, &
       normalize_species
  use rad_params_module, only : ngroups, nugroup, dlognu
  use radhydro_params_module, only : fspace_type, comoving
  use radhydro_nd_module, only : advect_in_fspace
  use fluxlimiter_module, only : Edd_factor
  use advection_module, only : normalize_species_fluxes

  implicit none

  integer nstep_fsp
  integer lo(2), hi(2)
  integer uin_l1,uin_l2,uin_h1,uin_h2
  integer uout_l1,uout_l2,uout_h1,uout_h2
  integer Erout_l1,Erout_l2,Erout_h1,Erout_h2
  integer Erin_l1,Erin_l2,Erin_h1,Erin_h2
  integer pgdx_l1,pgdx_l2,pgdx_h1,pgdx_h2
  integer pgdy_l1,pgdy_l2,pgdy_h1,pgdy_h2
  integer ergdx_l1,ergdx_l2,ergdx_h1,ergdx_h2
  integer ergdy_l1,ergdy_l2,ergdy_h1,ergdy_h2
  integer lmgdx_l1,lmgdx_l2,lmgdx_h1,lmgdx_h2
  integer lmgdy_l1,lmgdy_l2,lmgdy_h1,lmgdy_h2
  integer ugdx_l1,ugdx_l2,ugdx_h1,ugdx_h2
  integer ugdy_l1,ugdy_l2,ugdy_h1,ugdy_h2
  integer   src_l1,  src_l2,  src_h1,  src_h2
  integer    gv_l1,   gv_l2,   gv_h1,   gv_h2
  integer flux1_l1,flux1_l2,flux1_h1,flux1_h2
  integer flux2_l1,flux2_l2,flux2_h1,flux2_h2
  integer rflux1_l1,rflux1_l2,rflux1_h1,rflux1_h2
  integer rflux2_l1,rflux2_l2,rflux2_h1,rflux2_h2
  integer area1_l1,area1_l2,area1_h1,area1_h2
  integer area2_l1,area2_l2,area2_h1,area2_h2
  integer vol_l1,vol_l2,vol_h1,vol_h2

  double precision uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
  double precision  Erin( Erin_l1: Erin_h1, Erin_l2: Erin_h2,0:ngroups-1)
  double precision Erout(Erout_l1:Erout_h1,Erout_l2:Erout_h2,0:ngroups-1)
  double precision  pgdx( pgdx_l1: pgdx_h1, pgdx_l2: pgdx_h2)
  double precision  pgdy( pgdy_l1: pgdy_h1, pgdy_l2: pgdy_h2)
  double precision ergdx(ergdx_l1:ergdx_h1,ergdx_l2:ergdx_h2,0:ngroups-1)
  double precision ergdy(ergdy_l1:ergdy_h1,ergdy_l2:ergdy_h2,0:ngroups-1)
  double precision lmgdx(lmgdx_l1:lmgdx_h1,lmgdx_l2:lmgdx_h2,0:ngroups-1)
  double precision lmgdy(lmgdy_l1:lmgdy_h1,lmgdy_l2:lmgdy_h2,0:ngroups-1)
  double precision  ugdx( ugdx_l1: ugdx_h1, ugdx_l2: ugdx_h2)
  double precision  ugdy( ugdy_l1: ugdy_h1, ugdy_l2: ugdy_h2)
  double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,NVAR)
  double precision  grav(   gv_l1:   gv_h1,   gv_l2:   gv_h2,2)
  double precision  flux1( flux1_l1: flux1_h1, flux1_l2: flux1_h2,NVAR)
  double precision  flux2( flux2_l1: flux2_h1, flux2_l2: flux2_h2,NVAR)
  double precision rflux1(rflux1_l1:rflux1_h1,rflux1_l2:rflux1_h2,0:ngroups-1)
  double precision rflux2(rflux2_l1:rflux2_h1,rflux2_l2:rflux2_h2,0:ngroups-1)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
  double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
  double precision pdivu(lo(1):hi(1),lo(2):hi(2))
  double precision uy_xfc( ugdx_l1: ugdx_h1, ugdx_l2: ugdx_h2)
  double precision ux_yfc( ugdy_l1: ugdy_h1, ugdy_l2: ugdy_h2)
  double precision dx, dy, dt

  integer i, j, n, g

  double precision div1 
  double precision SrU, SrV
  double precision rho, Up, Vp, SrE

  double precision, dimension(0:ngroups-1) :: Erscale
  double precision, dimension(0:ngroups-1) :: ustar, af
  double precision :: Eddf, Eddfxm, Eddfxp, Eddfym, Eddfyp, f1, f2, f1xm, f1xp, f1ym, f1yp
  double precision :: Gf1E(2)
  double precision :: ux, uy, divu, lamc, Egdc
  double precision :: dudx(2), dudy(2), nhat(2), GnDotu(2), nnColonDotGu
  double precision :: dpdx, dprdx, dpdy, dprdy, ek1, ek2, dek

  if (ngroups .gt. 1) then
     if (fspace_type .eq. 1) then
        Erscale = dlognu
     else
        Erscale = nugroup*dlognu
     end if
  end if

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

  do g = 0, ngroups-1
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        div1 = .5d0*(div(i,j) + div(i,j+1))
        div1 = difmag*min(0.d0,div1)
        rflux1(i,j,g) = rflux1(i,j,g) &
             + dx*div1*(Erin(i,j,g) - Erin(i-1,j,g))
        rflux1(i,j,g) = area1(i,j)*rflux1(i,j,g)
     enddo
  enddo
  enddo

  do g = 0, ngroups-1
  do j = lo(2),hi(2)+1
     do i = lo(1),hi(1)
        div1 = .5d0*(div(i,j) + div(i+1,j))
        div1 = difmag*min(0.d0,div1)
        rflux2(i,j,g) = rflux2(i,j,g) &
             + dy*div1*(Erin(i,j,g) - Erin(i,j-1,g))
        rflux2(i,j,g) = area2(i,j)*rflux2(i,j,g)
     enddo
  enddo
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

  do g = 0, ngroups-1
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        Erout(i,j,g) = Erin(i,j,g) + dt * &
             ( rflux1(i,j,g) - rflux1(i+1,j,g) &
             +   rflux2(i,j,g) - rflux2(i,j+1,g) ) / vol(i,j) 
     enddo
  enddo
  end do

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

  ! Add gradp term to momentum equation
  do j = lo(2),hi(2)
     do i = lo(1),hi(1)
        dpdx  = ( pgdx(i+1,j)- pgdx(i,j))/ dx
        dpdy  = ( pgdy(i,j+1)- pgdy(i,j))/ dy

        dprdx = 0.d0
        dprdy = 0.d0
        do g=0,ngroups-1
           lamc = 0.25d0*(lmgdx(i,j,g)+lmgdx(i+1,j,g)+lmgdy(i,j,g)+lmgdy(i,j+1,g))
           dprdx = dprdx + lamc*(ergdx(i+1,j,g)-ergdx(i,j,g))/dx
           dprdy = dprdy + lamc*(ergdy(i,j+1,g)-ergdy(i,j,g))/dy
        end do

        uout(i,j,UMX) = uout(i,j,UMX) - dt * dpdx
        uout(i,j,UMY) = uout(i,j,UMY) - dt * dpdy
        ek1 = (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / (2.d0*uout(i,j,URHO))

        uout(i,j,UMX) = uout(i,j,UMX) - dt * dprdx
        uout(i,j,UMY) = uout(i,j,UMY) - dt * dprdy
        ek2 = (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / (2.d0*uout(i,j,URHO))
        dek = ek2 - ek1

        uout(i,j,UEDEN) = uout(i,j,UEDEN) + dek
        if (.not. comoving) then ! mixed-frame (single group only)
           Erout(i,j,0) = Erout(i,j,0) - dek
        end if
     enddo
  enddo

  ! Add radiation source term to rho*u, rhoE, and Er
  if (comoving) then
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           ux = 0.5d0*(ugdx(i,j) + ugdx(i+1,j))
           uy = 0.5d0*(ugdy(i,j) + ugdy(i,j+1))
           
           divu = (ugdx(i+1,j  )*area1(i+1,j  ) - ugdx(i,j)*area1(i,j) &
                +  ugdy(i  ,j+1)*area2(i  ,j+1) - ugdy(i,j)*area2(i,j) &
                & ) / vol(i,j) 

           dudx(1) = (ugdx(i+1,j)-ugdx(i,j))/dx 
           dudx(2) = (uy_xfc(i+1,j)-uy_xfc(i,j))/dx 

           dudy(1) = (ux_yfc(i,j+1)-ux_yfc(i,j))/dy 
           dudy(2) = (ugdy(i,j+1)-ugdy(i,j))/dy 

           ! Note that for single group, fspace_type is always 1
           do g=0, ngroups-1
              
              nhat(1) = (ergdx(i+1,j,g)-ergdx(i,j,g))/dx
              nhat(2) = (ergdy(i,j+1,g)-ergdy(i,j,g))/dy

              GnDotu(1) = dot_product(nhat, dudx)
              GnDotu(2) = dot_product(nhat, dudy)

              nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.d-50)

              lamc = 0.25d0*(lmgdx(i,j,g)+lmgdx(i+1,j,g)+lmgdy(i,j,g)+lmgdy(i,j+1,g))
              Eddf = Edd_factor(lamc)
              f1 = (1.d0-Eddf)*0.5d0
              f2 = (3.d0*Eddf-1.d0)*0.5d0
              af(g) = -(f1*divu + f2*nnColonDotGu)

              if (fspace_type .eq. 1) then
                 Eddfxp = Edd_factor(lmgdx(i+1,j  ,g))
                 Eddfxm = Edd_factor(lmgdx(i  ,j  ,g))
                 Eddfyp = Edd_factor(lmgdy(i  ,j+1,g))
                 Eddfym = Edd_factor(lmgdy(i  ,j  ,g)) 

                 f1xp = 0.5d0*(1.d0-Eddfxp)
                 f1xm = 0.5d0*(1.d0-Eddfxm)
                 f1yp = 0.5d0*(1.d0-Eddfyp)
                 f1ym = 0.5d0*(1.d0-Eddfym)

                 Gf1E(1) = (f1xp*ergdx(i+1,j,g) - f1xm*ergdx(i,j,g)) / dx
                 Gf1E(2) = (f1yp*ergdy(i,j+1,g) - f1ym*ergdy(i,j,g)) / dy
              
                 Egdc = 0.25d0*(ergdx(i,j,g)+ergdx(i+1,j,g)+ergdy(i,j,g)+ergdy(i,j+1,g))

                 Erout(i,j,g) = Erout(i,j,g) + dt*(ux*Gf1E(1)+uy*Gf1E(2)) &
                      - dt*f2*Egdc*nnColonDotGu
              end if
              
           end do
           
           if (ngroups.gt.1) then
              ustar = Erout(i,j,:) / Erscale
              call advect_in_fspace(ustar, af, dt, nstep_fsp)
              Erout(i,j,:) = ustar * Erscale
           end if
        end do
     end do
  end if


  do j = lo(2),hi(2)
     do i = lo(1),hi(1)+1
        flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR)
        flux1(i,j,   UMX) = flux1(i,j,UMX) + dt*area1(i,j)*pgdx(i,j)
     enddo
  enddo

  do g = 0, ngroups-1
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)+1
           rflux1(i,j,g) = dt * rflux1(i,j,g)
        enddo
     enddo
  end do

  do j = lo(2),hi(2)+1 
     do i = lo(1),hi(1)
        flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR)
        flux2(i,j,UMY) = flux2(i,j,UMY) + dt*area2(i,j)*pgdy(i,j)
     enddo
  enddo

  do g = 0, ngroups-1
     do j = lo(2),hi(2)+1 
        do i = lo(1),hi(1)
           rflux2(i,j,g) = dt * rflux2(i,j,g)
        enddo
     enddo
  end do  

end subroutine consup_rad


subroutine ppflaten(lof, hif, &
     flatn, q, q_l1,q_l2, q_h1,q_h2)
  use meth_params_module, only : QPRES, QU, QV
  use radhydro_params_module, only : flatten_pp_threshold, QRADVAR, qptot
  implicit none
  integer, intent(in) :: lof(2), hif(2), q_l1, q_h1, q_l2, q_h2
  double precision, intent(in) :: q(q_l1:q_h1,q_l2:q_h2,QRADVAR)
  double precision, intent(inout) :: flatn(q_l1:q_h1,q_l2:q_h2)

  integer :: i,j

  do j=lof(2),hif(2)
  do i=lof(1),hif(1)
     if (q(i-1,j,QU)+q(i,j-1,QV) > q(i+1,j,QU)+q(i,j+1,QV)) then
        if (q(i,j,QPRES) < flatten_pp_threshold* q(i,j,qptot)) then
           flatn(i,j) = 0.d0
        end if
     end if
  end do
  end do
end subroutine ppflaten

end module rad_advection_module
