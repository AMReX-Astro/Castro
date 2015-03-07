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
     src,src_l1,src_l2,src_h1,src_h2, &
     srcQ,srQ_l1,srQ_l2,srQ_h1,srQ_h2, &
     courno,dx,dy,dt,ngp,ngf,iflaten)

!     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
!     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
!     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
!     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
!     routine that computes flatn).  

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
       QVAR, QRHO, QU, QV, QGAME, QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
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
  integer srQ_l1,srQ_l2,srQ_h1,srQ_h2
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
  double precision :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,QVAR)
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

  q(:,:,QGAME) = 0.d0 ! QGAME is not used in radiation hydro. Setting it to 0 to mute valgrind.

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
  use meth_params_module, only : NVAR, ppm_type 
  use radhydro_params_module, only : QRADVAR
  use rad_params_module, only : ngroups
  use riemann_rad_module, only : cmpflx_rad
  use trace_ppm_rad_module, only : trace_ppm_rad
  
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
       URHO, UMX, UMY, UEDEN, UFA, UFS, UFX, nadv
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
       URHO, UMX, UMY, UEDEN, UFA, UFS, UFX, nadv
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
