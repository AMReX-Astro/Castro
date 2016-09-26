module rad_advection_module

  use bl_constants_module
  
  implicit none

  private

  public umeth2d_rad, ctoprim_rad, consup_rad, ppflaten

contains

! :::
! ::: ------------------------------------------------------------------
! :::
subroutine ctoprim_rad(lo,hi, &
                       uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                       Erin,Erin_l1,Erin_l2,Erin_h1,Erin_h2, &
                       lam,lam_l1,lam_l2,lam_h1,lam_h2, &
                       q,q_l1,q_l2,q_h1,q_h2, &
                       qaux,qa_l1,qa_l2,qa_h1,qa_h2, &
                       src,src_l1,src_l2,src_h1,src_h2, &
                       srcQ,srQ_l1,srQ_l2,srQ_h1,srQ_h2, &
                       dx,dy,dt,ngp,ngf)

  ! Will give primitive variables on lo-ngp:hi+ngp.  Declared
  ! dimensions of q,c,gamc,csml are given by DIMS(q).  This
  ! declared region is assumed to encompass lo-ngp:hi+ngp.

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, &
                                 QVAR, QRHO, QU, QV, QW, QGAME, QREINT, QPRES, &
                                 QTEMP, QFS, QFX, &
                                 NQAUX, QGAMC, QGAMCG, QC, QCG, QCSML, QDPDR, QDPDE,  &
                                 nadv, allow_negative_energy, small_temp, &
                                 npassive, upass_map, qpass_map
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, comoving
  use rad_params_module, only : ngroups
  use rad_util_module, only : compute_ptot_ctot

  implicit none

  double precision, parameter:: small = 1.d-8

  integer lo(2), hi(2)
  integer uin_l1,uin_l2,uin_h1,uin_h2
  integer Erin_l1,Erin_l2,Erin_h1,Erin_h2
  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer q_l1,q_l2,q_h1,q_h2
  integer qa_l1,qa_l2,qa_h1,qa_h2
  integer src_l1,src_l2,src_h1,src_h2
  integer srQ_l1,srQ_l2,srQ_h1,srQ_h2

  double precision :: uin(uin_l1:uin_h1,uin_l2:uin_h2,NVAR)
  double precision :: Erin(Erin_l1:Erin_h1,Erin_l2:Erin_h2,0:ngroups-1)
  double precision :: lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision :: q(q_l1:q_h1,q_l2:q_h2,QRADVAR)
  double precision :: qaux(qa_l1:qa_h1,qa_l2:qa_h2,QRADVAR)
  double precision :: src (src_l1:src_h1,src_l2:src_h2,NVAR)
  double precision :: srcQ(srQ_l1:srQ_h1,srQ_l2:srQ_h2,QVAR)
  double precision :: dx, dy, dt

  integer          :: i, j, g
  integer          :: ngp, ngf, loq(2), hiq(2)
  integer          :: n, nq, ipassive
  double precision :: courx, coury, courmx, courmy

  double precision :: csrad2, ptot, ctot, gamc_tot

  type(eos_t) :: eos_state

  do i=1,2
     loq(i) = lo(i)-ngp
     hiq(i) = hi(i)+ngp
  enddo

  ! Make q (all but p), except put e in slot for rho.e, fix after eos
  ! call.  The temperature is used as an initial guess for the eos
  ! call and will be overwritten
  do j = loq(2),hiq(2)
     do i = loq(1),hiq(1)

        if (uin(i,j,URHO) .le. ZERO) then
           print *,'   '
           print *,'>>> Error: Castro_2d::ctoprim ',i,j
           print *,'>>> ... negative density ', uin(i,j,URHO)
           print *,'    '
           call bl_error("Error:: Castro_2d.f90 :: ctoprim")
        end if

        q(i,j,QRHO) = uin(i,j,URHO)
        q(i,j,QU) = uin(i,j,UMX)/uin(i,j,URHO)
        q(i,j,QV) = uin(i,j,UMY)/uin(i,j,URHO)

        ! we should set this based on the kinetic energy, using dual_energy_eta1
        q(i,j,QREINT ) = uin(i,j,UEINT)/q(i,j,QRHO)

        q(i,j,QTEMP  ) = uin(i,j,UTEMP)
        q(i,j,qrad:qradhi) = Erin(i,j,:)
     enddo
  enddo

  ! Load advected / species / and auxiliary quatities, c, into q,
  ! assuming they arrived in uin as rho.c
  do ipassive = 1, npassive
     n  = upass_map(ipassive)
     nq = qpass_map(ipassive)
     do j = loq(2),hiq(2)
        do i = loq(1),hiq(1)
           q(i,j,nq) = uin(i,j,n)/q(i,j,QRHO)
        enddo
     enddo
  enddo


  ! Get gamc, p, T, c, csml using q state
  do j = loq(2), hiq(2)
     do i = loq(1), hiq(1)

        eos_state % rho = q(i,j,QRHO)
        eos_state % T   = q(i,j,QTEMP)
        eos_state % e   = q(i,j,QREINT)
        eos_state % xn  = q(i,j,QFS:QFS+nspec-1)
        eos_state % aux = q(i,j,QFX:QFX+naux-1)

        ! If necessary, reset the energy using small_temp -- TODO: we handle this 
        ! differently now
        if ((allow_negative_energy .eq. 0) .and. (q(i,j,QREINT) .lt. 0)) then
           q(i,j,QTEMP) = small_temp
           eos_state % T = q(i,j,QTEMP)
           call eos(eos_input_rt, eos_state)
           q(i,j,QPRES ) = eos_state % p
           q(i,j,QREINT) = eos_state % e
           if (q(i,j,QREINT) .lt. ZERO) then
              print *,'   '
              print *,'>>> Error: ctoprim ',i,j
              print *,'>>> ... new e from eos call is negative ',q(i,j,QREINT)
              print *,'    '
              call bl_error("Error:: ctoprim")
           end if
        end if

        call eos(eos_input_re, eos_state)

        q(i,j,QTEMP) = eos_state % T
        q(i,j,QREINT) = eos_state % e * q(i,j,QRHO)
        q(i,j,QPRES) = eos_state % p
        q(i,j,QGAME) = q(i,j,QPRES) / q(i,j,QREINT) + ONE

        qaux(i,j,QGAMCG)   = eos_state % gam1
        qaux(i,j,QCG)      = eos_state % cs
        qaux(i,j,QDPDR)    = eos_state % dpdr_e
        qaux(i,j,QDPDE)    = eos_state % dpde

        call compute_ptot_ctot(lam(i,j,:), q(i,j,:), qaux(i,j,QCG), &
                               ptot, ctot, gamc_tot)

        q(i,j,QPTOT) = ptot

        qaux(i,j,QC)    = ctot
        qaux(i,j,QGAMC) = gamc_tot

        qaux(i,j,QCSML) = max(small, small * qaux(i,j,QC))

        q(i,j,qreitot) = q(i,j,QREINT) + sum(q(i,j,qrad:qradhi))

     enddo
  enddo

  ! Compute sources in terms of Q
  do j = loq(2), hiq(2)
     do i = loq(1), hiq(1)

        srcQ(i,j,QRHO  ) = src(i,j,URHO)
        srcQ(i,j,QU    ) = (src(i,j,UMX) - q(i,j,QU) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
        srcQ(i,j,QV    ) = (src(i,j,UMY) - q(i,j,QV) * srcQ(i,j,QRHO)) / q(i,j,QRHO)
        srcQ(i,j,QREINT) = src(i,j,UEINT)
        srcQ(i,j,QPRES ) = qaux(i,j,QDPDE) * &
                (srcQ(i,j,QREINT) - q(i,j,QREINT)*srcQ(i,j,QRHO)/q(i,j,QRHO))/q(i,j,QRHO) + &
                qaux(i,j,QDPDR) * srcQ(i,j,QRHO)! + &
!                sum(dpdX_er(i,j,:)*(src(i,j,UFS:UFS+nspec-1) - &
!                    q(i,j,QFS:QFS+nspec-1)*srcQ(i,j,QRHO))) / q(i,j,QRHO)

        do ipassive = 1, npassive
           n  = upass_map(ipassive)
           nq = qpass_map(ipassive)
           srcQ(i,j,nq) = (src(i,j,n) - q(i,j,nq) * srcQ(i,j,QRHO))/q(i,j,QRHO)
        enddo

     end do
  end do

end subroutine ctoprim_rad


! ::: ---------------------------------------------------------------
! ::: :: UMETH2D     Compute hyperbolic fluxes using unsplit second
! ::: ::               order Godunov integrator.
! ::: ::
! ::: :: inputs/outputs
! ::: :: q           => (const)  input state, primitives
! ::: :: qaux        => (const)  auxillary hydro info
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

subroutine umeth2d_rad(q, qd_l1, qd_l2, qd_h1, qd_h2, &
                       qaux, qa_l1, qa_l2, qa_h1, qa_h2, &
                       lam, lam_l1, lam_l2, lam_h1, lam_h2, &
                       flatn, &
                       srcQ, src_l1, src_l2, src_h1, src_h2, &
                       ilo1, ilo2, ihi1, ihi2, dx, dy, dt, &
                       flux1, fd1_l1, fd1_l2, fd1_h1, fd1_h2, &
                       flux2, fd2_l1, fd2_l2, fd2_h1, fd2_h2, &
                       rflux1, rfd1_l1, rfd1_l2, rfd1_h1, rfd1_h2, &
                       rflux2, rfd2_l1, rfd2_l2, rfd2_h1, rfd2_h2, &
                       q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                       q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                       area1, area1_l1, area1_l2, area1_h1, area1_h2, &
                       area2, area2_l1, area2_l2, area2_h1, area2_h2, &
                       pdivu, vol, vol_l1, vol_l2, vol_h1, vol_h2, &
                       dloga, dloga_l1, dloga_l2, dloga_h1, dloga_h2, &
                       domlo, domhi)

  use network, only : nspec
  use meth_params_module, only : NVAR, QVAR, ppm_type, GDPRES, &
       GDU, GDV, GDERADS, GDLAMS, QC, QCG, QCSML, QGAMC, QGAMCG, NQAUX, NGDNV

  use radhydro_params_module, only : QRADVAR
  use rad_params_module, only : ngroups
  use riemann_module, only : cmpflx
  use trace_ppm_rad_module, only : trace_ppm_rad
  use transverse_module

  implicit none

  integer qd_l1, qd_l2, qd_h1, qd_h2
  integer qa_l1, qa_l2, qa_h1, qa_h2
  integer lam_l1,lam_l2,lam_h1,lam_h2
  integer ergdx_l1, ergdx_l2, ergdx_h1, ergdx_h2
  integer ergdy_l1, ergdy_l2, ergdy_h1, ergdy_h2
  integer lmgdx_l1, lmgdx_l2, lmgdx_h1, lmgdx_h2
  integer lmgdy_l1, lmgdy_l2, lmgdy_h1, lmgdy_h2
  integer rfd1_l1, rfd1_l2, rfd1_h1, rfd1_h2
  integer rfd2_l1, rfd2_l2, rfd2_h1, rfd2_h2
  integer q1_l1, q1_l2, q1_h1, q1_h2
  integer q2_l1, q2_l2, q2_h1, q2_h2
  integer dloga_l1, dloga_l2, dloga_h1, dloga_h2
  integer src_l1, src_l2, src_h1, src_h2
  integer fd1_l1, fd1_l2, fd1_h1, fd1_h2
  integer fd2_l1, fd2_l2, fd2_h1, fd2_h2
  integer area1_l1, area1_l2, area1_h1, area1_h2
  integer area2_l1, area2_l2, area2_h1, area2_h2
  integer vol_l1, vol_l2, vol_h1, vol_h2
  integer ilo1, ilo2, ihi1, ihi2
  integer domlo(2), domhi(2)
  
  double precision dx, dy, dt
  double precision     q(qd_l1:qd_h1,qd_l2:qd_h2,QRADVAR)
  double precision  qaux(qa_l1:qa_h1,qa_l2:qa_h2,NQAUX)
  double precision flatn(qd_l1:qd_h1,qd_l2:qd_h2)
  double precision  srcQ(src_l1:src_h1,src_l2:src_h2,QVAR)
  double precision dloga(dloga_l1:dloga_h1,dloga_l2:dloga_h2)
  double precision q1(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv)
  double precision q2(q2_l1:q2_h1,q2_l2:q2_h2,ngdnv)
  double precision flux1(fd1_l1:fd1_h1,fd1_l2:fd1_h2,NVAR)
  double precision flux2(fd2_l1:fd2_h1,fd2_l2:fd2_h2,NVAR)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
  double precision pdivu(ilo1:ihi1,ilo2:ihi2)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
  double precision lam(lam_l1:lam_h1,lam_l2:lam_h2,0:ngroups-1)
  double precision rflux1(rfd1_l1:rfd1_h1,rfd1_l2:rfd1_h2,0:ngroups-1)
  double precision rflux2(rfd2_l1:rfd2_h1,rfd2_l2:rfd2_h2,0:ngroups-1)

  ! Left and right state arrays (edge centered, cell centered)
  double precision, allocatable :: qm(:,:,:), qp(:,:,:)
  double precision, allocatable ::qxm(:,:,:), qym(:,:,:)
  double precision, allocatable ::qxp(:,:,:), qyp(:,:,:)

  ! Work arrays to hold riemann state and conservative fluxes
  double precision, allocatable::   fx(:,:,:),  fy(:,:,:)
  double precision, allocatable::  rfx(:,:,:), rfy(:,:,:)
  double precision, allocatable::   qgdxtmp(:,:,:)

  double precision, allocatable :: shk(:,:)

  ! Local scalar variables
  double precision :: dtdx
  double precision :: hdtdx, hdt, hdtdy
  integer          :: i,j
  double precision :: pggdx, pggdy

  allocate ( qgdxtmp(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv))

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

  allocate (shk(ilo1-1:ihi1+1,ilo2-1:ihi2+1))

  ! Local constants
  dtdx = dt/dx
  hdtdx = HALF*dtdx
  hdtdy = HALF*dt/dy
  hdt = HALF*dt

  ! multidimensional shock detection -- this will be used to do the
  ! hybrid Riemann solver
  !if (hybrid_riemann == 1) then
  !   call shock(q,qd_l1,qd_l2,qd_h1,qd_h2, &
  !              shk,ilo1-1,ilo2-1,ihi1+1,ihi2+1, &
  !              ilo1,ilo2,ihi1,ihi2,dx,dy)
  !else
  shk(:,:) = ZERO
  !endif


  ! NOTE: Geometry terms need to be punched through
  ! Trace to edges w/o transverse flux correction terms
  if (ppm_type .eq. 0) then
     call bl_error("ppm_type <=0 is not supported in umeth2d_rad")
  else
     call trace_ppm_rad(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
                        q,qaux(:,:,QC),qaux(:,:,QCG),flatn,qd_l1,qd_l2,qd_h1,qd_h2, &
                        dloga,dloga_l1,dloga_l2,dloga_h1,dloga_h2, &
                        qxm,qxp,qym,qyp,ilo1-1,ilo2-1,ihi1+2,ihi2+2, &
                        srcQ,src_l1,src_l2,src_h1,src_h2, &
                        qaux(:,:,QGAMC),qaux(:,:,QGAMCG),qa_l1,qa_l2,qa_h1,qa_h2, &
                        ilo1,ilo2,ihi1,ihi2,dx,dy,dt)
  end if

  call cmpflx(qxm, qxp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
              fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
              qgdxtmp, q1_l1, q1_l2, q1_h1, q1_h2, &
              lam, lam_l1, lam_l2, lam_h1, lam_h2, &
              rfx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
              qaux(:,:,QGAMCG), qaux(:,:,QGAMC), qaux(:,:,QCSML), qaux(:,:,QC), &
              qa_l1, qa_l2, qa_h1, qa_h2, &
              shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
              1, ilo1, ihi1, ilo2-1, ihi2+1, domlo, domhi)

  call cmpflx(qym, qyp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
              fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
              q2, q2_l1, q2_l2, q2_h1, q2_h2, &
              lam,lam_l1,lam_l2,lam_h1,lam_h2, &
              rfy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
              qaux(:,:,QGAMCG), qaux(:,:,QGAMC), qaux(:,:,QCSML), qaux(:,:,QC), &
              qa_l1, qa_l2, qa_h1, qa_h2, &
              shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
              2, ilo1-1, ihi1+1, ilo2, ihi2, domlo, domhi)
  
  call transy(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
              qxm, qm, qxp, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
              fy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
              rfy, ilo1-1, ilo2, ihi1+1, ihi2+1, &
              q2, q2_l1, q2_l2, q2_h1, q2_h2, &
              qaux(:,:,QGAMCG), qa_l1, qa_l2, qa_h1, qa_h2, &
              srcQ, src_l1, src_l2, src_h1, src_h2, &
              hdt, hdtdy, &
              ilo1-1, ihi1+1, ilo2, ihi2)

  call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
              flux1,  fd1_l1,  fd1_l2,  fd1_h1,  fd1_h2, &
              q1, q1_l1, q1_l2, q1_h1, q1_h2, &
              lam,lam_l1,lam_l2,lam_h1,lam_h2, &
              rflux1, rfd1_l1, rfd1_l2, rfd1_h1, rfd1_h2, &
              qaux(:,:,QGAMCG), qaux(:,:,QGAMC), qaux(:,:,QCSML), qaux(:,:,QC), &
              qa_l1, qa_l2, qa_h1, qa_h2, &
              shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
              1, ilo1, ihi1, ilo2, ihi2, domlo, domhi)

  call transx(lam,lam_l1,lam_l2,lam_h1,lam_h2, &
              qym, qm,qyp,qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
              fx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
              rfx, ilo1, ilo2-1, ihi1+1, ihi2+1, &
              qgdxtmp, q1_l1, q1_l2, q1_h1, q1_h2, &
              qaux(:,:,QGAMCG), qa_l1, qa_l2, qa_h1, qa_h2, &
              srcQ,  src_l1,  src_l2,  src_h1,  src_h2, &
              hdt, hdtdx, &
              area1, area1_l1, area1_l2, area1_h1, area1_h2, &
              vol, vol_l1, vol_l2, vol_h1, vol_h2, &
              ilo1, ihi1, ilo2-1, ihi2+1)

  call cmpflx(qm, qp, ilo1-1, ilo2-1, ihi1+2, ihi2+2, &
              flux2,  fd2_l1,  fd2_l2,  fd2_h1,  fd2_h2, &
              q2, q2_l1, q2_l2, q2_h1, q2_h2, &
              lam,lam_l1,lam_l2,lam_h1,lam_h2, &
              rflux2, rfd2_l1, rfd2_l2, rfd2_h1, rfd2_h2, &
              qaux(:,:,QGAMCG), qaux(:,:,QGAMC), qaux(:,:,QCSML), qaux(:,:,QC), &
              qa_l1, qa_l2, qa_h1, qa_h2, &
              shk, ilo1-1, ilo2-1, ihi1+1, ihi2+1, &
              2, ilo1, ihi1, ilo2, ihi2, domlo, domhi)


  ! Construct p div{U} -- this will be used as a source to the internal
  ! energy update.  Note we construct this using the interface states
  ! returned from the Riemann solver.
    do j = ilo2,ihi2
       do i = ilo1,ihi1
          pdivu(i,j) = HALF*( &
               (q1(i+1,j,GDPRES) + q1(i,j,GDPRES)) * &
               (q1(i+1,j,GDU)*area1(i+1,j) - q1(i,j,GDU)*area1(i,j)) + &
               (q2(i,j+1,GDPRES) + q2(i,j,GDPRES)) * &
               (q2(i,j+1,GDV)*area2(i,j+1) - q2(i,j,GDV)*area2(i,j)) ) / vol(i,j)
       end do
    end do

  deallocate(qm,qp,qxm,qxp,qym,qyp)
  deallocate(fx,fy)
  deallocate(rfx,rfy)
  deallocate(qgdxtmp)
  deallocate(shk)

end subroutine umeth2d_rad

! :::
! ::: ------------------------------------------------------------------
! :::

subroutine consup_rad(uin, uin_l1, uin_l2, uin_h1, uin_h2, &
                      uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                      Erin, Erin_l1, Erin_l2, Erin_h1, Erin_h2, &
                      Erout,Erout_l1,Erout_l2,Erout_h1,Erout_h2, &
                      q1, q1_l1, q1_l2, q1_h1, q1_h2, &
                      q2, q2_l1, q2_l2, q2_h1, q2_h2, &
                      src , src_l1, src_l2, src_h1, src_h2, &
                      flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                      flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                      rflux1,rflux1_l1,rflux1_l2,rflux1_h1,rflux1_h2, &
                      rflux2,rflux2_l1,rflux2_l2,rflux2_h1,rflux2_h2, &
                      area1,area1_l1,area1_l2,area1_h1,area1_h2, &
                      area2,area2_l1,area2_l2,area2_h1,area2_h2, &
                      vol,vol_l1,vol_l2,vol_h1,vol_h2, &
                      div,pdivu, &
                      lo,hi,dx,dy,dt, nstep_fsp)

  use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, &
       ngdnv, GDPRES, GDU, GDV, GDERADS, GDLAMS
  use rad_params_module, only : ngroups, nugroup, dlognu
  use radhydro_params_module, only : fspace_type, comoving
  use radhydro_nd_module, only : advect_in_fspace
  use fluxlimiter_module, only : Edd_factor
  use advection_util_2d_module, only : normalize_species_fluxes
  use prob_params_module, only : coord_type

  implicit none

  integer nstep_fsp
  integer lo(2), hi(2)
  integer uin_l1,uin_l2,uin_h1,uin_h2
  integer uout_l1,uout_l2,uout_h1,uout_h2
  integer Erout_l1,Erout_l2,Erout_h1,Erout_h2
  integer Erin_l1,Erin_l2,Erin_h1,Erin_h2
  integer q1_l1, q1_l2, q1_h1, q1_h2
  integer q2_l1, q2_l2, q2_h1, q2_h2
  integer   src_l1,  src_l2,  src_h1,  src_h2
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
  double precision q1(q1_l1:q1_h1,q1_l2:q1_h2,ngdnv)
  double precision q2(q2_l1:q2_h1,q2_l2:q2_h2,ngdnv)
  double precision   src(  src_l1:  src_h1,  src_l2:  src_h2,NVAR)
  double precision  flux1( flux1_l1: flux1_h1, flux1_l2: flux1_h2,NVAR)
  double precision  flux2( flux2_l1: flux2_h1, flux2_l2: flux2_h2,NVAR)
  double precision rflux1(rflux1_l1:rflux1_h1,rflux1_l2:rflux1_h2,0:ngroups-1)
  double precision rflux2(rflux2_l1:rflux2_h1,rflux2_l2:rflux2_h2,0:ngroups-1)
  double precision area1(area1_l1:area1_h1,area1_l2:area1_h2)
  double precision area2(area2_l1:area2_h1,area2_l2:area2_h2)
  double precision vol(vol_l1:vol_h1,vol_l2:vol_h2)
  double precision div(lo(1):hi(1)+1,lo(2):hi(2)+1)
  double precision pdivu(lo(1):hi(1),lo(2):hi(2))
  double precision dx, dy, dt

  integer i, j, n, g

  double precision div1

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
  call normalize_species_fluxes(flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                                flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                                lo,hi)

  ! correct the fluxes to include the effects of the artificial viscosity
  do n = 1, NVAR
     if (n ==  UTEMP) then
        flux1(:,:,n) = ZERO
        flux2(:,:,n) = ZERO
     else
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)+1
              div1 = HALF*(div(i,j) + div(i,j+1))
              div1 = difmag*min(ZERO,div1)

              flux1(i,j,n) = flux1(i,j,n) &
                   + dx*div1*(uin(i,j,n) - uin(i-1,j,n))

              flux1(i,j,n) = area1(i,j)*flux1(i,j,n)
           enddo
        enddo

        do j = lo(2), hi(2)+1
           do i = lo(1), hi(1)
              div1 = HALF*(div(i,j) + div(i+1,j))
              div1 = difmag*min(ZERO,div1)

              flux2(i,j,n) = flux2(i,j,n) &
                   + dy*div1*(uin(i,j,n) - uin(i,j-1,n))

              flux2(i,j,n) = area2(i,j)*flux2(i,j,n)
           enddo
        enddo
     endif
  enddo

  do g = 0, ngroups-1
     do j = lo(2),hi(2)
        do i = lo(1), hi(1)+1
           div1 = HALF*(div(i,j) + div(i,j+1))
           div1 = difmag*min(ZERO,div1)

           rflux1(i,j,g) = rflux1(i,j,g) &
                + dx*div1*(Erin(i,j,g) - Erin(i-1,j,g))

           rflux1(i,j,g) = area1(i,j)*rflux1(i,j,g)
        enddo
     enddo
  enddo

  do g = 0, ngroups-1
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
           div1 = HALF*(div(i,j) + div(i+1,j))
           div1 = difmag*min(ZERO,div1)

           rflux2(i,j,g) = rflux2(i,j,g) &
                + dy*div1*(Erin(i,j,g) - Erin(i,j-1,g))

           rflux2(i,j,g) = area2(i,j)*rflux2(i,j,g)
        enddo
     enddo
  enddo

  ! do the conservative update
  do n = 1, NVAR
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           uout(i,j,n) = uout(i,j,n) + dt * ( flux1(i,j,n) - flux1(i+1,j,n) + &
                                              flux2(i,j,n) - flux2(i,j+1,n) ) / vol(i,j)
        enddo
     enddo
  enddo

  do g = 0, ngroups-1
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           Erout(i,j,g) = Erin(i,j,g) + dt * &
                ( rflux1(i,j,g) - rflux1(i+1,j,g) + &
                  rflux2(i,j,g) - rflux2(i,j+1,g) ) / vol(i,j)
        enddo
     enddo
  end do

  ! Add source term to (rho e)
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        uout(i,j,UEINT) = uout(i,j,UEINT) - dt * pdivu(i,j)
     enddo
  enddo

  ! Add gradp term to momentum equation -- only for axisymmetry coords
  ! (and only for the radial flux)
  do j = lo(2), hi(2)
     do i = lo(1), hi(1)

        ! pgdnv from the Riemann solver is only the gas contribution,
        ! not the radiation contribution.  Note that we've already included
        ! the gas pressure in the momentum flux for all Cartesian coordinate
        ! directions
        if (coord_type == 1) then
           dpdx = ( q1(i+1,j,GDPRES) - q1(i,j,GDPRES))/ dx
        else
           dpdx = ZERO
        endif

        dpdy = ZERO

        ! radiation pressure contribution
        dprdx = ZERO
        dprdy = ZERO
        do g= 0, ngroups-1
           lamc = 0.25d0*(q1(i,j,GDLAMS+g) + q1(i+1,j,GDLAMS+g) + &
                          q2(i,j,GDLAMS+g) + q2(i,j+1,GDLAMS+g))
           dprdx = dprdx + lamc*(q1(i+1,j,GDERADS+g) - q1(i,j,GDERADS+g))/dx
           dprdy = dprdy + lamc*(q2(i,j+1,GDERADS+g) - q2(i,j,GDERADS+g))/dy
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
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           ux = HALF*(q1(i,j,GDU) + q1(i+1,j,GDU))
           uy = HALF*(q2(i,j,GDV) + q2(i,j+1,GDV))

           divu = (q1(i+1,j,GDU)*area1(i+1,j) - q1(i,j,GDU)*area1(i,j) + &
                   q2(i,j+1,GDV)*area2(i,j+1) - q2(i,j,GDV)*area2(i,j)) / vol(i,j)

           dudx(1) = (q1(i+1,j,GDU) - q1(i,j,GDU))/dx
           dudx(2) = (q1(i+1,j,GDV) - q1(i,j,GDV))/dx

           dudy(1) = (q2(i,j+1,GDU) - q2(i,j,GDU))/dy
           dudy(2) = (q2(i,j+1,GDV) - q2(i,j,GDV))/dy

           ! Note that for single group, fspace_type is always 1
           do g=0, ngroups-1

              nhat(1) = (q1(i+1,j,GDERADS+g) - q1(i,j,GDERADS+g))/dx
              nhat(2) = (q2(i,j+1,GDERADS+g) - q2(i,j,GDERADS+g))/dy

              GnDotu(1) = dot_product(nhat, dudx)
              GnDotu(2) = dot_product(nhat, dudy)

              nnColonDotGu = dot_product(nhat, GnDotu) / (dot_product(nhat,nhat)+1.d-50)

              lamc = 0.25d0*(q1(i,j,GDLAMS+g) + q1(i+1,j,GDLAMS+g) + &
                             q2(i,j,GDLAMS+g) + q2(i,j+1,GDLAMS+g))
              Eddf = Edd_factor(lamc)
              f1 = (ONE-Eddf)*HALF
              f2 = (3.d0*Eddf-ONE)*HALF
              af(g) = -(f1*divu + f2*nnColonDotGu)

              if (fspace_type .eq. 1) then
                 Eddfxp = Edd_factor(q1(i+1,j  ,GDLAMS+g))
                 Eddfxm = Edd_factor(q1(i  ,j  ,GDLAMS+g))
                 Eddfyp = Edd_factor(q2(i  ,j+1,GDLAMS+g))
                 Eddfym = Edd_factor(q2(i  ,j  ,GDLAMS+g))

                 f1xp = HALF*(ONE-Eddfxp)
                 f1xm = HALF*(ONE-Eddfxm)
                 f1yp = HALF*(ONE-Eddfyp)
                 f1ym = HALF*(ONE-Eddfym)

                 Gf1E(1) = (f1xp*q1(i+1,j,GDERADS+g) - f1xm*q1(i,j,GDERADS+g)) / dx
                 Gf1E(2) = (f1yp*q2(i,j+1,GDERADS+g) - f1ym*q2(i,j,GDERADS+g)) / dy

                 Egdc = 0.25d0*(q1(i,j,GDERADS+g) + q1(i+1,j,GDERADS+g) + &
                                q2(i,j,GDERADS+g) + q2(i,j+1,GDERADS+g))

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


  ! scale the fluxes (and correct the momentum flux with the grad p part)
  ! so we can use them in the flux correction at coarse-fine interfaces
  ! later.

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)+1
        flux1(i,j,1:NVAR) = dt * flux1(i,j,1:NVAR)
     enddo
  enddo

  do g = 0, ngroups-1
     do j = lo(2), hi(2)
        do i = lo(1),hi(1)+1
           rflux1(i,j,g) = dt * rflux1(i,j,g)
        enddo
     enddo
  end do

  do j = lo(2), hi(2)+1
     do i = lo(1), hi(1)
        flux2(i,j,1:NVAR) = dt * flux2(i,j,1:NVAR)
     enddo
  enddo

  do g = 0, ngroups-1
     do j = lo(2), hi(2)+1
        do i = lo(1), hi(1)
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
              flatn(i,j) = ZERO
           end if
        end if
     end do
  end do
end subroutine ppflaten

end module rad_advection_module
