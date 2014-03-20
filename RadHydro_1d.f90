module rad_advection_module

  implicit none

  private

  public umeth1d_rad, ctoprim_rad, consup_rad

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine ctoprim_rad(lo,hi,uin,uin_l1,uin_h1, &
     Erin, Erin_l1, Erin_h1, &
     lam, lam_l1, lam_h1, &
     q,c,cg,gamc,gamcg,csml,flatn,q_l1,q_h1,&
     src,srcQ,src_l1,src_h1, &
     courno,dx,dt,ngp,ngf,iflaten)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
       QVAR, QRHO, QU, QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
       nadv, small_temp, allow_negative_energy
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, comoving, &
       flatten_pp_threshold, first_order_hydro
  use rad_params_module, only : ngroups
  use flatten_module, only : uflaten
  use fluxlimiter_module, only : Edd_factor

  implicit none

  double precision, parameter:: small = 1.d-8

!     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
!     if iflaten=1.  Declared dimensions of q,c,gamc,csml,flatn are given
!     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
!     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
!     routine that computes flatn).  

  integer          :: lo(1), hi(1)
  integer          :: uin_l1,uin_h1, Erin_l1, Erin_h1, lam_l1, lam_h1 
  integer          :: q_l1,q_h1
  integer          ::  src_l1,src_h1
  double precision ::   uin(uin_l1:uin_h1,NVAR)
  double precision :: Erin(Erin_l1:Erin_h1, 0:ngroups-1)
  double precision :: lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision ::     q(  q_l1:  q_h1,QRADVAR)
  double precision ::     c(  q_l1:  q_h1)
  double precision ::    cg(  q_l1:  q_h1)
  double precision :: gamcg(  q_l1:  q_h1)
  double precision ::  gamc(  q_l1:  q_h1)
  double precision ::  csml(  q_l1:  q_h1)
  double precision :: flatn(  q_l1:  q_h1)
  double precision ::   src(src_l1:src_h1,NVAR)
  double precision ::  srcQ(src_l1:src_h1,QVAR)
  double precision :: dx, dt, courno
  integer iflaten

  integer          :: i, g
  integer          :: ngp, ngf, loq(1), hiq(1)
  integer          :: n, nq
  integer          :: iadv, ispec, iaux
  double precision :: courx, courmx

  double precision, allocatable :: dpdrho(:), dpde(:), flatg(:)

  double precision :: csrad2, prad, Eddf, gamr

  type(eos_t) :: eos_state

  loq(1) = lo(1)-ngp
  hiq(1) = hi(1)+ngp

  allocate(dpdrho(q_l1:q_h1))
  allocate(dpde  (q_l1:q_h1))
  allocate(flatg (q_l1:q_h1))

!     Make q (all but p), except put e in slot for rho.e, fix after eos call
!     The temperature is used as an initial guess for the eos call and will be overwritten
  do i = loq(1),hiq(1)

     if (uin(i,URHO) .le. 0.d0) then
        print *,'   '
        print *,'>>> Error: Castro_1d::ctoprim ',i
        print *,'>>> ... negative density ',uin(i,URHO)
        print *,'    '
        call bl_error("Error:: Castro_1d.f90 :: ctoprim")
     end if

     q(i,QRHO) = uin(i,URHO)
     q(i,QU) = uin(i,UMX)/uin(i,URHO)
     q(i,QREINT) = uin(i,UEINT)/q(i,QRHO)
     q(i,QTEMP ) = uin(i,UTEMP)
     q(i,qrad:qradhi) = Erin(i,0:ngroups-1)
  enddo

!     Load advected quatities, c, into q, assuming they arrived in uin as rho.c
  do iadv = 1, nadv
     n  = UFA + iadv - 1
     nq = QFA + iadv - 1
     q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
  enddo
  
!     Load species, c, into q, assuming they arrived in uin as rho.c
  do ispec = 1, nspec
     n  = UFS + ispec - 1
     nq = QFS + ispec - 1
     q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
  enddo

!     Load auxiliary variables which are needed in the EOS
  do iaux = 1, naux
     n  = UFX + iaux - 1
     nq = QFX + iaux - 1
     q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
  enddo

!     Get gamc, p, T, c, csml using q state
  do i = loq(1), hiq(1)

     eos_state % rho = q(i,QRHO)
     eos_state % T   = q(i,QTEMP)
     eos_state % e   = q(i,QREINT)
     eos_state % xn  = q(i,QFS:QFS+nspec-1)
     eos_state % aux = q(i,QFX:QFX+naux-1)

     call eos(eos_input_re, eos_state)

     q(i,QTEMP) = eos_state % T
     q(i,QPRES) = eos_state % p
     dpdrho(i)  = eos_state % dpdr_e
     dpde(i)    = eos_state % dpde
     gamcg(i)   = eos_state % gam1
     cg(i)      = eos_state % cs

     csrad2 = 0.d0
     prad = 0.d0
     do g=0, ngroups-1
        if (comoving) then
           Eddf = Edd_factor(lam(i,g))
           gamr = (3.d0-Eddf)/2.d0
        else
           gamr = lam(i,g) + 1.d0
        end if
        prad = prad + (lam(i,g)*q(i,qrad+g))
        csrad2 = csrad2 + gamr * (lam(i,g)*q(i,qrad+g)) / q(i,QRHO)
     end do

     q(i,qptot) = q(i,QPRES) + prad
     c(i) = cg(i)**2 + csrad2
     gamc(i) = c(i) * q(i,QRHO) / q(i,qptot)
     c(i) = sqrt(c(i))
     csml(i) = max(small, small * c(i))
  end do

!     Make this "rho e" instead of "e"
  do i = loq(1),hiq(1)
     q(i,QREINT ) = q(i,QREINT )*q(i,QRHO) 
     q(i,qreitot) = q(i,QREINT) + sum(q(i,qrad:qradhi))
  enddo

  ! compute srcQ terms
  do i = lo(1)-1, hi(1)+1
     srcQ(i,QRHO   ) = src(i,URHO)
     srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
     srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + 0.5d0 * q(i,QU)**2 * srcQ(i,QRHO)
     srcQ(i,QPRES  ) = dpde(i) * (srcQ(i,QREINT) - q(i,QREINT)*srcQ(i,QRHO)/q(i,QRHO)) / q(i,QRHO) + &
              dpdrho(i) * srcQ(i,QRHO)! - &
!              sum(dpdX_er(i,:)*(src(i,UFS:UFS+nspec-1) - q(i,QFS:QFS+nspec-1)*srcQ(i,QRHO)))/q(i,QRHO)

     do ispec=1,nspec
        srcQ(i,QFS+ispec-1) = ( src(i,UFS+ispec-1) - q(i,QFS+ispec-1) * srcQ(i,QRHO) ) / q(i,QRHO)
     end do
     
     do iaux=1,naux
        srcQ(i,QFX+iaux-1) = ( src(i,UFX+iaux-1) - q(i,QFX+iaux-1) * srcQ(i,QRHO) ) / q(i,QRHO)
     end do
     
     do iadv=1,nadv
        srcQ(i,QFA+iadv-1) = ( src(i,UFA+iadv-1) - q(i,QFA+iadv-1) * srcQ(i,QRHO) ) / q(i,QRHO)
     end do

  end do

  !     Compute running max of Courant number over grids
  courmx = courno
  do i = lo(1),hi(1)

     courx  = ( c(i)+abs(q(i,QU)) ) * dt/dx
     courmx = max( courmx, courx )

     if (courx .gt. 1.d0) then
        print *,'   '
        call bl_warning("Warning:: RadHydro_1d.f90 :: CFL violation in ctoprim")
        print *,'>>> ... (u+c) * dt / dx > 1 ', courx
        print *,'>>> ... at cell (i)       : ',i
        print *,'>>> ... u, c                ',q(i,QU), c(i)
        print *,'>>> ... density             ',q(i,QRHO)
     end if
  enddo
  courno = courmx

  !     Compute flattening coef for slope calculations
  if (first_order_hydro) then
     flatn = 0.d0
  else if (iflaten.eq.1) then
     loq(1)=lo(1)-ngf
     hiq(1)=hi(1)+ngf
     call uflaten(loq,hiq, &
          q(:,qptot), &
          q(:,QU), &
          flatn,q_l1,q_h1)
     call uflaten(loq,hiq, &
          q(:,qpres), &
          q(:,QU), &
          flatg,q_l1,q_h1)
     flatn = flatn * flatg

     if (flatten_pp_threshold > 0.d0) then
        call ppflaten(loq,hiq, &
             flatn, q, q_l1, q_h1)
     end if
  else
     flatn = 1.d0
  endif

  deallocate(dpdrho,dpde,flatg)

end subroutine ctoprim_rad

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

subroutine umeth1d_rad(lo,hi,domlo,domhi, &
     lam, lam_l1, lam_h1, &       
     q,c,cg,gamc,gamcg,csml,flatn,qd_l1,qd_h1, &
     srcQ,src_l1,src_h1, &
     grav, gv_l1, gv_h1, &
     ilo,ihi,dx,dt, &
     flux ,   fd_l1,   fd_h1, &
     rflux,  rfd_l1,  rfd_h1, &
     pgdnv,pgdnv_l1,pgdnv_h1, &
     ergdnv,ergdnv_l1,ergdnv_h1, &
     lamgdnv,lamgdnv_l1,lamgdnv_h1, &
     ugdnv,ugdnv_l1,ugdnv_h1, &
     dloga,dloga_l1,dloga_h1)

  use meth_params_module, only : QVAR, NVAR, ppm_type
  use radhydro_params_module, only : QRADVAR
  use rad_params_module, only : ngroups

  implicit none
  
  integer lo(1),hi(1)
  integer domlo(1),domhi(1)
  integer dloga_l1,dloga_h1
  integer lam_l1,lam_h1
  integer qd_l1,qd_h1
  integer src_l1,src_h1
  integer fd_l1,fd_h1
  integer rfd_l1,rfd_h1
  integer pgdnv_l1,pgdnv_h1
  integer ergdnv_l1,ergdnv_h1
  integer lamgdnv_l1,lamgdnv_h1
  integer ugdnv_l1,ugdnv_h1
  integer gv_l1,gv_h1
  integer ilo,ihi
  double precision dx, dt
  double precision lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision     q(   qd_l1:qd_h1,QRADVAR)
  double precision  gamc(   qd_l1:qd_h1)
  double precision gamcg(   qd_l1:qd_h1)
  double precision flatn(   qd_l1:qd_h1)
  double precision  csml(   qd_l1:qd_h1)
  double precision     c(   qd_l1:qd_h1)
  double precision    cg(   qd_l1:qd_h1)
  double precision  flux(fd_l1   :fd_h1,NVAR)
  double precision rflux(rfd_l1:rfd_h1, 0:ngroups-1)
  double precision  srcQ(src_l1  :src_h1,NVAR)
  double precision  grav(gv_l1   :gv_h1)
  double precision  pgdnv(pgdnv_l1:pgdnv_h1)
  double precision ergdnv(ergdnv_l1:ergdnv_h1, 0:ngroups-1)
  double precision lamgdnv(lamgdnv_l1:lamgdnv_h1, 0:ngroups-1)
  double precision ugdnv(ugdnv_l1:ugdnv_h1)
  double precision dloga(dloga_l1:dloga_h1)
  
  integer i

!     Left and right state arrays (edge centered, cell centered)
  double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

!     Work arrays to hold 3 planes of riemann state and conservative fluxes
  allocate ( dq(ilo-1:ihi+1,QRADVAR))
  allocate ( qm(ilo-1:ihi+1,QRADVAR))
  allocate ( qp(ilo-1:ihi+1,QRADVAR))
  
!     Trace to edges w/o transverse flux correction terms
  if (ppm_type .gt. 0) then
     call trace_ppm_rad(lam, lam_l1, lam_h1, &       
          q,dq,c,cg,flatn,qd_l1,qd_h1, &
          dloga,dloga_l1,dloga_h1, &
          srcQ,src_l1,src_h1, &
          grav,gv_l1,gv_h1, &
          qm,qp,ilo-1,ihi+1, &
          ilo,ihi,domlo,domhi,dx,dt)
  else
     call bl_error("ppm_type <=0 is not supported in umeth1d_rad")
     ! call trace(q,dq,c,flatn,qd_l1,qd_h1, &
     !      dloga,dloga_l1,dloga_h1, &
     !      srcQ,src_l1,src_h1, &
     !      grav,gv_l1,gv_h1, &
     !      qm,qp,ilo-1,ihi+1, &
     !      ilo,ihi,domlo,domhi,dx,dt)
  end if
  
  !     Solve Riemann problem, compute xflux from improved predicted states 
  call cmpflx_rad(lo, hi, domlo, domhi, &
       lam, lam_l1, lam_h1, & 
       qm, qp, ilo-1,ihi+1, &
       flux ,  fd_l1, fd_h1, &
       rflux, rfd_l1,rfd_h1, &
       pgdnv,pgdnv_l1,pgdnv_h1, &
       ergdnv,ergdnv_l1,ergdnv_h1, &
       lamgdnv,lamgdnv_l1,lamgdnv_h1, &
       ugdnv,ugdnv_l1,ugdnv_h1, &
       gamc,gamcg,csml,c,qd_l1,qd_h1,ilo,ihi)

  deallocate (dq,qm,qp)

end subroutine umeth1d_rad


! ::: 
! ::: ------------------------------------------------------------------
! ::: 
subroutine trace_ppm_rad(lam, lam_l1, lam_h1, &   
     q,dq,c,cg,flatn,qd_l1,qd_h1, &
     dloga,dloga_l1,dloga_h1, &
     srcQ,src_l1,src_h1,&
     grav,gv_l1,gv_h1, &
     qxm,qxp,qpd_l1,qpd_h1, &
     ilo,ihi,domlo,domhi,dx,dt)

  use network, only : nspec, naux
  use meth_params_module, only : QVAR, QRHO, QU, QREINT, QPRES, QFA, QFS, QFX, & 
       nadv, small_dens, ppm_type, fix_mass_flux
  use prob_params_module, only : physbc_lo, physbc_hi, Outflow
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot
  use rad_params_module, only : ngroups
  use ppm_module, only : ppm

  implicit none

  integer ilo,ihi
  integer domlo(1),domhi(1)
  integer lam_l1, lam_h1
  integer    qd_l1,   qd_h1
  integer dloga_l1,dloga_h1
  integer   qpd_l1,  qpd_h1
  integer   src_l1,  src_h1
  integer    gv_l1,   gv_h1
  double precision dx, dt
  double precision lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision     q( qd_l1: qd_h1,QRADVAR)
  double precision  srcQ(src_l1:src_h1,QVAR)
  double precision flatn(qd_l1:qd_h1)
  double precision     c(qd_l1:qd_h1)
  double precision    cg(qd_l1:qd_h1)
  double precision dloga(dloga_l1:dloga_h1)

  double precision   dq( qpd_l1: qpd_h1,QRADVAR)
  double precision  qxm( qpd_l1: qpd_h1,QRADVAR)
  double precision  qxp( qpd_l1: qpd_h1,QRADVAR)
  double precision grav(  gv_l1:  gv_h1)

  !     Local variables
  integer i, g
  integer n, iadv
  integer ns, ispec, iaux
  
  double precision hdt,dtdx

  double precision, dimension(0:ngroups-1) :: er,der,alphar,sourcer,qrtmp,hr
  double precision, dimension(0:ngroups-1) :: lam0, lamp, lamm
  double precision cc, csq, rho, u, p, ptot, rhoe, enth, cgassq
  double precision dum, dptotm
  double precision drho, drhoe, dptot
  double precision dup, dptotp
  
  double precision alpham, alphap, alpha0, alphae 
  double precision sourcr,sourcp,source,courn,eta,dlogatmp

  double precision rhoe_g, h_g, alphae_g, drhoe_g
  
  logical :: fix_mass_flux_lo, fix_mass_flux_hi
  
  double precision, allocatable :: Ip(:,:,:)
  double precision, allocatable :: Im(:,:,:)

  double precision :: er_foo

  fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) &
       .and. (ilo .eq. domlo(1))
  fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) &
       .and. (ihi .eq. domhi(1))

  if (ppm_type .eq. 0) then
     print *,'Oops -- shouldnt be in trace_ppm with ppm_type = 0'
     call bl_error("Error:: RadHydro_1d.f90 :: trace_ppm_rad")
  end if

  hdt = 0.5d0 * dt
  dtdx = dt/dx

  allocate(Ip(ilo-1:ihi+1,3,QRADVAR))
  allocate(Im(ilo-1:ihi+1,3,QRADVAR))

  ! Compute Ip and Im
  do n=1,QRADVAR
     call ppm(q(:,n),qd_l1,qd_h1,q(:,QU),c, flatn, &
          Ip(:,:,n),Im(:,:,n),ilo,ihi,dx,dt)
  end do

  ! Trace to left and right edges using upwind PPM
  do i = ilo-1, ihi+1

     do g=0, ngroups-1
        lam0(g) = lam(i,g)
        lamp(g) = lam(i,g)
        lamm(g) = lam(i,g)
     end do

     cgassq = cg(i)**2
     cc = c(i)
     csq = cc**2

     rho = q(i,QRHO)
     u = q(i,QU)
     p = q(i,QPRES)
     rhoe_g = q(i,QREINT)
     h_g = (p+rhoe_g) / rho
     er(:) = q(i,qrad:qradhi)
     hr(:) = (lam0+1.d0)*er/rho
     ptot = q(i,qptot)
     rhoe = q(i,qreitot)
     enth = ( (rhoe+ptot)/rho )/csq

     ! plus state on face i
     dum    = flatn(i)*(u    - Im(i,1,QU))
     dptotm = flatn(i)*(ptot - Im(i,1,qptot))
             
     drho  = flatn(i)*(rho  - Im(i,2,QRHO))
     drhoe_g = flatn(i)*(rhoe_g  - Im(i,2,QREINT))
     drhoe = flatn(i)*(rhoe - Im(i,2,qreitot))
     dptot = flatn(i)*(ptot - Im(i,2,qptot))
     der(:)= flatn(i)*(er(:)- Im(i,2,qrad:qradhi))
             
     dup    = flatn(i)*(u    - Im(i,3,QU))
     dptotp = flatn(i)*(ptot - Im(i,3,qptot))

     alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
     alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
     alpha0 = drho - dptot/csq
     alphae = drhoe - dptot*enth
     alphae_g = drhoe_g - dptot/csq*h_g
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
        alphae = 0.d0
        alphae_g = 0.d0
        alphar(:) = 0.d0
     else if (u .lt. 0.d0) then
        alpha0 = -alpha0
        alphae = -alphae
        alphae_g = -alphae_g
        alphar(:) = -alphar(:)
     else
        alpha0 = -0.5d0*alpha0
        alphae = -0.5d0*alphae
        alphae_g = -0.5d0*alphae_g
        alphar(:) = -0.5d0*alphar(:)
     endif

     if (i .ge. ilo) then
        qxp(i,QRHO) = rho + alphap + alpham + alpha0 
        qxp(i,QRHO) = max(small_dens,qxp(i,QRHO))
        qxp(i,QU) = u + (alphap - alpham)*cc/rho
        qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
        qxp(i,qrad:qradhi) = qrtmp
        qxp(i,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
        qxp(i,QPRES) = p + (alphap + alpham)*cgassq - sum(lamp(:)*alphar(:))
        qxp(i,qptot) = ptot + (alphap + alpham)*csq 
        qxp(i,qreitot) = qxp(i,QREINT) + sum(qrtmp)

        ! add non-gravitational source term
        qxp(i  ,QRHO  )  = qxp(i,QRHO   ) + hdt*srcQ(i,QRHO)
        qxp(i  ,QRHO  )  = max(small_dens,qxp(i,QRHO))
        qxp(i  ,QU    )  = qxp(i,QU     ) + hdt*srcQ(i,QU)
        qxp(i  ,QREINT)  = qxp(i,QREINT ) + hdt*srcQ(i,QREINT)
        qxp(i  ,QPRES )  = qxp(i,QPRES  ) + hdt*srcQ(i,QPRES)
        qxp(i  ,qptot )  = qxp(i,qptot  ) + hdt*srcQ(i,QPRES)
        qxp(i  ,qreitot) = qxp(i,qreitot) + hdt*srcQ(i,QREINT)

        ! add gravitational source term
        qxp(i  ,QU) = qxp(i,QU) + hdt*grav(i)

        do g=0, ngroups-1
           if (qxp(i,qrad+g) < 0.d0) then
              er_foo = - qxp(i,qrad+g)
              qxp(i,qrad+g) = 0.d0
              qxp(i,qptot) = qxp(i,qptot) + lamp(g) * er_foo
              qxp(i,qreitot) = qxp(i,qreitot) + er_foo
           end if
        end do

        if ( qxp(i,QPRES) < 0.d0 ) then
           qxp(i,QPRES) = p
        end if
     end if

     ! minus state on face i+1
     dum    = flatn(i)*(u    - Ip(i,1,QU))
     dptotm = flatn(i)*(ptot - Ip(i,1,qptot))

     drho  = flatn(i)*(rho  - Ip(i,2,QRHO))
     drhoe_g = flatn(i)*(rhoe_g - Ip(i,2,QREINT))
     drhoe = flatn(i)*(rhoe - Ip(i,2,qreitot))
     dptot = flatn(i)*(ptot - Ip(i,2,qptot))
     der(:)= flatn(i)*(er(:)- Ip(i,2,qrad:qradhi))

     dup    = flatn(i)*(u    - Ip(i,3,QU))
     dptotp = flatn(i)*(ptot - Ip(i,3,qptot))

     alpham = 0.5d0*(dptotm/(rho*cc) - dum)*rho/cc
     alphap = 0.5d0*(dptotp/(rho*cc) + dup)*rho/cc
     alpha0 = drho - dptot/csq
     alphae = drhoe - dptot*enth
     alphae_g = drhoe_g - dptot/csq*h_g
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
        alphae = -alphae
        alphae_g = -alphae_g
        alphar(:) = -alphar(:)
     else if (u .lt. 0.d0) then
        alpha0 = 0.d0
        alphae = 0.d0
        alphae_g = 0.d0
        alphar(:) = 0.d0
     else
        alpha0 = -0.5d0*alpha0
        alphae = -0.5d0*alphae
        alphae_g = -0.5d0*alphae_g
        alphar(:) = -0.5d0*alphar(:)
     endif

     if (i .le. ihi) then
        qxm(i+1,QRHO) = rho + alphap + alpham + alpha0 
        qxm(i+1,QRHO) = max(qxm(i+1,QRHO),small_dens)
        qxm(i+1,QU) = u + (alphap - alpham)*cc/rho  
        qrtmp = er(:) + (alphap + alpham)*hr + alphar(:)
        qxm(i+1,qrad:qradhi) = qrtmp
        qxm(i+1,QREINT) = rhoe_g + (alphap + alpham)*h_g + alphae_g
        qxm(i+1,QPRES) = p + (alphap + alpham)*cgassq - sum(lamm(:)*alphar(:))
        qxm(i+1,qptot) = ptot + (alphap + alpham)*csq
        qxm(i+1,qreitot) = qxm(i+1,QREINT) + sum(qrtmp)

        ! add non-gravitational source term
        qxm(i+1,QRHO   ) = qxm(i+1,QRHO   ) + hdt*srcQ(i,QRHO)
        qxm(i+1,QRHO   ) = max(small_dens, qxm(i+1,QRHO))
        qxm(i+1,QU     ) = qxm(i+1,QU     ) + hdt*srcQ(i,QU)
        qxm(i+1,QREINT ) = qxm(i+1,QREINT ) + hdt*srcQ(i,QREINT)
        qxm(i+1,QPRES  ) = qxm(i+1,QPRES  ) + hdt*srcQ(i,QPRES)
        qxm(i+1,qptot  ) = qxm(i+1,qptot  ) + hdt*srcQ(i,QPRES)
        qxm(i+1,qreitot) = qxm(i+1,qreitot) + hdt*srcQ(i,QREINT)

        ! add gravitational source term
        qxm(i+1,QU) = qxm(i+1,QU) + hdt*grav(i)

        do g=0, ngroups-1
           if (qxm(i+1,qrad+g) < 0.d0) then
              er_foo = - qxm(i+1,qrad+g)
              qxm(i+1,qrad+g) = 0.d0
              qxm(i+1,qptot) = qxm(i+1,qptot) + lamm(g) * er_foo
              qxm(i+1,qreitot) = qxm(i+1,qreitot) + er_foo
           end if
        end do

        if ( qxm(i+1,QPRES) < 0.d0 ) then
           qxm(i+1,QPRES) = p
        end if
     end if

     if(dloga(i).ne.0)then
        courn = dtdx*(cc+abs(u))
        eta = (1.d0-courn)/(cc*dt*abs(dloga(i)))
        dlogatmp = min(eta,1.d0)*dloga(i)
        sourcr = -0.5d0*dt*rho*dlogatmp*u
        sourcp = sourcr*cgassq
        source = sourcr*h_g
        sourcer(:) = -0.5d0*dt*dlogatmp*u*(lam0(:)+1.d0)*er(:)
        if (i .le. ihi) then
           qxm(i+1,QRHO  ) = qxm(i+1,QRHO  ) + sourcr
           qxm(i+1,QRHO  ) = max(small_dens, qxm(i+1,QRHO))
           qxm(i+1,QPRES ) = qxm(i+1,QPRES ) + sourcp
           qxm(i+1,QREINT) = qxm(i+1,QREINT) + source
           qxm(i+1,qrad:qradhi) = qxm(i+1,qrad:qradhi) + sourcer(:)
!           qxm(i+1,qptot ) = sum(lamm(:)*qxm(i+1,qrad:qradhi)) + qxm(i+1,QPRES)
           qxm(i+1,qptot) = qxm(i+1,qptot) + sum(lamm(:)*sourcer(:)) + sourcp
           qxm(i+1,qreitot) = sum(qxm(i+1,qrad:qradhi))  + qxm(i+1,QREINT)
        end if
        if (i .ge. ilo) then
           qxp(i  ,QRHO  ) = qxp(i  ,QRHO  ) + sourcr
           qxp(i  ,QRHO  ) = max(small_dens,qxp(i,QRHO))
           qxp(i  ,QPRES ) = qxp(i  ,QPRES ) + sourcp
           qxp(i  ,QREINT) = qxp(i  ,QREINT) + source
           qxp(i  ,qrad:qradhi) = qxp(i  ,qrad:qradhi) + sourcer(:)
!           qxp(i  ,qptot ) = sum(lamp(:)*qxp(i,qrad:qradhi)) + qxp(i,QPRES)
           qxp(i,qptot) = qxp(i,qptot) + sum(lamp(:)*sourcer(:)) + sourcp
           qxp(i  ,qreitot) = sum(qxp(i,qrad:qradhi))  + qxp(i,QREINT)
        end if
     endif

  end do

  ! Enforce constant mass flux rate if specified
  if (fix_mass_flux_lo) then
     qxm(ilo,QRHO   ) = q(domlo(1)-1,QRHO)
     qxm(ilo,QU     ) = q(domlo(1)-1,QU  )
     qxm(ilo,QPRES  ) = q(domlo(1)-1,QPRES)
     qxm(ilo,QREINT ) = q(domlo(1)-1,QREINT)
     qxm(ilo,qrad:qradhi) = q(domlo(1)-1,qrad:qradhi)
     qxm(ilo,qptot  ) = q(domlo(1)-1,qptot)
     qxm(ilo,qreitot) = q(domlo(1)-1,qreitot) 
  end if
   
  ! Enforce constant mass flux rate if specified
  if (fix_mass_flux_hi) then
     qxp(ihi+1,QRHO   ) = q(domhi(1)+1,QRHO)
     qxp(ihi+1,QU     ) = q(domhi(1)+1,QU  )
     qxp(ihi+1,QPRES  ) = q(domhi(1)+1,QPRES)
     qxp(ihi+1,QREINT ) = q(domhi(1)+1,QREINT)
     qxp(ihi+1,qrad:qradhi) = q(domhi(1)+1,qrad:qradhi)
     qxp(ihi+1,qptot  ) = q(domhi(1)+1,qptot)
     qxp(ihi+1,qreitot) = q(domhi(1)+1,qreitot)
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

end subroutine trace_ppm_rad

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine cmpflx_rad(lo,hi,domlo,domhi, &
     lam, lam_l1, lam_h1, &       
     qm,qp,qpd_l1,qpd_h1, &
     flx,flx_l1,flx_h1, &
     rflx,rflx_l1,rflx_h1, &
     pgdnv,pg_l1,pg_h1, &
     ergdnv,erg_l1,erg_h1, &
     lamgdnv,lg_l1,lg_h1, &
     ugdnv,ug_l1,ug_h1, &
     gamc,gamcg,csml,c,qd_l1,qd_h1,ilo,ihi)

  use meth_params_module, only : QVAR, NVAR
  use radhydro_params_module, only : QRADVAR
  use rad_params_module, only : ngroups

  implicit none

  integer lo(1),hi(1)
  integer domlo(1),domhi(1)
  integer ilo,ihi
  integer lam_l1, lam_h1
  integer qpd_l1,qpd_h1
  integer flx_l1, flx_h1
  integer rflx_l1, rflx_h1
  integer  pg_l1, pg_h1
  integer erg_l1,erg_h1
  integer  lg_l1, lg_h1
  integer  ug_l1, ug_h1
  integer  qd_l1,  qd_h1
  double precision lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision    qm(qpd_l1:qpd_h1, QRADVAR)
  double precision    qp(qpd_l1:qpd_h1, QRADVAR)
  double precision   flx(flx_l1:flx_h1, NVAR)
  double precision rflx(rflx_l1:rflx_h1, 0:ngroups-1)
  double precision  pgdnv( pg_l1: pg_h1)
  double precision ergdnv(erg_l1:erg_h1, 0:ngroups-1)
  double precision lamgdnv( lg_l1: lg_h1, 0:ngroups-1)
  double precision ugdnv( ug_l1: ug_h1)
  double precision  gamc( qd_l1: qd_h1)
  double precision gamcg( qd_l1: qd_h1)
  double precision     c( qd_l1: qd_h1)
  double precision  csml( qd_l1: qd_h1)

!     Local variables
  integer i
  double precision, allocatable :: smallc(:),cavg(:),gamcp(:), gamcm(:), gamcgp(:), gamcgm(:)
      
  allocate ( smallc(ilo:ihi+1) )
  allocate ( cavg(ilo:ihi+1) )
  allocate ( gamcp(ilo:ihi+1) )
  allocate ( gamcm(ilo:ihi+1) )
  allocate (gamcgp(ilo:ihi+1) )
  allocate (gamcgm(ilo:ihi+1) )
  
  do i = ilo, ihi+1 
     smallc(i) = max( csml(i), csml(i-1) )
     cavg(i) = 0.5d0*( c(i) + c(i-1) )
     gamcgm(i) = gamcg(i-1)
     gamcgp(i) = gamcg(i)
     gamcm (i) = gamc (i-1)
     gamcp (i) = gamc (i)
  enddo
  
  !     Solve Riemann problem (gdnv state passed back, but only (u,p,Er) saved)
  call riemannus_rad(lam, lam_l1, lam_h1, &       
       qm, qp,qpd_l1,qpd_h1, smallc, cavg, &
       gamcm, gamcp, gamcgm, gamcgp, &
       flx, flx_l1, flx_h1, &
       rflx, rflx_l1, rflx_h1, &
       pgdnv, pg_l1, pg_h1, &
       ergdnv,erg_l1,erg_h1, &
       lamgdnv, lg_l1, lg_h1, &
       ugdnv, ug_l1, ug_h1, ilo, ihi, domlo, domhi )
  
  deallocate (smallc,cavg,gamcp,gamcm,gamcgp,gamcgm)
  
end subroutine cmpflx_rad


! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine riemannus_rad(lam, lam_l1, lam_h1, & 
     ql,qr,qpd_l1,qpd_h1,smallc,cav, &
     gamcl,gamcr,gamcgl,gamcgr,uflx,uflx_l1,uflx_h1,&
     rflx,rflx_l1,rflx_h1, &
     pgdnv,pg_l1,pg_h1,  &
     ergdnv,erg_l1,erg_h1,  &
     lamgdnv,lg_l1,lg_h1,  &
     ugdnv,ug_l1,ug_h1, &
     ilo,ihi,domlo,domhi)

  use network, only : nspec, naux
  use meth_params_module, only : QVAR, NVAR, QRHO, QU, QPRES, QREINT, QFA, QFS, QFX, &
       URHO, UMX, UEDEN, UEINT, UFA, UFS, UFX, nadv, small_dens, small_pres, &
       fix_mass_flux
  use prob_params_module, only : physbc_lo, physbc_hi, Outflow, Symmetry
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, fspace_type
  use rad_params_module, only : ngroups
  use fluxlimiter_module, only : Edd_factor

  implicit none

  double precision, parameter:: small = 1.d-8

  integer ilo,ihi
  integer domlo(1),domhi(1)
  integer lam_l1, lam_h1
  integer  qpd_l1,  qpd_h1
  integer   pg_l1,   pg_h1
  integer  erg_l1,  erg_h1
  integer   lg_l1,   lg_h1
  integer   ug_l1,   ug_h1
  integer uflx_l1, uflx_h1
  integer rflx_l1, rflx_h1
  double precision lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision ql(qpd_l1:qpd_h1, QRADVAR)
  double precision qr(qpd_l1:qpd_h1, QRADVAR)
  double precision    cav(ilo:ihi+1),smallc(ilo:ihi+1)
  double precision  gamcl(ilo:ihi+1), gamcr(ilo:ihi+1)
  double precision gamcgl(ilo:ihi+1),gamcgr(ilo:ihi+1)
  double precision  uflx(uflx_l1:uflx_h1, NVAR)
  double precision  rflx(rflx_l1:rflx_h1, 0:ngroups-1)
  double precision  pgdnv( pg_l1: pg_h1)
  double precision ergdnv(erg_l1:erg_h1, 0:ngroups-1)
  double precision lamgdnv( lg_l1: lg_h1, 0:ngroups-1)
  double precision ugdnv( ug_l1: ug_h1)

  double precision rgdnv, ustar
  double precision, dimension(0:ngroups-1) :: erl, err
  double precision rl, ul, pl, rel, pl_g, rel_g, wl
  double precision rr, ur, pr, rer, pr_g, rer_g, wr
  double precision rstar, cstar, pstar
  double precision ro, uo, po, reo, co, gamco, reo_g, po_g, co_g, gamco_g
  double precision sgnm, spin, spout, ushock, frac
  double precision rhoetot, scr
  
  double precision wsmall, csmall
  integer iadv, n, nq
  integer k,ispec, iaux
  logical :: fix_mass_flux_lo, fix_mass_flux_hi

  double precision :: regdnv_g, pgdnv_g, pgdnv_t
  double precision :: drho, estar_g, pstar_g
  double precision, dimension(0:ngroups-1) :: lambda, reo_r, po_r, estar_r, regdnv_r
  double precision :: eddf, f1
  integer :: g

  !     Solve Riemann Problem

  fix_mass_flux_lo = (fix_mass_flux .eq. 1) .and. (physbc_lo(1) .eq. Outflow) .and. (ilo .eq. domlo(1))
  fix_mass_flux_hi = (fix_mass_flux .eq. 1) .and. (physbc_hi(1) .eq. Outflow) .and. (ihi .eq. domhi(1))

  do k = ilo, ihi+1

     rl  = ql(k,QRHO)
     ul  = ql(k,QU)
     pl  = ql(k,qptot)
     rel = ql(k,qreitot)
     erl(:) = ql(k,qrad:qradhi)
     pl_g = ql(k,QPRES)
     rel_g = ql(k,QREINT)

     rr  = qr(k,QRHO)
     ur  = qr(k,QU)
     pr  = qr(k,qptot)
     rer = qr(k,qreitot)
     err(:) = qr(k,qrad:qradhi)
     pr_g = qr(k,QPRES)
     rer_g = qr(k,QREINT) 

     csmall = smallc(k)
     wsmall = small_dens*csmall
     wl = max(wsmall,sqrt(abs(gamcl(k)*pl*rl)))
     wr = max(wsmall,sqrt(abs(gamcr(k)*pr*rr)))
     pstar = ((wr*pl + wl*pr) + wl*wr*(ul - ur))/(wl + wr)
     pstar = max(pstar,small_pres)
     ustar = ((wl*ul + wr*ur) + (pl - pr))/(wl + wr)

     if (ustar .gt. 0.d0) then
        lambda(:) = lam(k-1,:)
        ro = rl
        uo = ul
        po = pl
        po_g = pl_g
        po_r(:) = erl(:) * lambda(:)
        reo = rel
        reo_r(:) = erl(:)
        reo_g = rel_g
        gamco = gamcl(k)
        gamco_g = gamcgl(k)
     else if (ustar .lt. 0.d0) then
        lambda(:) = lam(k,:)
        ro = rr
        uo = ur
        po = pr
        po_g = pr_g
        po_r(:) = err(:) * lambda(:)
        reo = rer
        reo_r(:) = err(:)
        reo_g = rer_g
        gamco = gamcr(k)
        gamco_g = gamcgr(k)
     else

        do g=0, ngroups-1
           lambda(g) = 0.5d0*(lam(k-1,g)+lam(k,g))
        end do

        ro = 0.5d0*(rl+rr)
        uo = 0.5d0*(ul+ur)
        reo = 0.5d0*(rel+rer)
        reo_r(:) = 0.5d0*(erl(:)+err(:))
        reo_g = 0.5d0*(rel_g+rer_g)
        po = 0.5d0*(pl+pr)
        po_r(:) = lambda(:) * reo_r(:)
        po_g = 0.5*(pr_g+pl_g)
        gamco = 0.5d0*(gamcl(k)+gamcr(k))
        gamco_g = 0.5d0*(gamcgl(k)+gamcgr(k))
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
        scr = small*cav(k)  
     else
        scr = spout-spin
     endif
     frac = (1.d0 + (spout + spin)/scr)*0.5d0
     frac = max(0.d0,min(1.d0,frac))

     rgdnv = frac*rstar + (1.d0 - frac)*ro
     ugdnv(k) = frac*ustar + (1.d0 - frac)*uo
     pgdnv_t = frac*pstar + (1.d0 - frac)*po
     pgdnv_g = frac*pstar_g + (1.d0 - frac)*po_g
     regdnv_g = frac*estar_g + (1.d0 - frac)*reo_g
     regdnv_r(:) = frac*estar_r(:) + (1.d0 - frac)*reo_r(:)

     if (spout .lt. 0.d0) then
        rgdnv = ro
        ugdnv(k) = uo
        pgdnv_t = po
        pgdnv_g = po_g
        regdnv_g = reo_g
        regdnv_r(:) = reo_r(:)
     endif
     if (spin .ge. 0.d0) then
        rgdnv = rstar
        ugdnv(k) = ustar
        pgdnv_t = pstar
        pgdnv_g = pstar_g
        regdnv_g = estar_g
        regdnv_r(:) = estar_r(:)
     endif

     if (k.eq.0 .and. physbc_lo(1) .eq. Symmetry) ugdnv(k) = 0.d0

     if (fix_mass_flux_lo .and. k.eq.domlo(1) .and. ugdnv(k) .ge. 0.d0) then
        rgdnv    = ql(k,QRHO)
        ugdnv(k) = ql(k,QU)
        regdnv_g = rel_g
        regdnv_r(:) = erl(:)
     end if
     if (fix_mass_flux_hi .and. k.eq.domhi(1)+1 .and. ugdnv(k) .le. 0.d0) then
        rgdnv    = qr(k,QRHO)
        ugdnv(k) = qr(k,QU)
        regdnv_g = rer_g
        regdnv_r(:) = err(:)
     end if

     do g=0, ngroups-1
        ergdnv(k,g) = max(regdnv_r(g), 0.d0)
     end do

     pgdnv(k) = pgdnv_g

     lamgdnv(k,:) = lambda
     
     ! Compute fluxes, order as conserved state (not q)
     uflx(k,URHO) = rgdnv*ugdnv(k)
     uflx(k,UMX) = uflx(k,URHO)*ugdnv(k) 

     rhoetot = regdnv_g + 0.5d0*rgdnv*ugdnv(k)**2 
     uflx(k,UEDEN) = ugdnv(k)*(rhoetot + pgdnv_g)
     uflx(k,UEINT) = ugdnv(k)*regdnv_g

     if (fspace_type.eq.1) then
        do g=0,ngroups-1
           eddf = Edd_factor(lambda(g))
           f1 = 0.5d0*(1.d0-eddf)
           rflx(k,g) = (1.d0+f1) * ergdnv(k,g) * ugdnv(k)
        end do
     else ! type 2
        do g=0,ngroups-1
           rflx(k,g) = ergdnv(k,g) * ugdnv(k)
        end do
     end if

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
        n  = UFS + ispec - 1
        nq = QFS + ispec - 1
        if (ustar .ge. 0.d0) then
           uflx(k,n) = uflx(k,URHO)*ql(k,nq)
        else
           uflx(k,n) = uflx(k,URHO)*qr(k,nq)
        endif
     enddo

     do iaux = 1, naux
        n  = UFX + iaux - 1
        nq = QFX + iaux - 1
        if (ustar .ge. 0.d0) then
           uflx(k,n) = uflx(k,URHO)*ql(k,nq)
        else
           uflx(k,n) = uflx(k,URHO)*qr(k,nq)
        endif
     enddo
     
  enddo
end subroutine riemannus_rad


! ::: 
! ::: ------------------------------------------------------------------
! ::: 
     
subroutine consup_rad(uin,  uin_l1,  uin_h1, &
     uout, uout_l1 ,uout_h1, &
     Erin,Erin_l1,Erin_h1, &
     Erout,Erout_l1,Erout_h1, &
     pgdnv,pgdnv_l1,pgdnv_h1, &
     ergdnv,ergdnv_l1,ergdnv_h1, &
     lamgdnv,lamgdnv_l1,lamgdnv_h1, &
     ugdnv,ugdnv_l1,ugdnv_h1, &
     src,  src_l1,  src_h1, &
     grav, grav_l1, grav_h1, &
     flux, flux_l1, flux_h1, &
     rflux,rflux_l1,rflux_h1, &
     flat, flat_l1, flat_h1, &
     area,area_l1,area_h1, &
     vol,vol_l1,vol_h1, &
     div,pdivu,lo,hi,dx,dt, &
     nstep_fsp)

  use meth_params_module, only : difmag, NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
       normalize_species
  use rad_params_module, only : ngroups, nugroup, dlognu
  use radhydro_params_module, only : fspace_type, comoving, QRADVAR, QPTOT
  use radhydro_nd_module, only : advect_in_fspace
  use fluxlimiter_module, only : Edd_factor
  use advection_module, only : normalize_species_fluxes

  implicit none
  integer nstep_fsp
  integer lo(1), hi(1)
  integer   uin_l1,  uin_h1
  integer  uout_l1, uout_h1
  integer  Erin_l1, Erin_h1
  integer Erout_l1,Erout_h1
  integer pgdnv_l1,pgdnv_h1
  integer ugdnv_l1,ugdnv_h1
  integer ergdnv_l1,ergdnv_h1
  integer lamgdnv_l1,lamgdnv_h1
  integer   src_l1,  src_h1
  integer  grav_l1, grav_h1
  integer  flux_l1, flux_h1
  integer rflux_l1,rflux_h1
  integer  flat_l1, flat_h1
  integer  area_l1, area_h1
  integer   vol_l1,  vol_h1
  double precision   uin(uin_l1:uin_h1,NVAR)
  double precision  uout(uout_l1:uout_h1,NVAR)
  double precision  Erin( Erin_l1: Erin_h1, 0:ngroups-1)
  double precision Erout(Erout_l1:Erout_h1, 0:ngroups-1)
  double precision  pgdnv(pgdnv_l1:pgdnv_h1)
  double precision  ugdnv(ugdnv_l1:ugdnv_h1)
  double precision ergdnv(ergdnv_l1:ergdnv_h1, 0:ngroups-1)
  double precision lamgdnv(lamgdnv_l1:lamgdnv_h1, 0:ngroups-1)
  double precision   src(  src_l1:  src_h1,NVAR)
  double precision  grav( grav_l1: grav_h1)
  double precision  flux( flux_l1: flux_h1,NVAR)
  double precision rflux(rflux_l1:rflux_h1, 0:ngroups-1)
  double precision  flat( flat_l1: flat_h1)
  double precision  area( area_l1: area_h1)
  double precision    vol(vol_l1:vol_h1)
  double precision    div(lo(1):hi(1)+1)
  double precision  pdivu(lo(1):hi(1)  )
  double precision dx, dt

  integer          :: i, n, g, ng2, g0, g1
  double precision :: div1
  double precision :: SrU,Up,SrE

  double precision, dimension(0:ngroups-1) :: Erscale
  double precision, dimension(0:ngroups-1) :: ustar, af
  double precision :: Eddf, Eddflft, Eddfrgt, f1, f2, f1lft, f1rgt
  double precision :: ux, divu, dudx, Egdc, lamc
  double precision :: dpdx, dprdx, ek1, ek2, dek
  
  if (ngroups .gt. 1) then
     if (fspace_type .eq. 1) then
        Erscale = dlognu
     else
        Erscale = nugroup*dlognu
     end if
  end if

  ! Normalize the species fluxes
  if (normalize_species .eq. 1) &
       call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

  do n = 1, NVAR
     if ( n.eq.UTEMP ) then
        flux(:,n) = 0.d0
     else
        do i = lo(1),hi(1)+1
           div1 = difmag*min(0.d0,div(i))
           flux(i,n) = flux(i,n) &
                + dx*div1*(uin(i,n) - uin(i-1,n))
           flux(i,n) = area(i) * flux(i,n) * dt
        enddo
     endif
  enddo

  do g=0, ngroups-1
     do i = lo(1),hi(1)+1
        div1 = difmag*min(0.d0,div(i))
        rflux(i,g) = rflux(i,g) &
             + dx*div1*(Erin(i,g) - Erin(i-1,g))
        rflux(i,g) = area(i) * rflux(i,g) * dt
     enddo
  end do

  do n = 1, NVAR
     if ( n.eq.UTEMP) then
        do i = lo(1),hi(1)
           uout(i,n) = uin(i,n)
        enddo
     else
        do i = lo(1),hi(1)
           uout(i,n) = uin(i,n) &
                + ( flux(i,n) - flux(i+1,n) ) / vol(i) &
                + dt * src(i,n)
        enddo
     end if
  enddo

  do g=0, ngroups-1
     do i = lo(1),hi(1)
        Erout(i,g) = Erin(i,g) + (rflux(i,g) - rflux(i+1,g) ) / vol(i) 
     enddo
  end do

  ! Add source term to (rho e)
  do i = lo(1),hi(1)
     uout(i,UEINT) = uout(i,UEINT)  - dt * pdivu(i)
  enddo

  ! Add gravitational source terms to momentum and energy equations 
  do i = lo(1),hi(1)

     Up  = uin(i,UMX) / uin(i,URHO)
     SrU = uin(i,URHO) * grav(i)

     ! This doesn't work
     ! SrE = SrU*(Up + SrU*dt/(2.d0*rho))

     ! This works 
     ! SrE = SrU*Up 

     SrE = uin(i,UMX ) * grav(i)

     uout(i,UMX  ) = uout(i,UMX  ) + dt * SrU
     uout(i,UEDEN) = uout(i,UEDEN) + dt * SrE
     
  enddo
  
  ! Add gradp term to momentum equation
  do i = lo(1),hi(1)
     dpdx  = (  pgdnv(i+1)- pgdnv(i) ) / dx

     dprdx = 0.d0
     do g=0,ngroups-1
        lamc = 0.5d0*(lamgdnv(i,g)+lamgdnv(i+1,g))
        dprdx = dprdx + lamc*(ergdnv(i+1,g)-ergdnv(i,g))/dx
     end do

     uout(i,UMX) = uout(i,UMX) - dt * dpdx
     ek1 = uout(i,UMX)**2/(2.d0*uout(i,URHO))

     uout(i,UMX) = uout(i,UMX) - dt * dprdx
     ek2 = uout(i,UMX)**2/(2.d0*uout(i,URHO))

     dek = ek2-ek1

     uout(i,UEDEN) = uout(i,UEDEN) +dek
     if (.not. comoving) then ! mixed-frame (single group only)
        Erout(i,0) = Erout(i,0) - dek
     end if
  enddo

  ! Add radiation source term to rho*u, rhoE, and Er
  if (comoving) then 
     do i = lo(1),hi(1)

        ux = 0.5d0*(ugdnv(i) + ugdnv(i+1))

        divu = (ugdnv(i+1)*area(i+1)-ugdnv(i)*area(i))/vol(i)
        dudx = (ugdnv(i+1)-ugdnv(i))/dx

        ! Note that for single group, fspace_type is always 1
        do g=0, ngroups-1

           lamc = 0.5d0*(lamgdnv(i,g)+lamgdnv(i+1,g))
           Eddf = Edd_factor(lamc)
           f1 = (1.d0-Eddf)*0.5d0
           f2 = (3.d0*Eddf-1.d0)*0.5d0
           af(g) = -(f1*divu + f2*dudx)
           
           if (fspace_type .eq. 1) then
              Eddflft = Edd_factor(lamgdnv(i,g))
              f1lft = 0.5d0*(1.d0-Eddflft)
              Eddfrgt = Edd_factor(lamgdnv(i+1,g))
              f1rgt = 0.5d0*(1.d0-Eddfrgt)
           
              Egdc = 0.5d0*(ergdnv(i,g)+ergdnv(i+1,g))
              Erout(i,g) = Erout(i,g) + dt*ux*(f1rgt*ergdnv(i+1,g)-f1lft*ergdnv(i,g))/dx &
                   - dt*f2*Egdc*dudx
           end if

        end do

        if (ngroups.gt.1) then
           ustar = Erout(i,:) / Erscale
           call advect_in_fspace(ustar, af, dt, nstep_fsp)
           Erout(i,:) = ustar * Erscale
        end if 
     end do
  end if

  do i = lo(1),hi(1)+1
     flux(i,UMX) = flux(i,UMX) + dt*area(i)*pgdnv(i)
  enddo

end subroutine consup_rad


subroutine ppflaten(lof, hif, &
     flatn, q, q_l1, q_h1)
  use meth_params_module, only : QPRES, QU
  use radhydro_params_module, only : flatten_pp_threshold, QRADVAR, qptot
  implicit none
  integer, intent(in) :: lof(1), hif(1), q_l1, q_h1
  double precision, intent(in) :: q(q_l1:q_h1,QRADVAR)
  double precision, intent(inout) :: flatn(q_l1:q_h1)

  integer :: i

  do i=lof(1),hif(1)
     if (q(i-1,QU) > q(i+1,QU)) then
        if (q(i,QPRES) < flatten_pp_threshold* q(i,qptot)) then
           flatn(i) = 0.d0
        end if
     end if
  end do

end subroutine ppflaten

end module rad_advection_module
