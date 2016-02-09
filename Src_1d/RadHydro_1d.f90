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
     src,src_l1,src_h1, &
     srcQ,srQ_l1,srQ_h1, &
     courno,dx,dt,ngp,ngf,iflaten)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFA, UFS, UFX, &
       QVAR, QRHO, QU, QV, QW, QGAME, QREINT, QPRES, QTEMP, QFA, QFS, QFX, &
       npassive, upass_map, qpass_map, nadv, small_temp, allow_negative_energy
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
  integer          ::  srQ_l1,srQ_h1
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
  double precision ::  srcQ(srQ_l1:srQ_h1,QVAR)
  double precision :: dx, dt, courno
  integer iflaten

  integer          :: i, g
  integer          :: ngp, ngf, loq(1), hiq(1)
  integer          :: n, nq
  integer          :: iadv, ispec, iaux, ipassive
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

  do ipassive = 1, npassive
     n  = upass_map(ipassive)
     nq = qpass_map(ipassive)
     q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
  end do

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

     ! If necessary, reset the energy using small_temp
     if ((allow_negative_energy .eq. 0) .and. (q(i,QREINT) .lt. 0)) then
        q(i,QTEMP) = small_temp
        eos_state % T = q(i,QTEMP)
        call eos(eos_input_rt, eos_state)
        q(i,QPRES ) = eos_state % p
        q(i,QREINT) = eos_state % e
        if (q(i,QREINT) .lt. 0.d0) then
           print *,'   '
           print *,'>>> Error: ctoprim ',i
           print *,'>>> ... new e from eos call is negative ',q(i,QREINT)
           print *,'    '
           call bl_error("Error:: ctoprim")
        end if
     end if
     
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
  do i = loq(1), hiq(1)
     srcQ(i,QRHO   ) = src(i,URHO)
     srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
     srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + 0.5d0 * q(i,QU)**2 * srcQ(i,QRHO)
     srcQ(i,QPRES  ) = dpde(i) * (srcQ(i,QREINT) - q(i,QREINT)*srcQ(i,QRHO)/q(i,QRHO)) / q(i,QRHO) + &
              dpdrho(i) * srcQ(i,QRHO)! - &
!              sum(dpdX_er(i,:)*(src(i,UFS:UFS+nspec-1) - q(i,QFS:QFS+nspec-1)*srcQ(i,QRHO)))/q(i,QRHO)


     do ipassive = 1, npassive
        n  = upass_map(ipassive)
        nq = qpass_map(ipassive)
        srcQ(i,nq) = (src(i,n) - q(i,nq) * srcQ(i,QRHO))/q(i,QRHO)
     enddo

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

  ! Compute flattening coef for slope calculations
  if (first_order_hydro) then
     flatn = 0.d0

  else if (iflaten.eq.1) then
     loq(1)=lo(1)-ngf
     hiq(1)=hi(1)+ngf
     call uflaten([loq(1), 0, 0], [hiq(1), 0, 0], &
                  q(:,qptot), &
                  q(:,QU), q(:,QV), q(:,QW), &
                  flatn, [q_l1, 0, 0], [q_h1, 0, 0])
     call uflaten([loq(1), 0, 0], [hiq(1), 0, 0], &
                  q(:,qpres), &
                  q(:,QU), q(:,QV), q(:,QW), &
                  flatg, [q_l1, 0, 0], [q_h1, 0, 0])

     do i = loq(1), hiq(1)
        flatn(i) = flatn(i) * flatg(i)
     enddo

     if (flatten_pp_threshold > 0.d0) then
        call ppflaten(loq, hiq, flatn, q, q_l1, q_h1)
     end if

  else
     flatn = 1.d0
  endif

  deallocate(dpdrho,dpde,flatg)

  q(:,QGAME) = 0.d0 ! QGAME is not used in radiation hydro. Setting it to 0 to mute valgrind.
  
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
     ilo,ihi,dx,dt, &
     flux ,   fd_l1,   fd_h1, &
     rflux,  rfd_l1,  rfd_h1, &
     pgdnv,pgdnv_l1,pgdnv_h1, &
     ergdnv,ergdnv_l1,ergdnv_h1, &
     lamgdnv,lamgdnv_l1,lamgdnv_h1, &
     ugdnv,ugdnv_l1,ugdnv_h1, &
     dloga,dloga_l1,dloga_h1)

  use meth_params_module, only : NVAR, ppm_type
  use radhydro_params_module, only : QRADVAR
  use rad_params_module, only : ngroups
  use riemann_module, only : cmpflx
  use trace_ppm_rad_module, only : trace_ppm_rad
  
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
  double precision  pgdnv(pgdnv_l1:pgdnv_h1)
  double precision ergdnv(ergdnv_l1:ergdnv_h1, 0:ngroups-1)
  double precision lamgdnv(lamgdnv_l1:lamgdnv_h1, 0:ngroups-1)
  double precision ugdnv(ugdnv_l1:ugdnv_h1)
  double precision dloga(dloga_l1:dloga_h1)
  
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
          qm,qp,ilo-1,ihi+1, &
          ilo,ihi,domlo,domhi,dx,dt)
  else
     call bl_error("ppm_type <=0 is not supported in umeth1d_rad")
     ! call trace(q,dq,c,flatn,qd_l1,qd_h1, &
     !      dloga,dloga_l1,dloga_h1, &
     !      srcQ,src_l1,src_h1, &
     !      qm,qp,ilo-1,ihi+1, &
     !      ilo,ihi,domlo,domhi,dx,dt)
  end if
  
  !     Solve Riemann problem, compute xflux from improved predicted states 
  call cmpflx(lo, hi, domlo, domhi, &
              qm, qp, ilo-1,ihi+1, &
              flux ,  fd_l1, fd_h1, &
              pgdnv,pgdnv_l1,pgdnv_h1, &
              ugdnv,ugdnv_l1,ugdnv_h1, &
              lam, lam_l1, lam_h1, & 
              rflux, rfd_l1,rfd_h1, &
              ergdnv,ergdnv_l1,ergdnv_h1, &
              lamgdnv,lamgdnv_l1,lamgdnv_h1, &
              gamcg,gamc,csml,c,qd_l1,qd_h1,ilo,ihi)

  deallocate (dq,qm,qp)

end subroutine umeth1d_rad


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
  use radhydro_params_module, only : fspace_type, comoving
  use radhydro_nd_module, only : advect_in_fspace
  use fluxlimiter_module, only : Edd_factor
  use advection_util_module, only : normalize_species_fluxes

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
  double precision  flux( flux_l1: flux_h1,NVAR)
  double precision rflux(rflux_l1:rflux_h1, 0:ngroups-1)
  double precision  flat( flat_l1: flat_h1)
  double precision  area( area_l1: area_h1)
  double precision    vol(vol_l1:vol_h1)
  double precision    div(lo(1):hi(1)+1)
  double precision  pdivu(lo(1):hi(1)  )
  double precision dx, dt

  integer          :: i, n, g
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
                + ( flux(i,n) - flux(i+1,n) ) / vol(i)
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
