module rad_advection_module

  use bl_constants_module, only : ZERO, HALF, ONE

  implicit none

  private

  public umeth1d_rad, ctoprim_rad, consup_rad, ppflaten

contains

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

subroutine ctoprim_rad(lo,hi,uin,uin_l1,uin_h1, &
                       Erin, Erin_l1, Erin_h1, &
                       lam, lam_l1, lam_h1, &
                       q, q_l1, q_h1, &
                       qaux, qa_l1, qa_h1, &
                       src,src_l1,src_h1, &
                       srcQ,srQ_l1,srQ_h1, &
                       courno,dx,dt,ngp,ngf)

  use network, only : nspec, naux
  use eos_module
  use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
                                 QVAR, QRHO, QU, QV, QW, QGAME, QREINT, &
                                 QPRES, QTEMP, QFS, QFX, &
                                 NQAUX, QGAMC, QGAMCG, QC, QCG, QCSML, QDPDR, QDPDE,  &
                                 npassive, upass_map, qpass_map, nadv, small_temp, allow_negative_energy
  use radhydro_params_module, only : QRADVAR, qrad, qradhi, qptot, qreitot, comoving, &
       first_order_hydro
  use rad_params_module, only : ngroups
  use rad_util_module, only : compute_ptot_ctot
  
  implicit none

  double precision, parameter:: small = 1.d-8

  ! Will give primitive variables on lo-ngp:hi+ngp.  Declared dimensions of
  ! q,c,gamc,csml are given by DIMS(q).  This declared region is
  ! assumed to encompass lo-ngp:hi+ngp.  

  integer          :: lo(1), hi(1)
  integer          :: uin_l1,uin_h1, Erin_l1, Erin_h1, lam_l1, lam_h1 
  integer          :: q_l1, q_h1, qa_l1, qa_h1
  integer          ::  src_l1,src_h1
  integer          ::  srQ_l1,srQ_h1
  double precision ::   uin(uin_l1:uin_h1,NVAR)
  double precision :: Erin(Erin_l1:Erin_h1, 0:ngroups-1)
  double precision :: lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision ::     q(  q_l1:  q_h1,QRADVAR)
  double precision ::  qaux( qa_l1: qa_h1,NQAUX)
  double precision ::   src(src_l1:src_h1,NVAR)
  double precision ::  srcQ(srQ_l1:srQ_h1,QVAR)
  double precision :: dx, dt, courno

  integer          :: i, g
  integer          :: ngp, ngf, loq(1), hiq(1)
  integer          :: ipassive, n, nq
  double precision :: courx, courmx

  double precision, allocatable :: dpdrho(:), dpde(:)

  double precision :: ptot, ctot, gamc_tot

  type(eos_t) :: eos_state

  loq(1) = lo(1)-ngp
  hiq(1) = hi(1)+ngp

  allocate(dpdrho(q_l1:q_h1))
  allocate(dpde  (q_l1:q_h1))


  ! Make q (all but p), except put e in slot for rho.e, fix after eos
  ! call.  The temperature is used as an initial guess for the eos
  ! call and will be overwritten
  do i = loq(1), hiq(1)

     if (uin(i,URHO) .le. 0.d0) then
        print *,'   '
        print *,'>>> Error: Castro_1d::ctoprim ',i
        print *,'>>> ... negative density ', uin(i,URHO)
        print *,'    '
        call bl_error("Error:: Castro_1d.f90 :: ctoprim")
     end if

     q(i,QRHO) = uin(i,URHO)
     q(i,QU) = uin(i,UMX)/uin(i,URHO)

     ! we should set this based on the kinetical energy, using dual_energy_eta1
     q(i,QREINT) = uin(i,UEINT)/q(i,QRHO)

     q(i,QTEMP ) = uin(i,UTEMP)
     q(i,qrad:qradhi) = Erin(i,0:ngroups-1)
  enddo

  ! Load passive quantities, c, into q, assuming they arrived in uin
  ! as rho.c
  do ipassive = 1, npassive
     n  = upass_map(ipassive)
     nq = qpass_map(ipassive)
     q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
  enddo

  ! Get gamc, p, T, c, csml using q state
  do i = loq(1), hiq(1)

     eos_state % rho = q(i,QRHO)
     eos_state % T   = q(i,QTEMP)
     eos_state % e   = q(i,QREINT)
     eos_state % xn  = q(i,QFS:QFS+nspec-1)
     eos_state % aux = q(i,QFX:QFX+naux-1)

     ! If necessary, reset the energy using small_temp -- TODO: we handle this
     ! differently now
     if ((allow_negative_energy .eq. 0) .and. (q(i,QREINT) .lt. 0)) then
        q(i,QTEMP) = small_temp
        eos_state % T = q(i,QTEMP)
        call eos(eos_input_rt, eos_state)
        q(i,QPRES ) = eos_state % p
        q(i,QREINT) = eos_state % e
        if (q(i,QREINT) .lt. 0.d0) then
           print *,'   '
           print *,'>>> Error: ctoprim ',i
           print *,'>>> ... new e from eos call is negative ', q(i,QREINT)
           print *,'    '
           call bl_error("Error:: ctoprim")
        end if
     end if
     
     call eos(eos_input_re, eos_state)

     q(i,QTEMP) = eos_state % T
     q(i,QREINT) = eos_state % e * q(i,QRHO)
     q(i,QPRES) = eos_state % p
     q(i,QGAME) = q(i,QPRES) / q(i,QREINT) + ONE

     qaux(i,QGAMCG) = eos_state % gam1
     qaux(i,QCG)    = eos_state % cs
     qaux(i,QDPDR)  = eos_state % dpdr_e
     qaux(i,QDPDE)  = eos_state % dpde

     call compute_ptot_ctot(lam(i,:), q(i,:), qaux(i,QCG), ptot, ctot, gamc_tot)

     q(i,QPTOT) = ptot

     qaux(i,QC)    = ctot
     qaux(i,QGAMC) = gamc_tot

     qaux(i,QCSML) = max(small, small * ctot)
  end do

  ! Make this "rho e" instead of "e"
  do i = loq(1),hiq(1)
     q(i,qreitot) = q(i,QREINT) + sum(q(i,qrad:qradhi))
  enddo

  ! compute srcQ terms
  do i = loq(1), hiq(1)
     srcQ(i,QRHO   ) = src(i,URHO)
     srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
     srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + &
          HALF * q(i,QU)**2 * srcQ(i,QRHO)
     srcQ(i,QPRES  ) = dpde(i) * (srcQ(i,QREINT) - &
          q(i,QREINT)*srcQ(i,QRHO)/q(i,QRHO)) / q(i,QRHO) + &
          dpdrho(i) * srcQ(i,QRHO)! - &
!              sum(dpdX_er(i,:)*(src(i,UFS:UFS+nspec-1) - q(i,QFS:QFS+nspec-1)*srcQ(i,QRHO)))/q(i,QRHO)


     do ipassive = 1, npassive
        n  = upass_map(ipassive)
        nq = qpass_map(ipassive)
        srcQ(i,nq) = (src(i,n) - q(i,nq) * srcQ(i,QRHO))/q(i,QRHO)
     enddo

  end do

  ! Compute running max of Courant number over grids
  courmx = courno
  do i = lo(1),hi(1)

     courx  = ( qaux(i,QC) + abs(q(i,QU)) ) * dt/dx
     courmx = max( courmx, courx )

     if (courx .gt. ONE) then
        print *,'   '
        call bl_warning("Warning:: RadHydro_1d.f90 :: CFL violation in ctoprim")
        print *,'>>> ... (u+c) * dt / dx > 1 ', courx
        print *,'>>> ... at cell (i)       : ',i
        print *,'>>> ... u, c                ',q(i,QU), qaux(i,QC)
        print *,'>>> ... density             ',q(i,QRHO)
     end if
  enddo
  courno = courmx

  deallocate(dpdrho,dpde)
  
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
                       q, qd_l1, qd_h1, &
                       qaux, qa_l1, qa_h1, &
                       srcQ,src_l1,src_h1, &
                       ilo,ihi,dx,dt, &
                       flux ,   fd_l1,   fd_h1, &
                       rflux,  rfd_l1,  rfd_h1, &
                       q1, q1_l1, q1_h1, &
                       ergdnv,ergdnv_l1,ergdnv_h1, &
                       lamgdnv,lamgdnv_l1,lamgdnv_h1, &
                       dloga,dloga_l1,dloga_h1)

  use meth_params_module, only : NVAR, ppm_type, QC, QCG, QCSML, QGAMC, QGAMCG, NQAUX, NGDNV
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
  integer qa_l1,qa_h1
  integer src_l1,src_h1
  integer fd_l1,fd_h1
  integer rfd_l1,rfd_h1
  integer q1_l1, q1_h1
  integer ergdnv_l1,ergdnv_h1
  integer lamgdnv_l1,lamgdnv_h1
  integer ilo,ihi
  double precision dx, dt
  double precision lam(lam_l1:lam_h1, 0:ngroups-1)
  double precision     q(   qd_l1:qd_h1,QRADVAR)
  double precision  qaux(   qd_l1:qd_h1,NQAUX)
  double precision flatn(   qd_l1:qd_h1)
  double precision  flux(fd_l1   :fd_h1,NVAR)
  double precision rflux(rfd_l1:rfd_h1, 0:ngroups-1)
  double precision  srcQ(src_l1  :src_h1,NVAR)
  double precision    q1(   q1_l1:q1_h1, NGDNV)
  double precision ergdnv(ergdnv_l1:ergdnv_h1, 0:ngroups-1)
  double precision lamgdnv(lamgdnv_l1:lamgdnv_h1, 0:ngroups-1)
  double precision dloga(dloga_l1:dloga_h1)
  
  ! Left and right state arrays (edge centered, cell centered)
  double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

  ! Work arrays to hold 3 planes of riemann state and conservative fluxes
  allocate ( dq(ilo-1:ihi+1,QRADVAR))
  allocate ( qm(ilo-1:ihi+1,QRADVAR))
  allocate ( qp(ilo-1:ihi+1,QRADVAR))

  ! Trace to edges w/o transverse flux correction terms
  if (ppm_type .gt. 0) then
     call trace_ppm_rad(lam, lam_l1, lam_h1, &       
                        q,qaux(:,QC),qaux(:,QCG),qaux(:,QGAMC),qaux(:,QGAMCG),flatn,qd_l1,qd_h1, &
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
  
  ! Solve Riemann problem, compute xflux from improved predicted states 
  call cmpflx(lo, hi, domlo, domhi, &
              qm, qp, ilo-1,ihi+1, &
              flux ,  fd_l1, fd_h1, &
              q1, q1_l1, q1_h1, &
              lam, lam_l1, lam_h1, & 
              rflux, rfd_l1,rfd_h1, &
              ergdnv,ergdnv_l1,ergdnv_h1, &
              lamgdnv,lamgdnv_l1,lamgdnv_h1, &
              qaux(:,QGAMCG),qaux(:,QGAMC),qaux(:,QCSML),qaux(:,QC),qd_l1,qd_h1,ilo,ihi)

  deallocate (dq,qm,qp)

end subroutine umeth1d_rad


! ::: 
! ::: ------------------------------------------------------------------
! ::: 
     
subroutine consup_rad(uin,  uin_l1,  uin_h1, &
     uout, uout_l1 ,uout_h1, &
     Erin,Erin_l1,Erin_h1, &
     Erout,Erout_l1,Erout_h1, &
     q1,q1_l1,q1_h1, &
     ergdnv,ergdnv_l1,ergdnv_h1, &
     lamgdnv,lamgdnv_l1,lamgdnv_h1, &
     src,  src_l1,  src_h1, &
     flux, flux_l1, flux_h1, &
     rflux,rflux_l1,rflux_h1, &
     flat, flat_l1, flat_h1, &
     area,area_l1,area_h1, &
     vol,vol_l1,vol_h1, &
     div,pdivu,lo,hi,dx,dt, &
     nstep_fsp)

  use meth_params_module, only : difmag, NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, NGDNV, GDU, GDPRES
  use rad_params_module, only : ngroups, nugroup, dlognu
  use radhydro_params_module, only : fspace_type, comoving
  use radhydro_nd_module, only : advect_in_fspace
  use fluxlimiter_module, only : Edd_factor
  use advection_util_1d_module, only : normalize_species_fluxes
  use prob_params_module, only : coord_type

  implicit none
  integer nstep_fsp
  integer lo(1), hi(1)
  integer   uin_l1,  uin_h1
  integer  uout_l1, uout_h1
  integer  Erin_l1, Erin_h1
  integer Erout_l1,Erout_h1
  integer    q1_l1,   q1_h1
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
  double precision    q1(q1_l1:q1_h1, NGDNV)
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
  call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

  do n = 1, NVAR
     if ( n == UTEMP ) then
        flux(:,n) = ZERO
     else
        do i = lo(1),hi(1)+1
           div1 = difmag*min(ZERO,div(i))
           flux(i,n) = flux(i,n) &
                + dx*div1*(uin(i,n) - uin(i-1,n))
           flux(i,n) = area(i) * flux(i,n) * dt
        enddo
     endif
  enddo

  do g=0, ngroups-1
     do i = lo(1),hi(1)+1
        div1 = difmag*min(ZERO,div(i))
        rflux(i,g) = rflux(i,g) &
             + dx*div1*(Erin(i,g) - Erin(i-1,g))
        rflux(i,g) = area(i) * rflux(i,g) * dt
     enddo
  end do

  do n = 1, NVAR
     do i = lo(1),hi(1)
        uout(i,n) = uout(i,n) + ( flux(i,n) - flux(i+1,n) ) / vol(i)
     enddo
  enddo

  do g=0, ngroups-1
     do i = lo(1),hi(1)
        Erout(i,g) = Erin(i,g) + (rflux(i,g) - rflux(i+1,g) ) / vol(i) 
     enddo
  end do

  ! Add source term to (rho e)
  do i = lo(1),hi(1)
     uout(i,UEINT) = uout(i,UEINT) - dt * pdivu(i)
  enddo

  ! Add gradp term to momentum equation
  do i = lo(1),hi(1)
     dpdx  = (  q1(i+1,GDPRES)- q1(i,GDPRES) ) / dx

     dprdx = ZERO
     do g=0,ngroups-1
        lamc = HALF*(lamgdnv(i,g)+lamgdnv(i+1,g))
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

        ux = HALF*(q1(i,GDU) + q1(i+1,GDU))

        divu = (q1(i+1,GDU)*area(i+1)-q1(i,GDU)*area(i))/vol(i)
        dudx = (q1(i+1,GDU)-q1(i,GDU))/dx

        ! Note that for single group, fspace_type is always 1
        do g=0, ngroups-1

           lamc = HALF*(lamgdnv(i,g)+lamgdnv(i+1,g))
           Eddf = Edd_factor(lamc)
           f1 = (ONE-Eddf)*HALF
           f2 = (3.d0*Eddf-ONE)*HALF
           af(g) = -(f1*divu + f2*dudx)
           
           if (fspace_type .eq. 1) then
              Eddflft = Edd_factor(lamgdnv(i,g))
              f1lft = HALF*(ONE-Eddflft)
              Eddfrgt = Edd_factor(lamgdnv(i+1,g))
              f1rgt = HALF*(ONE-Eddfrgt)
           
              Egdc = HALF*(ergdnv(i,g)+ergdnv(i+1,g))
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

  if (coord_type .eq. 0) then
     do i = lo(1),hi(1)+1
        flux(i,UMX) = flux(i,UMX) + dt*area(i)*q1(i,GDPRES)
     enddo
  end if

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
           flatn(i) = ZERO
        end if
     end if
  end do

end subroutine ppflaten

end module rad_advection_module
