module advection_module

  implicit none

  private

  public umeth1d, ctoprim, consup, enforce_minimum_density, normalize_new_species, &
       normalize_species_fluxes

contains

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

  subroutine umeth1d(lo,hi,domlo,domhi, &
                     q,c,gamc,csml,flatn,qd_l1,qd_h1, &
                     srcQ,src_l1,src_h1, &
                     grav, gv_l1, gv_h1, &
                     ilo,ihi,dx,dt, &
                     flux ,   fd_l1,   fd_h1, &
                     pgdnv,pgdnv_l1,pgdnv_h1, &
                     ugdnv,ugdnv_l1,ugdnv_h1, &
                     dloga,dloga_l1,dloga_h1)

    use meth_params_module, only : QVAR, NVAR, ppm_type
    use riemann_module
    use trace_module, only : trace
    use trace_ppm_module, only : trace_ppm

    implicit none
    integer lo(1),hi(1)
    integer domlo(1),domhi(1)
    integer dloga_l1,dloga_h1
    integer qd_l1,qd_h1
    integer src_l1,src_h1
    integer fd_l1,fd_h1
    integer pgdnv_l1,pgdnv_h1
    integer ugdnv_l1,ugdnv_h1
    integer gv_l1,gv_h1
    integer ilo,ihi
    double precision dx, dt
    double precision     q(   qd_l1:qd_h1,QVAR)
    double precision  gamc(   qd_l1:qd_h1)
    double precision flatn(   qd_l1:qd_h1)
    double precision  csml(   qd_l1:qd_h1)
    double precision     c(   qd_l1:qd_h1)
    double precision  flux(fd_l1   :fd_h1,NVAR)
    double precision  srcQ(src_l1  :src_h1,NVAR)
    double precision  grav(gv_l1   :gv_h1)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision ugdnv(ugdnv_l1:ugdnv_h1)
    double precision dloga(dloga_l1:dloga_h1)
    
    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    allocate ( dq(ilo-1:ihi+1,QVAR))
    allocate ( qm(ilo-1:ihi+1,QVAR))
    allocate ( qp(ilo-1:ihi+1,QVAR))

    ! Trace to edges w/o transverse flux correction terms
    if (ppm_type .gt. 0) then
       call trace_ppm(q,dq,c,flatn,gamc,qd_l1,qd_h1, &
                      dloga,dloga_l1,dloga_h1, &
                      srcQ,src_l1,src_h1, &
                      grav,gv_l1,gv_h1, &
                      qm,qp,ilo-1,ihi+1, &
                      ilo,ihi,domlo,domhi,dx,dt)
    else
       call trace(q,dq,c,flatn,qd_l1,qd_h1, &
                  dloga,dloga_l1,dloga_h1, &
                  srcQ,src_l1,src_h1, &
                  grav,gv_l1,gv_h1, &
                  qm,qp,ilo-1,ihi+1, &
                  ilo,ihi,domlo,domhi,dx,dt)
    end if

    ! Solve Riemann problem, compute xflux from improved predicted states 
    call cmpflx(lo, hi, domlo, domhi, &
                qm, qp, ilo-1,ihi+1, &
                flux ,  fd_l1, fd_h1, &
                pgdnv,pgdnv_l1,pgdnv_h1, &
                ugdnv,ugdnv_l1,ugdnv_h1, &
                gamc, csml,c,qd_l1,qd_h1,ilo,ihi)

    deallocate (dq,qm,qp)

  end subroutine umeth1d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi,uin,uin_l1,uin_h1, &
                     q,c,gamc,csml,flatn,q_l1,q_h1,&
                     src,srcQ,src_l1,src_h1, &
                     courno,dx,dt,ngp,ngf)
    
    ! Will give primitive variables on lo-ngp:hi+ngp, and flatn on
    ! lo-ngf:hi+ngf if iflaten=1.  Declared dimensions of
    ! q,c,gamc,csml,flatn are given by DIMS(q).  This declared
    ! region is assumed to encompass lo-ngp:hi+ngp.  Also, uflaten
    ! call assumes ngp>=ngf+3 (ie, primitve data is used by the
    ! routine that computes flatn).

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
                                   UFA, UFS, UFX, &
                                   QVAR, QRHO, QU, QREINT, QPRES, QTEMP, &
                                   QFA, QFS, QFX, &
                                   nadv, small_temp, allow_negative_energy, use_flattening
    use flatten_module
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer          :: lo(1), hi(1)
    integer          :: uin_l1,uin_h1
    integer          :: q_l1,q_h1
    integer          ::  src_l1,src_h1
    double precision ::   uin(uin_l1:uin_h1,NVAR)
    double precision ::     q(  q_l1:  q_h1,QVAR)
    double precision ::     c(  q_l1:  q_h1)
    double precision ::  gamc(  q_l1:  q_h1)
    double precision ::  csml(  q_l1:  q_h1)
    double precision :: flatn(  q_l1:  q_h1)
    double precision ::   src(src_l1:src_h1,NVAR)
    double precision ::  srcQ(src_l1:src_h1,QVAR)
    double precision :: dx, dt, courno
    
    integer          :: i
    integer          :: pt_index(1)
    integer          :: ngp, ngf, loq(1), hiq(1)
    integer          :: n, nq
    integer          :: iadv, ispec, iaux
    double precision :: courx, courmx
    
    double precision, allocatable :: dpdrho(:), dpde(:) !, dpdX_er(:,:)
    
    type (eos_t) :: eos_state
    
    loq(1) = lo(1)-ngp
    hiq(1) = hi(1)+ngp
    
    allocate(dpdrho(q_l1:q_h1))
    allocate(dpde  (q_l1:q_h1))
    ! allocate(dpdX_er(q_l1:q_h1,nspec))

    ! Make q (all but p), except put e in slot for rho.e, fix after
    ! eos call The temperature is used as an initial guess for the eos
    ! call and will be overwritten
    do i = loq(1),hiq(1)

       if (uin(i,URHO) .le. ZERO) then
          print *,'   '
          print *,'>>> Error: Castro_1d::ctoprim ',i
          print *,'>>> ... negative density ',uin(i,URHO)
          print *,'    '
          call bl_error("Error:: Castro_1d.f90 :: ctoprim")
       end if
       
       q(i,QRHO) = uin(i,URHO)
       q(i,QU) = uin(i,UMX)/uin(i,URHO)
       !        eken = HALF*q(i,QU)**2
       !        q(i,QREINT) = uin(i,UEDEN)/q(i,QRHO) - eken
       q(i,QREINT) = uin(i,UEINT)/q(i,QRHO)
       q(i,QTEMP ) = uin(i,UTEMP)
    enddo
    
    ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
    do iadv = 1, nadv
       n  = UFA + iadv - 1
       nq = QFA + iadv - 1
       q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
    enddo
    
    ! Load species, c, into q, assuming they arrived in uin as rho.c
    do ispec = 1, nspec
       n  = UFS + ispec - 1
       nq = QFS + ispec - 1
       q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
    enddo
    
    ! Load auxiliary variables which are needed in the EOS
    do iaux = 1, naux
       n  = UFX + iaux - 1
       nq = QFX + iaux - 1
       q(loq(1):hiq(1),nq) = uin(loq(1):hiq(1),n)/q(loq(1):hiq(1),QRHO)
    enddo
    
    ! Get gamc, p, T, c, csml using q state
    do i = loq(1), hiq(1)
       
       pt_index(1) = i
       
       eos_state % T   = q(i,QTEMP)
       eos_state % rho = q(i,QRHO)
       eos_state % xn  = q(i,QFS:QFS+nspec-1)
       eos_state % aux = q(i,QFX:QFX+naux-1)
       
       ! If necessary, reset the energy using small_temp
       if (allow_negative_energy .eq. 0 .and. q(i,QREINT) .le. ZERO) then
          
          q(i,QTEMP) = small_temp
          eos_state % T = q(i,QTEMP)
          
          call eos(eos_input_rt, eos_state, pt_index = pt_index)
          q(i,QREINT) = eos_state % e
          
          if (q(i,QREINT) .lt. ZERO) then
             print *,'   '
             print *,'>>> Error: Castro_1d::ctoprim ',i
             print *,'>>> ... new e from eos (input_rt) call is negative ',q(i,QREINT)
             print *,'    '
             call bl_error("Error:: Castro_1d.f90 :: ctoprim")
          end if
       end if
       
       eos_state % e = q(i,QREINT)
       
       call eos(eos_input_re, eos_state, pt_index = pt_index)
       
       q(i,QTEMP)  = eos_state % T
       q(i,QREINT) = eos_state % e
       q(i,QPRES)  = eos_state % p
       
       dpdrho(i) = eos_state % dpdr_e
       dpde(i)   = eos_state % dpde
       c(i)      = eos_state % cs
       gamc(i)   = eos_state % gam1
       
       csml(i) = max(small, small * c(i))
    end do
    
    ! Make this "rho e" instead of "e"
    do i = loq(1),hiq(1)
       q(i,QREINT ) = q(i,QREINT )*q(i,QRHO)
    enddo
    
    ! compute srcQ terms
    do i = lo(1)-1, hi(1)+1
       srcQ(i,QRHO   ) = src(i,URHO)
       srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
       srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + HALF * q(i,QU)**2 * srcQ(i,QRHO)
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

    ! Compute running max of Courant number over grids
    courmx = courno
    do i = lo(1),hi(1)
       
       courx  = ( c(i)+abs(q(i,QU)) ) * dt/dx
       courmx = max( courmx, courx )
       
       if (courx .gt. ONE) then
          print *,'   '
          call bl_warning("Warning:: Castro_1d.f90 :: CFL violation in ctoprim")
          print *,'>>> ... (u+c) * dt / dx > 1 ', courx
          print *,'>>> ... at cell (i)       : ',i
          print *,'>>> ... u, c                ',q(i,QU), c(i)
          print *,'>>> ... density             ',q(i,QRHO)
       end if
    enddo
    courno = courmx
    
    ! Compute flattening coef for slope calculations
    if (use_flattening == 1) then
       loq(1)=lo(1)-ngf
       hiq(1)=hi(1)+ngf
       call uflaten(loq,hiq, &
                    q(q_l1,QPRES), &
                    q(q_l1,QU), &
                    flatn,q_l1,q_h1)
    else
       flatn = ONE
    endif
    
    deallocate(dpdrho,dpde)
    
  end subroutine ctoprim

! ::: 
! ::: ------------------------------------------------------------------
! ::: 
     
  subroutine consup( &
                    uin,  uin_l1,  uin_h1, &
                    uout, uout_l1 ,uout_h1, &
                    pgdnv,pgdnv_l1,pgdnv_h1, &
                    src,  src_l1,  src_h1, &
                    grav, grav_l1, grav_h1, &
                    flux, flux_l1, flux_h1, &
                    area,area_l1,area_h1, &
                    vol,vol_l1,vol_h1, &
                    div,pdivu,lo,hi,dx,dt)

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
                                   UFS, UFX, &
                                   normalize_species
    use bl_constants_module

    implicit none
    integer lo(1), hi(1)
    integer   uin_l1,  uin_h1
    integer  uout_l1, uout_h1
    integer pgdnv_l1,pgdnv_h1
    integer   src_l1,  src_h1
    integer  grav_l1, grav_h1
    integer  flux_l1, flux_h1
    integer  area_l1, area_h1
    integer   vol_l1,  vol_h1
    double precision   uin(uin_l1:uin_h1,NVAR)
    double precision  uout(uout_l1:uout_h1,NVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision   src(  src_l1:  src_h1,NVAR)
    double precision  grav( grav_l1: grav_h1)
    double precision  flux( flux_l1: flux_h1,NVAR)
    double precision  area( area_l1: area_h1)
    double precision    vol(vol_l1:vol_h1)
    double precision    div(lo(1):hi(1)+1)
    double precision  pdivu(lo(1):hi(1)  )
    double precision dx, dt
    
    integer          :: i, n
    double precision :: div1, dpdx
    double precision :: SrU,Up,SrE
    
    ! Normalize the species fluxes
    if (normalize_species .eq. 1) &
         call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)
    
    do n = 1, NVAR
       if ( n.eq.UTEMP ) then
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
    
    ! Add gradp term to momentum equation
    do i = lo(1),hi(1)
       
       !        dpdx = HALF * (area(i)+area(i+1)) * (pgdnv(i+1)-pgdnv(i)) / vol(i)
       dpdx = ( pgdnv(i+1)-pgdnv(i) ) / dx
       uout(i,UMX)   = uout(i,UMX) - dt * dpdx
       
    enddo

    ! Add source term to (rho e)
    do i = lo(1),hi(1)
       uout(i,UEINT) = uout(i,UEINT)  - dt * pdivu(i)
    enddo
    
    ! Add gravitational source terms to momentum and energy equations 
    do i = lo(1),hi(1)
       
       Up  = uin(i,UMX) / uin(i,URHO)
       SrU = uin(i,URHO) * grav(i)
       
       ! This doesn't work
       ! SrE = SrU*(Up + SrU*dt/(TWO*rho))
       
       ! This works 
       ! SrE = SrU*Up 
       
       SrE = uin(i,UMX ) * grav(i)
       
       uout(i,UMX  ) = uout(i,UMX  ) + dt * SrU
       uout(i,UEDEN) = uout(i,UEDEN) + dt * SrE
       
    enddo
    
    do i = lo(1),hi(1)+1
       flux(i,UMX) = flux(i,UMX) + dt*area(i)*pgdnv(i)
    enddo
    
  end subroutine consup
  
! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module

    implicit none
    
    integer          :: lo(1),hi(1)
    integer          :: flux_l1,flux_h1
    double precision :: flux(flux_l1:flux_h1,NVAR)
    
    ! Local variables
    integer          :: i,n
    double precision :: sum,fac
    
    do i = lo(1),hi(1)+1
       sum = ZERO
       do n = UFS, UFS+nspec-1
          sum = sum + flux(i,n)
       end do
       if (sum .ne. ZERO) then
          fac = flux(i,URHO) / sum
       else
          fac = ONE
       end if
       do n = UFS, UFS+nspec-1
          flux(i,n) = flux(i,n) * fac
       end do
    end do
    
  end subroutine normalize_species_fluxes

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine enforce_minimum_density( uin,  uin_l1, uin_h1, &
                                      uout,uout_l1,uout_h1,&
                                      lo, hi, mass_added, eint_added, eden_added, verbose)
    use network, only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, &
                                   UFS, UFX, UFA, small_dens, nadv
    use bl_constants_module

    implicit none
    
    integer          :: lo(1), hi(1), verbose
    integer          :: uin_l1 , uin_h1
    integer          :: uout_l1,uout_h1
    double precision :: uin(uin_l1:uin_h1,NVAR)
    double precision :: uout(uout_l1:uout_h1,NVAR)
    double precision :: mass_added, eint_added, eden_added
    
    ! Local variables
    integer          :: i,ii,n
    double precision :: min_dens
    double precision, allocatable :: fac(:)
    double precision :: initial_mass, final_mass
    double precision :: initial_eint, final_eint
    double precision :: initial_eden, final_eden
    
    allocate(fac(lo(1):hi(1)))
    
    initial_mass = ZERO
    final_mass = ZERO
    
    initial_eint = ZERO
    final_eint = ZERO
    
    initial_eden = ZERO
    final_eden = ZERO
    
    do i = lo(1),hi(1)
       
       initial_mass = initial_mass + uout(i,URHO)
       initial_eint = initial_eint + uout(i,UEINT)
       initial_eden = initial_eden + uout(i,UEDEN)
       
       if (uout(i,URHO) .eq. ZERO) then
          
          print *,'   '
          print *,'>>> Error: Castro_1d::enforce_minimum_density ',i
          print *,'>>> ... density exactly zero in grid ',lo(1),hi(1)
          print *,'    '
          call bl_error("Error:: Castro_1d.f90 :: enforce_minimum_density")
          
       else if (uout(i,URHO) < small_dens) then
          
          min_dens = uin(i,URHO)
          do ii = -1,1
             min_dens = min(min_dens,uin(i+ii,URHO))
             if (uout(i+ii,URHO) > small_dens) &
                  min_dens = min(min_dens,uout(i+ii,URHO))
          end do
          
          if (verbose .gt. 0) then
             if (uout(i,URHO) < ZERO) then
                print *,'   '
                print *,'>>> Warning: Castro_1d::enforce_minimum_density ',i
                print *,'>>> ... resetting negative density '
                print *,'>>> ... from ',uout(i,URHO),' to ',min_dens
                print *,'    '
             else
                print *,'   '
                print *,'>>> Warning: Castro_1d::enforce_minimum_density ',i
                print *,'>>> ... resetting small density '
                print *,'>>> ... from ',uout(i,URHO),' to ',min_dens
                print *,'    '
             end if
          end if
          
          fac(i) = min_dens / uout(i,URHO)
          
       end if

    enddo

    do i = lo(1),hi(1)

       if (uout(i,URHO) < small_dens) then

          uout(i,URHO ) = uout(i,URHO ) * fac(i)
          uout(i,UEDEN) = uout(i,UEDEN) * fac(i)
          uout(i,UEINT) = uout(i,UEINT) * fac(i)
          uout(i,UMX  ) = uout(i,UMX  ) * fac(i)
          
          do n = UFS, UFS+nspec-1
             uout(i,n) = uout(i,n) * fac(i)
          end do
          do n = UFX, UFX+naux-1
             uout(i,n) = uout(i,n) * fac(i)
          end do
          do n = UFA, UFA+nadv-1
             uout(i,n) = uout(i,n) * fac(i)
          end do
          
       end if
       
       final_mass = final_mass + uout(i,URHO)
       final_eint = final_eint + uout(i,UEINT)
       final_eden = final_eden + uout(i,UEDEN)
       
    enddo
    
    deallocate(fac)
    
    mass_added = mass_added + (final_mass - initial_mass)
    eint_added = eint_added + (final_eint - initial_eint)
    eden_added = eden_added + (final_eden - initial_eden)
    
  end subroutine enforce_minimum_density

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine normalize_new_species(u,u_l1,u_h1,lo,hi)

    use network, only : nspec
    use meth_params_module, only : NVAR, URHO, UFS
    use bl_constants_module    

    implicit none

    integer          :: lo(1), hi(1)
    integer          :: u_l1,u_h1
    double precision :: u(u_l1:u_h1,NVAR)
    
    ! Local variables
    integer          :: i,n
    double precision :: fac,sum
    
    do i = lo(1),hi(1)
       sum = ZERO
       do n = UFS, UFS+nspec-1
          sum = sum + u(i,n)
       end do
       if (sum .ne. ZERO) then
          fac = u(i,URHO) / sum
       else
          fac = ONE
       end if
       do n = UFS, UFS+nspec-1
          u(i,n) = u(i,n) * fac
       end do
    end do
    
  end subroutine normalize_new_species
  
end module advection_module
