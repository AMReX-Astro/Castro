module advection_module

  implicit none

  private

  public umeth1d, ctoprim, consup

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
    integer ilo,ihi
    double precision dx, dt
    double precision     q(   qd_l1:qd_h1,QVAR)
    double precision  gamc(   qd_l1:qd_h1)
    double precision flatn(   qd_l1:qd_h1)
    double precision  csml(   qd_l1:qd_h1)
    double precision     c(   qd_l1:qd_h1)
    double precision  flux(fd_l1   :fd_h1,NVAR)
    double precision  srcQ(src_l1  :src_h1,QVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision ugdnv(ugdnv_l1:ugdnv_h1)
    double precision dloga(dloga_l1:dloga_h1)
    
    ! Left and right state arrays (edge centered, cell centered)
    double precision, allocatable:: dq(:,:),  qm(:,:),   qp(:,:)

    ! Work arrays to hold 3 planes of riemann state and conservative fluxes
    allocate ( qm(ilo-1:ihi+1,QVAR))
    allocate ( qp(ilo-1:ihi+1,QVAR))

    ! Trace to edges w/o transverse flux correction terms
    if (ppm_type .gt. 0) then
       call trace_ppm(q,c,flatn,gamc,qd_l1,qd_h1, &
                      dloga,dloga_l1,dloga_h1, &
                      srcQ,src_l1,src_h1, &
                      qm,qp,ilo-1,ihi+1, &
                      ilo,ihi,domlo,domhi,dx,dt)
    else
       allocate ( dq(ilo-1:ihi+1,QVAR))

       call trace(q,dq,c,flatn,qd_l1,qd_h1, &
                  dloga,dloga_l1,dloga_h1, &
                  srcQ,src_l1,src_h1, &
                  qm,qp,ilo-1,ihi+1, &
                  ilo,ihi,domlo,domhi,dx,dt)

       deallocate(dq)
    end if

    ! Solve Riemann problem, compute xflux from improved predicted states 
    call cmpflx(lo, hi, domlo, domhi, &
                qm, qp, ilo-1,ihi+1, &
                flux ,  fd_l1, fd_h1, &
                pgdnv,pgdnv_l1,pgdnv_h1, &
                ugdnv,ugdnv_l1,ugdnv_h1, &
                gamc, csml,c,qd_l1,qd_h1,ilo,ihi)

    deallocate (qm,qp)

  end subroutine umeth1d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi,uin,uin_l1,uin_h1, &
                     q,c,gamc,csml,flatn,q_l1,q_h1,&
                     src,src_l1,src_h1, &
                     srcQ,srQ_l1,srQ_h1, &
                     courno,dx,dt,ngp,ngf)
    
    ! Will give primitive variables on lo-ngp:hi+ngp, and flatn on
    ! lo-ngf:hi+ngf if iflaten=1.  Declared dimensions of
    ! q,c,gamc,csml,flatn are given by DIMS(q).  This declared
    ! region is assumed to encompass lo-ngp:hi+ngp.  Also, uflaten
    ! call assumes ngp>=ngf+3 (ie, primitve data is used by the
    ! routine that computes flatn).

    use network, only : nspec, naux
    use eos_module
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS, UFX, &
                                   QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QTEMP, QGAME, &
                                   QFS, QFX, &
                                   npassive, upass_map, qpass_map, small_temp, allow_negative_energy, use_flattening, &
                                   dual_energy_eta1
    use flatten_module
    use bl_constants_module

    implicit none

    double precision, parameter:: small = 1.d-8

    integer          :: lo(1), hi(1)
    integer          :: uin_l1,uin_h1
    integer          :: q_l1,q_h1
    integer          ::  src_l1,src_h1
    integer          ::  srQ_l1,srQ_h1
    double precision ::   uin(uin_l1:uin_h1,NVAR)
    double precision ::     q(  q_l1:  q_h1,QVAR)
    double precision ::     c(  q_l1:  q_h1)
    double precision ::  gamc(  q_l1:  q_h1)
    double precision ::  csml(  q_l1:  q_h1)
    double precision :: flatn(  q_l1:  q_h1)
    double precision ::   src(src_l1:src_h1,NVAR)
    double precision ::  srcQ(srQ_l1:srQ_h1,QVAR)
    double precision :: dx, dt, courno
    
    integer          :: i
    integer          :: ngp, ngf, loq(1), hiq(1)
    integer          :: n, nq, ipassive
    double precision :: courx, courmx
    double precision :: kineng
    
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

       ! Load passively advected qunatities into q
       do ipassive = 1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)
          q(i,nq) = uin(i,n)/q(i,QRHO)
       enddo

       q(i,QU  ) = uin(i,UMX )/uin(i,URHO)

       ! Get the internal energy, which we'll use for determining the pressure.
       ! We use a dual energy formalism. If (E - K) < eta1 and eta1 is suitably small, 
       ! then we risk serious numerical truncation error in the internal energy.
       ! Therefore we'll use the result of the separately updated internal energy equation.
       ! Otherwise, we'll set e = E - K.

       kineng = HALF * q(i,QRHO) * sum(q(i,QU:QW)**2)

       if ( (uin(i,UEDEN) - kineng) / uin(i,UEDEN) .gt. dual_energy_eta1) then
          q(i,QREINT) = (uin(i,UEDEN) - kineng) / q(i,QRHO)
       else
          q(i,QREINT) = uin(i,UEINT) / q(i,QRHO)
       endif

       q(i,QTEMP ) = uin(i,UTEMP)

    enddo
    
    ! Get gamc, p, T, c, csml using q state
    do i = loq(1), hiq(1)
       
       eos_state % T   = q(i,QTEMP)
       eos_state % rho = q(i,QRHO)
       eos_state % xn  = q(i,QFS:QFS+nspec-1)
       eos_state % aux = q(i,QFX:QFX+naux-1)
       
       ! If necessary, reset the energy using small_temp
       if (allow_negative_energy .eq. 0 .and. q(i,QREINT) .le. ZERO) then
          
          q(i,QTEMP) = small_temp
          eos_state % T = q(i,QTEMP)
          
          call eos(eos_input_rt, eos_state)
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
       
       call eos(eos_input_re, eos_state)
       
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
       q(i,QGAME) = q(i,QPRES)/q(i,QREINT) + ONE
    enddo

    srcQ = ZERO

    ! compute srcQ terms
    do i = loq(1), hiq(1)
       srcQ(i,QRHO   ) = src(i,URHO)
       srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
       srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + HALF * q(i,QU)**2 * srcQ(i,QRHO)
       srcQ(i,QPRES  ) = dpde(i) * (srcQ(i,QREINT) - q(i,QREINT)*srcQ(i,QRHO)/q(i,QRHO)) / q(i,QRHO) + &
            dpdrho(i) * srcQ(i,QRHO)! - &
       !              sum(dpdX_er(i,:)*(src(i,UFS:UFS+nspec-1) - q(i,QFS:QFS+nspec-1)*srcQ(i,QRHO)))/q(i,QRHO)
       
       do ipassive=1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)
          srcQ(i,nq) = ( src(i,n) - q(i,nq) * srcQ(i,QRHO) ) / q(i,QRHO)
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
       call uflaten((/ loq(1), 0, 0 /), (/ hiq(1), 0, 0 /), &
                    q(q_l1,QPRES), &
                    q(q_l1,QU), &
                    q(q_l1,QV), &
                    q(q_l1,QW), &
                    flatn,(/ q_l1, 0, 0 /), (/ q_h1, 0, 0 /))
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
                    flux, flux_l1, flux_h1, &
                    area,area_l1,area_h1, &
                    vol,vol_l1,vol_h1, &
                    div,pdivu,lo,hi,dx,dt,E_added_flux, &
                    xmom_added_flux,ymom_added_flux,zmom_added_flux, &
                    verbose)

    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, UFS, UFX
    use bl_constants_module
    use advection_util_module, only: normalize_species_fluxes

    integer lo(1), hi(1)
    integer   uin_l1,  uin_h1
    integer  uout_l1, uout_h1
    integer pgdnv_l1,pgdnv_h1
    integer   src_l1,  src_h1
    integer  flux_l1, flux_h1
    integer  area_l1, area_h1
    integer   vol_l1,  vol_h1
    integer verbose
    double precision   uin(uin_l1:uin_h1,NVAR)
    double precision  uout(uout_l1:uout_h1,NVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision   src(  src_l1:  src_h1,NVAR)
    double precision  flux( flux_l1: flux_h1,NVAR)
    double precision  area( area_l1: area_h1)
    double precision    vol(vol_l1:vol_h1)
    double precision    div(lo(1):hi(1)+1)
    double precision  pdivu(lo(1):hi(1)  )
    double precision dx, dt
    double precision E_added_flux
    double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
    
    integer          :: i, n
    double precision :: div1, dpdx
    
    ! Normalize the species fluxes
    call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)
    
    do n = 1, NVAR
       if ( n.eq.UTEMP .or. n.eq.UMY .or. n.eq.UMZ ) then
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
       if ( n.eq.UTEMP ) then
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

    ! Add up some diagnostic quantities. Note that these are volumetric sums
    ! so we are not dividing by the cell volume.

    if (verbose .eq. 1) then

       do i = lo(1), hi(1)
                      
          xmom_added_flux = xmom_added_flux + ( flux(i,UMX) - flux(i+1,UMX) )
          ymom_added_flux = ymom_added_flux + ( flux(i,UMY) - flux(i+1,UMY) )
          zmom_added_flux = zmom_added_flux + ( flux(i,UMZ) - flux(i+1,UMZ) )
          E_added_flux = E_added_flux + ( flux(i,UEDEN) - flux(i+1,UEDEN) )

       enddo

    endif
    
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
    
    do i = lo(1),hi(1)+1
       flux(i,UMX) = flux(i,UMX) + dt*area(i)*pgdnv(i)
    enddo
    
  end subroutine consup
  
end module advection_module
