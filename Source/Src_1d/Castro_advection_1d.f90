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
                     q,flatn,qd_l1,qd_h1, &
                     srcQ,src_l1,src_h1, &
                     ilo,ihi,dx,dt, &
                     flux ,   fd_l1,   fd_h1, &
                     pgdnv,pgdnv_l1,pgdnv_h1, &
                     ugdnv,ugdnv_l1,ugdnv_h1, &
                     dloga,dloga_l1,dloga_h1)

    use meth_params_module, only : QVAR, NVAR, ppm_type, QC, QCSML, QGAMC
    use riemann_module, only : cmpflx
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
    double precision flatn(   qd_l1:qd_h1)
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
       call trace_ppm(q,q(:,QC),flatn,q(:,QGAMC),qd_l1,qd_h1, &
                      dloga,dloga_l1,dloga_h1, &
                      srcQ,src_l1,src_h1, &
                      qm,qp,ilo-1,ihi+1, &
                      ilo,ihi,domlo,domhi,dx,dt)
    else
       allocate ( dq(ilo-1:ihi+1,QVAR))

       call trace(q,dq,q(:,QC),flatn,qd_l1,qd_h1, &
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
                q(:,QGAMC),q(:,QCSML),q(:,QC),qd_l1,qd_h1,ilo,ihi)

    deallocate (qm,qp)

  end subroutine umeth1d

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

  subroutine ctoprim(lo,hi, &
                     uin,uin_l1,uin_h1, &
                     q,q_l1,q_h1,&
                     src,src_l1,src_h1, &
                     srcQ,srQ_l1,srQ_h1, &
                     dx,dt)

    use network, only : nspec, naux
    use eos_module, only : eos
    use eos_type_module, only : eos_t, eos_input_re
    use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, &
                                   QVAR, QRHO, QU, QV, QW, QREINT, QPRES, QTEMP, QGAME, &
                                   QFS, QFX, &
                                   QC, QCSML, QGAMC, QDPDR, QDPDE, &
                                   npassive, upass_map, qpass_map, allow_negative_energy, &
                                   dual_energy_eta1
    use bl_constants_module, only : ZERO, HALF, ONE

    implicit none

    double precision, parameter:: small = 1.d-8

    integer          :: lo(1), hi(1)
    integer          :: uin_l1,uin_h1
    integer          :: q_l1,q_h1
    integer          ::  src_l1,src_h1
    integer          ::  srQ_l1,srQ_h1
    double precision ::   uin(uin_l1:uin_h1,NVAR)
    double precision ::     q(  q_l1:  q_h1,QVAR)
    double precision ::   src(src_l1:src_h1,NVAR)
    double precision ::  srcQ(srQ_l1:srQ_h1,QVAR)
    double precision :: dx, dt

    integer          :: i
    integer          :: n, nq, ipassive
    double precision :: kineng

    type (eos_t) :: eos_state

    ! Make q (all but p), except put e in slot for rho.e, fix after
    ! eos call The temperature is used as an initial guess for the eos
    ! call and will be overwritten
    do i = lo(1),hi(1)

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

    if (allow_negative_energy .eq. 0) eos_state % reset = .true.

    ! Get gamc, p, T, c, csml using q state
    do i = lo(1), hi(1)
       eos_state % T   = q(i,QTEMP)
       eos_state % rho = q(i,QRHO)
       eos_state % e   = q(i,QREINT)
       eos_state % xn  = q(i,QFS:QFS+nspec-1)
       eos_state % aux = q(i,QFX:QFX+naux-1)

       call eos(eos_input_re, eos_state)

       q(i,QTEMP)  = eos_state % T
       q(i,QREINT) = eos_state % e
       q(i,QPRES)  = eos_state % p
       q(i,QDPDR)  = eos_state % dpdr_e
       q(i,QDPDE)  = eos_state % dpde
       q(i,QC)     = eos_state % cs
       q(i,QGAMC)  = eos_state % gam1
       q(i,QCSML)  = max(small, small * q(i,QC))
    end do

    ! Make this "rho e" instead of "e"
    do i = lo(1),hi(1)
       q(i,QREINT ) = q(i,QREINT )*q(i,QRHO)
       q(i,QGAME) = q(i,QPRES)/q(i,QREINT) + ONE
    enddo

    srcQ = ZERO

    ! compute srcQ terms
    do i = lo(1), hi(1)
       srcQ(i,QRHO   ) = src(i,URHO)
       srcQ(i,QU     ) = (src(i,UMX) - q(i,QU) * srcQ(i,QRHO)) / q(i,QRHO)
       srcQ(i,QREINT ) = src(i,UEDEN) - q(i,QU) * src(i,UMX) + HALF * q(i,QU)**2 * srcQ(i,QRHO)
       srcQ(i,QPRES  ) = q(i,QDPDE) * (srcQ(i,QREINT) - q(i,QREINT)*srcQ(i,QRHO)/q(i,QRHO)) / q(i,QRHO) + &
                         q(i,QDPDR) * srcQ(i,QRHO)

       do ipassive=1, npassive
          n  = upass_map(ipassive)
          nq = qpass_map(ipassive)
          srcQ(i,nq) = ( src(i,n) - q(i,nq) * srcQ(i,QRHO) ) / q(i,QRHO)
       end do

    end do

  end subroutine ctoprim

! :::
! ::: ------------------------------------------------------------------
! :::

  subroutine consup( &
                    uin,  uin_l1,  uin_h1, &
                    uout, uout_l1 ,uout_h1, &
                    update,updt_l1,updt_h1, &
                    pgdnv,pgdnv_l1,pgdnv_h1, &
                    src,  src_l1,  src_h1, &
                    flux, flux_l1, flux_h1, &
                    area,area_l1,area_h1, &
                    vol,vol_l1,vol_h1, &
                    div,pdivu,lo,hi,dx,dt,mass_added_flux,E_added_flux, &
                    xmom_added_flux,ymom_added_flux,zmom_added_flux, &
                    mass_lost,xmom_lost,ymom_lost,zmom_lost, &
                    eden_lost,xang_lost,yang_lost,zang_lost, &
                    verbose)

    use eos_module
    use meth_params_module, only : difmag, NVAR, URHO, UMX, UMY, UMZ, &
                                   UEDEN, UEINT, UTEMP, track_grid_losses
    use bl_constants_module
    use advection_util_1d_module, only: normalize_species_fluxes
    use prob_params_module, only : domlo_level, domhi_level, center
    use castro_util_module, only : position, linear_to_angular_momentum
    use amrinfo_module, only : amr_level

    integer lo(1), hi(1)
    integer   uin_l1,  uin_h1
    integer  uout_l1, uout_h1
    integer  updt_l1, updt_h1
    integer pgdnv_l1,pgdnv_h1
    integer   src_l1,  src_h1
    integer  flux_l1, flux_h1
    integer  area_l1, area_h1
    integer   vol_l1,  vol_h1
    integer verbose
    double precision   uin(uin_l1:uin_h1,NVAR)
    double precision  uout(uout_l1:uout_h1,NVAR)
    double precision update(updt_l1:updt_h1,NVAR)
    double precision pgdnv(pgdnv_l1:pgdnv_h1)
    double precision   src(  src_l1:  src_h1,NVAR)
    double precision  flux( flux_l1: flux_h1,NVAR)
    double precision  area( area_l1: area_h1)
    double precision    vol(vol_l1:vol_h1)
    double precision    div(lo(1):hi(1)+1)
    double precision  pdivu(lo(1):hi(1)  )
    double precision dx, dt
    double precision E_added_flux, mass_added_flux
    double precision xmom_added_flux, ymom_added_flux, zmom_added_flux
    double precision mass_lost, xmom_lost, ymom_lost, zmom_lost
    double precision eden_lost, xang_lost, yang_lost, zang_lost
    
    integer          :: i, j, k, n
    double precision :: div1
    integer          :: domlo(3), domhi(3)
    double precision :: loc(3), ang_mom(3)
    
    ! Normalize the species fluxes.

    call normalize_species_fluxes(flux,flux_l1,flux_h1,lo,hi)
    
    do n = 1, NVAR
       if ( n == UTEMP ) then
          flux(:,n) = ZERO
       else if ( n == UMY ) then
          flux(:,n) = ZERO
       else if ( n == UMZ ) then
          flux(:,n) = ZERO
       else
          do i = lo(1),hi(1)+1
             div1 = difmag*min(ZERO,div(i))
             flux(i,n) = flux(i,n) &
                  + dx*div1*(uin(i,n) - uin(i-1,n))
             flux(i,n) = area(i) * flux(i,n)
          enddo
       endif
    enddo

    ! Fill the update array.

    do n = 1, NVAR
       do i = lo(1), hi(1)

          update(i,n) = update(i,n) + ( flux(i,n) - flux(i+1,n) ) / vol(i)

          ! Add p div(u) source term to (rho e).

          if (n == UEINT) then

             update(i,n) = update(i,n) - pdivu(i)

          endif

       enddo
    enddo

    ! Add gradp term to momentum equation.

    do i = lo(1),hi(1)

       update(i,UMX) = update(i,UMX) - ( pgdnv(i+1) - pgdnv(i) ) / dx

    enddo

    ! Add up some diagnostic quantities. Note that we are not dividing by the cell volume.

    if (verbose .eq. 1) then

       do i = lo(1), hi(1)

          mass_added_flux = mass_added_flux + dt * ( flux(i,URHO ) - flux(i+1,URHO ) )
          xmom_added_flux = xmom_added_flux + dt * ( flux(i,UMX  ) - flux(i+1,UMX  ) )
          ymom_added_flux = ymom_added_flux + dt * ( flux(i,UMY  ) - flux(i+1,UMY  ) )
          zmom_added_flux = zmom_added_flux + dt * ( flux(i,UMZ  ) - flux(i+1,UMZ  ) )
          E_added_flux    = E_added_flux    + dt * ( flux(i,UEDEN) - flux(i+1,UEDEN) )

       enddo

    endif

    if (track_grid_losses .eq. 1) then

       domlo = domlo_level(:,amr_level)
       domhi = domhi_level(:,amr_level)

       j = 0
       k = 0

       if (lo(1) .le. domlo(1) .and. hi(1) .ge. domlo(1)) then

          i = domlo(1)

          loc = position(i,j,k,ccx=.false.)

          mass_lost = mass_lost - dt * flux(i,URHO)
          xmom_lost = xmom_lost - dt * flux(i,UMX)
          ymom_lost = ymom_lost - dt * flux(i,UMY)
          zmom_lost = zmom_lost - dt * flux(i,UMZ)
          eden_lost = eden_lost - dt * flux(i,UEDEN)

          ang_mom   = linear_to_angular_momentum(loc - center, dt * flux(i,UMX:UMZ))
          xang_lost = xang_lost - ang_mom(1)
          yang_lost = yang_lost - ang_mom(2)
          zang_lost = zang_lost - ang_mom(3)

       endif

       if (lo(1) .le. domhi(1) .and. hi(1) .ge. domhi(1)) then

          i = domhi(1) + 1

          loc = position(i,j,k,ccx=.false.)

          mass_lost = mass_lost + dt * flux(i,URHO)
          xmom_lost = xmom_lost + dt * flux(i,UMX)
          ymom_lost = ymom_lost + dt * flux(i,UMY)
          zmom_lost = zmom_lost + dt * flux(i,UMZ)
          eden_lost = eden_lost + dt * flux(i,UEDEN)

          ang_mom   = linear_to_angular_momentum(loc - center, dt * flux(i,UMX:UMZ))
          xang_lost = xang_lost + ang_mom(1)
          yang_lost = yang_lost + ang_mom(2)
          zang_lost = zang_lost + ang_mom(3)

       endif

    endif

    ! Scale the fluxes for the form we expect later in refluxing.

    do n = 1, NVAR
       do i = lo(1), hi(1)+1

          flux(i,n) = dt * flux(i,n)

          ! Correct the momentum flux with the grad p part.

          if (n == UMX) then
             flux(i,n) = flux(i,n) + dt * area(i) * pgdnv(i)
          endif

       enddo
    enddo

  end subroutine consup

end module advection_module
