module hybrid_advection_module

  implicit none

contains

  ! Takes the initial linear momentum data in a state and converts it
  ! to the hybrid momenta.

  subroutine init_hybrid_momentum(lo, hi, state, s_lo, s_hi) bind(C,name='init_hybrid_momentum')

    use meth_params_module, only: NVAR, UMR, UMP, UMX, UMZ
    use castro_util_module, only: position
    use prob_params_module, only: center

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer          :: i, j, k
    double precision :: loc(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             state(i,j,k,UMR:UMP) = linear_to_hybrid_momentum(loc, state(i,j,k,UMX:UMZ))

          enddo
       enddo
    enddo

  end subroutine init_hybrid_momentum



  ! Convert a linear momentum into a "hybrid" momentum that has
  ! an angular momentum component.

  function linear_to_hybrid_momentum(loc, mom_in) result(mom_out)

    implicit none

    double precision :: loc(3), mom_in(3), mom_out(3)

    double precision :: R, mom(3)

    R = sqrt( loc(1)**2 + loc(2)**2 )

    mom = mom_in
    
    mom_out(1) = mom(1) * (loc(1) / R) + mom(2) * (loc(2) / R)
    mom_out(2) = mom(2) * loc(1)       - mom(1) * loc(2)
    mom_out(3) = mom(3)

  end function linear_to_hybrid_momentum



  ! Convert a "hybrid" momentum into a linear momentum.

  function hybrid_to_linear_momentum(loc, mom_in) result(mom_out)

    implicit none

    double precision :: loc(3), mom_in(3), mom_out(3)

    double precision :: R, mom(3)

    mom = mom_in

    R = sqrt( loc(1)**2 + loc(2)**2 )
    
    mom_out(1) = mom(1) * (loc(1) / R)    - mom(2) * (loc(2) / R**2)
    mom_out(2) = mom(2) * (loc(1) / R**2) + mom(1) * (loc(2) / R)
    mom_out(3) = mom(3)
    
  end function hybrid_to_linear_momentum



  ! Update momentum to account for source term

  subroutine add_hybrid_momentum_source(loc, mom, source)

    implicit none

    double precision :: loc(3), mom(3), source(3)

    double precision :: R

    R = sqrt( loc(1)**2 + loc(2)**2 )

    mom(1) = mom(1) - source(1) * (loc(1) / R) - source(2) * (loc(2) / R)
    mom(2) = mom(2) + source(1) * loc(2) - source(2) * loc(1)
    mom(3) = mom(3) + source(3)

  end subroutine add_hybrid_momentum_source



  ! Update state to account for hybrid advection.

  subroutine hybrid_update(lo, hi, dx, dt, &
                           sold, sold_lo, sold_hi, &
                           snew, snew_lo, snew_hi, &
                           q1, q1_lo, q1_hi, &
                           q2, q2_lo, q2_hi, &
                           q3, q3_lo, q3_hi)

    use bl_constants_module, only: ZERO, HALF, ONE
    use meth_params_module, only: URHO, UMR, UML, UMP, UMX, UMZ, NVAR, QVAR, NGDNV, GDRHO, GDU, GDV, GDW, GDPRES, hybrid_hydro
    use prob_params_module, only: center
    use castro_util_module, only: position, area, volume
    use amrinfo_module, only: amr_level
    use prob_params_module, only: domlo_level, domhi_level, physbc_lo, physbc_hi, Symmetry, SlipWall, NoSlipWall
    
    implicit none

    integer :: lo(3), hi(3)
    integer :: sold_lo(3), sold_hi(3)
    integer :: snew_lo(3), snew_hi(3)
    integer :: q1_lo(3), q1_hi(3)
    integer :: q2_lo(3), q2_hi(3)
    integer :: q3_lo(3), q3_hi(3)

    double precision :: dx(3), dt
    double precision :: sold(sold_lo(1):sold_hi(1),sold_lo(2):sold_hi(2),sold_lo(3):sold_hi(3),NVAR)
    double precision :: snew(snew_lo(1):snew_hi(1),snew_lo(2):snew_hi(2),snew_lo(3):snew_hi(3),NVAR)
    double precision :: q1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),NGDNV)
    double precision :: q2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),NGDNV)
    double precision :: q3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),NGDNV)

    double precision :: flux1(q1_lo(1):q1_hi(1),q1_lo(2):q1_hi(2),q1_lo(3):q1_hi(3),3)
    double precision :: flux2(q2_lo(1):q2_hi(1),q2_lo(2):q2_hi(2),q2_lo(3):q2_hi(3),3)
    double precision :: flux3(q3_lo(1):q3_hi(1),q3_lo(2):q3_hi(2),q3_lo(3):q3_hi(3),3)

    integer          :: i, j, k

    double precision :: loc(3), R
    double precision :: hybrid_mom(3), linear_mom(3), rho_new, rho_old

    logical          :: special_bnd_lo, special_bnd_hi
    double precision :: bnd_fac
    integer          :: domlo(3), domhi(3)

    ! First, construct the fluxes.

    ! Account for boundaries where we want to explicitly zero the fluxes.

    domlo = domlo_level(:, amr_level)
    domhi = domhi_level(:, amr_level)

    special_bnd_lo = (physbc_lo(1) .eq. Symmetry &
         .or.         physbc_lo(1) .eq. SlipWall &
         .or.         physbc_lo(1) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(1) .eq. Symmetry &
         .or.         physbc_hi(1) .eq. SlipWall &
         .or.         physbc_hi(1) .eq. NoSlipWall)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)+1

             bnd_fac = ONE

             if ( i .eq. domlo(1)   .and. special_bnd_lo .or. &
                  i .eq. domhi(1)+1 .and. special_bnd_hi ) then
                bnd_fac = ZERO
             endif

             loc = position(i,j,k,ccx=.false.) - center

             linear_mom = q1(i,j,k,GDRHO) * q1(i,j,k,GDU:GDW)

             hybrid_mom = linear_to_hybrid_momentum(loc, linear_mom)

             flux1(i,j,k,1) = hybrid_mom(1) * q1(i,j,k,GDU)
             flux1(i,j,k,2) = hybrid_mom(2) * q1(i,j,k,GDU) + loc(2) * q1(i,j,k,GDPRES)
             flux1(i,j,k,3) = hybrid_mom(3) * q1(i,j,k,GDU)

             flux1(i,j,k,:) = flux1(i,j,k,:) * area(i,j,k,1) * dt * bnd_fac

          enddo
       enddo
    enddo

    special_bnd_lo = (physbc_lo(2) .eq. Symmetry &
         .or.         physbc_lo(2) .eq. SlipWall &
         .or.         physbc_lo(2) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(2) .eq. Symmetry &
         .or.         physbc_hi(2) .eq. SlipWall &
         .or.         physbc_hi(2) .eq. NoSlipWall)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)+1

          bnd_fac = ONE

          if ( j .eq. domlo(2)   .and. special_bnd_lo .or. &
               j .eq. domhi(2)+1 .and. special_bnd_hi ) then
             bnd_fac = ZERO
          endif

          do i = lo(1), hi(1)

             loc = position(i,j,k,ccy=.false.) - center

             linear_mom = q2(i,j,k,GDRHO) * q2(i,j,k,GDU:GDW)

             hybrid_mom = linear_to_hybrid_momentum(loc, linear_mom)

             flux2(i,j,k,1) = hybrid_mom(1) * q2(i,j,k,GDV)
             flux2(i,j,k,2) = hybrid_mom(2) * q2(i,j,k,GDV) - loc(1) * q2(i,j,k,GDPRES)
             flux2(i,j,k,3) = hybrid_mom(3) * q2(i,j,k,GDV)

             flux2(i,j,k,:) = flux2(i,j,k,:) * area(i,j,k,2) * dt * bnd_fac

          enddo
       enddo
    enddo

    special_bnd_lo = (physbc_lo(3) .eq. Symmetry &
         .or.         physbc_lo(3) .eq. SlipWall &
         .or.         physbc_lo(3) .eq. NoSlipWall)
    special_bnd_hi = (physbc_hi(3) .eq. Symmetry &
         .or.         physbc_hi(3) .eq. SlipWall &
         .or.         physbc_hi(3) .eq. NoSlipWall)

    do k = lo(3), hi(3)+1

       bnd_fac = ONE

       if ( k .eq. domlo(3)   .and. special_bnd_lo .or. &
            k .eq. domhi(3)+1 .and. special_bnd_hi ) then
          bnd_fac = ZERO
       endif

       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k,ccz=.false.) - center

             linear_mom = q3(i,j,k,GDRHO) * q3(i,j,k,GDU:GDW)

             hybrid_mom = linear_to_hybrid_momentum(loc, linear_mom)

             flux3(i,j,k,1) = hybrid_mom(1) * q3(i,j,k,GDW)
             flux3(i,j,k,2) = hybrid_mom(2) * q3(i,j,k,GDW)
             flux3(i,j,k,3) = hybrid_mom(3) * q3(i,j,k,GDW)

             flux3(i,j,k,:) = flux3(i,j,k,:) * area(i,j,k,3) * dt * bnd_fac

          enddo
       enddo
    enddo



    ! Now update state

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             rho_old = sold(i,j,k,URHO)
             rho_new = snew(i,j,k,URHO)

             snew(i,j,k,UMR:UMP) = sold(i,j,k,UMR:UMP) &
                                 + ( flux1(i,j,k,:) - flux1(i+1,j,k,:) &
                                 +   flux2(i,j,k,:) - flux2(i,j+1,k,:) &
                                 +   flux3(i,j,k,:) - flux3(i,j,k+1,:) ) / volume(i,j,k)

             ! Add the time-centered source term to the radial momentum

             snew(i,j,k,UMR) = snew(i,j,k,UMR) &
                             + dt * ( - (loc(1) / R) * (q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES)) / dx(1) &
                                      - (loc(2) / R) * (q2(i,j+1,k,GDPRES) - q2(i,j,k,GDPRES)) / dx(2) &
                                      + (HALF * (snew(i,j,k,UML) + sold(i,j,k,UML)))**2 / &
                                      ( (HALF * (rho_new + rho_old)) * R**3) )

             ! If we're doing the hybrid advection scheme, update the momenta accordingly.

             if (hybrid_hydro .eq. 1) then

                snew(i,j,k,UMX:UMZ) = hybrid_to_linear_momentum(loc, snew(i,j,k,UMR:UMP))

             endif

          enddo
       enddo
    enddo

  end subroutine hybrid_update

end module hybrid_advection_module

