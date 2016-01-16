module hybrid_advection_module

  implicit none

contains

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

  subroutine add_momentum_source(loc, mom, source)

    use meth_params_module, only: hybrid_hydro

    implicit none

    double precision :: loc(3), mom(3), source(3)

    double precision :: hybrid_mom(3), R

    R = sqrt( loc(1)**2 + loc(2)**2 )
    
    if (hybrid_hydro .eq. 1) then

       hybrid_mom = linear_to_hybrid_momentum(loc, mom)

       hybrid_mom(1) = hybrid_mom(1) - source(1) * (loc(1) / R) - source(2) * (loc(2) / R)
       hybrid_mom(2) = hybrid_mom(2) + source(1) * loc(2) - source(2) * loc(1)
       hybrid_mom(3) = hybrid_mom(3) + source(3)
                  
       mom = hybrid_to_linear_momentum(loc, hybrid_mom)    

    else
       
       mom = mom + source
       
    endif

  end subroutine add_momentum_source



  ! Update state to account for hybrid advection.

  subroutine hybrid_update(lo, hi, dx, dt, &
                           sold, sold_lo, sold_hi, &
                           snew, snew_lo, snew_hi, &
                           q1, q1_lo, q1_hi, &
                           q2, q2_lo, q2_hi, &
                           q3, q3_lo, q3_hi)

    use bl_constants_module, only: HALF
    use meth_params_module, only: URHO, UMX, UMZ, NVAR, QVAR, NGDNV, GDRHO, GDU, GDV, GDW, GDPRES
    use prob_params_module, only: center
    use castro_util_module, only: position, area, volume
    
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

    integer :: i, j, k    
    
    double precision :: loc(3), R
    double precision :: hybrid_mom_new(3), hybrid_mom_old(3), linear_mom(3), rho_new, rho_old

    ! First, construct the fluxes.
    
    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)+1

             loc = position(i,j,k,ccx=.false.) - center

             linear_mom = q1(i,j,k,GDRHO) * q1(i,j,k,GDU:GDW)

             hybrid_mom_old = linear_to_hybrid_momentum(loc, linear_mom)
                
             flux1(i,j,k,1) = hybrid_mom_old(1) * q1(i,j,k,GDU)
             flux1(i,j,k,2) = hybrid_mom_old(2) * q1(i,j,k,GDU) + loc(2) * q1(i,j,k,GDPRES)
             flux1(i,j,k,3) = hybrid_mom_old(3) * q1(i,j,k,GDU)

             flux1(i,j,k,:) = flux1(i,j,k,:) * area(i,j,k,1) * dt
             
          enddo
       enddo
    enddo

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)+1
          do i = lo(1),hi(1)

             loc = position(i,j,k,ccy=.false.) - center

             linear_mom = q2(i,j,k,GDRHO) * q2(i,j,k,GDU:GDW)

             hybrid_mom_old = linear_to_hybrid_momentum(loc, linear_mom)

             flux2(i,j,k,1) = hybrid_mom_old(1) * q2(i,j,k,GDV)
             flux2(i,j,k,2) = hybrid_mom_old(2) * q2(i,j,k,GDV) - loc(1) * q2(i,j,k,GDPRES)
             flux2(i,j,k,3) = hybrid_mom_old(3) * q2(i,j,k,GDV)

             flux2(i,j,k,:) = flux2(i,j,k,:) * area(i,j,k,2) * dt

          enddo
       enddo
    enddo

    do k = lo(3),hi(3)+1
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             loc = position(i,j,k,ccz=.false.) - center

             linear_mom = q3(i,j,k,GDRHO) * q3(i,j,k,GDU:GDW)

             hybrid_mom_old = linear_to_hybrid_momentum(loc, linear_mom)

             flux3(i,j,k,1) = hybrid_mom_old(1) * q3(i,j,k,GDW)
             flux3(i,j,k,2) = hybrid_mom_old(2) * q3(i,j,k,GDW)
             flux3(i,j,k,3) = hybrid_mom_old(3) * q3(i,j,k,GDW)

             flux3(i,j,k,:) = flux3(i,j,k,:) * area(i,j,k,3) * dt

          enddo
       enddo
    enddo



    ! Now update state
    
    do k = lo(3), hi(3)       
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)
             
             loc = position(i,j,k) - center

             R = sqrt( loc(1)**2 + loc(2)**2 )

             rho_old = sold(i,j,k,URHO)
             rho_new = snew(i,j,k,URHO)

             hybrid_mom_old(:) = linear_to_hybrid_momentum(loc, sold(i,j,k,UMX:UMZ))
             
             hybrid_mom_new(:) = hybrid_mom_old(:) &
                               + ( flux1(i,j,k,:) - flux1(i+1,j,k,:) &
                                 + flux2(i,j,k,:) - flux2(i,j+1,k,:) &
                                 + flux3(i,j,k,:) - flux3(i,j,k+1,:) ) / volume(i,j,k)

             ! Add the time-centered source term to the radial momentum

             hybrid_mom_new(1) = hybrid_mom_new(1) &
                               + dt * ( - (loc(1) / R) * (q1(i+1,j,k,GDPRES) - q1(i,j,k,GDPRES)) / dx(1) &
                                        - (loc(2) / R) * (q2(i,j+1,k,GDPRES) - q2(i,j,k,GDPRES)) / dx(2) &
                                        + (HALF * (hybrid_mom_new(2) + hybrid_mom_old(2)))**2 / &
                                        ( (HALF * (rho_new + rho_old)) * R**3) )

             snew(i,j,k,UMX:UMZ) = hybrid_to_linear_momentum(loc, hybrid_mom_new)

          enddo
       enddo
    enddo
    
  end subroutine hybrid_update
  
end module hybrid_advection_module

