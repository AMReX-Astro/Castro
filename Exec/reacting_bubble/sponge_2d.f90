module sponge_module
  
  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_l2,uout_h1,uout_h2,lo,hi, &
                    t,dt,dx,dy,domlo,domhi)

    use bl_constants_module, only : M_PI
    use meth_params_module , only : NVAR, URHO, UMX, UMY, UEDEN

    integer          :: lo(2), hi(2), domlo(2), domhi(2)
    integer          :: uout_l1,uout_l2,uout_h1,uout_h2
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision :: t,dt
    double precision :: dx, dy

    integer          :: i,j
    double precision :: rho, ke_old, ke_new, fac, sponge_mult, sponge_weighting

    sponge_weighting = 1.d5

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          rho = uout(i,j,URHO)
          if (rho .le. 1.e6) then

             if (rho < 1.e4) then
                fac =  1.d0
             else 
                fac =  0.5d0*(1.d0 - cos(M_PI*(rho - 100000.0)/(1000.0 - 100000.0)))
             end if

             sponge_mult = 1.d0 / (1.d0 + fac * sponge_weighting * dt)

             ke_old = 0.5d0 * ( uout(i,j,UMX)**2 + uout(i,j,UMY)**2 )/ rho

             uout(i,j,UMX  ) = uout(i,j,UMX  ) * sponge_mult
             uout(i,j,UMY  ) = uout(i,j,UMY  ) * sponge_mult
             
             ke_new = 0.5d0 * ( uout(i,j,UMX)**2 + uout(i,j,UMY)**2 )/ rho
             
             uout(i,j,UEDEN) = uout(i,j,UEDEN) + (ke_new-ke_old)
             
          end if
          
       enddo
    enddo
    
  end subroutine sponge

end module sponge_module

  
