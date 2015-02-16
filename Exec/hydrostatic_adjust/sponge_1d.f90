module sponge_module

  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_h1,lo,hi,t,dt,dx,domlo,domhi)

    use bl_constants_module, only : M_PI
    use meth_params_module , only : NVAR, URHO, UMX, UEDEN
    use probdata_module, only: sponge_start_density, sponge_width_factor, &
                               sponge_weighting

    implicit none
    integer          :: lo(1), hi(1), domlo(1), domhi(1)
    integer          :: uout_l1,uout_h1
    double precision :: uout(uout_l1:uout_h1,NVAR)
    double precision :: t,dt
    double precision :: dx
    
    integer          :: i
    double precision :: rho, ke_old, ke_new, fac, sponge_mult
    
    do i = lo(1),hi(1)
       
       rho = uout(i,URHO)
       if (rho <= sponge_start_density) then
          
          if (rho < sponge_start_density/sponge_width_factor) then
             fac =  1.d0
          else 
             fac =  0.5d0*(1.d0 - cos(M_PI*(sponge_start_density - rho)/ &
                  (sponge_start_density - sponge_start_density/sponge_width_factor)))
          end if
          
          sponge_mult = 1.d0 / (1.d0 + fac * sponge_weighting * dt)
          
          ke_old = 0.5d0 * uout(i,UMX)**2 / rho
          
          uout(i,UMX) = uout(i,UMX) * sponge_mult
          
          ke_new = 0.5d0 * uout(i,UMX)**2 / rho
          
          uout(i,UEDEN) = uout(i,UEDEN) + (ke_new-ke_old)
          
       end if
       
    enddo
    
  end subroutine sponge

end module sponge_module
