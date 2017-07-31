module sponge_module
  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_l2,&
                    uout_h1,uout_h2,lo,hi,t,dt, &
                    dx,dy,domlo,domhi, &
                    E_added,xmom_added,ymom_added) 

    use bl_constants_module, only : M_PI
    use meth_params_module , only : NVAR, URHO, UMX, UMY, UEDEN

    implicit none
    integer          :: lo(2), hi(2), domlo(2), domhi(2)
    integer          :: uout_l1,uout_l2,uout_h1,uout_h2
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision :: t,dt
    double precision :: dx,dy
    
    integer          :: i,j
    double precision :: rho, ke_old, ke_new, fac
    double precision :: sponge_mult, sponge_center_density, sponge_start_factor, sponge_kappa
    double precision :: sponge_start_density
    double precision :: E_added,xmom_added,ymom_added
    
    sponge_center_density = 4.d-8
    sponge_start_factor = 2.0d0
    sponge_kappa = 0.001d0
    
    sponge_start_density = sponge_center_density * sponge_start_factor
    
    do j = lo(2),hi(2)
       do i = lo(1),hi(1)
        
          rho = uout(i,j,URHO)
          
          if (rho <= sponge_start_density) then
             
             if (rho < sponge_center_density/sponge_start_factor) then
                fac = 1.d0
             else
                fac =  0.5d0*(1.d0 - cos(M_PI*(sponge_start_density - rho)/ &
                     (sponge_start_density - sponge_center_density/sponge_start_factor)))
             endif
             
             sponge_mult = 1.d0 / (1.d0 + fac * sponge_kappa * dt)
             
             ke_old = 0.5d0 * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / rho
             
             uout(i,j,UMX) = uout(i,j,UMX) * sponge_mult
             uout(i,j,UMY) = uout(i,j,UMY) * sponge_mult
             
             ke_new = 0.5d0 * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / rho
             
             uout(i,j,UEDEN) = uout(i,j,UEDEN) + (ke_new-ke_old)
             
          endif
          
       enddo
    enddo
    
  end subroutine sponge

end module sponge_module
