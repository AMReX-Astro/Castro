module sponge_module   

  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_h1,lo,hi,t,dt,dx,domlo,domhi, &
                    E_added, xmom_added)

    use bl_constants_module, only : M_PI
    use meth_params_module , only : NVAR, URHO, UMX, UEDEN

    integer          :: lo(1), hi(1), domlo(1), domhi(1)
    integer          :: uout_l1,uout_h1
    double precision :: uout(uout_l1:uout_h1,NVAR)
    double precision :: t,dt
    double precision :: dx
    double precision :: E_added, xmom_added

    integer          :: i
    double precision :: rho, ke_old, ke_new, fac, sponge_mult, sponge_weighting
    double precision :: umx_old

    sponge_weighting = 1000.d0

    xmom_added = 0.0d0
    E_added = 0.0d0

    do i = lo(1),hi(1)

       rho = uout(i,URHO)
       if (rho .le. 100.0) then

          if (rho < 1.e-4) then
             fac =  1.d0
          else 
             fac =  0.5d0*(1.d0 - cos(M_PI*(rho - 100.0)/(1.e-4 - 100.0)))
          end if

          sponge_mult = 1.d0 / (1.d0 + fac * sponge_weighting * dt)

          ke_old = 0.5d0 * uout(i,UMX)**2 / rho

          umx_old = uout(i,UMX)

          uout(i,UMX  ) = uout(i,UMX  ) * sponge_mult

          ke_new = 0.5d0 * uout(i,UMX)**2 / rho

          uout(i,UEDEN) = uout(i,UEDEN) + (ke_new-ke_old)

          xmom_added = xmom_added + (uout(i,UMX) - umx_old)
          E_added = E_added + ke_new-ke_old

       end if

    enddo

  end subroutine sponge

end module sponge_module
