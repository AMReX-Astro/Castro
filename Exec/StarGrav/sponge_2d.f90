module sponge_module

  implicit none

contains

  subroutine sponge(uout,uout_l1,uout_l2,&
                    uout_h1,uout_h2,lo,hi,t,dt, &
                    dx,dy,domlo,domhi,E_added,xmom_added,ymom_added)

    use bl_constants_module, only : M_PI
    use meth_params_module , only : NVAR, URHO, UMX, UMY, UEDEN

    integer          :: lo(2), hi(2), domlo(2), domhi(2)
    integer          :: uout_l1,uout_l2,uout_h1,uout_h2
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision :: t,dt
    double precision :: dx,dy
    double precision :: E_added, xmom_added, ymom_added

    integer          :: i,j
    double precision :: rho, ke_old, ke_new, fac, sponge_mult, sponge_weighting
    double precision :: umx_old, umy_old

    sponge_weighting = 1000.d0

    E_added = 0.0d0
    xmom_added = 0.0d0
    ymom_added = 0.0d0

    do j = lo(2),hi(2)
       do i = lo(1),hi(1)

          rho = uout(i,j,URHO)
          if (rho .le. 100.0) then

             if (rho < 1.e-4) then
                fac =  1.d0
             else 
                fac =  0.5d0*(1.d0 - cos(M_PI*(rho - 100.0)/(1.e-4 - 100.0)))
             end if

             sponge_mult = 1.d0 / (1.d0 + fac * sponge_weighting * dt)

             ke_old = 0.5d0 * ( uout(i,j,UMX)**2 + uout(i,j,UMY)**2 )/ rho

             umx_old = uout(i,j,UMX  )
             umy_old = uout(i,j,UMY  )

             uout(i,j,UMX  ) = uout(i,j,UMX  ) * sponge_mult
             uout(i,j,UMY  ) = uout(i,j,UMY  ) * sponge_mult

             ke_new = 0.5d0 * ( uout(i,j,UMX)**2 + uout(i,j,UMY)**2 )/ rho

             uout(i,j,UEDEN) = uout(i,j,UEDEN) + (ke_new-ke_old)

             E_added = E_added + (ke_new-ke_old)
             xmom_added = xmom_added + (uout(i,j,UMX  )-umx_old)
             ymom_added = ymom_added + (uout(i,j,UMY  )-umy_old)
          end if

       enddo
    enddo

  end subroutine sponge

end module sponge_module
