module sponge_module

  implicit none

contains
  
  subroutine sponge(uout,uout_l1,uout_l2,uout_l3,&
                    uout_h1,uout_h2,uout_h3,lo,hi,t,dt, &
                    dx,dy,dz,domlo,domhi, &
                    E_added,xmom_added,ymom_added,zmom_added)

    use bl_constants_module, only : M_PI
    use meth_params_module , only : NVAR, URHO, UMX, UMY, UMZ, UEDEN

    integer          :: lo(3), hi(3), domlo(3), domhi(3)
    integer          :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision :: t,dt
    double precision :: dx,dy,dz
    double precision :: E_added, xmom_added, ymom_added, zmom_added

    integer          :: i,j,k
    double precision :: rho, ke_old, ke_new, fac, sponge_mult, sponge_weighting
    double precision :: umx_old, umy_old, umz_old

    sponge_weighting = 1000.d0

    E_added = 0.0d0
    xmom_added = 0.0d0
    ymom_added = 0.0d0
    zmom_added = 0.0d0

    do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1),hi(1)

             rho = uout(i,j,k,URHO)
             if (rho .le. 100.0) then

                if (rho < 1.e-4) then
                   fac =  1.d0
                else 
                   fac =  0.5d0*(1.d0 - cos(M_PI*(rho - 100.0)/(1.e-4 - 100.0)))
                end if

                sponge_mult = 1.d0 / (1.d0 + fac * sponge_weighting * dt)

                ke_old = 0.5d0 * ( uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 &
                     +uout(i,j,k,UMZ)**2 )/ rho

                umx_old = uout(i,j,k,UMX  )
                umy_old = uout(i,j,k,UMY  )
                umz_old = uout(i,j,k,UMZ  )

                uout(i,j,k,UMX  ) = uout(i,j,k,UMX  ) * sponge_mult
                uout(i,j,k,UMY  ) = uout(i,j,k,UMY  ) * sponge_mult
                uout(i,j,k,UMZ  ) = uout(i,j,k,UMZ  ) * sponge_mult

                ke_new = 0.5d0 * ( uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 &
                     +uout(i,j,k,UMZ)**2 )/ rho

                uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + (ke_new-ke_old)

                E_added = E_added + (ke_new-ke_old)
                xmom_added = xmom_added + (uout(i,j,k,UMX  ) - umx_old)
                ymom_added = ymom_added + (uout(i,j,k,UMY  ) - umy_old)
                zmom_added = zmom_added + (uout(i,j,k,UMZ  ) - umz_old)

             end if

          enddo
       enddo
    enddo

  end subroutine sponge

end module sponge_module
