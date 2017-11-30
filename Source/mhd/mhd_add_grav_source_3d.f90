! :::
! ::: ------------------------------------------------------------------
! :::

    !===========================================================================
    ! This is called from within threaded loops in advance_gas_tile so *no* OMP here ...
    !===========================================================================
    subroutine add_grav_source(uin,uin_lo, uin_hi, &
                               uout,uout_lo,uout_hi, &
                               grav, gv_lo, gv_hi, &
                               lo,hi,dx,dy,dz,dt,e_added,ke_added)

      use amrex_fort_module, only : rt => amrex_real
      use eos_module
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, &
           UEDEN, grav_source_type

      implicit none

      integer lo(3), hi(3)
      integer uin_lo(3), uin_hi(3)
      integer  uout_lo(3), uout_hi(3)
      integer  gv_lo(3), gv_hi(3)

      real(rt) uin( uin_lo(1): uin_hi(1), uin_lo(2): uin_hi(2), uin_lo(3): uin_hi(3),NVAR)
      real(rt) uout(uout_lo(1):uout_hi(1),uout_lo(2):uout_hi(2),uout_lo(3):uout_hi(3),NVAR)
      real(rt) grav(gv_lo(1):  gv_hi(1),  gv_lo(2):  gv_hi(2),  gv_lo(3):  gv_hi(3),3)
      real(rt) dx, dy, dz, dt
      real(rt) e_added,ke_added

      !real(rt) :: a_half, a_oldsq, a_newsq, a_newsq_inv
      real(rt) :: rho
      real(rt) :: SrU, SrV, SrW, SrE
      real(rt) :: rhoInv, dt_a_new
      real(rt) :: old_rhoeint, new_rhoeint, old_ke, new_ke
      integer          :: i, j, k

      !a_half  = 0.5d0 * (a_old + a_new)
      !a_oldsq = a_old * a_old
      !a_newsq = a_new * a_new
      !a_newsq_inv = 1.d0 / a_newsq

      !dt_a_new    = dt / a_new

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 3: Puts all gravitational work into KE, not (rho e)

      ! Add gravitational source terms
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 
               old_rhoeint = uout(i,j,k,UEDEN) - old_ke
               ! ****   End Diagnostics ****

               rho    = uin(i,j,k,URHO)
               rhoInv = 1.0d0 / rho

               SrU = rho * grav(i,j,k,1)
               SrV = rho * grav(i,j,k,2)
               SrW = rho * grav(i,j,k,3)

               ! We use a_new here because we think of d/dt(a rho u) = ... + (rho g)
               uout(i,j,k,UMX)   = uout(i,j,k,UMX) + SrU 
               uout(i,j,k,UMY)   = uout(i,j,k,UMY) + SrV 
               uout(i,j,k,UMZ)   = uout(i,j,k,UMZ) + SrW

               if (grav_source_type .eq. 1) then

                   ! This does work (in 1-d)
                   ! Src = rho u dot g, evaluated with all quantities at t^n
                   SrE = uin(i,j,k,UMX) * grav(i,j,k,1) + &
                         uin(i,j,k,UMY) * grav(i,j,k,2) + &
                         uin(i,j,k,UMZ) * grav(i,j,k,3)
                   uout(i,j,k,UEDEN) = (uout(i,j,k,UEDEN) + SrE) 

               else if (grav_source_type .eq. 3) then

                   new_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                     uout(i,j,k,URHO) 
                   uout(i,j,k,UEDEN) = old_rhoeint + new_ke

               else 
                  call bl_error("Error:: Nyx_advection_3d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_ke = 0.5d0 * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 
               new_rhoeint = uout(i,j,k,UEDEN) - new_ke
 
                e_added =  e_added + (new_rhoeint - old_rhoeint)
               ke_added = ke_added + (new_ke      - old_ke     )
               ! ****   End Diagnostics ****

            enddo
         enddo
      enddo

      ! print *,' EADDED ',lo(1),lo(2),lo(3), e_added
      ! print *,'KEADDED ',lo(1),lo(2),lo(3),ke_added

      end subroutine add_grav_source
