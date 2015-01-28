module grav_sources_module

  implicit none

  private

  public add_grav_source

contains

! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine add_grav_source(uin,uin_l1,uin_l2,uin_h1,uin_h2,&
                               uout,uout_l1,uout_l2,uout_h1,uout_h2,&
                               grav, gv_l1, gv_l2, gv_h1, gv_h2, &
                               lo,hi,dt,E_added, &
                               xmom_added,ymom_added)

      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, grav_source_type
      use bl_constants_module

      implicit none

      integer lo(2), hi(2)
      integer uin_l1,uin_l2,uin_h1,uin_h2
      integer  uout_l1, uout_l2, uout_h1, uout_h2
      integer  gv_l1, gv_l2, gv_h1, gv_h2

      double precision  uin( uin_l1: uin_h1, uin_l2: uin_h2,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
      double precision grav(  gv_l1:  gv_h1,  gv_l2:  gv_h2,2)
      double precision dt,E_added
      double precision xmom_added,ymom_added

      double precision :: rho
      double precision :: SrU, SrV, SrE
      double precision :: rhoInv
      double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
      double precision :: old_xmom, old_ymom
      integer          :: i, j

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: same as original, except in correction, it uses updates U
      ! 3: Puts all gravitational work into KE, not (rho e)


      ! Add gravitational source terms
      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_re = uout(i,j,UEDEN)
               old_ke = HALF * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / &
                                uout(i,j,URHO) 
               old_rhoeint = uout(i,j,UEDEN) - old_ke
               old_xmom = uout(i,j,UMX)
               old_ymom = uout(i,j,UMY)
               ! ****   End Diagnostics ****

               rho    = uin(i,j,URHO)
               rhoInv = ONE / rho

               SrU = rho * grav(i,j,1)
               SrV = rho * grav(i,j,2)

               uout(i,j,UMX)   = uout(i,j,UMX) + SrU * dt
               uout(i,j,UMY)   = uout(i,j,UMY) + SrV * dt

               if (grav_source_type == 1 .or. grav_source_type == 2) then

                   ! Src = rho u dot g, evaluated with all quantities at t^n
                   SrE = uin(i,j,UMX) * grav(i,j,1) + &
                         uin(i,j,UMY) * grav(i,j,2) 
                   uout(i,j,UEDEN) = uout(i,j,UEDEN) + SrE * dt

               else if (grav_source_type .eq. 3) then

                   new_ke = HALF * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / &
                                    uout(i,j,URHO) 
                   uout(i,j,UEDEN) = old_rhoeint + new_ke

               else 
                  call bl_error("Error:: Castro_grav_sources_2d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               new_ke = HALF * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / &
                                uout(i,j,URHO) 

               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_rhoeint = uout(i,j,UEDEN) - new_ke
 
               E_added =  E_added + uout(i,j,UEDEN) - old_re
               xmom_added = xmom_added + uout(i,j,UMX) - old_xmom
               ymom_added = ymom_added + uout(i,j,UMY) - old_ymom
               ! ****   End Diagnostics ****

         enddo
      enddo

      end subroutine add_grav_source

end module grav_sources_module

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_corrgsrc(lo,hi, &
                             gold,gold_l1,gold_l2,gold_h1,gold_h2, &
                             gnew,gnew_l1,gnew_l2,gnew_h1,gnew_h2, &
                             uold,uold_l1,uold_l2,uold_h1,uold_h2, &
                             unew,unew_l1,unew_l2,unew_h1,unew_h2, &
                             dt,E_added)

      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, grav_source_type
      use bl_constants_module

      implicit none

      integer lo(2),hi(2)
      integer gold_l1,gold_l2,gold_h1,gold_h2
      integer gnew_l1,gnew_l2,gnew_h1,gnew_h2
      integer uold_l1,uold_l2,uold_h1,uold_h2
      integer unew_l1,unew_l2,unew_h1,unew_h2
      double precision   gold(gold_l1:gold_h1,gold_l2:gold_h2,2)
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,2)
      double precision  uold(uold_l1:uold_h1,uold_l2:uold_h2,NVAR)
      double precision  unew(unew_l1:unew_h1,unew_l2:unew_h2,NVAR)
      double precision  dt,E_added

      integer i,j

      double precision SrU_old, SrV_old
      double precision SrU_new, SrV_new
      double precision SrUcorr, SrVcorr, SrEcorr
      double precision rhoo, Upo, Vpo
      double precision rhon, Upn, Vpn

      double precision rhooinv, rhoninv
      double precision old_ke, old_rhoeint, old_re
      double precision new_ke, new_rhoeint

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: Modification of type 1 that updates the U before constructing SrEcorr
      ! 3: Puts all gravitational work into KE, not (rho e)

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_re = unew(i,j,UEDEN)
               old_ke = HALF * (unew(i,j,UMX)**2 + unew(i,j,UMY)**2) / &
                                unew(i,j,URHO) 
               old_rhoeint = unew(i,j,UEDEN) - old_ke
               ! ****   End Diagnostics ****

               rhoo    = uold(i,j,URHO)
               rhooinv = ONE / uold(i,j,URHO)
               Upo     = uold(i,j,UMX) * rhooinv
               Vpo     = uold(i,j,UMY) * rhooinv

               ! Define old source terms
               SrU_old = rhoo * gold(i,j,1)
               SrV_old = rhoo * gold(i,j,2)

               rhon    = unew(i,j,URHO)
               rhoninv = ONE / unew(i,j,URHO)
               Upn     = unew(i,j,UMX) * rhoninv
               Vpn     = unew(i,j,UMY) * rhoninv

               ! Define new source terms
               SrU_new = rhon * gnew(i,j,1)
               SrV_new = rhon * gnew(i,j,2)

               ! Define corrections to source terms
               SrUcorr = HALF*(SrU_new - SrU_old)
               SrVcorr = HALF*(SrV_new - SrV_old)

               ! This does work (in 1-d)
               if (grav_source_type .eq. 1) then
                   SrEcorr =  HALF * ( (SrU_new * Upn - SrU_old * Upo) + &
                                       (SrV_new * Vpn - SrV_old * Vpo) )
               end if

               ! Correct state with correction terms
               unew(i,j,UMX)   = unew(i,j,UMX)   + SrUcorr*dt
               unew(i,j,UMY)   = unew(i,j,UMY)   + SrVcorr*dt

               if (grav_source_type .eq. 1) then

                  ! Note SrEcorr was constructed before updating unew
                  unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr*dt
               else if (grav_source_type .eq. 2) then
                  ! Note SrEcorr  is constructed  after updating unew
                  Upn     = unew(i,j,UMX) * rhoninv
                  Vpn     = unew(i,j,UMY) * rhoninv

                  SrEcorr = HALF * ( (SrU_new * Upn - SrU_old * Upo) + &
                                     (SrV_new * Vpn - SrV_old * Vpo) )
                  unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr*dt
               else if (grav_source_type .eq. 3) then
                   new_ke = HALF * (unew(i,j,UMX)**2 + unew(i,j,UMY)**2) / &
                                     unew(i,j,URHO) 
                   unew(i,j,UEDEN) = old_rhoeint + new_ke
               else 
                  call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_ke = HALF * (unew(i,j,UMX)**2 + unew(i,j,UMY)**2) / &
                                unew(i,j,URHO) 
               new_rhoeint = unew(i,j,UEDEN) - new_ke
 
                E_added =  E_added + unew(i,j,UEDEN) - old_re
               ! ****   End Diagnostics ****

         enddo
      enddo

      end subroutine ca_corrgsrc

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_syncgsrc(lo,hi, &
                             gphi,gphi_l1,gphi_l2,gphi_h1,gphi_h2, &
                             gdphi,gdphi_l1,gdphi_l2,gdphi_h1,gdphi_h2, &
                             state,state_l1,state_l2,state_h1,state_h2, &
                             dstate,dstate_l1,dstate_l2,dstate_h1,dstate_h2, &
                             sync_src,src_l1,src_l2,src_h1,src_h2, &
                             dt)

      use meth_params_module, only : NVAR, URHO, UMX, UMY
      use bl_constants_module

      implicit none

      integer lo(2),hi(2)
      integer gphi_l1,gphi_l2,gphi_h1,gphi_h2
      integer gdphi_l1,gdphi_l2,gdphi_h1,gdphi_h2
      integer state_l1,state_l2,state_h1,state_h2
      integer dstate_l1,dstate_l2,dstate_h1,dstate_h2
      integer src_l1,src_l2,src_h1,src_h2
      double precision     gphi(gphi_l1:gphi_h1,gphi_l2:gphi_h2,2)
      double precision    gdphi(gdphi_l1:gdphi_h1,gdphi_l2:gdphi_h2,2)
      double precision    state(state_l1:state_h1,state_l2:state_h2,NVAR)
      double precision   dstate(dstate_l1:dstate_h1,dstate_l2:dstate_h2,2+1)
      double precision sync_src(src_l1:src_h1,src_l2:src_h2,2+1)
      double precision dt

!     Note that dstate is drho and drhoU, state is the entire state, and src
!     is S_rhoU and S_rhoE

      integer i,j
      double precision rho_pre, rhoU_pre, rhoV_pre
      double precision gx, gy,dgx, dgy, SrU, SrV, SrE

      do j = lo(2),hi(2)
         do i = lo(1),hi(1)
            
            rho_pre  = state(i,j,URHO) - dstate(i,j,1)
            rhoU_pre = state(i,j,UMX)  - dstate(i,j,2)
            rhoV_pre = state(i,j,UMY)  - dstate(i,j,3)
            
            gx  = gphi(i,j,1)
            gy  = gphi(i,j,2)

            dgx = gdphi(i,j,1)
            dgy = gdphi(i,j,2)

            SrU = dstate(i,j,1)*gx + rho_pre*dgx
            SrV = dstate(i,j,1)*gy + rho_pre*dgy

            SrE = SrU * (rhoU_pre + HALF*SrU*dt)/rho_pre &
                 +SrV * (rhoV_pre + HALF*SrV*dt)/rho_pre
            
            sync_src(i,j,1) = SrU
            sync_src(i,j,2) = SrV
            sync_src(i,j,3) = SrE

         enddo
      enddo
      end subroutine ca_syncgsrc
