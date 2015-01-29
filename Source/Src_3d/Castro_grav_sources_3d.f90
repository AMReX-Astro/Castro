module grav_sources_module

  implicit none

  private

  public add_grav_source

contains

! :::
! ::: ------------------------------------------------------------------
! :::

    subroutine add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               lo,hi,dt,E_added,xmom_added,ymom_added,zmom_added)

      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, grav_source_type
      use bl_constants_module

      implicit none

      integer lo(3), hi(3)
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer  uout_l1, uout_l2, uout_l3, uout_h1, uout_h2, uout_h3
      integer  gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3

      double precision  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
      double precision uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
      double precision grav(  gv_l1:  gv_h1,  gv_l2:  gv_h2,  gv_l3:  gv_h3,3)
      double precision dt
      double precision E_added, xmom_added, ymom_added, zmom_added

      double precision :: rho
      double precision :: SrU, SrV, SrW, SrE
      double precision :: rhoInv
      double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
      double precision :: old_xmom, old_ymom, old_zmom
      integer          :: i, j, k

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: same as original, except in correction, it uses updates U
      ! 3: Puts all gravitational work into KE, not (rho e)


      ! Add gravitational source terms
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_re = uout(i,j,k,UEDEN)
               old_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 
               old_rhoeint = uout(i,j,k,UEDEN) - old_ke
               old_xmom = uout(i,j,k,UMX)
               old_ymom = uout(i,j,k,UMY)
               old_zmom = uout(i,j,k,UMZ)
               ! ****   End Diagnostics ****

               rho    = uin(i,j,k,URHO)
               rhoInv = ONE / rho

               SrU = rho * grav(i,j,k,1)
               SrV = rho * grav(i,j,k,2)
               SrW = rho * grav(i,j,k,3)

               uout(i,j,k,UMX)   = uout(i,j,k,UMX) + SrU * dt
               uout(i,j,k,UMY)   = uout(i,j,k,UMY) + SrV * dt
               uout(i,j,k,UMZ)   = uout(i,j,k,UMZ) + SrW * dt

               if (grav_source_type == 1 .or. grav_source_type == 2) then

                   ! Src = rho u dot g, evaluated with all quantities at t^n
                   SrE = uin(i,j,k,UMX) * grav(i,j,k,1) + &
                         uin(i,j,k,UMY) * grav(i,j,k,2) + &
                         uin(i,j,k,UMZ) * grav(i,j,k,3)

                   uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + SrE * dt

               else if (grav_source_type .eq. 3) then

                   new_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                     uout(i,j,k,URHO) 
                   uout(i,j,k,UEDEN) = old_rhoeint + new_ke

               else if (grav_source_type .eq. 4) then

                  ! Do nothing here, for the conservative gravity option.

               else 
                  call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               new_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 

               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_rhoeint = uout(i,j,k,UEDEN) - new_ke
 
               E_added =  E_added + uout(i,j,k,UEDEN) - old_re

               xmom_added = xmom_added + uout(i,j,k,UMX) - old_xmom
               ymom_added = ymom_added + uout(i,j,k,UMY) - old_ymom
               zmom_added = zmom_added + uout(i,j,k,UMZ) - old_zmom
               ! ****   End Diagnostics ****

            enddo
         enddo
      enddo

      end subroutine add_grav_source

end module grav_sources_module

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_corrgsrc(lo,hi, &
                             gold,gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3, &
                             gnew,gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3, &
                             uold,uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3, &
                             unew,unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3, &
                             pold,pold_l1,pold_l2,pold_l3,pold_h1,pold_h2,pold_h3, &
                             pnew,pnew_l1,pnew_l2,pnew_l3,pnew_h1,pnew_h2,pnew_h3, &
                             flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                             flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                             flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                             dx,dt, &
                             vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                             xmom_added,ymom_added,zmom_added,E_added)

      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, grav_source_type
      use bl_constants_module
      use multifab_module
      use fundamental_constants_module, only: Gconst

      implicit none

      integer lo(3),hi(3)
      integer gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3
      integer gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3
      integer uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
      integer unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3
      integer pold_l1,pold_l2,pold_l3,pold_h1,pold_h2,pold_h3
      integer pnew_l1,pnew_l2,pnew_l3,pnew_h1,pnew_h2,pnew_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
      integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
      integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3
      double precision   gold(gold_l1:gold_h1,gold_l2:gold_h2,gold_l3:gold_h3,3)
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      double precision   grav(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      double precision  uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
      double precision  unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)
      double precision  pold(pold_l1:pold_h1,pold_l2:pold_h2,pold_l3:pold_h3)
      double precision  pnew(pnew_l1:pnew_h1,pnew_l2:pnew_h2,pnew_l3:pnew_h3)
      double precision   phi(pnew_l1:pnew_h1,pnew_l2:pnew_h2,pnew_l3:pnew_h3)
      double precision flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
      double precision flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
      double precision flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)
      double precision   vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)
      double precision  dx(3),dt,E_added

      integer i,j,k

      double precision SrU_old, SrV_old, SrW_old
      double precision SrU_new, SrV_new, SrW_new
      double precision SrUcorr, SrVcorr, SrWcorr, SrEcorr
      double precision rhoo, Upo, Vpo, Wpo
      double precision rhon, Upn, Vpn, Wpn

      double precision rhooinv, rhoninv
      double precision old_ke, old_rhoeint, old_re
      double precision new_ke, new_rhoeint
      double precision old_xmom, old_ymom, old_zmom, xmom_added, ymom_added, zmom_added
      double precision rho_E_added, flux_added

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: Modification of type 1 that updates the U before constructing SrEcorr
      ! 3: Puts all gravitational work into KE, not (rho e)
      ! 4: Conservative gravity approach from the AREPO code paper (Springel 2010).

      do k = lo(3)-1, hi(3)+1
         do j = lo(2)-1, hi(2)+1
            do i = lo(1)-1, hi(1)+1
               grav(i,j,k,:) = HALF * (gnew(i,j,k,:) + gold(i,j,k,:))
               phi(i,j,k)    = HALF * (pnew(i,j,k) + pold(i,j,k))               
            enddo
         enddo
      enddo

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_re = unew(i,j,k,UEDEN)
               old_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
               old_rhoeint = unew(i,j,k,UEDEN) - old_ke
               old_xmom = unew(i,j,k,UMX)
               old_ymom = unew(i,j,k,UMY)
               old_zmom = unew(i,j,k,UMZ)
               ! ****   End Diagnostics ****

               rhoo    = uold(i,j,k,URHO)
               rhooinv = ONE / uold(i,j,k,URHO)
               Upo     = uold(i,j,k,UMX) * rhooinv
               Vpo     = uold(i,j,k,UMY) * rhooinv
               Wpo     = uold(i,j,k,UMZ) * rhooinv

               ! Define old source terms
               SrU_old = rhoo * gold(i,j,k,1)
               SrV_old = rhoo * gold(i,j,k,2)
               SrW_old = rhoo * gold(i,j,k,3)

               rhon    = unew(i,j,k,URHO)
               rhoninv = ONE / unew(i,j,k,URHO)
               Upn     = unew(i,j,k,UMX) * rhoninv
               Vpn     = unew(i,j,k,UMY) * rhoninv
               Wpn     = unew(i,j,k,UMZ) * rhoninv

               ! Define new source terms
               SrU_new = rhon * gnew(i,j,k,1)
               SrV_new = rhon * gnew(i,j,k,2)
               SrW_new = rhon * gnew(i,j,k,3)

               ! Define corrections to source terms
               SrUcorr = HALF*(SrU_new - SrU_old)
               SrVcorr = HALF*(SrV_new - SrV_old)
               SrWcorr = HALF*(SrW_new - SrW_old)

               if (grav_source_type .eq. 1) then
                   SrEcorr =  HALF * ( (SrU_new * Upn - SrU_old * Upo) + &
                                       (SrV_new * Vpn - SrV_old * Vpo) + &
                                       (SrW_new * Wpn - SrW_old * Wpo) )
               end if

               ! Correct state with correction terms
               unew(i,j,k,UMX)   = unew(i,j,k,UMX)   + SrUcorr*dt
               unew(i,j,k,UMY)   = unew(i,j,k,UMY)   + SrVcorr*dt
               unew(i,j,k,UMZ)   = unew(i,j,k,UMZ)   + SrWcorr*dt

               if (grav_source_type .eq. 1) then

                   ! Note SrEcorr was constructed before updating unew
                   unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr*dt

               else if (grav_source_type .eq. 2) then

                   ! Note SrEcorr  is constructed  after updating unew
                   Upn     = unew(i,j,k,UMX) * rhoninv
                   Vpn     = unew(i,j,k,UMY) * rhoninv
                   Wpn     = unew(i,j,k,UMZ) * rhoninv
                   SrEcorr =  HALF * ( (SrU_new * Upn - SrU_old * Upo) + &
                                       (SrV_new * Vpn - SrV_old * Vpo) + &
                                       (SrW_new * Wpn - SrW_old * Wpo) )

                   unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr*dt

               else if (grav_source_type .eq. 3) then

                   new_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                    unew(i,j,k,URHO) 

                   unew(i,j,k,UEDEN) = old_rhoeint + new_ke

               else if (grav_source_type .eq. 4) then

                  SrEcorr = HALF * flux1(i  ,j,k,URHO) * (phi(i  ,j,k) - phi(i-1,j,k)) + &
                            HALF * flux1(i+1,j,k,URHO) * (phi(i+1,j,k) - phi(i  ,j,k)) + &
                            HALF * flux2(i,j  ,k,URHO) * (phi(i,j,  k) - phi(i,j-1,k)) + &
                            HALF * flux2(i,j+1,k,URHO) * (phi(i,j+1,k) - phi(i,j  ,k)) + &
                            HALF * flux3(i,j,k  ,URHO) * (phi(i,j,k  ) - phi(i,j,k-1)) + &
                            HALF * flux3(i,j,k+1,URHO) * (phi(i,j,k+1) - phi(i,j,k  ))

                  SrEcorr = SrEcorr / vol(i,j,k)

                  unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr

               else 
                  call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
               new_rhoeint = unew(i,j,k,UEDEN) - new_ke
               E_added =  E_added + unew(i,j,k,UEDEN) - old_re
               xmom_added = xmom_added + unew(i,j,k,UMX) - old_xmom
               ymom_added = ymom_added + unew(i,j,k,UMY) - old_ymom
               zmom_added = zmom_added + unew(i,j,k,UMZ) - old_zmom
               ! ****   End Diagnostics ****

            enddo
         enddo
      enddo

      end subroutine ca_corrgsrc

! :::
! ::: ------------------------------------------------------------------
! :::

     subroutine ca_syncgsrc(lo,hi, &
                            gphi,gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3, &
                            gdphi,gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3, &
                            state,state_l1,state_l2,state_l3,state_h1,state_h2,state_h3, &
                            dstate,dstate_l1,dstate_l2,dstate_l3, &
                            dstate_h1,dstate_h2,dstate_h3, &
                            sync_src,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3,dt)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ
     use bl_constants_module

     implicit none

     integer lo(3),hi(3)
     integer gphi_l1,gphi_l2,gphi_l3,gphi_h1,gphi_h2,gphi_h3
     integer gdphi_l1,gdphi_l2,gdphi_l3,gdphi_h1,gdphi_h2,gdphi_h3
     integer state_l1,state_l2,state_l3,state_h1,state_h2,state_h3
     integer dstate_l1,dstate_l2,dstate_l3,dstate_h1,dstate_h2,dstate_h3
     integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
     double precision   gphi(gphi_l1:gphi_h1,gphi_l2:gphi_h2,gphi_l3:gphi_h3,3)
     double precision  gdphi(gdphi_l1:gdphi_h1,gdphi_l2:gdphi_h2,gdphi_l3:gdphi_h3,3)
     double precision  state(state_l1:state_h1,state_l2:state_h2,state_l3:state_h3,NVAR)
     double precision dstate(dstate_l1:dstate_h1,dstate_l2:dstate_h2,dstate_l3:dstate_h3,3+1)
     double precision sync_src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3,3+1)
     double precision dt

     !    Note that dstate is drho and drhoU, state is the entire state, and src
     !    is S_rhoU and S_rhoE

     integer          :: i,j,k
     double precision :: rho_pre, rhoU_pre, rhoV_pre, rhoW_pre
     double precision :: gx, gy, gz, dgx, dgy, dgz, SrU, SrV, SrW, SrE

     do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rho_pre  = state(i,j,k,URHO) - dstate(i,j,k,1)
               rhoU_pre = state(i,j,k,UMX)  - dstate(i,j,k,2)
               rhoV_pre = state(i,j,k,UMY)  - dstate(i,j,k,3)
               rhoW_pre = state(i,j,k,UMZ)  - dstate(i,j,k,4)

               gx  = gphi(i,j,k,1)
               gy  = gphi(i,j,k,2)
               gz  = gphi(i,j,k,3)

               dgx = gdphi(i,j,k,1)
               dgy = gdphi(i,j,k,2)
               dgz = gdphi(i,j,k,3)

               SrU = dstate(i,j,k,1)*gx + rho_pre*dgx
               SrV = dstate(i,j,k,1)*gy + rho_pre*dgy
               SrW = dstate(i,j,k,1)*gz + rho_pre*dgz

               SrE = ( SrU * (rhoU_pre + (HALF*dt)*SrU) + &
                       SrV * (rhoV_pre + (HALF*dt)*SrV) + &
                       SrW * (rhoW_pre + (HALF*dt)*SrW) ) / rho_pre

               sync_src(i,j,k,1) = SrU
               sync_src(i,j,k,2) = SrV
               sync_src(i,j,k,3) = SrW
               sync_src(i,j,k,4) = SrE

            enddo
         enddo
     enddo

     end subroutine ca_syncgsrc
