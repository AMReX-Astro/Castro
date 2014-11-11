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
                               lo,hi,dt,E_added)

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
      double precision E_added

      double precision :: rho
      double precision :: SrU, SrV, SrW, SrE
      double precision :: rhoInv
      double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
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

               else 
                  call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               new_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                 uout(i,j,k,URHO) 

               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_rhoeint = uout(i,j,k,UEDEN) - new_ke
 
               E_added =  E_added + uout(i,j,k,UEDEN) - old_re
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
                             dt,E_added)

      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, grav_source_type
      use bl_constants_module

      implicit none

      integer lo(3),hi(3)
      integer gold_l1,gold_l2,gold_l3,gold_h1,gold_h2,gold_h3
      integer gnew_l1,gnew_l2,gnew_l3,gnew_h1,gnew_h2,gnew_h3
      integer uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
      integer unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3
      double precision   gold(gold_l1:gold_h1,gold_l2:gold_h2,gold_l3:gold_h3,3)
      double precision   gnew(gnew_l1:gnew_h1,gnew_l2:gnew_h2,gnew_l3:gnew_h3,3)
      double precision  uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
      double precision  unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)
      double precision  dt,E_added

      integer i,j,k

      double precision SrU_old, SrV_old, SrW_old
      double precision SrU_new, SrV_new, SrW_new
      double precision SrUcorr, SrVcorr, SrWcorr, SrEcorr
      double precision rhoo, Upo, Vpo, Wpo
      double precision rhon, Upn, Vpn, Wpn

      double precision rhooinv, rhoninv
      double precision old_ke, old_rhoeint, old_re
      double precision new_ke, new_rhoeint

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: Modification of type 1 that updates the U before constructing SrEcorr
      ! 3: Puts all gravitational work into KE, not (rho e)


      !$OMP PARALLEL DO PRIVATE(i,j,k,rhoo,Upo,Vpo,Wpo,SrU_old,SrV_old,SrW_old,rhon,Upn,Vpn,Wpn,SrU_new) &
      !$OMP PRIVATE(SrV_new,SrW_new,SrUcorr,SrVcorr,SrWcorr,SrEcorr,rhooinv,rhoninv) &
      !$OMP PRIVATE(old_ke,new_ke,old_rhoeint,new_rhoeint) reduction(+:E_added) 
      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               ! **** Start Diagnostics ****
               old_re = unew(i,j,k,UEDEN)
               old_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
               old_rhoeint = unew(i,j,k,UEDEN) - old_ke
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
               else 
                  call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
               end if

               ! **** Start Diagnostics ****
               ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
               new_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
               new_rhoeint = unew(i,j,k,UEDEN) - new_ke
 
                E_added =  E_added + unew(i,j,k,UEDEN) - old_re
               ! ****   End Diagnostics ****
            enddo
         enddo
      enddo
      !$OMP END PARALLEL DO

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

     !$OMP PARALLEL DO PRIVATE(i,j,k,rho_pre,rhoU_pre,rhoV_pre,rhoW_pre,gx,gy,gz,dgx,dgy,dgz,SrU,SrV,SrW,SrE)
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
     !$OMP END PARALLEL DO

     end subroutine ca_syncgsrc
