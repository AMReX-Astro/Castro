! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_corrgsrc(lo,hi, &
           gold,gold_l1,gold_h1, &
           gnew,gnew_l1,gnew_h1, &
           uold,uold_l1,uold_h1, &
           unew,unew_l1,unew_h1, &
           dx,dt, &
           vol,vol_l1,vol_h1, &
           xmom_added,E_added)

      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, grav_source_type
      use bl_constants_module

      implicit none

      integer lo(1),hi(1)
      integer gold_l1,gold_h1
      integer gnew_l1,gnew_h1
      integer uold_l1,uold_h1
      integer unew_l1,unew_h1
      integer vol_l1,vol_h1
      double precision   gold(gold_l1:gold_h1)
      double precision   gnew(gnew_l1:gnew_h1)
      double precision   grav(gnew_l1:gnew_h1)
      double precision  uold(uold_l1:uold_h1,NVAR)
      double precision  unew(unew_l1:unew_h1,NVAR)
      double precision   vol(vol_l1:vol_h1)
      double precision dx(1), dt
      double precision E_added, xmom_added

      integer i
      double precision :: rhoo, rhooinv, Upo
      double precision :: rhon, rhoninv, Upn
      double precision :: SrU_new, SrU_old
      double precision :: SrUcorr,SrEcorr

      double precision :: old_ke, old_rhoeint, old_re
      double precision :: new_ke, new_rhoeint
      double precision :: old_xmom

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: Modification of type 1 that updates the U before constructing SrEcorr
      ! 3: Puts all gravitational work into KE, not (rho e)

      if (grav_source_type .eq. 4) then
         call bl_error("Error:: grav_source_type == 4 not enabled in one dimension.")
      endif

      do i = lo(1),hi(1)

         ! **** Start Diagnostics ****
         old_re = unew(i,UEDEN)
         old_ke = HALF * (unew(i,UMX)**2) / unew(i,URHO) 
         old_rhoeint = unew(i,UEDEN) - old_ke
         old_xmom = unew(i,UMX)
         ! ****   End Diagnostics ****

         ! Define old source term
         rhoo = uold(i,URHO)
         rhooinv = ONE / uold(i,URHO)
         Upo = uold(i,UMX) * rhooinv

         SrU_old = rhoo * gold(i)

         ! Define new source term            
         rhon = unew(i,URHO)
         rhoninv = ONE / unew(i,URHO)
         Upn  = unew(i,UMX) * rhoninv
         
         SrU_new = rhon * gnew(i)
         
         ! Define corrections to source terms
         SrUcorr = HALF*(SrU_new - SrU_old)

         if (grav_source_type .eq. 1) then
             SrEcorr =  HALF * ( SrU_new * Upn - SrU_old * Upo )
         end if

         ! Correct state with correction terms
         unew(i,UMX) = unew(i,UMX) + SrUcorr * dt

         if (grav_source_type .eq. 1) then

            ! Note SrEcorr was constructed before updating unew
            unew(i,UEDEN) = unew(i,UEDEN) + SrEcorr*dt
         else if (grav_source_type .eq. 2) then
            ! Note SrEcorr  is constructed  after updating unew
            Upn     = unew(i,UMX) * rhoninv

            SrEcorr = HALF * ( SrU_new * Upn - SrU_old * Upo ) 

            unew(i,UEDEN) = unew(i,UEDEN) + SrEcorr*dt
         else if (grav_source_type .eq. 3) then
             new_ke = HALF * (unew(i,UMX)**2) / unew(i,URHO) 
             unew(i,UEDEN) = old_rhoeint + new_ke
         else 
            call bl_error("Error:: Castro_grav_sources_1d.f90 :: bogus grav_source_type")
         end if

         ! **** Start Diagnostics ****
         ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
         new_ke = HALF * (unew(i,UMX)**2) / unew(i,URHO) 
         new_rhoeint = unew(i,UEDEN) - new_ke

         E_added =  E_added + unew(i,UEDEN) - old_re
         xmom_added = xmom_added + unew(i,UMX) - old_xmom
         ! ****   End Diagnostics ****

      enddo

      end subroutine ca_corrgsrc

! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_syncgsrc(lo,hi, &
           gphi,gphi_l1,gphi_h1, &
           gdphi,gdphi_l1,gdphi_h1, &
           state,state_l1,state_h1, &
           dstate,dstate_l1,dstate_h1, &
           sync_src,src_l1,src_h1,dt)

      use meth_params_module, only : NVAR, URHO, UMX
      use bl_constants_module

      implicit none

      integer lo(1),hi(1)
      integer gphi_l1,gphi_h1
      integer gdphi_l1,gdphi_h1
      integer state_l1,state_h1
      integer dstate_l1,dstate_h1
      integer src_l1,src_h1
      double precision     gphi(gphi_l1:gphi_h1)
      double precision    gdphi(gdphi_l1:gdphi_h1)
      double precision    state(state_l1:state_h1,NVAR)
      double precision   dstate(dstate_l1:dstate_h1,1+1)
      double precision sync_src(src_l1:src_h1,1+1)
      double precision dt

!     Note that dstate is drho and drhoU, state is the entire state, and src
!     is S_rhoU and S_rhoE

      integer i
      double precision rho_pre, rhoU_pre
      double precision gx, dgx, SrU, SrE

      do i = lo(1),hi(1)
            
         rho_pre  = state(i,URHO) - dstate(i,1)
         rhoU_pre = state(i,UMX)  - dstate(i,2)
         
         gx  = gphi(i)
         dgx = gdphi(i)

         SrU = dstate(i,1)*gx + rho_pre*dgx

         SrE = SrU * (rhoU_pre + HALF*SrU*dt)/rho_pre
         
         sync_src(i,1) = SrU
         sync_src(i,2) = SrE

      enddo
      end subroutine ca_syncgsrc

