! ::: 
! ::: ------------------------------------------------------------------
! ::: 

      subroutine ca_corrgsrc(lo,hi, &
           gold,gold_l1,gold_h1, &
           gnew,gnew_l1,gnew_h1, &
           uold,uold_l1,uold_h1, &
           unew,unew_l1,unew_h1,dt)

      use meth_params_module, only : NVAR, URHO, UMX, UEDEN

      implicit none

      integer lo(1),hi(1)
      integer gold_l1,gold_h1
      integer gnew_l1,gnew_h1
      integer uold_l1,uold_h1
      integer unew_l1,unew_h1
      double precision   gold(gold_l1:gold_h1)
      double precision   gnew(gnew_l1:gnew_h1)
      double precision  uold(uold_l1:uold_h1,NVAR)
      double precision  unew(unew_l1:unew_h1,NVAR)
      double precision dt

      integer i
      double precision :: rhon, Upn
      double precision :: SrU_new, SrU_old
      double precision :: SrUcorr,SrEcorr
      double precision :: Upo

      do i = lo(1),hi(1)

         ! Define old source term
         SrU_old = uold(i,URHO) * gold(i)
            
         rhon = unew(i,URHO)
         Upn  = unew(i,UMX) / rhon
         Upo  = uold(i,UMX) / uold(i,URHO)
         
         ! Define new source term
         SrU_new = rhon * gnew(i)
         
         ! Define corrections to source terms
         SrUcorr = 0.5d0*(SrU_new - SrU_old)

         unew(i,UMX) = unew(i,UMX) + dt*SrUcorr

         ! This doesn't work
         ! SrEcorr = SrUcorr*(Upn + SrUcorr*dt/(2*rhon))

         ! This works
         SrEcorr = 0.5d0*(SrU_new * Upn - SrU_old * Upo)

         unew(i,UEDEN) = unew(i,UEDEN) + SrEcorr*dt

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

         SrE = SrU * (rhoU_pre + 0.5*SrU*dt)/rho_pre
         
         sync_src(i,1) = SrU
         sync_src(i,2) = SrE

      enddo
      end subroutine ca_syncgsrc

