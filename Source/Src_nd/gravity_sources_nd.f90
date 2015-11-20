    subroutine ca_gsrc(lo,hi,phi,phi_lo,phi_hi,grav,grav_lo,grav_hi, &
                       uold,uold_lo,uold_hi,unew,unew_lo,unew_hi,dx,dt,time,E_added,mom_added)

      use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, grav_source_type
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: phi_lo(3), phi_hi(3)
      integer          :: grav_lo(3), grav_hi(3)
      integer          :: uold_lo(3), uold_hi(3)
      integer          :: unew_lo(3), unew_hi(3)

      double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
      double precision :: grav(grav_lo(1):grav_hi(1),grav_lo(2):grav_hi(2),grav_lo(3):grav_hi(3),3)
      double precision :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
      double precision :: unew(unew_lo(1):unew_hi(1),unew_lo(2):unew_hi(2),unew_lo(3):unew_hi(3),NVAR)
      double precision :: dx(3), dt, time
      double precision :: E_added, mom_added(3)

      double precision :: rho, rhoInv
      double precision :: Sr(3), SrE
      double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
      double precision :: old_mom(3)
      integer          :: i, j, k

      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: same as original, except in correction, it uses updates U
      ! 3: Puts all gravitational work into KE, not (rho e)
      ! 4: Conservative energy formulation

      ! Add gravitational source terms
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)
               rho    = uold(i,j,k,URHO)
               rhoInv = ONE / rho

               ! **** Start Diagnostics ****
               old_re = unew(i,j,k,UEDEN)
               old_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
               old_rhoeint = unew(i,j,k,UEDEN) - old_ke
               old_mom = unew(i,j,k,UMX:UMZ)
               ! ****   End Diagnostics ****

               Sr = rho * grav(i,j,k,:) * dt

               unew(i,j,k,UMX:UMZ) = unew(i,j,k,UMX:UMZ) + Sr

               if (grav_source_type == 1 .or. grav_source_type == 2) then

                   ! Src = rho u dot g, evaluated with all quantities at t^n

                   SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

               else if (grav_source_type .eq. 3) then

                   new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
                   SrE = new_ke - old_ke

               else if (grav_source_type .eq. 4) then

                  ! Do nothing here, for the conservative gravity option.

                  SrE = ZERO

               else 
                  call bl_error("Error:: gravity_sources_nd.f90 :: invalid grav_source_type")
               end if

               unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrE

               ! **** Start Diagnostics ****
               new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
               new_rhoeint = unew(i,j,k,UEDEN) - new_ke
               E_added =  E_added + unew(i,j,k,UEDEN) - old_re
               mom_added = mom_added + unew(i,j,k,UMX:UMZ) - old_mom
               ! ****   End Diagnostics ****
            enddo
         enddo
      enddo

      end subroutine ca_gsrc

! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ca_corrgsrc(lo,hi, &
                             pold,po_lo,po_hi, &
                             pnew,pn_lo,pn_hi, &
                             gold,go_lo,go_hi, &
                             gnew,gn_lo,gn_hi, &
                             uold,uo_lo,uo_hi, &
                             unew,un_lo,un_hi, &
                             flux1,f1_lo,f1_hi, &
                             flux2,f2_lo,f2_hi, &
                             flux3,f3_lo,f3_hi, &
                             dx,dt,time, &
                             vol,vol_lo,vol_hi, &
                             E_added,mom_added)

      use mempool_module, only : bl_allocate, bl_deallocate
      use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, grav_source_type, gravity_type, get_g_from_phi
      use prob_params_module, only : dg
      use bl_constants_module
      use multifab_module
      use fundamental_constants_module, only: Gconst

      implicit none

      integer          :: lo(3), hi(3)

      integer          :: po_lo(3), po_hi(3)
      integer          :: pn_lo(3), pn_hi(3)
      integer          :: go_lo(3), go_hi(3)
      integer          :: gn_lo(3), gn_hi(3)
      integer          :: uo_lo(3), uo_hi(3)
      integer          :: un_lo(3), un_hi(3)
      integer          :: f1_lo(3), f1_hi(3)
      integer          :: f2_lo(3), f2_hi(3)
      integer          :: f3_lo(3), f3_hi(3)
      integer          :: vol_lo(3), vol_hi(3)

      ! Old and new time gravitational potential

      double precision :: pold(po_lo(1):po_hi(1),po_lo(2):po_hi(2),po_lo(3):po_hi(3))
      double precision :: pnew(pn_lo(1):pn_hi(1),pn_lo(2):pn_hi(2),pn_lo(3):pn_hi(3))

      ! Old and new time gravitational acceleration

      double precision :: gold(go_lo(1):go_hi(1),go_lo(2):go_hi(2),go_lo(3):go_hi(3),3)
      double precision :: gnew(gn_lo(1):gn_hi(1),gn_lo(2):gn_hi(2),gn_lo(3):gn_hi(3),3)

      ! Old and new time state data

      double precision :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
      double precision :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)

      ! Hydrodynamics fluxes

      double precision :: flux1(f1_lo(1):f1_hi(1),f1_lo(2):f1_hi(2),f1_lo(3):f1_hi(3),NVAR)
      double precision :: flux2(f2_lo(1):f2_hi(1),f2_lo(2):f2_hi(2),f2_lo(3):f2_hi(3),NVAR)
      double precision :: flux3(f3_lo(1):f3_hi(1),f3_lo(2):f3_hi(2),f3_lo(3):f3_hi(3),NVAR)

      double precision :: vol(vol_lo(1):vol_hi(1),vol_lo(2):vol_hi(2),vol_lo(3):vol_hi(3))
      double precision :: dx(3), dt, time
      double precision :: E_added, mom_added(3)

      integer          :: i, j, k

      double precision :: Sr_old(3), Sr_new(3), Srcorr(3)
      double precision :: vold(3), vnew(3)
      double precision :: SrE_old, SrE_new, SrEcorr
      double precision :: rhoo, rhooinv, rhon, rhoninv

      double precision :: old_ke, old_rhoeint, old_re
      double precision :: new_ke, new_rhoeint
      double precision :: old_mom(3)

      double precision, pointer :: phi(:,:,:)
      double precision, pointer :: drho1(:,:,:)
      double precision, pointer :: drho2(:,:,:)
      double precision, pointer :: drho3(:,:,:)
      double precision, pointer :: grav(:,:,:,:)
      double precision, pointer :: gravx(:,:,:)
      double precision, pointer :: gravy(:,:,:)
      double precision, pointer :: gravz(:,:,:)
      
      ! Gravitational source options for how to add the work to (rho E):
      ! grav_source_type = 
      ! 1: Original version ("does work")
      ! 2: Modification of type 1 that updates the U before constructing SrEcorr
      ! 3: Puts all gravitational work into KE, not (rho e)
      ! 4: Conservative gravity approach (discussed in first white dwarf merger paper).

      if (grav_source_type .eq. 4) then
         call bl_allocate(phi,   lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
         call bl_allocate(drho1, lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))
         call bl_allocate(drho2, lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))
         call bl_allocate(drho3, lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)
         call bl_allocate(grav,  lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1,1,3)
         call bl_allocate(gravx, lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))
         call bl_allocate(gravy, lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))
         call bl_allocate(gravz, lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)
         
         ! For our purposes, we want the time-level n+1/2 phi because we are 
         ! using fluxes evaluated at that time. To second order we can 
         ! average the new and old potentials.

         ! We will also negate the answer so that phi is negative,
         ! the usual physics convention, which will make the energy 
         ! update more easy to understand.

         phi = ZERO

         do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
            do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
               do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                  phi(i,j,k) = - HALF * (pnew(i,j,k) + pold(i,j,k))
                  grav(i,j,k,:) = HALF * (gnew(i,j,k,:) + gold(i,j,k,:))
               enddo
            enddo
         enddo

         ! We need to perform the following hack to deal with the fact that
         ! the potential is defined on cell edges, not cell centers, for ghost
         ! zones. We redefine the boundary zone values as equal to the adjacent
         ! cell minus the original value. Then later when we do the adjacent zone
         ! minus the boundary zone, we'll get the boundary value, which is what we want.
         
         if (dg(3) > 0) then
            
            k = lo(3) - 1 * dg(3)         
            do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
               do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                  phi(i,j,k) = phi(i,j,k+1) - phi(i,j,k)
               enddo
            enddo

            k = hi(3) + 1 * dg(3)         
            do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
               do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                  phi(i,j,k) = phi(i,j,k-1) - phi(i,j,k)
               enddo
            enddo            

         endif

         if (dg(2) > 0) then

            j = lo(2) - 1 * dg(2)
            do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
               do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                  phi(i,j,k) = phi(i,j+1,k) - phi(i,j,k)
               enddo
            enddo

            j = hi(2) + 1 * dg(2)
            do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
               do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                  phi(i,j,k) = phi(i,j-1,k) - phi(i,j,k)
               enddo
            enddo            

         endif

         if (dg(1) > 0) then

            i = lo(1) - 1 * dg(1)
            do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
               do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
                  phi(i,j,k) = phi(i+1,j,k) - phi(i,j,k)
               enddo
            enddo

            i = hi(1) + 1 * dg(1)
            do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
               do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
                  phi(i,j,k) = phi(i-1,j,k) - phi(i,j,k)
               enddo
            enddo

         endif
            
         ! Construct the mass changes using the density flux from the hydro step. 
         ! Note that in the hydrodynamics step, these fluxes were already 
         ! multiplied by dA and dt, so dividing by the cell volume is enough to 
         ! get the density change (flux * dt / dx). This will be fine in the usual 
         ! case where the volume is the same in every cell, but may need to be 
         ! generalized when this assumption does not hold.

         ! Also while we're doing this, construct the time-averaged
         ! edge-centered gravity.
         
         drho1 = ZERO

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)+1*dg(1)
                  drho1(i,j,k) = flux1(i,j,k,URHO) / vol(i,j,k)
                  gravx(i,j,k) = HALF * (grav(i,j,k,1) + grav(i-1,j,k,1))
               enddo
            enddo
         enddo

         drho2 = ZERO

         do k = lo(3), hi(3)
            do j = lo(2), hi(2)+1*dg(2)
               do i = lo(1), hi(1)
                  drho2(i,j,k) = flux2(i,j,k,URHO) / vol(i,j,k)
                  gravy(i,j,k) = HALF * (grav(i,j,k,2) + grav(i,j-1,k,2))
               enddo
            enddo
         enddo

         drho3 = ZERO

         do k = lo(3), hi(3)+1*dg(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  drho3(i,j,k) = flux3(i,j,k,URHO) / vol(i,j,k)
                  gravz(i,j,k) = HALF * (grav(i,j,k,3) + grav(i,j,k-1,3))
               enddo
            enddo
         enddo

      end if

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               rhoo    = uold(i,j,k,URHO)
               rhooinv = ONE / uold(i,j,k,URHO)

               rhon    = unew(i,j,k,URHO)
               rhoninv = ONE / unew(i,j,k,URHO)

               ! **** Start Diagnostics ****
               old_re = unew(i,j,k,UEDEN)
               old_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
               old_rhoeint = unew(i,j,k,UEDEN) - old_ke
               old_mom = unew(i,j,k,UMX:UMZ)
               ! ****   End Diagnostics ****

               ! Define old source terms
               
               vold = uold(i,j,k,UMX:UMZ) * rhooinv

               Sr_old = rhoo * gold(i,j,k,:) * dt
               SrE_old = dot_product(vold, Sr_old)

               ! Define new source terms

               vnew = unew(i,j,k,UMX:UMZ) * rhoninv

               Sr_new = rhon * gnew(i,j,k,:) * dt
               SrE_new = dot_product(vnew, Sr_new)

               ! Define corrections to source terms

               Srcorr = HALF * (Sr_new - Sr_old)

               ! Correct momenta

               unew(i,j,k,UMX:UMZ) = unew(i,j,k,UMX:UMZ) + Srcorr

               ! Correct energy

               if (grav_source_type .eq. 1) then

                  ! If grav_source_type == 1, then calculate SrEcorr before updating the velocities.

                  SrEcorr = HALF * (SrE_new - SrE_old)

               else if (grav_source_type .eq. 2) then

                  ! For this source type, we first update the momenta
                  ! before we calculate the energy source term.

                  vnew = unew(i,j,k,UMX:UMZ) * rhoninv
                  SrE_new = dot_product(vnew, Sr_new)

                  SrEcorr = HALF * (SrE_new - SrE_old)

               else if (grav_source_type .eq. 3) then

                  ! Instead of calculating the energy source term explicitly,
                  ! we simply update the kinetic energy.

                   new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
                   SrEcorr = new_ke - old_ke

               else if (grav_source_type .eq. 4) then

                  ! The change in the gas energy is equal in magnitude to, and opposite in sign to,
                  ! the change in the gravitational potential energy, rho * phi.
                  ! This must be true for the total energy, rho * E_g + rho * phi, to be conserved.
                  ! Consider as an example the zone interface i+1/2 in between zones i and i + 1.
                  ! There is an amount of mass drho_{i+1/2} leaving the zone. From this zone's perspective
                  ! it starts with a potential phi_i and leaves the zone with potential phi_{i+1/2} =
                  ! (1/2) * (phi_{i-1}+phi_{i}). Therefore the new rotational energy is equal to the mass
                  ! change multiplied by the difference between these two potentials.
                  ! This is a generalization of the cell-centered approach implemented in 
                  ! the other source options, which effectively are equal to 
                  ! SrEcorr = - drho(i,j,k) * phi(i,j,k),
                  ! where drho(i,j,k) = HALF * (unew(i,j,k,URHO) - uold(i,j,k,URHO)).

                  if (gravity_type == "PoissonGrav" .or. (gravity_type == "MonopoleGrav" &
                      .and. get_g_from_phi) ) then
                  
                     SrEcorr = - HALF * ( drho1(i  ,j,k) * (phi(i,j,k) - phi(i-1,j,k)) - &
                                          drho1(i+1,j,k) * (phi(i,j,k) - phi(i+1,j,k)) + &
                                          drho2(i,j  ,k) * (phi(i,j,k) - phi(i,j-1,k)) - &
                                          drho2(i,j+1,k) * (phi(i,j,k) - phi(i,j+1,k)) + &
                                          drho3(i,j,k  ) * (phi(i,j,k) - phi(i,j,k-1)) - &
                                          drho3(i,j,k+1) * (phi(i,j,k) - phi(i,j,k+1)) )

                  ! However, at present phi is only actually filled for Poisson gravity.
                  ! Here's an alternate version that only requires the use of the
                  ! gravitational acceleration. It relies on the concept that, to second order,
                  ! g_{i+1/2} = -( phi_{i+1} - phi_{i} ) / dx.
                     
                  else

                     SrEcorr = HALF * ( drho1(i  ,j,k) * gravx(i  ,j,k) * dx(1) + &
                                        drho1(i+1,j,k) * gravx(i+1,j,k) * dx(1) + &
                                        drho2(i,j  ,k) * gravy(i,j  ,k) * dx(2) + &
                                        drho2(i,j+1,k) * gravy(i,j+1,k) * dx(2) + &
                                        drho3(i,j,k  ) * gravz(i,j,k  ) * dx(3) + &
                                        drho3(i,j,k+1) * gravz(i,j,k+1) * dx(3) )
                  endif
                     

               else 
                  call bl_error("Error:: gravity_sources_nd.f90 :: invalid grav_source_type")
               end if

               unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr

               ! **** Start Diagnostics ****
               new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
               new_rhoeint = unew(i,j,k,UEDEN) - new_ke
               E_added =  E_added + unew(i,j,k,UEDEN) - old_re
               mom_added = mom_added + unew(i,j,k,UMX:UMZ) - old_mom
               ! ****   End Diagnostics ****

            enddo
         enddo
      enddo
      
      if (grav_source_type .eq. 4) then
         call bl_deallocate(phi)
         call bl_deallocate(drho1)
         call bl_deallocate(drho2)
         call bl_deallocate(drho3)
         call bl_deallocate(grav)
         call bl_deallocate(gravx)
         call bl_deallocate(gravy)
         call bl_deallocate(gravz)
      endif

      end subroutine ca_corrgsrc

! :::
! ::: ------------------------------------------------------------------
! :::

     subroutine ca_syncgsrc(lo,hi, &
                            gphi,gphi_lo,gphi_hi, &
                            gdphi,gdphi_lo,gdphi_hi, &
                            state,state_lo,state_hi, &
                            dstate,dstate_lo,dstate_hi, &
                            sync_src,src_lo,src_hi,dt)

     use meth_params_module, only : NVAR, URHO, UMX, UMZ
     use bl_constants_module

     implicit none

     integer          :: lo(3), hi(3)
     integer          :: gphi_lo(3), gphi_hi(3)
     integer          :: gdphi_lo(3), gdphi_hi(3)
     integer          :: state_lo(3), state_hi(3)
     integer          :: dstate_lo(3), dstate_hi(3)
     integer          :: src_lo(3), src_hi(3)
     double precision :: gphi(gphi_lo(1):gphi_hi(1),gphi_lo(2):gphi_hi(2),gphi_lo(3):gphi_hi(3),3)
     double precision :: gdphi(gdphi_lo(1):gdphi_hi(1),gdphi_lo(2):gdphi_hi(2),gdphi_lo(3):gdphi_hi(3),3)
     double precision :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
     double precision :: dstate(dstate_lo(1):dstate_hi(1),dstate_lo(2):dstate_hi(2),dstate_lo(3):dstate_hi(3),3+1)
     double precision :: sync_src(src_lo(1):src_hi(1),src_lo(2):src_hi(2),src_lo(3):src_hi(3),3+1)
     double precision :: dt

     !    Note that dstate is drho and drhoU, state is the entire state, and src
     !    is S_rhoU and S_rhoE

     integer          :: i, j, k
     double precision :: rho_pre, rhoU_pre(3)
     double precision :: grav(3), dgrav(3), Sr(3), SrE

     do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               rho_pre  = state(i,j,k,URHO) - dstate(i,j,k,1)
               rhoU_pre = state(i,j,k,UMX:UMZ)  - dstate(i,j,k,2:4)

               grav = gphi(i,j,k,:)
               dgrav = gdphi(i,j,k,:)

               Sr = dstate(i,j,k,1) * grav + rho_pre * dgrav
               SrE = ( dot_product(Sr, rhoU_pre) + HALF * dt * dot_product(Sr,Sr) ) / rho_pre

               sync_src(i,j,k,1:3) = Sr
               sync_src(i,j,k,4) = SrE

            enddo
         enddo
     enddo

     end subroutine ca_syncgsrc
