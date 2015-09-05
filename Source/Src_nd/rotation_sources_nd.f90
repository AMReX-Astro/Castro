  subroutine ca_rsrc(lo,hi,phi,phi_lo,phi_hi,rot,rot_lo,rot_hi, &
                     uold,uold_lo,uold_hi,unew,unew_lo,unew_hi,dx,dt,time,E_added,mom_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type
    use prob_params_module, only: coord_type, problo, center
    use bl_constants_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: phi_lo(3), phi_hi(3)
    integer         , intent(in   ) :: rot_lo(3), rot_hi(3)
    integer         , intent(in   ) :: uold_lo(3), uold_hi(3)
    integer         , intent(in   ) :: unew_lo(3), unew_hi(3)

    double precision, intent(in   ) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3))
    double precision, intent(in   ) :: rot(rot_lo(1):rot_hi(1),rot_lo(2):rot_hi(2),rot_lo(3):rot_hi(3),3)
    double precision, intent(in   ) :: uold(uold_lo(1):uold_hi(1),uold_lo(2):uold_hi(2),uold_lo(3):uold_hi(3),NVAR)
    double precision, intent(inout) :: unew(unew_lo(1):unew_hi(1),unew_lo(2):unew_hi(2),unew_lo(3):unew_hi(3),NVAR)
    double precision, intent(in   ) :: dx(3), dt, time

    integer          :: i, j ,k
    double precision :: Sr(3),SrE
    double precision :: rho, rhoInv

    double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
    double precision :: old_mom(3)
    double precision :: E_added, mom_added(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             rho = uold(i,j,k,URHO)
             rhoInv = ONE / rho

             ! **** Start Diagnostics ****
             old_re = unew(i,j,k,UEDEN)
             old_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
             old_rhoeint = unew(i,j,k,UEDEN) - old_ke
             old_mom = unew(i,j,k,UMX:UMZ)
             ! ****   End Diagnostics ****

             Sr = rho * rot(i,j,k,:) * dt

             unew(i,j,k,UMX:UMZ) = unew(i,j,k,UMX:UMZ) + Sr

             ! Kinetic energy source: this is v . the momentum source.
             ! We don't apply in the case of the conservative energy
             ! formulation.

             if (rot_source_type == 1 .or. rot_source_type == 2) then

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else if (rot_source_type .eq. 3) then

                new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
                SrE = new_ke - old_ke

             else if (rot_source_type .eq. 4) then

                ! Do nothing here, for the conservative rotation option.

                SrE = ZERO

             else 
                call bl_error("Error:: rotation_sources_nd.f90 :: invalid rot_source_type")
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

  end subroutine ca_rsrc



  subroutine ca_corrrsrc(lo,hi, &
                         pold,po_lo,po_hi, &
                         pnew,pn_lo,pn_hi, &
                         rold,ro_lo,ro_hi, &
                         rnew,rn_lo,rn_hi, &
                         uold,uo_lo,uo_hi, &
                         unew,un_lo,un_hi, &
                         flux1,f1_lo,f1_hi, &
                         flux2,f2_lo,f2_hi, &
                         flux3,f3_lo,f3_hi, &
                         dx,dt,time, &
                         vol,vol_lo,vol_hi, &
                         E_added,mom_added)

    ! Corrector step for the rotation source terms. This is applied after the hydrodynamics 
    ! update to fix the time-level n prediction and add the time-level n+1 data.
    ! This subroutine exists outside of the Fortran module above because it needs to be called 
    ! directly from C++.

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type, rot_axis
    use prob_params_module, only: coord_type, problo, center, dg
    use bl_constants_module
    use rotation_module, only: cross_product, get_omega, get_domegadt, rotational_acceleration

    implicit none

    integer          :: lo(3), hi(3)

    integer          :: po_lo(3),po_hi(3)
    integer          :: pn_lo(3),pn_hi(3)
    integer          :: ro_lo(3),ro_hi(3)
    integer          :: rn_lo(3),rn_hi(3)
    integer          :: uo_lo(3),uo_hi(3)
    integer          :: un_lo(3),un_hi(3)
    integer          :: f1_lo(3),f1_hi(3)
    integer          :: f2_lo(3),f2_hi(3)
    integer          :: f3_lo(3),f3_hi(3)
    integer          :: vol_lo(3),vol_hi(3)

    ! Old and new time rotational potential

    double precision :: pold(po_lo(1):po_hi(1),po_lo(2):po_hi(2),po_lo(3):po_hi(3))
    double precision :: pnew(pn_lo(1):pn_hi(1),pn_lo(2):pn_hi(2),pn_lo(3):pn_hi(3))

    ! Old and new time rotational acceleration

    double precision :: rold(ro_lo(1):ro_hi(1),ro_lo(2):ro_hi(2),ro_lo(3):ro_hi(3),3)
    double precision :: rnew(rn_lo(1):rn_hi(1),rn_lo(2):rn_hi(2),rn_lo(3):rn_hi(3),3)

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

    integer          :: i,j,k
    double precision :: r(3)
    double precision :: vnew(3),vold(3),omega_old(3),omega_new(3),domegadt_old(3),domegadt_new(3)
    double precision :: Sr_old(3), Sr_new(3), Srcorr(3), SrEcorr, SrE_old, SrE_new
    double precision :: rhoo, rhon, rhooinv, rhoninv

    double precision :: old_ke, old_rhoeint, old_re, new_ke, new_rhoeint
    double precision :: old_mom(3)

    double precision, pointer :: phi(:,:,:)
    double precision, pointer :: drho1(:,:,:), drho2(:,:,:), drho3(:,:,:)

    double precision :: mom1, mom2

    integer :: idir1, idir2, midx1, midx2

    ! Note that the time passed to this function
    ! is the new time at time-level n+1.
    
    omega_old = get_omega(time-dt)
    omega_new = get_omega(time   )

    domegadt_old = get_domegadt(time-dt)
    domegadt_new = get_domegadt(time   )

    if (rot_source_type == 4) then

       call bl_allocate(phi,   lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)
       call bl_allocate(drho1, lo(1),hi(1)+1,lo(2),hi(2),lo(3),hi(3))
       call bl_allocate(drho2, lo(1),hi(1),lo(2),hi(2)+1,lo(3),hi(3))
       call bl_allocate(drho3, lo(1),hi(1),lo(2),hi(2),lo(3),hi(3)+1)

       phi = ZERO

       do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
          do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
             do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                phi(i,j,k) = HALF * (pold(i,j,k) + pnew(i,j,k))
             enddo
          enddo
       enddo

       ! Construct the mass changes using the density flux from the hydro step. 
       ! Note that in the hydrodynamics step, these fluxes were already 
       ! multiplied by dA and dt, so dividing by the cell volume is enough to 
       ! get the density change (flux * dt / dx). This will be fine in the usual 
       ! case where the volume is the same in every cell, but may need to be 
       ! generalized when this assumption does not hold.
       
       drho1 = ZERO

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)+1*dg(1)
                drho1(i,j,k) = flux1(i,j,k,URHO) / vol(i,j,k)
             enddo
          enddo
       enddo

       drho2 = ZERO

       do k = lo(3), hi(3)
          do j = lo(2), hi(2)+1*dg(2)
             do i = lo(1), hi(1)
                drho2(i,j,k) = flux2(i,j,k,URHO) / vol(i,j,k)
             enddo
          enddo
       enddo

       drho3 = ZERO

       do k = lo(3), hi(3)+1*dg(3)
          do j = lo(2), hi(2)
             do i = lo(1), hi(1)
                drho3(i,j,k) = flux3(i,j,k,URHO) / vol(i,j,k)
             enddo
          enddo
       enddo

    endif

    do k = lo(3), hi(3)
       r(3) = problo(3) + dx(3)*(dble(k)+HALF) - center(3)
       do j = lo(2), hi(2)
          r(2) = problo(2) + dx(2)*(dble(j)+HALF) - center(2)
          do i = lo(1), hi(1)
             r(1) = problo(1) + dx(1)*(dble(i)+HALF) - center(1)

             rhoo = uold(i,j,k,URHO)
             rhooinv = ONE / uold(i,j,k,URHO)

             rhon = unew(i,j,k,URHO)
             rhoninv = ONE / unew(i,j,k,URHO)

             ! **** Start Diagnostics ****
             old_re = unew(i,j,k,UEDEN)
             old_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
             old_rhoeint = unew(i,j,k,UEDEN) - old_ke
             old_mom = unew(i,j,k,UMX:UMZ)
             ! ****   End Diagnostics ****

             ! Define old source terms

             vold = uold(i,j,k,UMX:UMZ) * rhooinv

             Sr_old = rhoo * rold(i,j,k,:) * dt
             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             vnew = unew(i,j,k,UMX:UMZ) * rhoninv

             Sr_new = rhon * rnew(i,j,k,:) * dt
             SrE_new = dot_product(vnew, Sr_new)

             ! Define correction terms

             Srcorr = HALF * (Sr_new - Sr_old)

             ! Correct momenta

             unew(i,j,k,UMX:UMZ) = unew(i,j,k,UMX:UMZ) + Srcorr

             ! Correct energy

             if (rot_source_type == 1) then

               ! If rot_source_type == 1, then calculate SrEcorr before updating the velocities.

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (rot_source_type == 2) then

                ! For this source type, we first update the momenta
                ! before we calculate the energy source term.

                vnew = unew(i,j,k,UMX:UMZ) * rhoninv
                Sr_new = rhon * rotational_acceleration(r, vnew, time) * dt
                SrE_new = dot_product(vnew, Sr_new)

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (rot_source_type == 3) then

                ! Instead of calculating the energy source term explicitly,
                ! we simply update the kinetic energy.

                new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
                SrEcorr = new_ke - old_ke

             else if (rot_source_type == 4) then

                ! Coupled momentum update.
                ! See Section 2.4 in the first wdmerger paper.

                ! Figure out which directions are updated, and then determine the right 
                ! array index relative to UMX (this works because UMX, UMY, UMZ are consecutive
                ! in the state array).

                idir1 = 1 + MOD(rot_axis    , 3)
                idir2 = 1 + MOD(rot_axis + 1, 3)

                midx1 = UMX + idir1 - 1
                midx2 = UMX + idir2 - 1

                mom1 = unew(i,j,k,midx1)
                mom2 = unew(i,j,k,midx2)

                ! Now do the implicit solve for the time-level n+1 Coriolis term. 
                ! It would be nice if this all could be generalized so that we don't 
                ! have to break it up by coordinate axis (in case the user wants to 
                ! rotate about multiple axes).

                unew(i,j,k,midx1) = (mom1 + dt * omega_new(rot_axis) * mom2) / (ONE + (dt * omega_new(rot_axis))**2)
                unew(i,j,k,midx2) = (mom2 - dt * omega_new(rot_axis) * mom1) / (ONE + (dt * omega_new(rot_axis))**2)

                ! Do the full corrector step with the centrifugal force (add 1/2 the new term, subtract 1/2 the old term)
                ! and do the remaining part of the corrector step for the Coriolis term (subtract 1/2 the old term). 

                unew(i,j,k,midx1) = unew(i,j,k,midx1) - dt * omega_old(rot_axis) * uold(i,j,k,midx2) &
                                  + HALF * dt * r(idir1) * omega_new(rot_axis)**2 * rhon &
                                  - HALF * dt * r(idir1) * omega_old(rot_axis)**2 * rhoo
                unew(i,j,k,midx2) = unew(i,j,k,midx2) + dt * omega_old(rot_axis) * uold(i,j,k,midx1) &
                                  + HALF * dt * r(idir2) * omega_new(rot_axis)**2 * rhon &
                                  - HALF * dt * r(idir2) * omega_old(rot_axis)**2 * rhoo

                ! The change in the gas energy is equal in magnitude to, and opposite in sign to,
                ! the change in the rotational potential energy, rho * phi.
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

                SrEcorr = - HALF * ( drho1(i  ,j,k) * (phi(i,j,k) - phi(i-1,j,k)) - &
                                     drho1(i+1,j,k) * (phi(i,j,k) - phi(i+1,j,k)) + &
                                     drho2(i,j  ,k) * (phi(i,j,k) - phi(i,j-1,k)) - &
                                     drho2(i,j+1,k) * (phi(i,j,k) - phi(i,j+1,k)) + &
                                     drho3(i,j,k  ) * (phi(i,j,k) - phi(i,j,k-1)) - &
                                     drho3(i,j,k+1) * (phi(i,j,k) - phi(i,j,k+1)) )

                ! Correct for the time rate of change of the potential, which acts 
                ! purely as a source term. For the velocities this is a corrector step
                ! and for the energy we add the full source term.

                Sr_old = - rhoo * cross_product(domegadt_old, r)
                Sr_new = - rhon * cross_product(domegadt_new, r)
               
                unew(i,j,k,UMX:UMZ) = unew(i,j,k,UMX:UMZ) + HALF * (Sr_new - Sr_old) * dt

                vnew = unew(i,j,k,UMX:UMZ) / rhon
                                
                SrEcorr = SrEcorr + HALF * (dot_product(vold, Sr_old) + dot_product(vnew, Sr_new)) * dt
                
             else 
                call bl_error("Error:: rotation_sources_nd.f90 :: invalid rot_source_type")
             end if

             unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr

             ! **** Start Diagnostics ****
             ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
             new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
             new_rhoeint = unew(i,j,k,UEDEN) - new_ke
             E_added =  E_added + unew(i,j,k,UEDEN) - old_re
             mom_added = mom_added + unew(i,j,k,UMX:UMZ) - old_mom
             ! ****   End Diagnostics ****

          enddo
       enddo
    enddo

    if (rot_source_type .eq. 4) then
       call bl_deallocate(phi)
       call bl_deallocate(drho1)
       call bl_deallocate(drho2)
       call bl_deallocate(drho3)
    endif

    end subroutine ca_corrrsrc



