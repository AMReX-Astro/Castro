  subroutine ca_rsrc(lo,hi,domlo,domhi,phi,phi_lo,phi_hi,rot,rot_lo,rot_hi, &
                     uold,uold_lo,uold_hi,unew,unew_lo,unew_hi,dx,dt,time,E_added,mom_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, rot_period, rot_source_type
    use prob_params_module, only: coord_type, problo, center
    use bl_constants_module
    use castro_util_module, only: position
    use hybrid_advection_module, only: add_momentum_source

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: domlo(3), domhi(3)
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
    double precision :: old_mom(3), loc(3)
    double precision :: E_added, mom_added(3)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center
               
             rho = uold(i,j,k,URHO)
             rhoInv = ONE / rho

             ! **** Start Diagnostics ****
             old_re = unew(i,j,k,UEDEN)
             old_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
             old_rhoeint = unew(i,j,k,UEDEN) - old_ke
             old_mom = unew(i,j,k,UMX:UMZ)
             ! ****   End Diagnostics ****

             Sr = rho * rot(i,j,k,:) * dt

             call add_momentum_source(loc, unew(i,j,k,UMX:UMZ), Sr)

             ! Kinetic energy source: this is v . the momentum source.
             ! We don't apply in the case of the conservative energy
             ! formulation.

             if (rot_source_type == 1 .or. rot_source_type == 2) then

                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)

             else if (rot_source_type .eq. 3) then

                new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoInv
                SrE = new_ke - old_ke

             else if (rot_source_type .eq. 4) then

                ! Add a predictor here; we'll remove this later.
                
                SrE = dot_product(uold(i,j,k,UMX:UMZ) * rhoInv, Sr)
                
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
                         domlo,domhi, &
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

    ! Corrector step for the rotation source terms. This is applied
    ! after the hydrodynamics update to fix the time-level n
    ! prediction and add the time-level n+1 data.  This subroutine
    ! exists outside of the Fortran module above because it needs to
    ! be called directly from C++.

    use mempool_module, only : bl_allocate, bl_deallocate
    use meth_params_module, only: NVAR, URHO, UMX, UMZ, UEDEN, rot_period, rot_source_type
    use prob_params_module, only: coord_type, problo, center, dg
    use bl_constants_module
    use math_module, only: cross_product
    use rotation_module, only: get_omega, get_domegadt, rotational_acceleration
    use castro_util_module, only: position
    use hybrid_advection_module, only: add_momentum_source

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: domlo(3), domhi(3)

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
    double precision :: loc(3)
    double precision :: vnew(3),vold(3),omega_old(3),omega_new(3),domegadt_old(3),domegadt_new(3)
    double precision :: Sr_old(3), Sr_new(3), Srcorr(3), SrEcorr, SrE_old, SrE_new
    double precision :: rhoo, rhon, rhooinv, rhoninv

    double precision :: old_ke, old_rhoeint, old_re, new_ke, new_rhoeint
    double precision :: old_mom(3), dt_omega_matrix(3,3), dt_omega(3), new_mom(3)

    double precision, pointer :: phi(:,:,:)

    double precision :: mom1, mom2

    integer :: idir1, idir2, midx1, midx2

    ! Rotation source options for how to add the work to (rho E):
    ! rot_source_type = 
    ! 1: Standard version ("does work")
    ! 2: Modification of type 1 that updates the momentum before constructing the energy corrector
    ! 3: Puts all work into KE, not (rho e)
    ! 4: Conservative rotation approach (discussed in first white dwarf merger paper)
    
    ! Note that the time passed to this function
    ! is the new time at time-level n+1.
    
    omega_old = get_omega(time-dt)
    omega_new = get_omega(time   )

    domegadt_old = get_domegadt(time-dt)
    domegadt_new = get_domegadt(time   )

    if (rot_source_type == 4) then

       call bl_allocate(phi,lo(1)-1,hi(1)+1,lo(2)-1,hi(2)+1,lo(3)-1,hi(3)+1)

       phi = ZERO

       do k = lo(3)-1*dg(3), hi(3)+1*dg(3)
          do j = lo(2)-1*dg(2), hi(2)+1*dg(2)
             do i = lo(1)-1*dg(1), hi(1)+1*dg(1)
                phi(i,j,k) = HALF * (pold(i,j,k) + pnew(i,j,k))
             enddo
          enddo
       enddo

       dt_omega = dt * omega_new
       
       dt_omega_matrix(1,1) = ONE + dt_omega(1)**2
       dt_omega_matrix(1,2) = dt_omega(1) * dt_omega(2) + dt_omega(3)
       dt_omega_matrix(1,3) = dt_omega(1) * dt_omega(3) - dt_omega(2)

       dt_omega_matrix(2,1) = dt_omega(2) * dt_omega(1) - dt_omega(3)
       dt_omega_matrix(2,2) = ONE + dt_omega(2)**2
       dt_omega_matrix(2,3) = dt_omega(2) * dt_omega(3) + dt_omega(1)

       dt_omega_matrix(3,1) = dt_omega(3) * dt_omega(1) + dt_omega(2)
       dt_omega_matrix(3,2) = dt_omega(3) * dt_omega(2) - dt_omega(1)
       dt_omega_matrix(3,3) = ONE + dt_omega(3)**2

       dt_omega_matrix = dt_omega_matrix / (ONE + dt_omega(1)**2 + dt_omega(2)**2 + dt_omega(3)**2)
       
    endif

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             loc = position(i,j,k) - center
             
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

             call add_momentum_source(loc, unew(i,j,k,UMX:UMZ), Srcorr)
 
             ! Correct energy

             if (rot_source_type == 1) then

               ! If rot_source_type == 1, then calculate SrEcorr before updating the velocities.

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (rot_source_type == 2) then

                ! For this source type, we first update the momenta
                ! before we calculate the energy source term.

                vnew = unew(i,j,k,UMX:UMZ) * rhoninv
                Sr_new = rhon * rotational_acceleration(loc, vnew, time) * dt
                SrE_new = dot_product(vnew, Sr_new)

                SrEcorr = HALF * (SrE_new - SrE_old)

             else if (rot_source_type == 3) then

                ! Instead of calculating the energy source term explicitly,
                ! we simply update the kinetic energy.

                new_ke = HALF * sum(unew(i,j,k,UMX:UMZ)**2) * rhoninv
                SrEcorr = new_ke - old_ke

             else if (rot_source_type == 4) then

                ! Coupled/implicit momentum update (wdmerger paper I; Section 2.4)

                ! Do the full corrector step with the centrifugal force (add 1/2 the new term, subtract 1/2 the old term)
                ! and do the time-level n part of the corrector step for the Coriolis term (subtract 1/2 the old term). 

                new_mom = unew(i,j,k,UMX:UMZ) - dt * cross_product(omega_old, uold(i,j,k,UMX:UMZ)) &
                                              + HALF * dt * loc * omega_new**2 * rhon &
                                              - HALF * dt * loc * omega_old**2 * rhoo                
                
                ! The following is the general solution to the 3D coupled system,
                ! assuming that the rotation vector has components along all three
                ! axes, obtained using Cramer's rule (the coefficient matrix is
                ! defined above). In practice the user will probably only be using
                ! one axis for rotation; if it's the z-axis, then this reduces to
                ! Equations 25 and 26 in the wdmerger paper.

                new_mom = matmul(dt_omega_matrix, new_mom)
                
                Srcorr = new_mom - unew(i,j,k,UMX:UMZ)
                
                call add_momentum_source(loc, unew(i,j,k,UMX:UMZ), Srcorr)
                
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

                ! Note that in the hydrodynamics step, the fluxes used here were already 
                ! multiplied by dA and dt, so dividing by the cell volume at the end is enough to 
                ! get the density change (flux * dt * dA / dV).
                
                SrEcorr = - HALF * ( flux1(i        ,j,k,URHO) * (phi(i,j,k) - phi(i-1,j,k)) - &
                                     flux1(i+1*dg(1),j,k,URHO) * (phi(i,j,k) - phi(i+1,j,k)) + &
                                     flux2(i,j        ,k,URHO) * (phi(i,j,k) - phi(i,j-1,k)) - &
                                     flux2(i,j+1*dg(2),k,URHO) * (phi(i,j,k) - phi(i,j+1,k)) + &
                                     flux3(i,j,k        ,URHO) * (phi(i,j,k) - phi(i,j,k-1)) - &
                                     flux3(i,j,k+1*dg(3),URHO) * (phi(i,j,k) - phi(i,j,k+1)) )

                ! Now normalize by the volume of this cell to get the specific energy change.                
                
                SrEcorr = SrEcorr / vol(i,j,k)
                
                ! Correct for the time rate of change of the potential, which acts 
                ! purely as a source term. For the velocities this is a corrector step
                ! and for the energy we add the full source term.

                Sr_old = - rhoo * cross_product(domegadt_old, loc)
                Sr_new = - rhon * cross_product(domegadt_new, loc)
               
                unew(i,j,k,UMX:UMZ) = unew(i,j,k,UMX:UMZ) + HALF * (Sr_new - Sr_old) * dt

                vnew = unew(i,j,k,UMX:UMZ) / rhon
                                
                SrEcorr = SrEcorr + HALF * (dot_product(vold, Sr_old) + dot_product(vnew, Sr_new)) * dt

                ! Finally, remove the predictor step we applied earlier.
                
                SrEcorr = SrEcorr - dot_product(uold(i,j,k,UMX:UMZ) * rhooinv, Sr_old)
                
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
    endif

    end subroutine ca_corrrsrc



