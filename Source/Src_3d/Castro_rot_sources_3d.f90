module rot_sources_module

  implicit none

  private

  public add_rot_source, cross_product, fill_rotation_field, get_omega

contains

  function get_omega() result(omega)

    use prob_params_module, only: coord_type
    use meth_params_module, only: rot_period, rot_axis
    use bl_constants_module, only: ZERO, TWO, M_PI

    implicit none

    double precision :: omega(3)

    if (coord_type == 0) then
       ! If rot_period is less than zero, that means rotation is disabled, and so we should effectively
       ! shut off the source term by setting omega = 0.

       omega = (/ ZERO, ZERO, ZERO /)

       if (rot_period > ZERO) then
          omega(rot_axis) = TWO * M_PI / rot_period
       endif
    else
       call bl_error("Error:: Rotate_3d.f90 :: unknown coord_type")
    endif

  end function



  subroutine fill_rotation_field(rot,rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3, &
                                 q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,lo,hi,dx)

    ! fill_rotation_field returns the sources to the velocity
    ! equations (not the conserved momentum equations) that are used
    ! in predicting the interface states
    use meth_params_module, only: QVAR, QU, QV, QW
    use prob_params_module, only: problo, center
    use bl_constants_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: rot_l1,rot_l2,rot_l3,rot_h1,rot_h2,rot_h3
    integer         , intent(in   ) :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3

    double precision, intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,rot_l3:rot_h3,3)
    double precision, intent(in   ) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
    double precision, intent(in   ) :: dx(3)

    integer          :: i,j,k
    double precision :: x,y,z,r(3)
    double precision :: v(3),omega(3)
    double precision :: omegadotr,omegacrossr(3),omegacrossomegacrossr(3),omegacrossv(3)

    omega = get_omega()

    do k = lo(3)-1, hi(3)+1
       z = problo(3) + dx(3)*(float(k)+HALF) - center(3)

       do j = lo(2)-1, hi(2)+1
          y = problo(2) + dx(2)*(float(j)+HALF) - center(2)

          do i = lo(1)-1, hi(1)+1
             x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

             r = (/ x, y, z /)

             omegacrossr = cross_product(omega,r)
             omegacrossomegacrossr = cross_product(omega,omegacrossr)

             v = (/ q(i,j,k,QU), q(i,j,k,QV), q(i,j,k,QW) /)

             omegacrossv = cross_product(omega,v)

             rot(i,j,k,:) = -TWO * omegacrossv - omegacrossomegacrossr
             
          enddo
       enddo
    enddo

  end subroutine fill_rotation_field


  subroutine add_rot_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                            uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                            lo,hi,dx,dt,E_added,xmom_added,ymom_added,zmom_added)

    ! Here we add dt * S_rot^n -- the time-level n rotation source to
    ! the momentum equation.  Note that uin here is the state at time
    ! n -- no sources should have been added to it (including
    ! gravity).
    
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type
    use prob_params_module, only: coord_type, problo, center
    use bl_constants_module

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    integer         , intent(in   ) :: uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
    integer         , intent(in   ) :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3

    double precision, intent(in   ) ::  uin( uin_l1: uin_h1, uin_l2: uin_h2, uin_l3: uin_h3,NVAR)
    double precision, intent(inout) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision, intent(in   ) :: dx(3), dt

    integer          :: i,j,k
    double precision :: x,y,z,r(3)
    double precision :: v(3),omega(3)
    double precision :: Sr(3),SrU,SrV,SrW,SrE
    double precision :: dens
    double precision :: omegadotr,omegacrossr(3),omegacrossv(3),omega2

    double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
    double precision :: old_xmom, old_ymom, old_zmom
    double precision :: E_added, xmom_added, ymom_added, zmom_added

    omega = get_omega()

    omega2 = dot_product(omega,omega)
    
    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(float(k)+HALF) - center(3)

       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(float(j)+HALF) - center(2)

          do i = lo(1), hi(1)
             x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

             ! **** Start Diagnostics ****
             old_re = uout(i,j,k,UEDEN)
             old_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                               uout(i,j,k,URHO) 
             old_rhoeint = uout(i,j,k,UEDEN) - old_ke
             old_xmom = uout(i,j,k,UMX)
             old_ymom = uout(i,j,k,UMY)
             old_zmom = uout(i,j,k,UMZ)
             ! ****   End Diagnostics ****

             r = (/ x, y, z /)

             dens = uin(i,j,k,URHO)

             v = (/ uin(i,j,k,UMX)/dens, &
                    uin(i,j,k,UMY)/dens, &
                    uin(i,j,k,UMZ)/dens /)

             omegacrossr = cross_product(omega,r)
             omegacrossv = cross_product(omega,v)
             omegadotr   = dot_product(omega,r)

             ! momentum sources: this is the Coriolis force
             ! (-2 rho omega x v) and the centrifugal force
             ! (-rho omega x ( omega x r))

             Sr = -TWO * dens * omegacrossv(:) - dens * (omegadotr * omega(:) - omega2 * r(:))

             SrU = Sr(1)
             SrV = Sr(2)
             SrW = Sr(3)

             uout(i,j,k,UMX) = uout(i,j,k,UMX) + SrU * dt
             uout(i,j,k,UMY) = uout(i,j,k,UMY) + SrV * dt
             uout(i,j,k,UMZ) = uout(i,j,k,UMZ) + SrW * dt

             ! Kinetic energy source: this is v . the momentum source.
             ! We don't apply in the case of the conservative energy
             ! formulation.

             if (rot_source_type == 1 .or. rot_source_type == 2) then

               SrE = dot_product(v, Sr)

               uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + SrE * dt

             else if (rot_source_type .eq. 3) then

                new_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                                  uout(i,j,k,URHO) 
                uout(i,j,k,UEDEN) = old_rhoeint + new_ke

             else if (rot_source_type .eq. 4) then

                ! Do nothing here, for the conservative rotation option.

             else 
                call bl_error("Error:: Castro_rot_sources_3d.f90 :: bogus rot_source_type")
             end if

             ! **** Start Diagnostics ****
             new_ke = HALF * (uout(i,j,k,UMX)**2 + uout(i,j,k,UMY)**2 + uout(i,j,k,UMZ)**2) / &
                               uout(i,j,k,URHO) 

             new_rhoeint = uout(i,j,k,UEDEN) - new_ke
 
             E_added =  E_added + uout(i,j,k,UEDEN) - old_re

             xmom_added = xmom_added + uout(i,j,k,UMX) - old_xmom
             ymom_added = ymom_added + uout(i,j,k,UMY) - old_ymom
             zmom_added = zmom_added + uout(i,j,k,UMZ) - old_zmom
             ! ****   End Diagnostics ****

          enddo
       enddo
    enddo

  end subroutine add_rot_source

  function cross_product(x,y) result(r)

    implicit none

    double precision :: x(3), y(3)
    double precision :: r(3)

    r(1) = x(2)*y(3) - x(3)*y(2)
    r(2) = x(3)*y(1) - x(1)*y(3)
    r(3) = x(1)*y(2) - x(2)*y(1)

  end function cross_product

end module rot_sources_module



  subroutine ca_corrrsrc(lo,hi, &
                         uold,uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3, &
                         unew,unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3, &
                         flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
                         flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
                         flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
                         dx,dt, &
                         vol,vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3, &
                         xmom_added,ymom_added,zmom_added,E_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type, rot_axis
    use prob_params_module, only: coord_type, problo, center
    use bl_constants_module
    use rot_sources_module, only: cross_product, get_omega

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3), dt

    integer :: uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
    integer :: unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3

    integer :: flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
    integer :: flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
    integer :: flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3

    integer :: vol_l1,vol_l2,vol_l3,vol_h1,vol_h2,vol_h3

    double precision :: uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
    double precision :: unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)

    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,flux1_l3:flux1_h3,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,flux2_l3:flux2_h3,NVAR)
    double precision :: flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2,flux3_l3:flux3_h3,NVAR)

    double precision :: vol(vol_l1:vol_h1,vol_l2:vol_h2,vol_l3:vol_h3)

    integer          :: i,j,k
    double precision :: x,y,z,r(3)
    double precision :: vnew(3),vold(3),omega(3)
    double precision :: Sr_old(3), Sr_new(3), SrUcorr, SrVcorr, SrWcorr, SrEcorr, SrE_old, SrE_new
    double precision :: rhoo, rhon, rhooinv, rhoninv
    double precision :: omegadotr,omegacrossr(3),omega2
    double precision :: omegacrossvold(3),omegacrossvnew(3)

    double precision :: old_ke, old_rhoeint, old_re, new_ke, new_rhoeint
    double precision :: old_xmom, old_ymom, old_zmom
    double precision :: E_added, xmom_added, ymom_added, zmom_added

    double precision :: phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    double precision :: mom1, mom2

    integer :: idir1, idir2, midx1, midx2

    omega = get_omega()

    omega2 = dot_product(omega,omega)

    if (rot_source_type == 4) then

       do k = lo(3)-1, hi(3)+1
          z = problo(3) + dx(3)*(float(k)+HALF) - center(3)
          do j = lo(2)-1, hi(2)+1
             y = problo(2) + dx(2)*(float(j)+HALF) - center(2)
             do i = lo(1)-1, hi(1)+1
                x = problo(1) + dx(1)*(float(i)+HALF) - center(1)
                
                r = (/ x, y, z /)
                omegacrossr = cross_product(omega,r)

                phi(i,j,k) = - HALF * dot_product(omegacrossr,omegacrossr)

             enddo
          enddo
       enddo

    endif

    do k = lo(3), hi(3)
       z = problo(3) + dx(3)*(float(k)+HALF) - center(3)
       do j = lo(2), hi(2)
          y = problo(2) + dx(2)*(float(j)+HALF) - center(2)
          do i = lo(1), hi(1)
             x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

             ! **** Start Diagnostics ****
             old_re = unew(i,j,k,UEDEN)
             old_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                               unew(i,j,k,URHO) 
             old_rhoeint = unew(i,j,k,UEDEN) - old_ke
             old_xmom = unew(i,j,k,UMX)
             old_ymom = unew(i,j,k,UMY)
             old_zmom = unew(i,j,k,UMZ)
             ! ****   End Diagnostics ****

             r = (/ x, y, z /)

             omegadotr = dot_product(omega,r)

             ! Define old source terms

             rhoo = uold(i,j,k,URHO)
             rhooinv = ONE / uold(i,j,k,URHO)

             vold = (/ uold(i,j,k,UMX) * rhooinv, &
                       uold(i,j,k,UMY) * rhooinv, &
                       uold(i,j,k,UMZ) * rhooinv /)

             omegacrossvold = cross_product(omega,vold)

             Sr_old = -TWO * rhoo * omegacrossvold(:) - rhoo * (omegadotr * omega(:) - omega2 * r(:))
             SrE_old = dot_product(vold, Sr_old)

             ! Define new source terms

             rhon = unew(i,j,k,URHO)
             rhoninv = ONE / unew(i,j,k,URHO)
             
             vnew = (/ unew(i,j,k,UMX) * rhoninv, &
                       unew(i,j,k,UMY) * rhoninv, &
                       unew(i,j,k,UMZ) * rhoninv /)

             omegacrossvnew = cross_product(omega,vnew)

             Sr_new = -TWO * rhon * omegacrossvnew(:) - rhon * (omegadotr * omega(:) - omega2 * r(:))

             ! Define correction terms

             SrUcorr = HALF * (Sr_new(1) - Sr_old(1))
             SrVcorr = HALF * (Sr_new(2) - Sr_old(2))
             SrWcorr = HALF * (Sr_new(3) - Sr_old(3))

             if (rot_source_type .le. 3) then

                ! Correct state with correction terms

                unew(i,j,k,UMX) = unew(i,j,k,UMX) + SrUcorr * dt
                unew(i,j,k,UMY) = unew(i,j,k,UMY) + SrVcorr * dt
                unew(i,j,k,UMZ) = unew(i,j,k,UMZ) + SrWcorr * dt

             endif

             if (rot_source_type == 1) then

               ! If rot_source_type == 1, then calculate SrEcorr before updating the velocities.

                SrE_new = dot_product(vnew, Sr_new)
                SrEcorr = HALF * (SrE_new - SrE_old)

                unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr * dt

             else if (rot_source_type == 2) then

                ! For this source type, we first update the momenta
                ! before we calculate the energy source term.

                vnew(1) = unew(i,j,k,UMX) * rhoninv
                vnew(2) = unew(i,j,k,UMY) * rhoninv
                vnew(3) = unew(i,j,k,UMZ) * rhoninv

                omegacrossvnew = cross_product(omega,vnew)

                Sr_new = -TWO * rhon * omegacrossvnew(:) - rhon * (omegadotr * omega(:) - omega2 * r(:))
                SrE_new = dot_product(vnew, Sr_new)

                SrEcorr = HALF * (SrE_new - SrE_old)

                unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr * dt

             else if (rot_source_type == 3) then

                ! Instead of calculating the energy source term explicitly,
                ! we simply set the total energy equal to the old internal
                ! energy plus the new kinetic energy.

                new_ke = HALF * (unew(i,j,k,UMX)**2 + unew(i,j,k,UMY)**2 + unew(i,j,k,UMZ)**2) / &
                                 unew(i,j,k,URHO) 
 
                unew(i,j,k,UEDEN) = old_rhoeint + new_ke

             else if (rot_source_type == 4) then

                ! Coupled momentum update.

                ! Figure out which directions are updated, and then determine the right 
                ! array index relative to UMX (this works because UMX, UMY, UMZ are consecutive
                ! in the state array).

                idir1 = 1 + MOD(rot_axis    , 3)
                idir2 = 1 + MOD(rot_axis + 1, 3)

                midx1 = UMX + idir1 - 1
                midx2 = UMX + idir2 - 1

                ! Update with the centrifugal term and the time-level n Coriolis term.

                unew(i,j,k,midx1) = unew(i,j,k,midx1) + HALF * dt * omega(rot_axis)**2 * r(idir1) * (rhon - rhoo) &
                                + dt * omega(rot_axis) * uold(i,j,k,midx2)
                unew(i,j,k,midx2) = unew(i,j,k,midx2) + HALF * dt * omega(rot_axis)**2 * r(idir2) * (rhon - rhoo) & 
                                - dt * omega(rot_axis) * uold(i,j,k,midx1)

                ! Now do the implicit solve for the time-level n+1 Coriolis term.

                mom1 = unew(i,j,k,midx1)
                mom2 = unew(i,j,k,midx2)

                unew(i,j,k,midx1) = (mom1 + dt * omega(rot_axis) * mom2) / (ONE + (dt * omega(rot_axis))**2)
                unew(i,j,k,midx2) = (mom2 - dt * omega(rot_axis) * mom1) / (ONE + (dt * omega(rot_axis))**2)

                ! Conservative energy formulation.
                ! note that the fluxes here have already been normalized by A (the cell face)
                ! and dt
                SrEcorr = HALF * flux1(i  ,j,k,URHO) * (phi(i  ,j,k) - phi(i-1,j,k)) + &
                          HALF * flux1(i+1,j,k,URHO) * (phi(i+1,j,k) - phi(i  ,j,k)) + &
                          HALF * flux2(i,j  ,k,URHO) * (phi(i,j,  k) - phi(i,j-1,k)) + &
                          HALF * flux2(i,j+1,k,URHO) * (phi(i,j+1,k) - phi(i,j  ,k)) + &
                          HALF * flux3(i,j,k  ,URHO) * (phi(i,j,k  ) - phi(i,j,k-1)) + &
                          HALF * flux3(i,j,k+1,URHO) * (phi(i,j,k+1) - phi(i,j,k  ))

                SrEcorr = SrEcorr / vol(i,j,k)

                unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr

             else 
                call bl_error("Error:: Castro_rot_sources_3d.f90 :: bogus rot_source_type")
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

    end subroutine ca_corrrsrc



