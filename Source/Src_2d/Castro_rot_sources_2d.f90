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

    ! If rot_period is less than zero, that means rotation is disabled, and so we should effectively
    ! shut off the source term by setting omega = 0.

    omega = (/ ZERO, ZERO, ZERO /)

    if (coord_type == 0) then
       if (rot_period > ZERO) then
          omega(3) = TWO * M_PI / rot_period
       endif
    elseif (coord_type == 1) then
       if (rot_period > ZERO) then
          omega(2) = TWO * M_PI / rot_period
       endif
    else
       call bl_error("Error:: Rotate_2d.f90 :: unknown coord_type")
    endif

  end function



  subroutine fill_rotation_field(rot,rot_l1,rot_l2,rot_h1,rot_h2, &
                                 q,q_l1,q_l2,q_h1,q_h2,lo,hi,dx)

    use meth_params_module, only: QVAR, QU, QV, QW
    use prob_params_module, only: problo, center
    use bl_constants_module

    implicit none

    integer         , intent(in   ) :: lo(2), hi(2)
    integer         , intent(in   ) :: rot_l1,rot_l2,rot_h1,rot_h2
    integer         , intent(in   ) :: q_l1,q_l2,q_h1,q_h2

    double precision, intent(inout) :: rot(rot_l1:rot_h1,rot_l2:rot_h2,2)
    double precision, intent(in   ) :: q(q_l1:q_h1,q_l2:q_h2,QVAR)
    double precision, intent(in   ) :: dx(2)

    integer          :: i,j
    double precision :: x,y,r(3)
    double precision :: v(3),omega(3)
    double precision :: omegadotr,omegacrossr(3),omegacrossomegacrossr(3),omegacrossv(3)

    omega = get_omega()

    do j = lo(2)-1, hi(2)+1
       y = problo(2) + dx(2)*(float(j)+HALF) - center(2)
       do i = lo(1)-1, hi(1)+1
          x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

          r = (/ x, y, ZERO /)

          omegacrossr = cross_product(omega,r)
          omegacrossomegacrossr = cross_product(omega,omegacrossr)
          

          v = (/ q(i,j,QU), q(i,j,QV), ZERO /)

          omegacrossv = cross_product(omega,v)

          rot(i,j,1:2) = -TWO * omegacrossv(1:2) - omegacrossomegacrossr(1:2)
             
       enddo
    enddo

  end subroutine fill_rotation_field



  subroutine add_rot_source(uin,uin_l1,uin_l2,uin_h1,uin_h2, &
                            uout,uout_l1,uout_l2,uout_h1,uout_h2, &
                            lo,hi,dx,dt,E_added,xmom_added,ymom_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type
    use prob_params_module, only: coord_type, problo, center
    use bl_constants_module

    implicit none

    integer         , intent(in   ) :: lo(2), hi(2)
    integer         , intent(in   ) :: uin_l1,uin_l2,uin_h1,uin_h2
    integer         , intent(in   ) :: uout_l1,uout_l2,uout_h1,uout_h2

    double precision, intent(in   ) ::  uin( uin_l1: uin_h1, uin_l2: uin_h2,NVAR)
    double precision, intent(inout) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,NVAR)
    double precision, intent(in   ) :: dx(2), dt

    integer          :: i,j
    double precision :: x,y,r(3)
    double precision :: v(3),omega(3)
    double precision :: Sr(3),SrU,SrV,SrW,SrE
    double precision :: dens
    double precision :: omegadotr,omegacrossr(3),omegacrossv(3),omega2

    double precision :: old_rhoeint, new_rhoeint, old_ke, new_ke, old_re
    double precision :: old_xmom, old_ymom, old_zmom
    double precision :: E_added, xmom_added, ymom_added, zmom_added

    omega = get_omega()

    omega2 = dot_product(omega,omega)
    
    do j = lo(2), hi(2)
       y = problo(2) + dx(2)*(float(j)+HALF) - center(2)
       do i = lo(1), hi(1)
          x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

          ! **** Start Diagnostics ****
          old_re = uout(i,j,UEDEN)
          old_ke = HALF * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / &
                            uout(i,j,URHO) 
          old_rhoeint = uout(i,j,UEDEN) - old_ke
          old_xmom = uout(i,j,UMX)
          old_ymom = uout(i,j,UMY)
          ! ****   End Diagnostics ****

          r = (/ x, y, ZERO /)

          dens = uin(i,j,URHO)

          v = (/ uin(i,j,UMX)/dens, &
                 uin(i,j,UMY)/dens, &
                 ZERO /)

          omegacrossr = cross_product(omega,r)
          omegacrossv = cross_product(omega,v)
          omegadotr   = dot_product(omega,r)

          ! momentum sources: this is the Coriolis force
          ! (-2 rho omega x v) and the centrifugal force
          ! (-rho omega x ( omega x r))

          Sr = -TWO * dens * omegacrossv(:) - dens * (omegadotr * omega(:) - omega2 * r(:))

          SrU = Sr(1)
          SrV = Sr(2)

          uout(i,j,UMX) = uout(i,j,UMX) + SrU * dt
          uout(i,j,UMY) = uout(i,j,UMY) + SrV * dt

          ! Kinetic energy source: this is v . the momentum source.
          ! We don't apply in the case of the conservative energy
          ! formulation.

          if (rot_source_type == 1 .or. rot_source_type == 2) then

            SrE = dot_product(v, Sr)

            uout(i,j,UEDEN) = uout(i,j,UEDEN) + SrE * dt

          else if (rot_source_type .eq. 3) then

             new_ke = HALF * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / &
                               uout(i,j,URHO) 
             uout(i,j,UEDEN) = old_rhoeint + new_ke

          else if (rot_source_type .eq. 4) then

             ! Do nothing here, for the conservative gravity option.

          else 
             call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
          end if

          ! **** Start Diagnostics ****
          new_ke = HALF * (uout(i,j,UMX)**2 + uout(i,j,UMY)**2) / &
                            uout(i,j,URHO) 

          new_rhoeint = uout(i,j,UEDEN) - new_ke

          E_added =  E_added + uout(i,j,UEDEN) - old_re

          xmom_added = xmom_added + uout(i,j,UMX) - old_xmom
          ymom_added = ymom_added + uout(i,j,UMY) - old_ymom
          ! ****   End Diagnostics ****

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
                         uold,uold_l1,uold_l2,uold_h1,uold_h2, &
                         unew,unew_l1,unew_l2,unew_h1,unew_h2, &
                         flux1,flux1_l1,flux1_l2,flux1_h1,flux1_h2, &
                         flux2,flux2_l1,flux2_l2,flux2_h1,flux2_h2, &
                         dx,dt, &
                         vol,vol_l1,vol_l2,vol_h1,vol_h2, &
                         xmom_added,ymom_added,E_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, rot_period, rot_source_type
    use prob_params_module, only: coord_type, problo, center
    use bl_constants_module
    use rot_sources_module, only: cross_product, get_omega

    implicit none

    integer         , intent(in   ) :: lo(2), hi(2)
    double precision, intent(in   ) :: dx(2), dt

    integer :: uold_l1,uold_l2,uold_h1,uold_h2
    integer :: unew_l1,unew_l2,unew_h1,unew_h2

    integer :: flux1_l1,flux1_l2,flux1_h1,flux1_h2
    integer :: flux2_l1,flux2_l2,flux2_h1,flux2_h2

    integer :: vol_l1,vol_l2,vol_h1,vol_h2

    double precision :: uold(uold_l1:uold_h1,uold_l2:uold_h2,NVAR)
    double precision :: unew(unew_l1:unew_h1,unew_l2:unew_h2,NVAR)

    double precision :: flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2,NVAR)
    double precision :: flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2,NVAR)

    double precision :: vol(vol_l1:vol_h1,vol_l2:vol_h2)

    integer          :: i,j
    double precision :: x,y,r(3)
    double precision :: vnew(3),vold(3),omega(3)
    double precision :: Sr_old(3), Sr_new(3), SrUcorr, SrVcorr, SrWcorr, SrEcorr, SrE_old, SrE_new
    double precision :: rhoo, rhon, rhooinv, rhoninv
    double precision :: omegadotr,omegacrossr(3),omega2
    double precision :: omegacrossvold(3),omegacrossvnew(3)

    double precision :: old_ke, old_rhoeint, old_re, new_ke, new_rhoeint
    double precision :: old_xmom, old_ymom, old_zmom
    double precision :: E_added, xmom_added, ymom_added

    double precision :: phi(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1)

    omega = get_omega()

    omega2 = dot_product(omega,omega)

    if (rot_source_type == 4) then

       do j = lo(2)-1, hi(2)+1
          y = problo(2) + dx(2)*(float(j)+HALF) - center(2)
          do i = lo(1)-1, hi(1)+1
             x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

             r = (/ x, y, ZERO /)
             omegacrossr = cross_product(omega,r)

             phi(i,j) = - HALF * dot_product(omegacrossr,omegacrossr)

          enddo
       enddo

    endif

    do j = lo(2), hi(2)
       y = problo(2) + dx(2)*(float(j)+HALF) - center(2)
       do i = lo(1), hi(1)
          x = problo(1) + dx(1)*(float(i)+HALF) - center(1)

          ! **** Start Diagnostics ****
          old_re = unew(i,j,UEDEN)
          old_ke = HALF * (unew(i,j,UMX)**2 + unew(i,j,UMY)**2) / &
                            unew(i,j,URHO) 
          old_rhoeint = unew(i,j,UEDEN) - old_ke
          old_xmom = unew(i,j,UMX)
          old_ymom = unew(i,j,UMY)
          ! ****   End Diagnostics ****

          r = (/ x, y, ZERO /)

          omegadotr = dot_product(omega,r)

          ! Define old source terms

          rhoo = uold(i,j,URHO)
          rhooinv = ONE / uold(i,j,URHO)

          vold = (/ uold(i,j,UMX) * rhooinv, &
                    uold(i,j,UMY) * rhooinv, &
                    ZERO /)

          omegacrossvold = cross_product(omega,vold)

          Sr_old = -TWO * rhoo * omegacrossvold(:) - rhoo * (omegadotr * omega(:) - omega2 * r(:))
          SrE_old = dot_product(vold, Sr_old)

          ! Define new source terms

          rhon = unew(i,j,URHO)
          rhoninv = ONE / unew(i,j,URHO)

          vnew = (/ unew(i,j,UMX) * rhoninv, &
                    unew(i,j,UMY) * rhoninv, &
                    ZERO /)

          omegacrossvnew = cross_product(omega,vnew)

          Sr_new = -TWO * rhon * omegacrossvnew(:) - rhon * (omegadotr * omega(:) - omega2 * r(:))

          ! Define correction terms

          SrUcorr = HALF * (Sr_new(1) - Sr_old(1))
          SrVcorr = HALF * (Sr_new(2) - Sr_old(2))

          ! Correct state with correction terms

          unew(i,j,UMX) = unew(i,j,UMX) + SrUcorr * dt
          unew(i,j,UMY) = unew(i,j,UMY) + SrVcorr * dt

          if (rot_source_type == 1) then

            ! If rot_source_type == 1, then calculate SrEcorr before updating the velocities.

             SrE_new = dot_product(vnew, Sr_new)
             SrEcorr = HALF * (SrE_new - SrE_old)

             unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr * dt

          else if (rot_source_type == 2) then

             ! For this source type, we first update the momenta
             ! before we calculate the energy source term.

             vnew(1) = unew(i,j,UMX) * rhoninv
             vnew(2) = unew(i,j,UMY) * rhoninv
             vnew(3) = ZERO

             omegacrossvnew = cross_product(omega,vnew)

             Sr_new = -TWO * rhon * omegacrossvnew(:) - rhon * (omegadotr * omega(:) - omega2 * r(:))
             SrE_new = dot_product(vnew, Sr_new)

             SrEcorr = HALF * (SrE_new - SrE_old)

             unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr * dt

          else if (rot_source_type == 3) then

             ! Instead of calculating the energy source term explicitly,
             ! we simply set the total energy equal to the old internal
             ! energy plus the new kinetic energy.

             new_ke = HALF * (unew(i,j,UMX)**2 + unew(i,j,UMY)**2) / &
                              unew(i,j,URHO) 

             unew(i,j,UEDEN) = old_rhoeint + new_ke

          else if (rot_source_type == 4) then

             ! Conservative energy formulation.

             SrEcorr = HALF * flux1(i  ,j,URHO) * (phi(i  ,j) - phi(i-1,j)) + &
                       HALF * flux1(i+1,j,URHO) * (phi(i+1,j) - phi(i  ,j)) + &
                       HALF * flux2(i,j  ,URHO) * (phi(i,j  ) - phi(i,j-1)) + &
                       HALF * flux2(i,j+1,URHO) * (phi(i,j+1) - phi(i,j  )) 

             SrEcorr = SrEcorr / vol(i,j)

             unew(i,j,UEDEN) = unew(i,j,UEDEN) + SrEcorr

          else 
             call bl_error("Error:: Castro_grav_sources_3d.f90 :: bogus grav_source_type")
          end if

          ! **** Start Diagnostics ****
          ! This is the new (rho e) as stored in (rho E) after the gravitational work is added
          new_ke = HALF * (unew(i,j,UMX)**2 + unew(i,j,UMY)**2) / &
                            unew(i,j,URHO) 
          new_rhoeint = unew(i,j,UEDEN) - new_ke
          E_added =  E_added + unew(i,j,UEDEN) - old_re
          xmom_added = xmom_added + unew(i,j,UMX) - old_xmom
          ymom_added = ymom_added + unew(i,j,UMY) - old_ymom
          ! ****   End Diagnostics ****

       enddo
    enddo

    end subroutine ca_corrrsrc
