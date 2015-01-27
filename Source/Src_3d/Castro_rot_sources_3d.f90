module rot_sources_module

  implicit none

  private

  public add_rot_source, cross_product

contains

  subroutine add_rot_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                            uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                            lo,hi,dx,dt,E_added,xmom_added,ymom_added,zmom_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type
    use probdata_module, only: center
    use prob_params_module, only: coord_type, xmin, ymin, zmin
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

    double precision :: E_added, xmom_added, ymom_added, zmom_added

    if (coord_type == 0) then
       ! If rot_period is zero, that means rotation is disabled, and so we should effectively
       ! shut off the source term by setting omega = 0.

       if (rot_period > ZERO) then
          omega = (/ ZERO, ZERO, TWO * M_PI / rot_period /)
       else
          omega = (/ ZERO, ZERO, ZERO /)
       endif
    else
       call bl_error("Error:: Rotate_3d.f90 :: unknown coord_type")
    endif

    omega2 = dot_product(omega,omega)
    
    do k = lo(3), hi(3)
       z = zmin + dx(3)*(float(k)+HALF) - center(3)
       do j = lo(2), hi(2)
          y = ymin + dx(2)*(float(j)+HALF) - center(2)
          do i = lo(1), hi(1)
             x = xmin + dx(1)*(float(i)+HALF) - center(1)

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

             ! kinetic energy source: this is v . the momentum
             ! force -- note that the Coriolis term drops out

             SrE = dot_product(v, Sr)

             uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + SrE * dt

          enddo
       enddo
    enddo

  end subroutine add_rot_source

  function cross_product(x,y) result(r)

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
                         dx,dt,E_added,xmom_added,ymom_added,zmom_added)

    use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, rot_period, rot_source_type
    use probdata_module, only: center
    use prob_params_module, only: coord_type, xmin, ymin, zmin
    use bl_constants_module
    use rot_sources_module, only: cross_product

    implicit none

    integer         , intent(in   ) :: lo(3), hi(3)
    double precision, intent(in   ) :: dx(3), dt

    integer :: uold_l1,uold_l2,uold_l3,uold_h1,uold_h2,uold_h3
    integer :: unew_l1,unew_l2,unew_l3,unew_h1,unew_h2,unew_h3

    double precision :: uold(uold_l1:uold_h1,uold_l2:uold_h2,uold_l3:uold_h3,NVAR)
    double precision :: unew(unew_l1:unew_h1,unew_l2:unew_h2,unew_l3:unew_h3,NVAR)

    integer          :: i,j,k
    double precision :: x,y,z,r(3)
    double precision :: vnew(3),vold(3),omega(3)
    double precision :: Sr_old(3), Sr_new(3), SrUcorr, SrVcorr, SrWcorr, SrEcorr, SrE_old, SrE_new
    double precision :: rhoo, rhon
    double precision :: omegadotr,omegacrossr(3),omega2
    double precision :: omegacrossvold(3),omegacrossvnew(3)

    double precision :: E_added, xmom_added, ymom_added, zmom_added

    if (coord_type == 0) then
       ! If rot_period is zero, that means rotation is disabled, and so we should effectively
       ! shut off the source term by setting omega = 0.

       if (rot_period > ZERO) then
          omega = (/ ZERO, ZERO, TWO * M_PI / rot_period /)
       else
          omega = (/ ZERO, ZERO, ZERO /)
       endif
    else
       call bl_error("Error:: Rotate_3d.f90 :: unknown coord_type")
    endif

    omega2 = dot_product(omega,omega)
    
    do k = lo(3), hi(3)
       z = zmin + dx(3)*(float(k)+HALF) - center(3)
       do j = lo(2), hi(2)
          y = ymin + dx(2)*(float(j)+HALF) - center(2)
          do i = lo(1), hi(1)
             x = xmin + dx(1)*(float(i)+HALF) - center(1)

             r = (/ x, y, z /)

             omegacrossr = cross_product(omega,r)
             omegadotr   = dot_product(omega,r)

             ! Define old source terms

             rhoo = uold(i,j,k,URHO)

             vold = (/ uold(i,j,k,UMX)/rhoo, &
                       uold(i,j,k,UMY)/rhoo, &
                       uold(i,j,k,UMZ)/rhoo /)

             omegacrossvold = cross_product(omega,vold)

             Sr_old = -TWO * rhoo * omegacrossvold(:) - rhoo * (omegadotr * omega(:) - omega2 * r(:))

             ! Define new source terms

             rhon = unew(i,j,k,URHO)
             
             vnew = (/ unew(i,j,k,UMX)/rhon, &
                       unew(i,j,k,UMY)/rhon, &
                       unew(i,j,k,UMZ)/rhon /)

             omegacrossvnew = cross_product(omega,vnew)

             Sr_new = -TWO * rhon * omegacrossvnew(:) - rhon * (omegadotr * omega(:) - omega2 * r(:))

             ! Define correction terms

             SrUcorr = HALF * (Sr_new(1) - Sr_old(1))
             SrVcorr = HALF * (Sr_new(2) - Sr_old(2))
             SrWcorr = HALF * (Sr_new(3) - Sr_old(3))

             ! Correct state with correction terms

             unew(i,j,k,UMX) = unew(i,j,k,UMX) + SrUcorr * dt
             unew(i,j,k,UMY) = unew(i,j,k,UMY) + SrVcorr * dt
             unew(i,j,k,UMZ) = unew(i,j,k,UMZ) + SrWcorr * dt             

             SrE_old = dot_product(vold, Sr_old)
             SrE_new = dot_product(vnew, Sr_new)

             SrEcorr = HALF * (SrE_new - SrE_old)

             unew(i,j,k,UEDEN) = unew(i,j,k,UEDEN) + SrEcorr * dt

          enddo
       enddo
    enddo

    end subroutine ca_corrrsrc



