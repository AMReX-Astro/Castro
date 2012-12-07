module sponge_module

  implicit none

  private

  public sponge

contains

  subroutine sponge(uout,uout_l1,uout_l2,uout_l3,&
       uout_h1,uout_h2,uout_h3,lo,hi,t,dt, &
       dx,dy,dz,domlo,domhi)

    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN
    use probdata_module

    implicit none
    integer          :: lo(3),hi(3),domlo(3),domhi(3)
    integer          :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision :: t,dt,dx,dy,dz

    integer :: i,j,k,n

    double precision :: sponge_kappa, smdamp, sponge_mult, x, y, z, radius, r_sp, r_tp
    double precision :: r, r_guess, dr, ke_old, ke_new

    logical :: converged

    double precision, parameter :: TOL = 1.d-6
    integer,          parameter :: max_iter = 25

    sponge_kappa = 1000.d0

    ! we first compute the radius of the edge of the dust

    ! find the analytic solution via Newton iteration
    converged = .false.
    r_guess = 0.95*r_old

    do n = 1, max_iter
       dr = -f(r_guess,t)/dfdr(r_guess,t)

       ! we have lots of sqrt(1 - r/r_0), so r(t) has to be less than r_0
       if (r_guess + dr > 6.5d8) then
          r_guess = 0.5d0*(r_guess + 6.5d8)
       else
          r_guess = r_guess + dr
       endif

       if (abs(dr/r_guess) < TOL) then
          converged = .true.
          r = r_guess
          r_old = r_guess
          exit
       endif

    enddo

    if (.not. converged) then
       print *, "ERROR: Newton iterations failed to converge"
       stop
    endif

    ! set the starting point to be 2.5d7 behind
    ! set the end point to be 5.0d7 behind
    r_sp = r + 2.5d7
    r_tp = r + 5.0d7

    do k = lo(3),hi(3)
       z = (dble(k)+0.5d0)*dz - center(3)

       do j = lo(2),hi(2)
          y = (dble(j)+0.5d0)*dy - center(2)

          do i = lo(1),hi(1)
             x = (dble(i)+0.5d0)*dx - center(1)

             ! compute radius
             radius = sqrt(x**2 + y**2 + z**2)

             ! apply sponge
             if (radius .ge. r_sp) then
                if (radius .lt. r_tp) then
                   smdamp = 0.5d0*(1.d0 - &
                        cos(3.14159265358979323846d0*(radius - r_sp)/(r_tp - r_sp)))
                else
                   smdamp = 1.d0
                endif
                sponge_mult = 1.d0 / (1.d0 + dt * smdamp * sponge_kappa)

                ke_old = 0.5d0 * ( uout(i,j,k,UMX)**2+uout(i,j,k,UMY)**2+uout(i,j,k,UMZ)**2) &
                     / uout(i,j,k,URHO)

                uout(i,j,k,UMX  ) = uout(i,j,k,UMX  ) * sponge_mult
                uout(i,j,k,UMY  ) = uout(i,j,k,UMY  ) * sponge_mult
                uout(i,j,k,UMZ  ) = uout(i,j,k,UMZ  ) * sponge_mult

                ke_new = 0.5d0 * ( uout(i,j,k,UMX)**2+uout(i,j,k,UMY)**2+uout(i,j,k,UMZ)**2) &
                     / uout(i,j,k,URHO)

                uout(i,j,k,UEDEN) = uout(i,j,k,UEDEN) + (ke_new-ke_old)
             endif

          enddo
       enddo

    enddo

  end subroutine sponge

  function f(r,t) result (func)

    implicit none

    double precision, intent(in) :: r, t
    double precision :: func

    ! the analytic solution from Colgate and White, Eq. 5 (with everything moved
    ! onto the LHS).  We use 1.d-20 to prevent against zeros in some places.
    func = sqrt(8.d0*3.14159265358979323846d0*6.67428d-8*1.d9/3.d0)*t - &
         sqrt(1.d0 - r/6.5d8)*sqrt(r/6.5d8) - asin(sqrt(1.0 - r/6.5d8 + 1.d-20))

    return
  end function f



  function dfdr(r,t) result (dfuncdr)

    implicit none

    double precision, intent(in) :: r, t
    double precision :: dfuncdr

    ! the derivative (wrt x) of the analytic function f, defined above.
    ! We use 1.d-20 to prevent against zeros in some places.
    dfuncdr = (0.5d0/6.5d8)*(-sqrt(1.d0 - r/6.5d8 + 1.d-20)/sqrt(r/6.5d8) + &
         sqrt(r/6.5d8)/sqrt(1 - r/6.5d8 + 1.d-20) + &
         1.d0/(sqrt(r/6.5d8)*sqrt(1.d0 - r/6.5d8 + 1.d-20)))


    return
  end function dfdr

end module sponge_module
