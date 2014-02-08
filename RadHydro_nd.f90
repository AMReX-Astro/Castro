module radhydro_nd_module

  implicit none

  integer, parameter :: rk_order = 3
  logical, parameter :: use_WENO = .false.

  double precision, parameter :: cfl = 0.5d0
  
  double precision, parameter :: onethird=1.d0/3.d0, twothirds=2.d0/3.d0, onesixth=1.d0/6.d0
  
  ! RK5
  double precision, parameter :: B1=0.5d0, B2=1.d0/16.d0, B3=0.5d0, B4=9.d0/16.d0, &
       B5=8.d0/7.d0, B6=7.d0/90.d0
  double precision, parameter :: C20=5.d0/8.d0 , C21=3.d0/8.d0 
  double precision, parameter :: C40=17.d0/8.d0 , C41=9.d0/8.d0 , C42=-3.d0 , C43=0.75d0  
  double precision, parameter :: C50=-5.d0/21.d0 , C51=2.d0/7.d0 , C52=0.d0 , &
       C53=4.d0 , C54=-64.d0/21.d0 
  double precision, parameter :: C60=-8.d0/27.d0 , C61=-1.d0/5.d0 , C62=32.d0/45.d0 , &
       C63=-32.d0/45.d0, C64=32.d0/27.d0 , C65=14.d0/45.d0 
  
  contains

    subroutine advect_in_fspace(ustar, af, dt, nstep_fsp) 
      use rad_params_module, only : ngroups, nnuspec, ng0, ng1, dlognu
      double precision, intent(inout) :: ustar(0:ngroups-1)
      double precision, intent(in) :: af(0:ngroups-1)
      double precision, intent(in) :: dt
      integer, intent(inout) :: nstep_fsp
      integer :: ng2

      if (nnuspec .eq. 0) then
              
         call update_one_species(ngroups, ustar, af, dlognu, dt, nstep_fsp)
              
      else
              
         call update_one_species(ng0, ustar(0:ng0-1), &
              &                          af(0:ng0-1), &
              &                      dlognu(0:ng0-1), &
              dt, nstep_fsp)
              
         if (nnuspec >= 2) then
            call update_one_species(ng1, ustar(ng0:ng0+ng1-1), &
                 &                          af(ng0:ng0+ng1-1), &
                 &                      dlognu(ng0:ng0+ng1-1), &
                 dt, nstep_fsp)
         end if
              
         if (nnuspec == 3) then
            ng2 = ngroups-ng0-ng1
            call update_one_species(ng2, ustar(ng0+ng1:), &
                 &                          af(ng0+ng1:), &
                 &                      dlognu(ng0+ng1:), &
                 dt, nstep_fsp)
         end if

      end if ! end of if nnuspec

    end subroutine advect_in_fspace

    subroutine update_one_species(n, u, a, dx, tend, nstepmax)
      integer, intent(in) :: n
      double precision, intent(inout) :: u(0:n-1)
      double precision, intent(in) :: a(0:n-1), dx(0:n-1)
      double precision, intent(in) :: tend
      integer, intent(inout) :: nstepmax

      double precision :: dt, acfl
      double precision :: f(0:n), u1(0:n-1), u2(0:n-1), u3(0:n-1), u4(0:n-1), u5(0:n-1)
      integer :: i, istep, nstep

      dt = 1.d50
      do i=0, n-1
         acfl = 1.d-50 + abs(a(i))
         dt = min(dt, dx(i)/acfl*cfl)
      end do

      if (dt >= tend) then
         nstep = 1
         dt = tend
      else
         nstep = ceiling(tend/dt)
         dt = tend / dble(nstep)
      end if

      do istep = 1, nstep
         if (rk_order .eq. 5) then
            ! RK5
            u1 = u + B1 * dt * dudt(u,a,dx,n)
            u2 = (C20*u + C21*u1) + B2 * dt * dudt(u1,a,dx,n)
            u3 = u + B3 * dt * dudt(u2,a,dx,n)
            u4 = (C40*u + C41*u1 + C42*u2 + C43*u3) + B4 * dt * dudt(u3,a,dx,n)
            u5 = (C50*u + C51*u1 + C52*u2 + C53*u3 + C54*u4) + B5 * dt * dudt(u4,a,dx,n)
            u  = (C60*u + C61*u1 + C62*u2 + C63*u3 + C64*u4 + C65*u5) & 
                 + B6 * dt * dudt(u5,a,dx,n)
         else if (rk_order .eq. 4) then
            ! RK4
            u1 = u + 0.5d0*dt*dudt(u,a,dx,n)
            u2 = u + 0.5d0*dt*dudt(u1,a,dx,n)
            u3 = u + dt*dudt(u2,a,dx,n)
            u = onethird*(u1+2.d0*u2+u3-u) + onesixth*dt*dudt(u3,a,dx,n)
         else if (rk_order .eq. 3) then
            ! RK3
            u1 = u + dt * dudt(u,a,dx,n)
            u1 = 0.75d0*u + 0.25d0*(u1 + dt*dudt(u1,a,dx,n))
            u = onethird*u + twothirds*(u1 + dt*dudt(u1,a,dx,n))
         else
            ! first-order
            u = u + dt * dudt(u,a,dx,n)
         end if
      end do

      nstepmax = max(nstepmax, nstep)

    end subroutine update_one_species

    function dudt(u,a,dx,n)
      integer, intent(in) :: n
      double precision, intent(in) :: u(0:n-1), a(0:n-1), dx(0:n-1)
      double precision :: dudt(0:n-1)

      integer :: i
      double precision :: f(0:n), ag(-2:n+1), ug(-2:n+1)
      double precision :: ul, ur, al, ar, fl, fr, r, a_plus, a_minus
      double precision :: fg(-2:n+1), fp(5), fm(5), fpw, fmw, alpha
        
      if (use_WENO) then

         ag(-2) = -a(1)
         ag(-1) = -a(0)
         ag(0:n-1) = a(0:n-1)
         ag(n) = -ag(n-1)
         ag(n+1) = -ag(n-2)
     
         ug(-2) = u(1)
         ug(-1) = u(0)
         ug(0:n-1) = u(0:n-1)
         ug(n) = ug(n-1)
         ug(n+1) = ug(n-2)

         fg = ag*ug
         ag = abs(ag)

         f(0) = 0.d0
         do i=1,n-1
            alpha = maxval(ag(i-3:i+2))
            fp = 0.5d0 * (fg(i-3:i+1) + alpha*ug(i-3:i+1))
            fm = 0.5d0 * (fg(i-2:i+2) - alpha*ug(i-2:i+2))
            call weno5(fp(1),fp(2),fp(3),fp(4),fp(5),fpw)
            call weno5(fm(5),fm(4),fm(3),fm(2),fm(1),fmw)
            f(i) = fpw + fmw
         end do
         f(n) = 0.d0

      else

         ag(-1) = -a(0)
         ag(0:n-1) = a(0:n-1)
         ag(n) = -a(n-1)
     
         ug(-1) = u(0)
         ug(0:n-1) = u(0:n-1)
         ug(n) = u(n-1)

         f(0) = 0.d0
         do i=1,n-1
            r = (ug(i-1)-ug(i-2)) / (ug(i)-ug(i-1) + 1.d-50)
            ul = ug(i-1) + 0.5d0 * (ug(i)-ug(i-1)) * MC(r)
            
            r = (ag(i-1)-ag(i-2)) / (ag(i)-ag(i-1) + 1.d-50)
            al = ag(i-1) + 0.5d0 * (ag(i)-ag(i-1)) * MC(r)
            
            fl = al*ul 
            
            r = (ug(i) - ug(i-1)) / (ug(i+1) - ug(i) + 1.d-50)
            ur = ug(i) - 0.5d0 * (ug(i+1) - ug(i)) * MC(r)
            
            r = (ag(i) - ag(i-1)) / (ag(i+1) - ag(i) + 1.d-50)
            ar = ag(i) - 0.5d0 * (ag(i+1) - ag(i)) * MC(r)
            
            fr = ar*ur 
            
            a_plus = max(0.d0, al, ar)
            a_minus = max(0.d0, -al, -ar)
            f(i) = (a_plus*fl + a_minus*fr - a_plus*a_minus*(ur-ul)) &
                 / (a_plus + a_minus + 1.d-50)
         end do
         f(n) = 0.d0

      end if
         
      do i = 0, n-1
         dudt(i) = (f(i) - f(i+1)) / dx(i)
      end do

    end function dudt

    ! function dm3(x,y,z)
    !   double precision :: x, y, z, dm3
    !   dm3 = 0.25d0 * (sign(1.0d0,x)+sign(1.0d0,y)) * abs(sign(1.0d0,x)+sign(1.0d0,z)) &
    !        * min(abs(x),abs(y),abs(z)) 
    ! end function dm3

    ! function superbee(r)
    !   double precision, intent(in) :: r
    !   double precision :: superbee
    !   superbee = max(0.d0, min(2.d0*r,1.d0), min(r,2.d0))
    ! end function superbee

    function MC(r)
      double precision, intent(in) :: r
      double precision :: MC
      MC = max(0.d0, min(2.d0*r, 0.5d0*(1.d0+r), 2.d0))
    end function MC

    ! function minmod(r)
    !   double precision, intent(in) :: r
    !   double precision :: minmod
    !   minmod = max(0.d0, min(1.d0, r))
    ! end function MINMOD

    subroutine weno5(vm2, vm1, v, vp1, vp2, v_weno5)

      implicit none

      double precision, intent(in) :: vm2, vm1, v, vp1, vp2
      double precision, intent(out) :: v_weno5 

      double precision, parameter :: epsw=1.0d-6, b1=13.d0/12.d0, b2=1.d0/6.d0 

      double precision :: djm1, ejm1, dj, ej, djp1, ejp1, dis0, dis1, dis2, &
           q30, q31, q32, d01, d02, a1ba0, a2ba0, w0, w1, w2

      djm1 = vm2 - 2.d0*vm1 + v
      ejm1 = vm2 - 4.d0*vm1 + 3.d0*v
      dj   = vm1 - 2.d0*v + vp1
      ej   = vm1 - vp1
      djp1 = v - 2.d0*vp1 + vp2
      ejp1 = 3.d0*v - 4.d0*vp1 + vp2

      dis0 = b1*djm1*djm1 + 0.25d0*ejm1*ejm1 + epsw
      dis1 = b1*dj*dj     + 0.25d0*ej*ej     + epsw
      dis2 = b1*djp1*djp1 + 0.25d0*ejp1*ejp1 + epsw

      q30 = 2.d0*vm2 - 7.d0*vm1 + 11.d0*v
      q31 = -vm1 + 5.d0*v + 2.d0*vp1
      q32 = 2.d0*v + 5.d0*vp1 - vp2

      d01 = dis0 / dis1
      d02 = dis0 / dis2
      a1ba0 = 6.d0 * d01 * d01
      a2ba0 = 3.d0 * d02 * d02
      w0 = 1.d0 / (1.d0 + a1ba0 + a2ba0)
      w1 = a1ba0 * w0
      w2 = 1.d0 - w0 - w1 
    
      if (w0.lt.1.0d-10) w0 = 0.d0
      if (w1.lt.1.0d-10) w1 = 0.d0
      if (w2.lt.1.0d-10) w2 = 0.d0

      v_weno5 = b2*(w0*q30 + w1*q31 + w2*q32)

      return 

    end subroutine weno5

end module radhydro_nd_module
