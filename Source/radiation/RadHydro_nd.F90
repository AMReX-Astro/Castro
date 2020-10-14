module radhydro_nd_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, parameter :: rk_order = 3
  logical, parameter :: use_WENO = .false.

  real(rt), parameter :: cfl = 0.5e0_rt

  real(rt), parameter :: onethird=1.e0_rt/3.e0_rt, twothirds=2.e0_rt/3.e0_rt, onesixth=1.e0_rt/6.e0_rt

  ! RK5
  real(rt), parameter :: B1=0.5e0_rt, B2=1.e0_rt/16.e0_rt, B3=0.5e0_rt, B4=9.e0_rt/16.e0_rt, &
                         B5=8.e0_rt/7.e0_rt, B6=7.e0_rt/90.e0_rt
  real(rt), parameter :: C20=5.e0_rt/8.e0_rt , C21=3.e0_rt/8.e0_rt
  real(rt), parameter :: C40=17.e0_rt/8.e0_rt , C41=9.e0_rt/8.e0_rt , C42=-3.e0_rt , C43=0.75e0_rt
  real(rt), parameter :: C50=-5.e0_rt/21.e0_rt , C51=2.e0_rt/7.e0_rt , C52=0.e0_rt , &
                         C53=4.e0_rt , C54=-64.e0_rt/21.e0_rt
  real(rt), parameter :: C60=-8.e0_rt/27.e0_rt , C61=-1.e0_rt/5.e0_rt , C62=32.e0_rt/45.e0_rt , &
                         C63=-32.e0_rt/45.e0_rt, C64=32.e0_rt/27.e0_rt , C65=14.e0_rt/45.e0_rt

contains

  subroutine advect_in_fspace(ustar, af, dt, nstep_fsp)

    use rad_params_module, only: ngroups, ng0, ng1, dlognu
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(inout) :: ustar(0:ngroups-1)
    real(rt), intent(in) :: af(0:ngroups-1)
    real(rt), intent(in) :: dt
    integer, intent(inout) :: nstep_fsp
    integer :: ng2

    !$gpu

    call update_one_species(ngroups, ustar, af, dlognu, dt, nstep_fsp)

  end subroutine advect_in_fspace


  subroutine update_one_species(n, u, a, dx, tend, nstepmax) bind(C, name="update_one_species")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in), value :: n
    real(rt), intent(inout) :: u(0:n-1)
    real(rt), intent(in) :: a(0:n-1), dx(0:n-1)
    real(rt), intent(in), value :: tend
    integer, intent(inout) :: nstepmax

    real(rt) :: dt, acfl
    real(rt) :: f(0:n), u1(0:n-1), u2(0:n-1), u3(0:n-1), u4(0:n-1), u5(0:n-1)
    integer :: i, istep, nstep

    !$gpu

    dt = 1.e50_rt
    do i=0, n-1
       acfl = 1.e-50_rt + abs(a(i))
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
          u1 = u + 0.5e0_rt*dt*dudt(u,a,dx,n)
          u2 = u + 0.5e0_rt*dt*dudt(u1,a,dx,n)
          u3 = u + dt*dudt(u2,a,dx,n)
          u = onethird*(u1+2.e0_rt*u2+u3-u) + onesixth*dt*dudt(u3,a,dx,n)
       else if (rk_order .eq. 3) then
          ! RK3
          u1 = u + dt * dudt(u,a,dx,n)
          u1 = 0.75e0_rt*u + 0.25e0_rt*(u1 + dt*dudt(u1,a,dx,n))
          u = onethird*u + twothirds*(u1 + dt*dudt(u1,a,dx,n))
       else
          ! first-order
          u = u + dt * dudt(u,a,dx,n)
       end if
    end do

    nstepmax = max(nstepmax, nstep)

  end subroutine update_one_species


  function dudt(u,a,dx,n)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer, intent(in) :: n
    real(rt), intent(in) :: u(0:n-1), a(0:n-1), dx(0:n-1)
    real(rt) :: dudt(0:n-1)

    integer :: i
    real(rt) :: f(0:n), ag(-2:n+1), ug(-2:n+1)
    real(rt) :: ul, ur, al, ar, fl, fr, r, a_plus, a_minus
    real(rt) :: fg(-2:n+1), fp(5), fm(5), fpw, fmw, alpha

    !$gpu

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

       f(0) = 0.e0_rt
       do i=1,n-1
          alpha = maxval(ag(i-3:i+2))
          fp = 0.5e0_rt * (fg(i-3:i+1) + alpha*ug(i-3:i+1))
          fm = 0.5e0_rt * (fg(i-2:i+2) - alpha*ug(i-2:i+2))
          call weno5(fp(1),fp(2),fp(3),fp(4),fp(5),fpw)
          call weno5(fm(5),fm(4),fm(3),fm(2),fm(1),fmw)
          f(i) = fpw + fmw
       end do
       f(n) = 0.e0_rt

    else

       ag(-1) = -a(0)
       ag(0:n-1) = a(0:n-1)
       ag(n) = -a(n-1)

       ug(-1) = u(0)
       ug(0:n-1) = u(0:n-1)
       ug(n) = u(n-1)

       f(0) = 0.e0_rt
       do i=1,n-1
          r = (ug(i-1)-ug(i-2)) / (ug(i)-ug(i-1) + 1.e-50_rt)
          ul = ug(i-1) + 0.5e0_rt * (ug(i)-ug(i-1)) * MC(r)

          r = (ag(i-1)-ag(i-2)) / (ag(i)-ag(i-1) + 1.e-50_rt)
          al = ag(i-1) + 0.5e0_rt * (ag(i)-ag(i-1)) * MC(r)

          fl = al*ul

          r = (ug(i) - ug(i-1)) / (ug(i+1) - ug(i) + 1.e-50_rt)
          ur = ug(i) - 0.5e0_rt * (ug(i+1) - ug(i)) * MC(r)

          r = (ag(i) - ag(i-1)) / (ag(i+1) - ag(i) + 1.e-50_rt)
          ar = ag(i) - 0.5e0_rt * (ag(i+1) - ag(i)) * MC(r)

          fr = ar*ur

          a_plus = max(0.e0_rt, al, ar)
          a_minus = max(0.e0_rt, -al, -ar)
          f(i) = (a_plus*fl + a_minus*fr - a_plus*a_minus*(ur-ul)) &
                 / (a_plus + a_minus + 1.e-50_rt)
       end do
       f(n) = 0.e0_rt

    end if

    do i = 0, n-1
       dudt(i) = (f(i) - f(i+1)) / dx(i)
    end do

  end function dudt


  function MC(r) result(MCr)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: r
    real(rt) :: MCr

    !$gpu

    MCr = max(0.e0_rt, min(2.e0_rt*r, 0.5e0_rt*(1.e0_rt+r), 2.e0_rt))

  end function MC


  subroutine weno5(vm2, vm1, v, vp1, vp2, v_weno5)

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: vm2, vm1, v, vp1, vp2
    real(rt), intent(out) :: v_weno5

    real(rt), parameter :: epsw=1.0e-6_rt, b1=13.e0_rt/12.e0_rt, b2=1.e0_rt/6.e0_rt

    real(rt) :: djm1, ejm1, dj, ej, djp1, ejp1, dis0, dis1, dis2, &
                q30, q31, q32, d01, d02, a1ba0, a2ba0, w0, w1, w2

    !$gpu

    djm1 = vm2 - 2.e0_rt*vm1 + v
    ejm1 = vm2 - 4.e0_rt*vm1 + 3.e0_rt*v
    dj   = vm1 - 2.e0_rt*v + vp1
    ej   = vm1 - vp1
    djp1 = v - 2.e0_rt*vp1 + vp2
    ejp1 = 3.e0_rt*v - 4.e0_rt*vp1 + vp2

    dis0 = b1*djm1*djm1 + 0.25e0_rt*ejm1*ejm1 + epsw
    dis1 = b1*dj*dj     + 0.25e0_rt*ej*ej     + epsw
    dis2 = b1*djp1*djp1 + 0.25e0_rt*ejp1*ejp1 + epsw

    q30 = 2.e0_rt*vm2 - 7.e0_rt*vm1 + 11.e0_rt*v
    q31 = -vm1 + 5.e0_rt*v + 2.e0_rt*vp1
    q32 = 2.e0_rt*v + 5.e0_rt*vp1 - vp2

    d01 = dis0 / dis1
    d02 = dis0 / dis2
    a1ba0 = 6.e0_rt * d01 * d01
    a2ba0 = 3.e0_rt * d02 * d02
    w0 = 1.e0_rt / (1.e0_rt + a1ba0 + a2ba0)
    w1 = a1ba0 * w0
    w2 = 1.e0_rt - w0 - w1

    if (w0.lt.1.0e-10_rt) w0 = 0.e0_rt
    if (w1.lt.1.0e-10_rt) w1 = 0.e0_rt
    if (w2.lt.1.0e-10_rt) w2 = 0.e0_rt

    v_weno5 = b2*(w0*q30 + w1*q31 + w2*q32)

  end subroutine weno5


  subroutine inelastic_scatter(temp, u, ks, dt, pt_index)
    ! reference: Larsen, Levermore, Pomraning, and Sanderson, 1985, JCP, 61, 359

    use rad_params_module, only: ngroups, xnu, nugroup, dlognu
    use fundamental_constants_module, only: k_B, m_e, c_light, hplanck
    use lapack_module, only: dgtsv
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    real(rt), intent(in) :: temp, ks, dt
    real(rt), intent(inout) :: u(ngroups)
    integer, intent(in), optional :: pt_index(:)

    integer :: i
    real(rt) :: theta, sigmadt
    integer :: N, NRHS, LDB, INFO
    real(rt) :: DL(ngroups-1), D(ngroups), DU(ngroups-1)
    real(rt) :: B(ngroups,1)
    real(rt) :: ah(2:ngroups), uxh(2:ngroups), bh(2:ngroups), cc(ngroups)

    real(rt), save :: tfac, gamma
    real(rt), allocatable, save :: x(:), xh(:), dlognuinv(:)
    logical, save :: first_call = .true.
!$omp threadprivate(tfac,gamma,x,xh,dlognuinv,first_call)

    if (first_call) then
       first_call = .false.
       tfac = k_B/(m_e*c_light**2)
       gamma = hplanck**2/(8.e0_rt*3.141596565968186e0_rt*(m_e*c_light)**3)
       allocate(x(ngroups))
       allocate(xh(2:ngroups))
       allocate(dlognuinv(ngroups))
       x = nugroup(0:ngroups-1) * (hplanck/(m_e*c_light**2))
       xh = xnu(1:ngroups-1) * (hplanck/(m_e*c_light**2))
       dlognuinv = 1.e0_rt/dlognu
    end if

    theta = temp*tfac
    sigmadt = ks*c_light*dt

    do i = 2, ngroups
       uxh(i) = 0.5e0_rt*(u(i-1)/x(i-1)+u(i)/x(i))
       bh(i) = exp(min(150.e0_rt,(x(i)-x(i-1))/theta))
       ah(i) = sigmadt*(xh(i)**2+gamma*uxh(i))**2 / (bh(i)-1.e0_rt)
    end do

    do i = 1, ngroups
       cc(i) = 1.e0_rt / (x(i)**3+gamma*u(i))
    end do

    B(:,1) = u
    i = 1
    D (i  ) = 1.e0_rt + dlognuinv(i)*cc(i  )*ah(i+1)
    DU(i  ) =         - dlognuinv(i)*cc(i+1)*ah(i+1)*bh(i+1)
    do i = 2, ngroups-1
       DL(i-1) =         - dlognuinv(i)*cc(i-1)*ah(i)
       D (i  ) = 1.e0_rt + dlognuinv(i)*cc(i  )*(ah(i)*bh(i)+ah(i+1))
       DU(i  ) =         - dlognuinv(i)*cc(i+1)*ah(i+1)*bh(i+1)
    end do
    i = ngroups
    DL(i-1) =         - dlognuinv(i)*cc(i-1)*ah(i)
    D (i  ) = 1.e0_rt + dlognuinv(i)*cc(i  )*ah(i)*bh(i)

    N = ngroups
    NRHS = 1
    LDB = ngroups
    call dgtsv(n,nrhs,DL,D,DU,B,LDB,INFO)

    if (INFO .eq. 0) then
       u = B(:,1)
    else
       stop 'inelastic_scatter failed'
    end if

  end subroutine inelastic_scatter

end module radhydro_nd_module
