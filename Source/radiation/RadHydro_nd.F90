module radhydro_nd_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

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
