!-----------------------------------------------------------------------

subroutine derpi(p, p_lo, p_hi, ncomp_p, &
                 u, u_lo, u_hi, ncomp_u, &
                 lo, hi, domlo, domhi, &
                 dx, time) &
                 bind(C, name="derpi")

  ! derive the dynamic pressure

  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use actual_eos_module, only : gamma_const
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
       const_grav, T_guess
  use probdata_module, only: pres_base, dens_base, do_isentropic
  use prob_params_module, only: center
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  real(rt) :: e, T
  real(rt) :: rhoInv
  integer  :: i, j, k, npts_1d
  real(rt) :: H, z, xn(1), const
  real(rt), allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(AMREX_SPACEDIM) + 1.e-8_rt) / dx(AMREX_SPACEDIM)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = ONE

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j = 0, npts_1d-1

     ! initial guess
     temp(j) = T_guess

     if (do_isentropic) then
        z = dble(j) * dx(AMREX_SPACEDIM)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
             (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
     else
        z = (dble(j)+HALF) * dx(AMREX_SPACEDIM)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             dx(AMREX_SPACEDIM) * HALF * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%rho = density(j)
     eos_state%T = temp(j)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  ! Compute pressure from the EOS
  do k = lo(3), hi(3)
     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           rhoInv = ONE/u(i,j,k,URHO)
           T = u(i,j,k,UTEMP)

           eos_state%rho = u(i,j,k,URHO)
           eos_state%T = u(i,j,k,UTEMP)
           eos_state%e = u(i,j,k,UEINT)*rhoInv
           eos_state%xn(:) = u(i,j,k,UFS:UFS-1+nspec)/u(i,j,k,URHO)
           eos_state%aux(:) = u(i,j,k,UFX:UFX-1+naux)/u(i,j,k,URHO)

           if (e <= ZERO) then
              call eos(eos_input_rt, eos_state)
              p(i,j,k,1) = eos_state%p

           else
              call eos(eos_input_re, eos_state)
              p(i,j,k,1) = eos_state%p

           end if
#if AMREX_SPACEDIM == 1
           p(i,j,k,1) = p(i,j,k,1) - pressure(i)
#elif AMREX_SPACEDIM == 2
           p(i,j,k,1) = p(i,j,k,1) - pressure(j)
#else
           p(i,j,k,1) = p(i,j,k,1) - pressure(k)
#endif

        end do
     end do
  end do

end subroutine derpi

!-----------------------------------------------------------------------

subroutine derpioverp0(p, p_lo, p_hi, ncomp_p, &
                       u, u_lo, u_hi, ncomp_u, &
                       lo, hi, domlo, domhi, &
                       dx, time) &
                       bind(C, name="derpioverp0")

  ! derive the ratio of the dynamic pressure to thermodynamic pressure

  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use actual_eos_module, only : gamma_const
  use meth_params_module, only : URHO, UEINT, UTEMP, UFS, UFX, &
       const_grav, T_guess
  use probdata_module, only: pres_base, dens_base, do_isentropic
  use prob_params_module, only: center

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  real(rt) :: e, T
  real(rt) :: rhoInv
  integer :: i, j, k, npts_1d
  real(rt) :: H, z, xn(1), const
  real(rt), allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(AMREX_SPACEDIM) + 1.e-8_rt) / dx(AMREX_SPACEDIM)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = ONE

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j = 0, npts_1d-1

     ! initial guess
     temp(j) = T_guess

     if (do_isentropic) then
        z = dble(j) * dx(AMREX_SPACEDIM)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
             (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
     else
        z = (dble(j)+HALF) * dx(AMREX_SPACEDIM)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             dx(AMREX_SPACEDIM) * HALF * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%rho = density(j)
     eos_state%T = temp(j)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  ! Compute pressure from the EOS
  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rhoInv = ONE/u(i,j,k,URHO)
           e = u(i,j,k,UEINT)*rhoInv
           T = u(i,j,k,UTEMP)

           eos_state%rho = u(i,j,k,URHO)
           eos_state%T = T
           eos_state%xn(:) = u(i,j,k,UFS:UFS-1+nspec)*rhoInv
           eos_state%aux(:) = u(i,j,k,UFX:UFX-1+naux)*rhoInv
           eos_state%e = e

           if (e <= ZERO) then
              call eos(eos_input_rt, eos_state)
              p(i,j,k,1) = eos_state%p

           else
              call eos(eos_input_re, eos_state)
              p(i,j,k,1) = eos_state%p

           end if

#if AMREX_SPACEDIM == 1
           p(i,j,k,1) = (p(i,j,k,1) - pressure(i)) / pressure(i)
#elif AMREX_SPACEDIM == 2
           p(i,j,k,1) = (p(i,j,k,1) - pressure(j)) / pressure(j)
#else
           p(i,j,k,1) = (p(i,j,k,1) - pressure(k)) / pressure(k)
#endif

        end do
     end do
  end do

end subroutine derpioverp0

!-----------------------------------------------------------------------

subroutine derrhopert(p, p_lo, p_hi, ncomp_p, &
                      u, u_lo, u_hi, ncomp_u, &
                      lo, hi, domlo, domhi, &
                      dx, time) &
                      bind(C, name="derrhopert")

  ! derive the perturbational density

  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use network, only : nspec, naux
  use meth_params_module, only : URHO, const_grav, T_guess
  use actual_eos_module, only: gamma_const
  use probdata_module
  use prob_params_module, only : center
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  ! local
  integer :: i, j, k, npts_1d

  real(rt) :: z, dens, H
  real(rt), allocatable :: density(:)

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(AMREX_SPACEDIM) + 1.e-8_rt) / dx(AMREX_SPACEDIM)

  allocate(density (0:npts_1d-1))

  do j = 0, npts_1d-1

     if (do_isentropic) then
        z = dble(j) * dx(AMREX_SPACEDIM)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
             (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
     else
        z = (dble(j)+HALF)*dx(AMREX_SPACEDIM)
        density(j) = dens_base * exp(-z/H)
     end if
  end do

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 1
           p(i,j,k,1) = u(i,j,k,URHO) - density(i)
#elif AMREX_SPACEDIM == 2
           p(i,j,k,1) = u(i,j,k,URHO) - density(j)
#else
           p(i,j,k,1) = u(i,j,k,URHO) - density(k)
#endif
        end do
     end do
  end do

end subroutine derrhopert

!-----------------------------------------------------------------------

subroutine dertpert(p, p_lo, p_hi, ncomp_p, &
                    u, u_lo, u_hi, ncomp_u, &
                    lo, hi, domlo, domhi, &
                    dx, time) &
                    bind(C, name="dertpert")
  ! derive the temperature perturbation

  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use network, only : nspec, naux
  use eos_module
  use eos_type_module
  use actual_eos_module, only : gamma_const
  use meth_params_module, only : UTEMP, const_grav, T_guess
  use probdata_module, only: pres_base, dens_base, do_isentropic
  use prob_params_module, only: center

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)

  integer,  intent(in   ), value :: ncomp_p
  integer,  intent(in   ), value :: ncomp_u
  real(rt), intent(in   ), value :: time

  real(rt) :: e, T
  real(rt) :: rhoInv
  integer  :: i, j, k, npts_1d
  real(rt) :: H, z, xn(1), const
  real(rt), allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(AMREX_SPACEDIM) + 1.e-8_rt) / dx(AMREX_SPACEDIM)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = ONE

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j = 0, npts_1d-1

     ! initial guess
     temp(j) = T_guess

     if (do_isentropic) then
        z = dble(j) * dx(AMREX_SPACEDIM)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - ONE)*z/ &
             (gamma_const*pres_base) + ONE)**(ONE/(gamma_const - ONE))
     else
        z = (dble(j)+HALF) * dx(AMREX_SPACEDIM)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j .gt. 0) then
        pressure(j) = pressure(j-1) - &
             dx(AMREX_SPACEDIM) * HALF * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%rho = density(j)
     eos_state%T = temp(j)
     eos_state%xn(:) = xn(:)
     eos_state%p = pressure(j)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

#if AMREX_SPACEDIM == 1
           p(i,j,k,1) = u(i,j,k,UTEMP) - temp(i)
#elif AMREX_SPACEDIM == 2
           p(i,j,k,1) = u(i,j,k,UTEMP) - temp(j)
#else
           p(i,j,k,1) = u(i,j,k,UTEMP) - temp(k)
#endif

        end do
     end do
  end do

end subroutine dertpert
