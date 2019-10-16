subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use castro_error_module
  use probdata_module
  use prob_params_module, only: center

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

  ! set center variable in prob_params_module
  center(1) = frac*(problo(1)+probhi(1))
  center(2) = frac*(problo(2)+probhi(2))
  center(3) = frac*(problo(3)+probhi(3))

end subroutine amrex_probinit


subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       delta, xlo, xhi)

  use amrex_constants_module, only : ZERO, HALF, ONE, TWO
  use probdata_module
  use prob_params_module, only: center
  use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP, const_grav, T_guess
  use eos_module
  use eos_type_module
  use actual_eos_module, only : gamma_const

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, delta(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  integer :: i, j, k, npts_1d
  real(rt) :: H, z, xn(1), x, y, x1, y1, z1, r1, const, T_height
  real(rt), allocatable :: pressure(:), density(:), temp(:), eint(:)

  type (eos_t) :: eos_state

  ! first make a 1D initial model for the entire domain
  npts_1d = (2.e0_rt*center(AMREX_SPACEDIM) + 1.e-8_rt) / delta(AMREX_SPACEDIM)

  allocate(pressure(0:npts_1d-1))
  allocate(density (0:npts_1d-1))
  allocate(temp    (0:npts_1d-1))
  allocate(eint    (0:npts_1d-1))

  const = pres_base/dens_base**gamma_const

  pressure(0) = pres_base
  density(0)  = dens_base

  ! only initialize the first species
  xn(1) = 1.e0_rt

  ! compute the pressure scale height (for an isothermal, ideal-gas
  ! atmosphere)
  H = pres_base / dens_base / abs(const_grav)

  do j = 0, npts_1d-1

     ! initial guess
     temp(j) = T_guess

     if (do_isentropic) then
        z = dble(j) * delta(AMREX_SPACEDIM)
        density(j) = dens_base*(const_grav*dens_base*(gamma_const - 1.0)*z/ &
             (gamma_const*pres_base) + 1.e0_rt)**(1.e0_rt/(gamma_const - 1.e0_rt))
     else
        z = (dble(j)+0.5e0_rt) * delta(AMREX_SPACEDIM)
        density(j) = dens_base * exp(-z/H)
     end if

     if (j > 0) then
        pressure(j) = pressure(j-1) - &
             delta(AMREX_SPACEDIM) * 0.5e0_rt * (density(j)+density(j-1)) * abs(const_grav)
     end if

     eos_state%p = pressure(j)
     eos_state%T = temp(j)
     eos_state%rho = density(j)
     eos_state%xn(:) = xn(:)

     call eos(eos_input_rp, eos_state)

     eint(j) = eos_state%e
     temp(j) = eos_state%T

  end do


  ! add an isobaric perturbation
  x1 = center(1)
  y1 = y_pert_center
  z1 = center(3)

  do k = lo(3), hi(3)
     z = (dble(k)+HALF)*delta(3)

     do j = lo(2), hi(2)
        y = (dble(j)+HALF)*delta(2)

        do i=lo(1),hi(1)
           x = (dble(i)+HALF)*delta(1)

           r1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2) / pert_width

#if AMREX_SPACEDIM == 1
           T_height = temp(i)
#elif AMREX_SPACEDIM == 2
           T_height = temp(j)
#else
           T_height = temp(k)
#endif

           state(i,j,k,UTEMP) = T_height * (ONE + (pert_factor * (ONE + tanh(TWO-r1))))
           state(i,j,k,UFS) = ONE

           eos_state%T = state(i,j,k,UTEMP)
           eos_state%rho = state(i,j,k,URHO)
           eos_state%p = pressure(j)
           eos_state%xn(:) = xn(:)

           call eos(eos_input_tp, eos_state)

           state(i,j,k,URHO) = eos_state%rho
           state(i,j,k,UEINT) = eos_state%e

           ! make state conservative
           state(i,j,k,UFS) = state(i,j,k,UFS)*state(i,j,k,URHO)
           state(i,j,k,UEINT) = state(i,j,k,UEINT)*state(i,j,k,URHO)

           ! assumes ke=0
           state(i,j,k,UEDEN) = state(i,j,k,UEINT)

           state(i,j,k,UMX:UMZ) = ZERO

        end do
     end do
  end do

  deallocate(pressure,density,temp,eint)

end subroutine ca_initdata
