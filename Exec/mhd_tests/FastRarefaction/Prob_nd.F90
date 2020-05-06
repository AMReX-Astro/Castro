subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use eos_module
  use eos_type_module
  use castro_error_module
  use network
  use probdata_module
  use amrex_constants_module, only : M_SQRT_PI
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(3), probhi(3)
  real(rt)         xn(nspec)

  type (eos_t) :: eos_state

  call probdata_init(name, namlen)

  split(1) = frac*(problo(1)+probhi(1))
  split(2) = frac*(problo(2)+probhi(2))
  split(3) = frac*(problo(3)+probhi(3))

  ! compute the internal energy (erg/cc) for the left and right state
  xn(:) = 0.0e0_rt
  xn(1) = 1.0e0_rt

  eos_state%rho = rho_l
  eos_state%p = p_l
  eos_state%T = 100000.e0_rt  ! initial guess
  eos_state%xn(:) = xn(:)

  call eos(eos_input_rp, eos_state)

  rhoe_l = rho_l*eos_state%e
  T_l = eos_state%T

  eos_state%rho = rho_r
  eos_state%p = p_r
  eos_state%T = 100000.e0_rt  ! initial guess
  eos_state%xn(:) = xn(:)

  call eos(eos_input_rp, eos_state)

  rhoe_r = rho_r*eos_state%e
  T_r = eos_state%T

end subroutine amrex_probinit


subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec
  use probdata_module
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, UFS
  use prob_params_module, only : problo

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: x, y, z
  integer  :: i, j, k

  do k = lo(3), hi(3)
     z = problo(3) + dx(3)*(dble(k) + 0.5e0_rt)

     do j = lo(2), hi(2)
        y = problo(2) + dx(2)*(dble(j) + 0.5e0_rt)

        do i = lo(1), hi(1)
           x = problo(1) + dx(1)*(dble(i) + 0.5e0_rt)

           if (idir == 1) then
              if (x <= split(1)) then
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = rho_l*u_l_x
                 state(i,j,k,UMY) = rho_l*u_l_y
                 state(i,j,k,UMZ) = rho_l*u_l_z
                 state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*(u_l_x**2+u_l_y**2+u_l_z**2) + 0.5e0_rt * (B_x_l**2 + B_y_l**2 + B_z_l**2)
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = rho_r*u_r_x
                 state(i,j,k,UMY) = rho_r*u_r_y
                 state(i,j,k,UMZ) = rho_r*u_r_z
                 state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*(u_r_x**2+u_r_y**2+u_r_z**2) + 0.5e0_rt *(B_x_r**2 + B_y_r**2 + B_z_r**2)
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              endif

           else if (idir == 2) then
              if (y <= split(2)) then
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = rho_l*u_l_y
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l_y*u_l_y + 0.5e0_rt * (B_x_l**2 + B_y_l**2 + B_z_l**2)
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = rho_r*u_r_y
                 state(i,j,k,UMZ) = 0.e0_rt
                 state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r_y*u_r_y + 0.5e0_rt *(B_x_r**2 + B_y_r**2 + B_z_r**2)
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              endif

           else if (idir == 3) then
              if (z <= split(3)) then
                 state(i,j,k,URHO) = rho_l
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = rho_l*u_l_z
                 state(i,j,k,UEDEN) = rhoe_l + 0.5*rho_l*u_l_z*u_l_z + 0.5e0_rt * (B_x_l**2 + B_y_l**2 + B_z_l**2)
                 state(i,j,k,UEINT) = rhoe_l
                 state(i,j,k,UTEMP) = T_l
              else
                 state(i,j,k,URHO) = rho_r
                 state(i,j,k,UMX) = 0.e0_rt
                 state(i,j,k,UMY) = 0.e0_rt
                 state(i,j,k,UMZ) = rho_r*u_r_z
                 state(i,j,k,UEDEN) = rhoe_r + 0.5*rho_r*u_r_z*u_r_z + 0.5e0_rt *(B_x_r**2 + B_y_r**2 + B_z_r**2)
                 state(i,j,k,UEINT) = rhoe_r
                 state(i,j,k,UTEMP) = T_r
              endif

           else
              call castro_error('invalid idir')
           endif

           state(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
           state(i,j,k,UFS  ) = state(i,j,k,URHO)


        enddo
     enddo
  enddo

end subroutine ca_initdata
