subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use eos_module
  use eos_type_module
  use castro_error_module
  use network
  use probdata_module
  use amrex_constants_module, only : M_SQRT_PI, TWO
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer init, namlen
  integer name(namlen)
  real(rt)         problo(3), probhi(3)
  real(rt)         xn(nspec)

  type (eos_t) :: eos_state

  call probdata_init(name, namlen)

  !temporarily give B_0 value here
  B_0 = 1.0 / (TWO*M_SQRT_PI)

  ! compute the internal energy (erg/cc) for the left and right state
  xn(:) = 0.0e0_rt
  xn(1) = 1.0e0_rt

  eos_state%rho = rho_0
  eos_state%p = p_0
  eos_state%T = 100000.e0_rt  ! initial guess
  eos_state%xn(:) = xn(:)

  call eos(eos_input_rp, eos_state)

  rhoe_0 = rho_0*eos_state%e
  T_0 = eos_state%T

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
  use amrex_constants_module, only : M_PI, TWO, FOUR

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

           state(i,j,k,URHO) = rho_0
           state(i,j,k,UMX) = -rho_0*u_0*sin(TWO*M_PI*y)
           state(i,j,k,UMY) = rho_0*u_0*sin(TWO*M_PI*x)
           state(i,j,k,UMZ) = 0.0e0_rt
           state(i,j,k,UEDEN) = rhoe_0 + 0.5*rho_0*((u_0*sin(TWO*M_PI*y))**2+(u_0*sin(TWO*M_PI*x))**2) 
           state(i,j,k,UEINT) = rhoe_0
           state(i,j,k,UTEMP) = T_0
          
           state(i,j,k,UFS:UFS-1+nspec) = 0.0e0_rt
           state(i,j,k,UFS  ) = state(i,j,k,URHO)


        enddo
     enddo
  enddo

end subroutine ca_initdata
