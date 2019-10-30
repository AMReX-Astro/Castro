subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use prob_params_module, only: center
  use probdata_module, only: rho_in, e_in, T_in
  use amrex_fort_module, only: rt => amrex_real
  use network, only : nspec
  use eos_type_module, only : eos_t, eos_input_re
  use eos_module, only : eos
  use meth_params_module, only : T_guess

  implicit none

  integer,  intent(in   ) :: init, namlen
  integer,  intent(in   ) :: name(namlen)
  real(rt), intent(in   ) :: problo(3), probhi(3)

  real(rt) :: xn(nspec)
  type(eos_t) :: eos_state

  call probdata_init(name, namlen)

  ! compute T_in from the input data
  xn(:) = 0.0_rt
  xn(1) = 1.0_rt

  eos_state % rho = rho_in
  eos_state % e = e_in
  eos_state % T = T_guess
  eos_state % xn(:) = xn(:)

  call eos(eos_input_re, eos_state)

  T_in = eos_state % T

  ! set center, domain extrema
  center(:) = (problo(:) + probhi(:)) / 2.e0_rt

end subroutine amrex_probinit

subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use amrex_constants_module, only: ZERO, HALF, ONE
  use probdata_module, only: rho_in, e_in, u_in, T_in
  use network, only: nspec
  use meth_params_module, only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only: center, problo
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  integer  :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           state(i,j,k,URHO) = rho_in

           state(i,j,k,UMX) = rho_in * u_in
           state(i,j,k,UMY) = ZERO
           state(i,j,k,UMZ) = ZERO

           state(i,j,k,UEDEN) = rho_in * e_in + 0.5_rt * rho_in * u_in**2
           state(i,j,k,UEINT) = rho_in * e_in
           state(i,j,k,UTEMP) = T_in

           state(i,j,k,UFS:UFS-1+nspec) = 0.0_rt
           state(i,j,k,UFS) = rho_in

        end do
     end do
  end do

end subroutine ca_initdata

