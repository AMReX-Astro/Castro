subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use amrex_constants_module, only: HALF
  use probdata_module, only: rho_ambient, rho_peak, t_ambient, sigma
  use prob_params_module, only: center
  use castro_error_module, only: castro_error
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: init, namlen
  integer,  intent(in   ) :: name(namlen)
  real(rt), intent(in   ) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

  center(:) = HALF * (problo(:) + probhi(:))

end subroutine amrex_probinit

subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use amrex_constants_module, only: ZERO, HALF
  use probdata_module, only: rho_ambient, rho_peak, t_ambient, sigma
  use eos_type_module, only: eos_input_rt, eos_t
  use eos_module, only: eos
  use network, only: nspec
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS, UTEMP
  use prob_params_module, only: problo, center
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xx, yy, zz
  integer  :: i, j, k

  type(eos_t) :: eos_state

  do k = lo(3), hi(3)
     zz = problo(3) + dx(3) * dble(k + HALF) - center(3)

     do j = lo(2), hi(2)
        yy = problo(2) + dx(2) * dble(j + HALF) - center(2)

        do i = lo(1), hi(1)
           xx = problo(1) + dx(1) * dble(i + HALF) - center(1)

           state(i,j,k,URHO) = rho_ambient + (rho_peak - rho_ambient) * exp(-(xx**2 + yy**2 + zz**2) / sigma**2)
           state(i,j,k,UMX:UMZ) = ZERO

           state(i,j,k,UFS) = state(i,j,k,URHO)

           state(i,j,k,UTEMP) = t_ambient

           eos_state % rho = state(i,j,k,URHO)
           eos_state % T = state(i,j,k,UTEMP)
           eos_state % xn(:) = state(i,j,k,UFS:UFS-1+nspec) / eos_state % rho

           call eos(eos_input_rt, eos_state)

           ! assuming no KE
           state(i,j,k,UEDEN) = eos_state % rho * eos_state % e
           state(i,j,k,UEINT) = eos_state % rho * eos_state % e

        end do
     end do
  end do

end subroutine ca_initdata
