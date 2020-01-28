subroutine amrex_probinit (init,name,namlen,problo,probhi) bind(c)

  use amrex_constants_module, only: ZERO, HALF, ONE
  use probdata_module, only: rho0, drho0, xn_zone
  use prob_params_module, only: center
  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in) :: init, namlen
  integer,  intent(in) :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

  ! set explosion center
  center(:) = HALF * (problo(:) + probhi(:))

  xn_zone(:) = ZERO
  xn_zone(1) = ONE

end subroutine amrex_probinit


subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use probdata_module, only: rho0, drho0, xn_zone
  use amrex_constants_module, only: M_PI, FOUR3RD, ZERO, HALF, ONE
  use meth_params_module , only: NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only: center, coord_type, problo
  use amrex_fort_module, only: rt => amrex_real
  use network, only: nspec
  use extern_probin_module, only: eos_gamma

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: xx, yy, zz
  real(rt) :: dist, p, eint
  integer  :: i, j, k

  do k = lo(3), hi(3)
     zz = problo(3) + dx(3) * (dble(k) + HALF)
     do j = lo(2), hi(2)
        yy = problo(2) + dx(2) * (dble(j) + HALF)
        do i = lo(1), hi(1)
           xx = problo(1) + dx(1) * (dble(i) + HALF)

           dist = sqrt((center(1)-xx)**2 + (center(2)-yy)**2 + (center(3)-zz)**2)

           if (dist <= HALF) then
              state(i,j,k,URHO) = rho0 + drho0 * exp(-16.e0_rt*dist**2) * cos(M_PI*dist)**6
           else
              state(i,j,k,URHO) = rho0
           endif

           state(i,j,k,UMX:UMZ) = 0.e0_rt

           ! we are isentropic, so p = (dens/rho0)**Gamma_1
           p = (state(i,j,k,URHO)/rho0)**eos_gamma
           eint = p / (eos_gamma - ONE)

           state(i,j,k,UEDEN) = eint
           state(i,j,k,UEINT) = eint

           state(i,j,k,UFS:UFS-1+nspec) = state(i,j,k,URHO) * xn_zone(:)

        end do
     end do
  end do

end subroutine ca_initdata
