subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use meth_params_module, only: small_temp, small_pres, small_dens, small_ener
  use amrex_fort_module, only : rt => amrex_real
  use probdata_module
  use eos_type_module, only : eos_t, eos_input_rt
  use eos_module, only : eos
  use network, only : nspec

  implicit none

  integer, intent(in)  :: init, namlen
  integer, intent(in)  :: name(namlen)
  real(rt), intent(in) :: problo(3), probhi(3)

  type (eos_t) :: eos_state

  call probdata_init(name, namlen)

  ! Given the inputs of small_dens and small_temp, figure out small_pres.

  if (small_dens > 0.0e0_rt .and. small_temp > 0.0e0_rt) then
     eos_state % rho = small_dens
     eos_state % T   = small_temp
     eos_state % xn  = 1.0e0_rt / nspec

     call eos(eos_input_rt, eos_state)

     small_pres = eos_state % p
     small_ener = eos_state % e
  endif

end subroutine amrex_probinit



subroutine ca_initdata(lo, hi, &
                       state, state_lo, state_hi, &
                       dx, problo) bind(C, name='ca_initdata')

  use probdata_module
  use eos_module, only: eos
  use eos_type_module, only: eos_t, eos_input_re
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UTEMP, &
                                 UEDEN, UEINT, UFS
  use network, only : nspec
  use amrex_constants_module, only: ZERO, HALF, ONE, TWO, M_PI
  use fundamental_constants_module, only: Gconst, M_solar
  use prob_params_module, only: center
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: state_lo(3), state_hi(3)
  real(rt), intent(inout) :: state(state_lo(1):state_hi(1),state_lo(2):state_hi(2),state_lo(3):state_hi(3),NVAR)
  real(rt), intent(in   ) :: dx(3), problo(3)

  real(rt) :: loc(3)
  real(rt) :: radius

  type (eos_t) :: zone_state

  integer :: i,j,k

  !$gpu

  ! Loop through the zones and set the zone state depending on whether we are
  ! inside the sphere or if we are in an ambient zone.

  !$OMP PARALLEL DO PRIVATE(i, j, k, loc, radius, zone_state)
  do k = lo(3), hi(3)
     loc(3) = problo(3) + dx(3) * (dble(k)+HALF)

     do j = lo(2), hi(2)
        loc(2) = problo(2) + dx(2) * (dble(j)+HALF)

        do i = lo(1), hi(1)
           loc(1) = problo(1) + dx(1) * (dble(i)+HALF)

           radius = sum( (loc - center)**2 )**HALF

           if (radius <= sphere_radius) then
              zone_state % rho = (sphere_mass * M_solar) / (TWO * M_PI * sphere_radius**2 * radius)
           else
              zone_state % rho = ambient_density
           endif

           zone_state % e     = 0.05 * Gconst * sphere_mass * m_solar / radius
           zone_state % xn(:) = ONE / nspec

           call eos(eos_input_re, zone_state)

           state(i,j,k,URHO)  = zone_state % rho
           state(i,j,k,UTEMP) = zone_state % T
           state(i,j,k,UEINT) = zone_state % e * zone_state % rho
           state(i,j,k,UFS:UFS+nspec-1) = zone_state % xn(:) * zone_state % rho

           state(i,j,k,UMX) = state(i,j,k,URHO) * smallu
           state(i,j,k,UMY) = state(i,j,k,URHO) * smallu
           state(i,j,k,UMZ) = state(i,j,k,URHO) * smallu

           state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                ( state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2 ) / &
                ( TWO * state(i,j,k,URHO) )

        enddo
     enddo
  enddo
  !$OMP END PARALLEL DO

end subroutine ca_initdata
