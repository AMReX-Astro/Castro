! This sets up the Gresho vortex problem as described in 
! Miczek, Roeple, and Edelmann 2015
!
! By choosing the reference pressure, p0, we can specify the
! Mach number

subroutine amrex_probinit(init, name, namlen, problo, probhi) bind(c)

  use probdata_module, only: p0, rho0, t_r, nsub, x_r, q_r
  use prob_params_module, only: center
  use amrex_constants_module, only: M_pi
  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: init, namlen
  integer,  intent(in   ) :: name(namlen)
  real(rt), intent(in   ) :: problo(3), probhi(3)

  call probdata_init(name, namlen)

  ! problem center
  center(:) = (problo(:) + probhi(:)) / 2.e0_rt

  ! characteristic scales
  x_r = probhi(1) - problo(1)
  q_r = 0.4_rt*M_pi*x_r/t_r

end subroutine amrex_probinit

subroutine ca_initdata(level, time, lo, hi, nscal, &
                       state, s_lo, s_hi, &
                       dx, xlo, xhi)

  use probdata_module, only: p0, rho0, t_r, nsub, x_r, q_r
  use actual_eos_module, only: gamma_const
  use amrex_constants_module, only: ZERO, HALF, ONE, TWO, FOUR, FIVE
  use meth_params_module , only: NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UFS
  use prob_params_module, only: center, problo, dim
  use network, only: nspec
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer,  intent(in   ) :: level, nscal
  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: s_lo(3), s_hi(3)
  real(rt), intent(in   ) :: xlo(3), xhi(3), time, dx(3)
  real(rt), intent(inout) :: state(s_lo(1):s_hi(1), s_lo(2):s_hi(2), s_lo(3):s_hi(3), NVAR)

  real(rt) :: x, y, z, xl, yl, zl, xx, yy, zz
  real(rt) :: r
  real(rt) :: reint, p, u_phi, u_tot

  integer  :: i, j, k, ii, jj, kk

  do k = lo(3), hi(3)
     zl = problo(3) + dx(3) * dble(k)
     z = problo(3) + dx(3) * dble(k + HALF)

     do j = lo(2), hi(2)
        yl = problo(2) + dx(2) * dble(j)
        y = problo(2) + dx(2) * dble(j + HALF)

        do i = lo(1), hi(1)
           xl = problo(1) + dx(1) * dble(i)
           x = problo(1) + dx(1) * dble(i + HALF)

           reint = ZERO
           u_tot = ZERO

           do kk = 0, nsub - 1
              zz = zl + dx(3) * dble(kk + HALF) / nsub

              do jj = 0, nsub - 1
                 yy = yl + dx(2) * dble(jj + HALF) / nsub

                 do ii = 0, nsub - 1
                    xx = xl + dx(1) * dble(ii + HALF) / nsub

                    r = sqrt((xx - center(1))**2 + (yy - center(2))**2)

                    if (r < 0.2_rt) then
                       u_phi = FIVE*r
                       p = p0 + 12.5_rt*r**2

                    else if (r < 0.4_rt) then
                       u_phi = TWO - FIVE*r
                       p = p0 + 12.5_rt*r**2 + FOUR*(ONE - FIVE*r - log(0.2_rt) + log(r))

                    else
                       u_phi = ZERO
                       p = p0 - TWO + FOUR*log(TWO)
                    end if

                    u_tot = u_tot + u_phi
                    reint = reint + p/(gamma_const - ONE)

                 end do
              end do
           end do

           u_phi = u_tot / (nsub**3)
           reint = reint / (nsub**3)

           state(i,j,k,URHO) = rho0

           ! phi unit vector: \hat{\phi} = -sin(phi) \hat{x} + cos(phi) \hat{y}
           ! with cos(phi) = x/r; sin(phi) = y/r
           r = sqrt((x - center(1))**2 + (y - center(2))**2)
           state(i,j,k,UMX) = -rho0 * q_r * u_phi * ((y - center(2)) / r) ! -sin(phi) = y/r
           state(i,j,k,UMY) =  rho0 * q_r * u_phi * ((x - center(1)) / r) !  cos(phi) = x/r
           state(i,j,k,UMZ) = ZERO

           state(i,j,k,UEDEN) = reint + &
                                0.5e0_rt*(state(i,j,k,UMX)**2 / state(i,j,k,URHO) + &
                                          state(i,j,k,UMY)**2 / state(i,j,k,URHO) + &
                                          state(i,j,k,UMZ)**2 / state(i,j,k,URHO))

           state(i,j,k,UEINT) = reint

           state(i,j,k,UFS:UFS-1+nspec) = ZERO
           state(i,j,k,UFS) = state(i,j,k,URHO)

        end do
     end do
  end do


end subroutine ca_initdata
