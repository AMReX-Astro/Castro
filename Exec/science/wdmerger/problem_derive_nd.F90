! Derive momentum, given the grid momenta.

subroutine derinertialmomentumx(p, p_lo, p_hi, ncomp_p, &
                                u, u_lo, u_hi, ncomp_u, &
                                lo, hi, domlo, domhi, &
                                dx, time) &
                                bind(C, name='derinertialmomentumx')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_p ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), vel(3), mom(3), rho

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(1)

        enddo
     enddo
  enddo

end subroutine derinertialmomentumx



! Derive momentum, given the grid momenta.

subroutine derinertialmomentumy(p, p_lo, p_hi, ncomp_p, &
                                u, u_lo, u_hi, ncomp_u, &
                                lo, hi, domlo, domhi, &
                                dx, time) &
                                bind(C, name='derinertialmomentumy')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_p ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), vel(3), mom(3), rho

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(2)

        enddo
     enddo
  enddo

end subroutine derinertialmomentumy



! Derive momentum, given the grid momenta.

subroutine derinertialmomentumz(p, p_lo, p_hi, ncomp_p, &
                                u, u_lo, u_hi, ncomp_u, &
                                lo, hi, domlo, domhi, &
                                dx, time) &
                                bind(C, name='derinertialmomentumz')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_p ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), vel(3), mom(3), rho

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(3)

        enddo
     enddo
  enddo

end subroutine derinertialmomentumz



! Derive angular momentum, given the grid momenta.

subroutine derinertialangmomx(L, L_lo, L_hi, ncomp_L, &
                              u, u_lo, u_hi, ncomp_u, &
                              lo, hi, domlo, domhi, &
                              dx, time) &
                              bind(C, name='derinertialangmomx')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use math_module, only: cross_product ! function
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: L_lo(3), L_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_L ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), vel(3), mom(3), ang_mom(3), rho

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)

           ang_mom = cross_product(loc, mom)

           L(i,j,k,1) = ang_mom(1)

        enddo
     enddo
  enddo

end subroutine derinertialangmomx



subroutine derinertialangmomy(L, L_lo, L_hi, ncomp_L, &
                              u, u_lo, u_hi, ncomp_u, &
                              lo, hi, domlo, domhi, &
                              dx, time) &
                              bind(C, name='derinertialangmomy')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use math_module, only: cross_product ! function
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: L_lo(3), L_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_L ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), vel(3), mom(3), ang_mom(3), rho

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           ang_mom = cross_product(loc, mom)

           L(i,j,k,1) = ang_mom(2)

        enddo
     enddo
  enddo

end subroutine derinertialangmomy



subroutine derinertialangmomz(L, L_lo, L_hi, ncomp_L, &
                              u, u_lo, u_hi, ncomp_u, &
                              lo, hi, domlo, domhi, &
                              dx, time) &
                              bind(C, name='derinertialangmomz')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF
  use math_module, only: cross_product ! function
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: L_lo(3), L_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_L ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), vel(3), mom(3), ang_mom(3), rho

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           ang_mom = cross_product(loc, mom)

           L(i,j,k,1) = ang_mom(3)

        enddo
     enddo
  enddo

end subroutine derinertialangmomz



! Derive radial momentum, given the grid momenta.

subroutine derinertialradmomx(R, R_lo, R_hi, ncomp_R, &
                              u, u_lo, u_hi, ncomp_u, &
                              lo, hi, domlo, domhi, &
                              dx, time) &
                              bind(C, name='derinertialradmomx')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF, ONE
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: R_lo(3), R_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: R(R_lo(1):R_hi(1),R_lo(2):R_hi(2),R_lo(3):R_hi(3),ncomp_R)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_R ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), mom(3), radInv

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           mom = u(i,j,k,2:4)
           mom = inertial_velocity(loc, mom, time)

           radInv = ONE / sqrt( loc(2)**2 + loc(3)**2 )

           R(i,j,k,1) = loc(2) * radInv * mom(2) + loc(3) * radInv * mom(3)

        enddo
     enddo
  enddo

end subroutine derinertialradmomx



subroutine derinertialradmomy(R, R_lo, R_hi, ncomp_R, &
                              u, u_lo, u_hi, ncomp_u, &
                              lo, hi, domlo, domhi, &
                              dx, time) &
                              bind(C, name='derinertialradmomy')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF, ONE
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: R_lo(3), R_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: R(R_lo(1):R_hi(1),R_lo(2):R_hi(2),R_lo(3):R_hi(3),ncomp_R)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_R ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), mom(3), radInv

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           mom = u(i,j,k,2:4)
           mom = inertial_velocity(loc, mom, time)

           radInv = ONE / sqrt( loc(1)**2 + loc(3)**2 )

           R(i,j,k,1) = loc(1) * radInv * mom(1) + loc(3) * radInv * mom(3)

        enddo
     enddo
  enddo

end subroutine derinertialradmomy



subroutine derinertialradmomz(R, R_lo, R_hi, ncomp_R, &
                              u, u_lo, u_hi, ncomp_u, &
                              lo, hi, domlo, domhi, &
                              dx, time) &
                              bind(C, name='derinertialradmomz')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF, ONE
  use wdmerger_util_module, only: inertial_velocity ! function
  use prob_params_module, only: problo, center

  implicit none

  integer,  intent(in   ) :: R_lo(3), R_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: R(R_lo(1):R_hi(1),R_lo(2):R_hi(2),R_lo(3):R_hi(3),ncomp_R)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_R ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 4
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), mom(3), radInv

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

           mom = u(i,j,k,2:4)
           mom = inertial_velocity(loc, mom, time)

           radInv = ONE / sqrt( loc(1)**2 + loc(2)**2 )

           R(i,j,k,1) = loc(1) * radInv * mom(1) + loc(2) * radInv * mom(2)

        enddo
     enddo
  enddo

end subroutine derinertialradmomz



! Derive the effective potential phiEff = phiGrav + phiRot

subroutine derphieff(phi, phi_lo, phi_hi, ncomp_phi, &
                     u, u_lo, u_hi, ncomp_u, &
                     lo, hi, domlo, domhi, &
                     dx, time) bind(C, name='derphieff')

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_phi ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 2
  real(rt), intent(in   ), value :: time

  integer :: i, j, k

  !$gpu

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           phi(i,j,k,1) = u(i,j,k,1) + u(i,j,k,2)

        enddo
     enddo
  enddo

end subroutine derphieff



! Derive an approximation to the effective potential of the primary only,
! by treating it as a point-mass at its center of mass.
! The u array contains the rotational potential, so we only need to calculate
! the gravitational potential from the point-mass.

subroutine derphieffpm_p(phi, phi_lo, phi_hi, ncomp_phi, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo, hi, domlo, domhi, &
                         dx, time) bind(C,name='derphieffpm_p')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF
  use wdmerger_util_module, only: mass_P, com_P
  use fundamental_constants_module, only: Gconst
  use prob_params_module, only: problo

  implicit none

  integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_phi ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 1
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), r

  !$gpu

  phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  ! Don't do anything here if the star no longer exists

  if (mass_P == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

           r = sqrt( sum( (loc - com_P)**2 ) )

           phi(i,j,k,1) = -Gconst * mass_P / r + u(i,j,k,1)

        enddo
     enddo
  enddo

end subroutine derphieffpm_p



! Same as above, but for the secondary.

subroutine derphieffpm_s(phi, phi_lo, phi_hi, ncomp_phi, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo, hi, domlo, domhi, &
                         dx, time) bind(C,name='derphieffpm_s')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF
  use wdmerger_util_module, only: mass_S, com_S
  use fundamental_constants_module, only: Gconst
  use prob_params_module, only: problo

  implicit none

  integer,  intent(in   ) :: phi_lo(3), phi_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_phi ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 1
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), r

  !$gpu

  phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  ! Don't do anything here if the star no longer exists

  if (mass_S == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

           r = sqrt( sum( (loc - com_S)**2 ) )

           phi(i,j,k,1) = -Gconst * mass_S / r + u(i,j,k,1)

        enddo
     enddo
  enddo

end subroutine derphieffpm_s



subroutine derrhophiGrav(rhophi, p_lo, p_hi, nk, &
                         dat, d_lo, d_hi, nc, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name="derrhophiGrav")

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: d_lo(3), d_hi(3)
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(inout) :: rhophi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),nk)
  real(rt), intent(in   ) ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: nk, nc
  real(rt), intent(in   ), value :: time

  integer :: i, j, k

  !$gpu

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhophi(i,j,k,1) = dat(i,j,k,1) * dat(i,j,k,2)
        end do
     end do
  end do

end subroutine derrhophiGrav



subroutine derrhophiRot(rhophi, p_lo, p_hi, nk, &
                        dat, d_lo, d_hi, nc, &
                        lo, hi, domlo, domhi, &
                        dx, time) &
                        bind(C, name="derrhophiRot")

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer,  intent(in   ) :: lo(3), hi(3)
  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: d_lo(3), d_hi(3)
  integer,  intent(in   ) :: domlo(3), domhi(3)
  real(rt), intent(inout) :: rhophi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),nk)
  real(rt), intent(in   ) ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: nk, nc
  real(rt), intent(in   ), value :: time

  integer :: i, j, k

  !$gpu

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhophi(i,j,k,1) = dat(i,j,k,1) * dat(i,j,k,2)
        end do
     end do
  end do

end subroutine derrhophiRot



! Create a mask for all zones considered to be within the primary star.
! It uses the same prescription as above for the effective potential of the
! star, and uses the stellar density threshold input parameter to determine
! what parts of the domain should be considered stellar material.
! The convention will be that the mask is positive (1) for zones inside the
! star and negative (-1) for zones outside the star.

subroutine derprimarymask(mask, mask_lo, mask_hi, ncomp_mask, &
                          u, u_lo, u_hi, ncomp_u, &
                          lo, hi, domlo, domhi, &
                          dx, time) bind(C, name='derprimarymask')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF, ONE
  use wdmerger_util_module, only: mass_P, com_P, mass_S, com_S, stellar_density_threshold
  use fundamental_constants_module, only: Gconst
  use prob_params_module, only: problo

  implicit none

  integer,  intent(in   ) :: mask_lo(3), mask_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3),ncomp_mask)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_mask ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 2 (density, rotational potential)
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), r_P, r_S, phi_P, phi_S

  !$gpu

  ! By default, assume we're not inside the star.

  mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = -ONE

  ! Don't do anything here if the star no longer exists

  if (mass_P == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

           ! Ignore zones whose density is too low.

           if (u(i,j,k,1) < stellar_density_threshold) cycle

           r_P = sqrt( sum( (loc - com_P)**2 ) )
           r_S = sqrt( sum( (loc - com_S)**2 ) )

           phi_p = -Gconst * mass_P / r_P + u(i,j,k,2)
           phi_s = -Gconst * mass_S / r_S + u(i,j,k,2)

           if (phi_p < ZERO .and. phi_p < phi_s) then
              mask(i,j,k,1) = ONE
           endif

        enddo
     enddo
  enddo

end subroutine derprimarymask



! Same as above, but for the secondary.

subroutine dersecondarymask(mask, mask_lo, mask_hi, ncomp_mask, &
                            u, u_lo, u_hi, ncomp_u, &
                            lo, hi, domlo, domhi, &
                            dx, time) bind(C, name='dersecondarymask')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: ZERO, HALF, ONE
  use wdmerger_util_module, only: mass_P, com_P, mass_S, com_S, stellar_density_threshold
  use fundamental_constants_module, only: Gconst
  use prob_params_module, only: problo

  implicit none

  integer,  intent(in   ) :: mask_lo(3), mask_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3),ncomp_mask)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_mask ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 2 (density, rotational potential)
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), r_P, r_S, phi_P, phi_S

  !$gpu

  ! By default, assume we're not inside the star.

  mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = -ONE

  ! Don't do anything here if the star no longer exists

  if (mass_S == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

           ! Ignore zones whose density is too low.

           if (u(i,j,k,1) < stellar_density_threshold) cycle

           r_P = sqrt( sum( (loc - com_P)**2 ) )
           r_S = sqrt( sum( (loc - com_S)**2 ) )

           phi_p = -Gconst * mass_P / r_P + u(i,j,k,2)
           phi_s = -Gconst * mass_S / r_S + u(i,j,k,2)

           if (phi_s < ZERO .and. phi_s < phi_p) then
              mask(i,j,k,1) = ONE
           endif

        enddo
     enddo
  enddo

end subroutine dersecondarymask
