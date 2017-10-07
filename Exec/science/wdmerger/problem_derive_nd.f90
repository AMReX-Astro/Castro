! Derive momentum, given input vector of the grid momenta.

subroutine ca_derinertialmomentumx(p,p_lo,p_hi,ncomp_p, &
                                   u,u_lo,u_hi,ncomp_u, &
                                   lo,hi,domlo,domhi, &
                                   dx,xlo,time,dt,bc,level,grid_no) &
                                   bind(C,name='ca_derinertialmomentumx')

  use bl_constants_module, only: HALF
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: p_lo(3), p_hi(3), ncomp_p ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(1)

        enddo
     enddo
  enddo

end subroutine ca_derinertialmomentumx



! Derive momentum, given input vector of the grid momenta.

subroutine ca_derinertialmomentumy(p,p_lo,p_hi,ncomp_p, &
                                   u,u_lo,u_hi,ncomp_u, &
                                   lo,hi,domlo,domhi, &
                                   dx,xlo,time,dt,bc,level,grid_no) &
                                   bind(C,name='ca_derinertialmomentumy')

  use bl_constants_module, only: HALF
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: p_lo(3), p_hi(3), ncomp_p ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(2)

        enddo
     enddo
  enddo

end subroutine ca_derinertialmomentumy



! Derive momentum, given input vector of the grid momenta.

subroutine ca_derinertialmomentumz(p,p_lo,p_hi,ncomp_p, &
                                   u,u_lo,u_hi,ncomp_u, &
                                   lo,hi,domlo,domhi, &
                                   dx,xlo,time,dt,bc,level,grid_no) &
                                   bind(C,name='ca_derinertialmomentumz')

  use bl_constants_module, only: HALF
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: p_lo(3), p_hi(3), ncomp_p ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           mom = rho * inertial_velocity(loc, vel, time)
           p(i,j,k,1) = mom(3)

        enddo
     enddo
  enddo

end subroutine ca_derinertialmomentumz



! Derive angular momentum, given input vector of the grid momenta.

subroutine ca_derinertialangmomx(L,L_lo,L_hi,ncomp_L, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialangmomx')

  use bl_constants_module, only: HALF
  use math_module, only: cross_product
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), ang_mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           ang_mom = cross_product(loc, rho * inertial_velocity(loc, vel, time))

           L(i,j,k,1) = ang_mom(1)

        enddo
     enddo
  enddo

end subroutine ca_derinertialangmomx



subroutine ca_derinertialangmomy(L,L_lo,L_hi,ncomp_L, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialangmomy')

  use bl_constants_module, only: HALF
  use math_module, only: cross_product
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), ang_mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           ang_mom = cross_product(loc, rho * inertial_velocity(loc, vel, time))          

           L(i,j,k,1) = ang_mom(2)

        enddo
     enddo
  enddo

end subroutine ca_derinertialangmomy



subroutine ca_derinertialangmomz(L,L_lo,L_hi,ncomp_L, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialangmomz')

  use bl_constants_module, only: HALF
  use math_module, only: cross_product
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: L_lo(3), L_hi(3), ncomp_L ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: L(L_lo(1):L_hi(1),L_lo(2):L_hi(2),L_lo(3):L_hi(3),ncomp_L)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), vel(3), ang_mom(3), rho

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           rho = u(i,j,k,1)
           vel = u(i,j,k,2:4) / rho
           ang_mom = cross_product(loc, rho * inertial_velocity(loc, vel, time))

           L(i,j,k,1) = ang_mom(3)

        enddo
     enddo
  enddo

end subroutine ca_derinertialangmomz



! Derive radial momentum, given input vector of the grid momenta.

subroutine ca_derinertialradmomx(R,R_lo,R_hi,ncomp_R, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialradmomx')

  use bl_constants_module, only: HALF, ONE
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: R_lo(3), R_hi(3), ncomp_R ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: R(R_lo(1):R_hi(1),R_lo(2):R_hi(2),R_lo(3):R_hi(3),ncomp_R)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), mom(3), radInv

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           mom = inertial_velocity(loc, u(i,j,k,2:4), time)

           radInv = ONE / sqrt( loc(2)**2 + loc(3)**2 )

           R(i,j,k,1) = loc(2) * radInv * mom(2) + loc(3) * radInv * mom(3)

        enddo
     enddo
  enddo

end subroutine ca_derinertialradmomx



subroutine ca_derinertialradmomy(R,R_lo,R_hi,ncomp_R, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialradmomy')

  use bl_constants_module, only: HALF, ONE
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: R_lo(3), R_hi(3), ncomp_R ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: R(R_lo(1):R_hi(1),R_lo(2):R_hi(2),R_lo(3):R_hi(3),ncomp_R)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), mom(3), radInv

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           mom = inertial_velocity(loc, u(i,j,k,2:4), time)

           radInv = ONE / sqrt( loc(1)**2 + loc(3)**2 )

           R(i,j,k,1) = loc(1) * radInv * mom(1) + loc(3) * radInv * mom(3)

        enddo
     enddo
  enddo

end subroutine ca_derinertialradmomy



subroutine ca_derinertialradmomz(R,R_lo,R_hi,ncomp_R, &
                                 u,u_lo,u_hi,ncomp_u, &
                                 lo,hi,domlo,domhi, &
                                 dx,xlo,time,dt,bc,level,grid_no) &
                                 bind(C,name='ca_derinertialradmomz')

  use bl_constants_module, only: HALF, ONE
  use wdmerger_util_module, only: inertial_velocity
  use prob_params_module, only: center

  implicit none

  integer          :: R_lo(3), R_hi(3), ncomp_R ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 4
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: R(R_lo(1):R_hi(1),R_lo(2):R_hi(2),R_lo(3):R_hi(3),ncomp_R)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), mom(3), radInv

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3) - center(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2) - center(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1) - center(1)

           mom = inertial_velocity(loc, u(i,j,k,2:4), time)

           radInv = ONE / sqrt( loc(1)**2 + loc(2)**2 )

           R(i,j,k,1) = loc(1) * radInv * mom(1) + loc(2) * radInv * mom(2)

        enddo
     enddo
  enddo

end subroutine ca_derinertialradmomz



! Derive the effective potential phiEff = phiGrav + phiRot

subroutine ca_derphieff(phi,phi_lo,phi_hi,ncomp_phi, &
                        u,u_lo,u_hi,ncomp_u, &
                        lo,hi,domlo,domhi, &
                        dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derphieff')

  implicit none

  integer          :: phi_lo(3), phi_hi(3), ncomp_phi ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 2
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           phi(i,j,k,1) = u(i,j,k,1) + u(i,j,k,2)

        enddo
     enddo
  enddo

end subroutine ca_derphieff



! Derive an approximation to the effective potential of the primary only,
! by treating it as a point-mass at its center of mass.
! The u array contains the rotational potential, so we only need to calculate
! the gravitational potential from the point-mass.

subroutine ca_derphieffpm_p(phi,phi_lo,phi_hi,ncomp_phi, &
                            u,u_lo,u_hi,ncomp_u, &
                            lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derphieffpm_p')

  use bl_constants_module, only: ZERO, HALF
  use probdata_module, only: mass_P, com_P
  use fundamental_constants_module, only: Gconst

  implicit none

  integer          :: phi_lo(3), phi_hi(3), ncomp_phi ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 1
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), r

  phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  ! Don't do anything here if the star no longer exists

  if (mass_P == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           r = sqrt( sum( (loc - com_P)**2 ) )

           phi(i,j,k,1) = -Gconst * mass_P / r + u(i,j,k,1)

        enddo
     enddo
  enddo

end subroutine ca_derphieffpm_p



! Same as above, but for the secondary.

subroutine ca_derphieffpm_s(phi,phi_lo,phi_hi,ncomp_phi, &
                            u,u_lo,u_hi,ncomp_u, &
                            lo,hi,domlo,domhi, &
                            dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derphieffpm_s')

  use bl_constants_module, only: ZERO, HALF
  use probdata_module, only: mass_S, com_S
  use fundamental_constants_module, only: Gconst

  implicit none

  integer          :: phi_lo(3), phi_hi(3), ncomp_phi ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 1
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: phi(phi_lo(1):phi_hi(1),phi_lo(2):phi_hi(2),phi_lo(3):phi_hi(3),ncomp_phi)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), r

  phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = ZERO

  ! Don't do anything here if the star no longer exists

  if (mass_S == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

           r = sqrt( sum( (loc - com_S)**2 ) )

           phi(i,j,k,1) = -Gconst * mass_S / r + u(i,j,k,1)

        enddo
     enddo
  enddo

end subroutine ca_derphieffpm_s



subroutine ca_derrhophiGrav(rhophi,p_lo,p_hi,nk, &
                            dat,d_lo,d_hi,nc, &
                            lo,hi,domlo,domhi,delta, &
                            xlo,time,dt,bc,level,grid_no) &
                            bind(C, name="ca_derrhophiGrav")

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: p_lo(3), p_hi(3), nk
  integer          :: d_lo(3), d_hi(3), nc
  integer          :: domlo(3), domhi(3)
  integer          :: bc(3,2,nc)
  double precision :: delta(3), xlo(3), time, dt
  double precision :: rhophi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),nk)
  double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
  integer          :: level, grid_no

  integer          :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhophi(i,j,k,1) = dat(i,j,k,1) * dat(i,j,k,2)
        end do
     end do
  end do

end subroutine ca_derrhophiGrav



subroutine ca_derrhophiRot(rhophi,p_lo,p_hi,nk, &
                           dat,d_lo,d_hi,nc, &
                           lo,hi,domlo,domhi,delta, &
                           xlo,time,dt,bc,level,grid_no) &
                           bind(C, name="ca_derrhophiRot")

  implicit none

  integer          :: lo(3), hi(3)
  integer          :: p_lo(3), p_hi(3), nk
  integer          :: d_lo(3), d_hi(3), nc
  integer          :: domlo(3), domhi(3)
  integer          :: bc(3,2,nc)
  double precision :: delta(3), xlo(3), time, dt
  double precision :: rhophi(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),nk)
  double precision ::    dat(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
  integer          :: level, grid_no

  integer          :: i, j, k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           rhophi(i,j,k,1) = dat(i,j,k,1) * dat(i,j,k,2)
        end do
     end do
  end do

end subroutine ca_derrhophiRot



! Create a mask for all zones considered to be within the primary star.
! It uses the same prescription as above for the effective potential of the
! star, and uses the stellar density threshold input parameter to determine
! what parts of the domain should be considered stellar material.
! The convention will be that the mask is positive (1) for zones inside the
! star and negative (-1) for zones outside the star.

subroutine ca_derprimarymask(mask,mask_lo,mask_hi,ncomp_mask, &
                             u,u_lo,u_hi,ncomp_u, &
                             lo,hi,domlo,domhi, &
                             dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derprimarymask')

  use bl_constants_module, only: ZERO, HALF, ONE
  use probdata_module, only: mass_P, com_P, mass_S, com_S, stellar_density_threshold
  use fundamental_constants_module, only: Gconst

  implicit none

  integer          :: mask_lo(3), mask_hi(3), ncomp_mask ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 2 (density, rotational potential)
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3),ncomp_mask)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), r_P, r_S, phi_P, phi_S

  ! By default, assume we're not inside the star.

  mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = -ONE

  ! Don't do anything here if the star no longer exists

  if (mass_P == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

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

end subroutine ca_derprimarymask



! Same as above, but for the secondary.

subroutine ca_dersecondarymask(mask,mask_lo,mask_hi,ncomp_mask, &
                               u,u_lo,u_hi,ncomp_u, &
                               lo,hi,domlo,domhi, &
                               dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_dersecondarymask')

  use bl_constants_module, only: ZERO, HALF, ONE
  use probdata_module, only: mass_P, com_P, mass_S, com_S, stellar_density_threshold
  use fundamental_constants_module, only: Gconst

  implicit none

  integer          :: mask_lo(3), mask_hi(3), ncomp_mask ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == 2 (density, rotational potential)
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: mask(mask_lo(1):mask_hi(1),mask_lo(2):mask_hi(2),mask_lo(3):mask_hi(3),ncomp_mask)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k
  double precision :: loc(3), r_P, r_S, phi_P, phi_S

  ! By default, assume we're not inside the star.

  mask(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = -ONE

  ! Don't do anything here if the star no longer exists

  if (mass_S == ZERO) return

  do k = lo(3), hi(3)
     loc(3) = xlo(3) + (dble(k - lo(3)) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = xlo(2) + (dble(j - lo(2)) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = xlo(1) + (dble(i - lo(1)) + HALF) * dx(1)

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

end subroutine ca_dersecondarymask



subroutine ca_derignitionradius(ignition_radius,ir_lo,ir_hi,ncomp_ir, &
                                u,u_lo,u_hi,ncomp_u, &
                                lo,hi,domlo,domhi, &
                                dx,xlo,time,dt,bc,level,grid_no) bind(C,name='ca_derignitionradius')

  use bl_constants_module, only: ZERO, HALF, ONE
  use meth_params_module, only: URHO, UTEMP, UFS
  use network, only: network_species_index
  use prob_params_module, only: dg

  implicit none

  integer          :: ir_lo(3), ir_hi(3), ncomp_ir ! == 1
  integer          :: u_lo(3), u_hi(3), ncomp_u ! == NVAR
  integer          :: lo(3), hi(3), domlo(3), domhi(3)
  double precision :: ignition_radius(ir_lo(1):ir_hi(1),ir_lo(2):ir_hi(2),ir_lo(3):ir_hi(3),ncomp_ir)
  double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  double precision :: dx(3), xlo(3), time, dt
  integer          :: bc(3,2,ncomp_u), level, grid_no

  integer          :: i, j, k

  double precision, parameter :: beta = 0.018d0
  double precision, parameter :: T_0 = 1.5d9

  integer          :: iC12
  double precision :: alpha, theta, deltaT, T_m, T_b, T_a
  double precision :: X_C, rho, T

  ! The following critical T gradient calculation comes from
  ! Garg and Chang (2017). It is appropriate for C+O material
  ! with X(C) ~ 0.5.

  iC12 = network_species_index("carbon-12")

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           T = u(i,j,k,UTEMP)
           rho = u(i,j,k,URHO)
           X_C = u(i,j,k,UFS+ic12-1) / u(i,j,k,URHO)

           ! Search all adjacent zones for the largest temperature gradient.

           T_m    = T
           deltaT = ZERO

           T_a = u(i-1*dg(1),j,k,UTEMP)

           if (abs(T_a - T) > deltaT) then
              T_b = min(T_a, T)
              T_m = max(T_a, T)
              deltaT = T_m - T_b
           end if

           T_a = u(i+1*dg(1),j,k,UTEMP)

           if (abs(T_a - T) > deltaT) then
              T_b = min(T_a, T)
              T_m = max(T_a, T)
              deltaT = T_m - T_b
           end if

           T_a = u(i,j-1*dg(2),k,UTEMP)

           if (abs(T_a - T) > deltaT) then
              T_b = min(T_a, T)
              T_m = max(T_a, T)
              deltaT = T_m - T_b
           end if

           T_a = u(i,j+1*dg(2),k,UTEMP)

           if (abs(T_a - T) > deltaT) then
              T_b = min(T_a, T)
              T_m = max(T_a, T)
              deltaT = T_m - T_b
           end if

           T_a = u(i,j,k-1*dg(3),UTEMP)

           if (abs(T_a - T) > deltaT) then
              T_b = min(T_a, T)
              T_m = min(T_a, T)
              deltaT = T_m - T_b
           end if

           T_a = u(i,j,k+1*dg(3),UTEMP)

           if (abs(T_a - T) > deltaT) then
              T_b = min(T_a, T)
              T_m = max(T_a, T)
              deltaT = T_m - T_b
           end if

           ! If we're below or above the temperature range for which
           ! we have a fit, skip this zone by setting it to zero.

           if (T_m < 1.6d9 .or. T_m > 2.8d9) then

              ignition_radius(i,j,k,1) = ZERO
              cycle

           endif

           ! Estimate alpha and theta by linearly interpolating
           ! between the values given in Table 2.

           if (T_m >= 1.6d9 .and. T_m < 1.8d9) then
              theta = 8.9d-3 + (8.1d-4 - 8.9d-3) / (1.8d9 - 1.6d9) * (T_m - 1.6d9)
              alpha = 22.0d0 + (22.0d0 - 21.0d0) / (1.8d9 - 1.6d9) * (T_m - 1.6d9)
           else if (T_m >= 1.8d9 .and. T_m < 2.0d9) then
              theta = 8.1d-4 + (1.0d-4 - 8.1d-4) / (2.0d9 - 1.8d9) * (T_m - 1.8d9)
              alpha = 21.0d0 + (20.3d0 - 21.0d0) / (2.0d9 - 1.8d9) * (T_m - 1.8d9)
           else if (T_m >= 2.0d9 .and. T_m < 2.4d9) then
              theta = 1.0d-4 + (3.4d-6 - 1.0d-4) / (2.4d9 - 2.0d9) * (T_m - 2.0d9)
              alpha = 20.3d0 + (18.8d0 - 20.3d0) / (2.4d9 - 2.0d9) * (T_m - 2.0d9)
           else if (T_m >= 2.4d9 .and. T_m <= 2.8d9) then
              theta = 3.4d-6 + (2.7d-7 - 3.4d-6) / (2.8d9 - 2.4d9) * (T_m - 2.4d9)
              alpha = 18.8d0 + (17.7d0 - 18.8d0) / (2.8d9 - 2.4d9) * (T_m - 2.4d9)
           end if

           ignition_radius(i,j,k,1) = 4.6d5 * (beta / 0.018d0) * ((alpha - ONE) / 20.0d0) * &
                                              (theta / 1.0d-3) * (deltaT / T_m) * (X_C / HALF)**(-2) * &
                                              (1.d7 / rho) * (HALF * (ONE + T_0 / T_m))**(-alpha)

        enddo
     enddo
  enddo

end subroutine ca_derignitionradius
