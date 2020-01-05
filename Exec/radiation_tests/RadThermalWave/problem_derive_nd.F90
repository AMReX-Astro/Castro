! Analytic solution

subroutine dertexact(p, p_lo, p_hi, ncomp_p, &
                     u, u_lo, u_hi, ncomp_u, &
                     lo, hi, domlo, domhi, &
                     dx, time) &
                     bind(C, name='dertexact')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF, M_PI
  use fundamental_constants_module, only: sigma_SB
  use extern_probin_module, only: const_kappa_r, kappa_r_exp_n
  use prob_params_module, only: problo
  use probdata_module, only: Eexp, rhocv

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_p ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 1
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), r2

  real(rt) :: a, pe, pfac, xi0, xf, Tbar, Tc, Q

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

           Q = Eexp / rhocv

           a = (16.0_rt * sigma_SB) / (3.0_rt * const_kappa_r) / rhocv
           pe = kappa_r_exp_n + 3.0_rt

           pfac = exp(log_gamma(2.5_rt + 1.0_rt / pe) - log_gamma(1.0_rt + 1.0_rt / pe) - log_gamma(1.5_rt))

           xi0 = ((3.0_rt * pe + 2.0_rt) / &
                  (2.0_rt**(pe - 1.0_rt) * pe * M_PI**pe) &
                 )**(1.0_rt / (3.0_rt * pe + 2.0_rt)) * &
                 pfac**(pe / (3.0_rt * pe + 2.0_rt))

           xf = xi0 * (a * Q**pe * max(time, 1.d-50))**(1.0_rt / (3.0_rt * pe + 2.0_rt))

           Tbar = Q / xf**3

           Tc = xi0**3 * (pe * xi0 * xi0 / (6.0_rt * pe + 4.0_rt))**(1.0_rt / pe) * Tbar

           r2 = loc(1)**2 + loc(2)**2 + loc(3)**2
           p(i,j,k,1) = Tc * max((1.e0_rt - r2 / xf**2), 0.0_rt)**(1.0_rt / pe)

        end do
     end do
  end do

end subroutine dertexact



! Error compared to analytic solution

subroutine derterror(p, p_lo, p_hi, ncomp_p, &
                     u, u_lo, u_hi, ncomp_u, &
                     lo, hi, domlo, domhi, &
                     dx, time) &
                     bind(C, name='derterror')

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: HALF, M_PI
  use fundamental_constants_module, only: sigma_SB
  use extern_probin_module, only: const_kappa_r, kappa_r_exp_n
  use prob_params_module, only: problo
  use probdata_module, only: Eexp, rhocv

  implicit none

  integer,  intent(in   ) :: p_lo(3), p_hi(3)
  integer,  intent(in   ) :: u_lo(3), u_hi(3)
  integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
  real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),ncomp_p)
  real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),ncomp_u)
  real(rt), intent(in   ) :: dx(3)
  integer,  intent(in   ), value :: ncomp_p ! == 1
  integer,  intent(in   ), value :: ncomp_u ! == 1
  real(rt), intent(in   ), value :: time

  integer  :: i, j, k
  real(rt) :: loc(3), r2

  real(rt) :: a, pe, pfac, xi0, xf, Tbar, Tc, Q

  !$gpu

  do k = lo(3), hi(3)
     loc(3) = problo(3) + (dble(k) + HALF) * dx(3)
     do j = lo(2), hi(2)
        loc(2) = problo(2) + (dble(j) + HALF) * dx(2)
        do i = lo(1), hi(1)
           loc(1) = problo(1) + (dble(i) + HALF) * dx(1)

           Q = Eexp / rhocv

           a = (16.0_rt * sigma_SB) / (3.0_rt * const_kappa_r) / rhocv
           pe = kappa_r_exp_n + 3.0_rt

           pfac = exp(log_gamma(2.5_rt + 1.0_rt / pe) - log_gamma(1.0_rt + 1.0_rt / pe) - log_gamma(1.5_rt))

           xi0 = ((3.0_rt * pe + 2.0_rt) / &
                  (2.0_rt**(pe - 1.0_rt) * pe * M_PI**pe) &
                 )**(1.0_rt / (3.0_rt * pe + 2.0_rt)) * &
                 pfac**(pe / (3.0_rt * pe + 2.0_rt))

           xf = xi0 * (a * Q**pe * max(time, 1.d-50))**(1.0_rt / (3.0_rt * pe + 2.0_rt))

           Tbar = Q / xf**3

           Tc = xi0**3 * (pe * xi0 * xi0 / (6.0_rt * pe + 4.0_rt))**(1.0_rt / pe) * Tbar

           r2 = loc(1)**2 + loc(2)**2 + loc(3)**2
           p(i,j,k,1) = u(i,j,k,1) - Tc * max((1.e0_rt - r2 / xf**2), 0.0_rt)**(1.0_rt / pe)

        end do
     end do
  end do

end subroutine derterror
