! Compute the analytic solution for the multigroup radiating sphere
! problem.  We follow the problem outline from Swesty and Myra (2009).
!
! The output from this is the radiation energy density as a function
! of energy group at a specified distance (r_obs) from the radiating
! sphere, at a time of t_obs.

module problem_derive_module

  use amrex_fort_module, only: rt => amrex_real
  use amrex_constants_module, only: M_PI
  use fundamental_constants_module, only: hplanck, k_B, c_light, ev2erg

  implicit none

  ! problem parameters

  ! we will specify the opacity as kappa = kapp_0 (nu_0 / nu)**3
  ! where nu_0 is a reference frequency
  ! 
  ! Swesty and Myra (2009) only say that the "opacity of the material
  ! is proportional to 1/E**3", but they don't give the proportionality
  ! constant.  The output depends critically on the constant.  Eric
  ! says that they used kappa = 1.e13 * (1.5e-6 MeV / E)**3
  !
  ! Converting 1.5e-6 MeV to a frequency gives nu_0 = 3.63e14 Hz
  real(rt), parameter :: nu_0 = 3.6e14  ! ref. freq (Hz)
  real(rt), parameter :: kappa_0 = 1.e13  ! scattering opacity (1/cm)

  ! geometry parameters
  real(rt), parameter :: R_sphere = 0.02_rt ! sphere radius (cm)
  real(rt), parameter :: r_obs = 0.06_rt  ! observer location (cm)

  ! physical parameters
  real(rt), parameter :: T_0 = 50.0_rt * ev2erg / k_B ! ambient temp (K)
  real(rt), parameter :: T_sphere = 1500.0_rt * ev2erg / k_B ! sphere temp (K)

contains

  ! to guard against over/underflows, we define some 'safe' functions
  ! (they come from Eric)

  function sferfc(arg1)

    implicit none

    real(rt), intent(in) :: arg1

    real(rt) :: sferfc

    sferfc = erfc(max(-26.0_rt, min(26.0_rt, arg1)))

  end function sferfc



  function sfexp(arg1)

    implicit none

    real(rt), intent(in) :: arg1

    real(rt) :: sfexp

    sfexp = exp(max(-650.0_rt, min(650.0_rt, arg1)))

  end function sfexp



  function planck(nu, T) result (B)

    ! the Planck function for a Blackbody (actually, energy density,
    ! B = (4 pi / c) I, where I is the normal Planck function
    ! 
    ! nu = frequency (Hz)
    ! T  = temperature (K)

    ! see Swesty and Myra (2009), eq. 23, but note that we are working
    ! in Hz here, so we have units of erg / cm^3 / Hz, where they
    ! have units of erg / cm^3 / MeV.  As a result, we have one
    ! less factor of hplanck.

    implicit none

    real(rt), intent(in) :: nu, T
    real(rt) :: B

    B = (8.0_rt * M_PI * hplanck * nu**3 / c_light**3) / (exp(hplanck * nu / (k_B * T)) - 1.0_rt)

    B = min(max(B, 1.e-50_rt), 1.e200_rt)

  end function planck



  function F_radsphere(r, t, nu) result (F)

    ! the function F(r,t) as defined in Swesty and Myra
    !
    ! r      = position
    ! t      = time (s)
    ! nu     = frequency (Hz)

    implicit none

    real(rt), intent(in) :: r, t, nu

    real(rt) :: kappa, erfc_term1, erfc_term2, F

    kappa = kappa_0 * (nu_0 / nu)**3

    erfc_term1 = sferfc(sqrt(3.0_rt * kappa / (4.0_rt * c_light * max(t, 1.e-50_rt))) * (r - R_sphere) - &
                        sqrt(c_light * t * kappa))

    erfc_term2 = sferfc(sqrt(3.0_rt * kappa / (4.0_rt * c_light * max(t, 1.e-50_rt))) * (r - R_sphere) + &
                        sqrt(c_light * t * kappa))

    F = 0.5_rt * (sfexp(-sqrt(3.0_rt) * kappa * (r - R_sphere)) * erfc_term1 + &
                  sfexp( sqrt(3.0_rt) * kappa * (r - R_sphere)) * erfc_term2)

    F = min(max(F, 1.e-50_rt), 1.e200_rt)

  end function F_radsphere




  ! Compute the analytic solution to the radiating sphere, as given
  ! by Swesty and Myra (2009), Eq. 76.

  subroutine deranalytic(p, p_lo, p_hi, ncomp_p, &
                         u, u_lo, u_hi, ncomp_u, &
                         lo, hi, domlo, domhi, &
                         dx, time) &
                         bind(C, name='deranalytic')

    use amrex_fort_module, only: rt => amrex_real
    use amrex_constants_module, only: HALF, M_PI
    use fundamental_constants_module, only: sigma_SB
    use prob_params_module, only: problo, center
    use rad_params_module, only: ngroups, nugroup, dnugroup

    implicit none

    integer,  intent(in   ) :: p_lo(3), p_hi(3)
    integer,  intent(in   ) :: u_lo(3), u_hi(3)
    integer,  intent(in   ) :: lo(3), hi(3), domlo(3), domhi(3)
    real(rt), intent(inout) :: p(p_lo(1):p_hi(1),p_lo(2):p_hi(2),p_lo(3):p_hi(3),0:ncomp_p-1)
    real(rt), intent(in   ) :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),0:ncomp_u-1)
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: ncomp_p ! == nGroups
    integer,  intent(in   ), value :: ncomp_u ! == nGroups
    real(rt), intent(in   ), value :: time

    integer  :: i, j, k
    real(rt) :: loc(3), r

    real(rt) :: E, F
    integer :: n

    do k = lo(3), hi(3)
       loc(3) = problo(3) + (dble(k) + HALF) * dx(3) - center(3)
       do j = lo(2), hi(2)
          loc(2) = problo(2) + (dble(j) + HALF) * dx(2) - center(2)
          do i = lo(1), hi(1)
             loc(1) = problo(1) + (dble(i) + HALF) * dx(1) - center(1)

             r = sqrt(loc(1)**2 + loc(2)**2 + loc(3)**2)

             ! compute the analytic radiating sphere solution for each point
             do n = 0, ngroups-1
                F = F_radsphere(r, time, nugroup(n))
                E = planck(nugroup(n), T_0) + &
                     (R_sphere / r) * (planck(nugroup(n), T_sphere) - &
                     planck(nugroup(n), T_0)) * F
                p(i,j,k,n) = E * dnugroup(n)
             end do

          end do
       end do
    end do

  end subroutine deranalytic

end module problem_derive_module
