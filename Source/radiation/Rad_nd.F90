
module rad_nd_module

  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), parameter, private :: tiny = 1.e-50_rt
  real(rt), parameter, private :: BIGKR = 1.e25_rt

contains

  subroutine gcv(lo, hi, &
                 cv, c_lo, c_hi, &
                 temp, t_lo, t_hi, &
                 const, em, en, tf, &
                 state, s_lo, s_hi) bind(C, name="gcv")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains temp on input
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: const, em, en, tf

    real(rt) :: alpha, teff, frhoal
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (em == 0.e0_rt) then
                alpha = const
             else
                alpha = const * state(i,j,k,URHO)**em
             end if

             frhoal = state(i,j,k,URHO) * alpha + tiny

             if (en == 0.e0_rt) then
                cv(i,j,k) = alpha
             else
                teff = max(temp(i,j,k), tiny)
                teff = teff + tf * exp(-teff / (tf + tiny))
                cv(i,j,k) = alpha * teff**(-en)
             end if

          end do
       end do
    end do

  end subroutine gcv



  subroutine ca_compute_c_v(lo, hi, &
                            cv, c_lo, c_hi, &
                            temp, t_lo, t_hi, &
                            state, s_lo, s_hi) &
                            bind(C, name="ca_compute_c_v")

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: c_lo(3), c_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: cv(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer     :: i, j, k
    real(rt)    :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             call eos(eos_input_rt, eos_state)

             cv(i,j,k) = eos_state % cv

          end do
       end do
    end do

  end subroutine ca_compute_c_v



  subroutine ca_get_rhoe(lo, hi, &
                         rhoe, r_lo, r_hi, &
                         temp, t_lo, t_hi, &
                         state, s_lo, s_hi) &
                         bind(C, name="ca_get_rhoe")

    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use network, only: nspec, naux
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: rhoe(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer     :: i, j, k
    real(rt)    :: rhoInv
    type(eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   =  temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             eos_state % aux = state(i,j,k,UFX:UFX+naux -1) * rhoInv

             call eos(eos_input_rt, eos_state)

             rhoe(i,j,k) = eos_state % rho * eos_state % e

          end do
       end do
    end do

  end subroutine ca_get_rhoe



  subroutine gtemp(lo, hi, &
                   temp, t_lo, t_hi, &
                   const, em, en, &
                   state, s_lo, s_hi) bind(C, name="gtemp")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))  ! temp contains frhoe on input
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value  :: const, em, en

    real(rt) :: alpha, teff, ex, frhoal
    integer  :: i, j, k

    !$gpu

#ifndef AMREX_USE_GPU
    if (en >= 1.e0_rt) then
       call castro_error("Bad exponent for cv calculation")
    end if
#endif

    ex = 1.e0_rt / (1.e0_rt - en)
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (em == 0.e0_rt) then
                alpha = const
             else
                alpha = const * state(i,j,k,URHO)**em
             end if

             frhoal = state(i,j,k,URHO) * alpha + tiny

             if (en == 0.e0_rt) then
                temp(i,j,k) = temp(i,j,k) / frhoal
             else
                teff = max(temp(i,j,k), tiny)
                temp(i,j,k) = ((1.e0_rt - en) * teff / frhoal)**ex
             end if

          end do
       end do
    end do

  end subroutine gtemp



  subroutine ca_compute_temp_given_rhoe(lo, hi, &
                                        temp, t_lo, t_hi, &
                                        state, s_lo, s_hi) &
                                        bind(C, name="ca_compute_temp_given_rhoe")

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: NVAR, URHO, UTEMP, UFS, UFX, small_temp
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains rhoe as input

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (temp(i,j,k) .le. 0.e0_rt) then

                temp(i,j,k) = small_temp

             else

                rhoInv = 1.e0_rt / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO)
                eos_state % T   = state(i,j,k,UTEMP)
                eos_state % e   =  temp(i,j,k)*rhoInv 
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
                eos_state % aux = state(i,j,k,UFX:UFX+naux -1) * rhoInv

                call eos(eos_input_re, eos_state)

                temp(i,j,k) = eos_state % T

             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_rhoe



  subroutine ca_compute_temp_given_cv(lo, hi, &
                                      temp, t_lo, t_hi, &
                                      state, s_lo, s_hi, &
                                      const_c_v, c_v_exp_m, c_v_exp_n) &
                                      bind(C, name="ca_compute_temp_given_cv")

    use meth_params_module, only: NVAR, URHO
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3)) ! temp contains rhoe as input
    real(rt), intent(in   ), value :: const_c_v, c_v_exp_m, c_v_exp_n

    integer  :: i, j, k
    real(rt) :: ex, alpha, rhoal, teff

    !$gpu

    ex = 1.e0_rt / (1.e0_rt - c_v_exp_n)

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (c_v_exp_m .eq. 0.e0_rt) then
                alpha = const_c_v
             else
                alpha = const_c_v * state(i,j,k,URHO)**c_v_exp_m
             endif

             rhoal = state(i,j,k,URHO) * alpha + 1.e-50_rt

             if (c_v_exp_n .eq. 0.e0_rt) then
                temp(i,j,k) = temp(i,j,k) / rhoal
             else
                teff = max(temp(i,j,k), 1.e-50_rt)
                temp(i,j,k) = ((1.e0_rt - c_v_exp_n) * teff / rhoal)**ex

             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_cv



  subroutine cfrhoe(lo, hi, &
                    frhoe, f_lo, f_hi, &
                    state, s_lo, s_hi) &
                    bind(C, name='cfrhoe')

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, UEINT

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: f_lo(3), f_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: frhoe(f_lo(1):f_hi(1),f_lo(2):f_hi(2),f_lo(3):f_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
             ! kin = 0.5e0_rt * (state(i,j,k,XMOM)   ** 2 +
             !                   state(i,j,k,XMOM+1) ** 2 +
             !                   state(i,j,k,XMOM+2) ** 2) /
             !                   state(i,j,k,DEN)
             ! frhoe(i,j,k) = state(i,j,k,EDEN) - kin
             frhoe(i,j,k) = state(i,j,k,UEINT)
          end do
       end do
    end do

  end subroutine cfrhoe



  subroutine rosse1(lo, hi, &
                    const, em, en, &
                    ep, nu, &
                    tf, kfloor, &
                    state, s_lo, s_hi, &
                    kappar, k_lo, k_hi) bind(C, name="rosse1")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UTEMP

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(inout) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ), value :: const, em, en, ep, nu, tf, kfloor

    real(rt) :: kf, teff
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             teff = max(state(i,j,k,UTEMP), tiny)
             teff = teff + tf * exp(-teff / (tf + tiny))

             kf = const * &
                  (state(i,j,k,URHO) ** em) * &
                  (teff ** (-en)) * &
                  (nu ** (ep))

             kappar(i,j,k) = max(kf, kfloor)

          end do
       end do
    end do

  end subroutine rosse1



  subroutine rosse1s(lo, hi, &
                     const, em, en, &
                     ep, sconst, &
                     sem, sen, &
                     sep, nu, &
                     tf, kfloor, &
                     state, s_lo, s_hi, &
                     kappar, k_lo, k_hi) bind(C, name="rosse1s")

    use amrex_fort_module, only: rt => amrex_real
    use meth_params_module, only: NVAR, URHO, UTEMP

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    integer,  intent(in   ) :: k_lo(3), k_hi(3)
    real(rt), intent(inout) :: kappar(k_lo(1):k_hi(1),k_lo(2):k_hi(2),k_lo(3):k_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    real(rt), intent(in   ), value :: const, em, en, ep, sconst, sem, sen, sep, nu, tf, kfloor

    real(rt) :: kf, teff, sct
    integer  :: i, j, k

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             teff = max(state(i,j,k,UTEMP), tiny)
             teff = teff + tf * exp(-teff / (tf + tiny))

             kf = const * &
                  (state(i,j,k,URHO) ** em) * &
                  (teff ** (-en)) * &
                  (nu ** (ep))

             sct = sconst * &
                  (state(i,j,k,URHO) ** sem) * &
                  (teff ** (-sen)) * &
                  (nu ** (sep))

             kappar(i,j,k) = max(kf + sct, kfloor)

          end do
       end do
    end do

  end subroutine rosse1s



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The following routined are used by neutrinos only.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ca_compute_temp_given_reye(lo, hi, &
                                        temp, t_lo, t_hi, &
                                        rhoe, r_lo, r_hi, &
                                        ye, y_lo, y_hi, &
                                        state, s_lo, s_hi) &
                                        bind(C, name='ca_compute_temp_given_reye')

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_re
    use meth_params_module, only: NVAR, URHO, UMX, UMY, UFS, UFX, small_temp
#ifndef AMREX_USE_GPU
    use castro_error_module, only: castro_error
#endif
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: r_lo(3), r_hi(3)
    integer,  intent(in   ) :: y_lo(3), y_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: rhoe(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3))
    real(rt), intent(in   ) :: ye(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             if (rhoe(i,j,k) .le. 0.e0_rt) then

                temp(i,j,k) = small_temp

             else

                rhoInv = 1.e0_rt / state(i,j,k,URHO)
                eos_state % rho = state(i,j,k,URHO)
                ! set initial guess of temperature
                eos_state % T = temp(i,j,k)
                eos_state % e = rhoe(i,j,k)*rhoInv 
                eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
                if (naux > 0) then
                   eos_state % aux = ye(i,j,k)
                end if

                call eos(eos_input_re, eos_state)

                temp(i,j,k) = eos_state % T

#ifndef AMREX_USE_GPU
                if (temp(i,j,k) .lt. 0.e0_rt) then
                   print *, 'negative temp in compute_temp_given_reye ', temp(i,j,k)
                   call castro_error("Error:: ca_compute_temp_given_reye")
                endif
#endif

             end if

          end do
       end do
    end do

  end subroutine ca_compute_temp_given_reye



  subroutine ca_compute_reye_given_ty(lo, hi, &
                                      rhoe, re_lo, re_hi, &
                                      rhoY, rY_lo, rY_hi, &
                                      temp, t_lo, t_hi, &
                                      ye, y_lo, y_hi, &
                                      state, s_lo, s_hi) &
                                      bind(C, name='ca_compute_reye_given_ty')

    use network, only: nspec, naux
    use eos_module, only: eos
    use eos_type_module, only: eos_t, eos_input_rt
    use meth_params_module, only: NVAR, URHO, UFS, UFX
    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: re_lo(3), re_hi(3)
    integer,  intent(in   ) :: rY_lo(3), ry_hi(3)
    integer,  intent(in   ) :: t_lo(3), t_hi(3)
    integer,  intent(in   ) :: y_lo(3), y_hi(3)
    integer,  intent(in   ) :: s_lo(3), s_hi(3)
    real(rt), intent(inout) :: rhoe(re_lo(1):re_hi(1),re_lo(2):re_hi(2),re_lo(3):re_hi(3))
    real(rt), intent(inout) :: rhoY(rY_lo(1):rY_hi(1),rY_lo(2):rY_hi(2),rY_lo(3):rY_hi(3))
    real(rt), intent(in   ) :: temp(t_lo(1):t_hi(1),t_lo(2):t_hi(2),t_lo(3):t_hi(3))
    real(rt), intent(in   ) :: ye(y_lo(1):y_hi(1),y_lo(2):y_hi(2),y_lo(3):y_hi(3))
    real(rt), intent(in   ) :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)

    integer      :: i, j, k
    real(rt)     :: rhoInv
    type (eos_t) :: eos_state

    !$gpu

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = 1.e0_rt / state(i,j,k,URHO)
             eos_state % rho = state(i,j,k,URHO)
             eos_state % T   = temp(i,j,k)
             eos_state % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv

             if (naux > 0) then
                eos_state % aux = ye(i,j,k)
                rhoY(i,j,k) = state(i,j,k,URHO)*ye(i,j,k)        
             end if

             call eos(eos_input_rt, eos_state)

             rhoe(i,j,k) = eos_state % rho * eos_state % e

          end do
       end do
    end do

  end subroutine ca_compute_reye_given_ty

end module rad_nd_module



!! -----------------------------------------------------------
!> @brief This routine is called at problem setup time and is used
!! to initialize values of physical constants used by the
!! radiation package.
!! -----------------------------------------------------------
subroutine ca_initradconstants(p, c, h, k, s, a, m, J_is_used) bind(C, name="ca_initradconstants")

  use fundamental_constants_module, only: c_fcm=>c_light, h_fcm=>hplanck, &
                                          k_fcm=>k_B, s_fcm=>sigma_SB, a_fcm=>n_A, ev2erg_fcm=>ev2erg
  use rad_params_module, only: pi, clight, hplanck
  use rad_params_module, only: kboltz, stefbol, arad, avogadro
  use rad_params_module, only: Hz2MeV, mev2erg, tiny
  use rad_params_module, only: radtoE  !, radtoJ, Etorad, radfluxtoF
  use rad_params_module, only: etafactor
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  real(rt) :: p, c, h, k, s, a, m
  integer :: J_is_used

  c = c_fcm
  h = h_fcm
  k = k_fcm
  s = s_fcm
  a = a_fcm
  m = 1.e6_rt * ev2erg_fcm

  pi       = p
  clight   = c
  hplanck  = h
  kboltz   = k
  stefbol  = s
  arad     = 4.*stefbol/clight
  avogadro = a
  mev2erg  = m
  Hz2MeV   = h / m
  tiny     = 1.e-50_rt

  if (J_is_used > 0) then
     radtoE = 4.e0_rt*pi/clight
     !           radtoJ = 1.0e0_rt
     !           Etorad = 1.e0_rt/radtoE
     !           radfluxtoF = 4.e0_rt*pi
     etafactor = 1.e0_rt
  else
     radtoE = 1.0e0_rt
     !           radtoJ = clight/(4.e0_rt*pi)
     !           Etorad = 1.0e0_rt
     !           radfluxtoF = 1.e0_rt
     etafactor = 4.e0_rt*pi/clight
  end if

end subroutine ca_initradconstants

! For single group, let set ngroups to 1.
subroutine ca_initsinglegroup(ngr) bind(C, name="ca_initsinglegroup")

  use rad_params_module, only: ngroups, nugroup, dnugroup, ng0, ng1
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer :: ngr

  ! Local variables
  integer :: i

  ng0 = 0
  ng1 = 0

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = 1.e0_rt  ! dummy
     dnugroup(i) = 1.e0_rt
  enddo

end subroutine ca_initsinglegroup

!! -----------------------------------------------------------
!> @brief This routine is called at problem setup time and is used
!! to initialize the arrays nugroup and dnugroup in
!! probdata with the neutrino group energies and widths.
!!
!! The widths are used to derive neutrino spectrum for plot files
!! -----------------------------------------------------------
subroutine ca_initgroups(nugr, dnugr, ngr, ngr0, ngr1)

  use rad_params_module, only: ngroups, ng0, ng1, nugroup, dnugroup, &
                               current_group, nnuspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt) :: nugr(0:ngr-1), dnugr(0:ngr-1)
  integer :: ngr, ngr0, ngr1

  ! Local variables
  integer :: i

  allocate(current_group, ng0, ng1, nnuspec)

  ng0     = ngr0
  ng1     = ngr1

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))

  do i = 0, ngroups-1
     nugroup(i)  = nugr(i)
     dnugroup(i) = dnugr(i)
  enddo

end subroutine ca_initgroups

subroutine ca_initgroups2(nugr, dnugr, xnugr, ngr)

  use rad_params_module, only: ngroups, nugroup, dnugroup, xnu, dlognu, lognugroup, &
                               current_group, ng0, ng1, nnuspec
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), xnugr(0:ngr)
  integer :: ngr

  ! Local variables
  integer   :: i

  allocate(current_group, ng0, ng1, nnuspec)

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups))
  allocate(dlognu(0:ngroups-1))
  allocate(lognugroup(0:ngroups-1))

  nugroup(:) = nugr(:)
  dnugroup(:) = dnugr(:)
  xnu(:) = xnugr(:)
  lognugroup(:) = log(nugroup)

  dlognu(0:ngroups-1) = log(xnu(1:ngroups)) - log(xnu(0:ngroups-1))

end subroutine ca_initgroups2

subroutine ca_initgroups3(nugr, dnugr, dlognugr, xnugr, ngr, ngr0, ngr1)
  ! used by MGFLDSolver

  use rad_params_module, only: ngroups, ng0, ng1, nnuspec, nradspec, nugroup, dnugroup, &
                               xnu, dlognu, lognugroup, erg2rhoYe, avogadro, hplanck, &
                               current_group
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  real(rt), intent(in) :: nugr(0:ngr-1), dnugr(0:ngr-1), dlognugr(0:ngr-1), xnugr(0:ngr+2)
  integer :: ngr, ngr0, ngr1

  ! Local variables
  integer :: i

  allocate(current_group, ng0, ng1, nnuspec)

  ng0     = ngr0
  ng1     = ngr1

  if (ng0 > 0) then
     if (ng1 .eq. 0) then
        nnuspec = 1  ! one neutrino species
     else if (ngroups .eq. ng0+ng1) then
        nnuspec = 2  ! two neutrino species
     else
        nnuspec = 3  ! three neutrino species
     end if
  else
     nnuspec = 0
  end if

  nradspec = max(nnuspec, 1)

  allocate(nugroup( 0:ngroups-1))
  allocate(dnugroup(0:ngroups-1))
  allocate(xnu(0:ngroups+2))
  allocate(dlognu(0:ngroups-1))
  allocate(erg2rhoYe(0:ngroups-1))
  allocate(lognugroup( 0:ngroups-1))

  nugroup(:) = nugr(:)
  dnugroup(:) = dnugr(:)
  xnu(:) = xnugr(:)
  dlognu(:) = dlognugr(:)
  lognugroup(:) = log(nugroup)

  erg2rhoYe = 0.e0_rt
  if (ng0 > 0) then
     erg2rhoYe(0:ng0-1) = 1.e0_rt / (avogadro*hplanck*nugroup(0:ng0-1))
     if (ng1 > 0) then
        erg2rhoYe(ng0:ng0+ng1-1) = -1.e0_rt / (avogadro*hplanck*nugroup(ng0:ng0+ng1-1))
     end if
  end if

end subroutine ca_initgroups3

!! -----------------------------------------------------------

subroutine ca_setgroup(igroup)

  use rad_params_module, only: current_group
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer :: igroup

  current_group = igroup

end subroutine ca_setgroup

!! -----------------------------------------------------------

subroutine ca_inelastic_sct(lo, hi, &
                            uu,uu_lo,uu_hi, &
                            Er,Er_lo,Er_hi, &
                            ks,ks_lo,ks_hi, &
                            dt) bind(C)

  use meth_params_module, only: NVAR, UEDEN, UEINT, UTEMP
  use rad_params_module, only: ngroups, nugroup, dlognu
  use radhydro_nd_module, only: inelastic_scatter
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: uu_lo(3),uu_hi(3)
  integer, intent(in) :: Er_lo(3),Er_hi(3)
  integer, intent(in) :: ks_lo(3),ks_hi(3)
  real(rt), intent(inout) :: uu(uu_lo(1):uu_hi(1),uu_lo(2):uu_hi(2),uu_lo(3):uu_hi(3),NVAR)
  real(rt), intent(inout) :: Er(Er_lo(1):Er_hi(1),Er_lo(2):Er_hi(2),Er_lo(3):Er_hi(3),0:ngroups-1)
  real(rt), intent(in   ) :: ks(ks_lo(1):ks_hi(1),ks_lo(2):ks_hi(2),ks_lo(3):ks_hi(3))
  real(rt), intent(in) :: dt

  integer :: i, j, k
  real(rt) :: Ertotold, Ertmp(0:ngroups-1), dEr
  real(rt) :: Erscale(0:ngroups-1)

  Erscale = nugroup*dlognu

  do       k = lo(3), hi(3)
     do    j = lo(2), hi(2)
        do i = lo(1), hi(1)
           Ertmp = Er(i,j,k,:)
           Ertotold = sum(Ertmp)
           Ertmp = Ertmp / Erscale

           call inelastic_scatter(uu(i,j,k,UTEMP), Ertmp, ks(i,j,k), dt, (/i,j,k/))

           Ertmp = Ertmp * Erscale
           dEr = sum(Ertmp) - Ertotold
           Er(i,j,k,:) = Ertmp
           uu(i,j,k,UEINT) = uu(i,j,k,UEINT) - dEr
           uu(i,j,k,UEDEN) = uu(i,j,k,UEDEN) - dEr
        end do
     end do
  end do

end subroutine ca_inelastic_sct

!! -----------------------------------------------------------

subroutine ca_compute_scattering(lo, hi, &
                                 kps,kps_lo,kps_hi, &
                                 sta,sta_lo,sta_hi)

  use rad_params_module, only : ngroups, nugroup
  use opacity_table_module, only : get_opacities
  use network, only : naux
  use meth_params_module, only : NVAR, URHO, UTEMP, UFX
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: kps_lo(3),kps_hi(3)
  integer, intent(in) :: sta_lo(3),sta_hi(3)
  real(rt), intent(inout) :: kps(kps_lo(1):kps_hi(1),kps_lo(2):kps_hi(2),kps_lo(3):kps_hi(3))
  real(rt), intent(in   ) :: sta(sta_lo(1):sta_hi(2),sta_lo(2):sta_hi(2),sta_lo(3):sta_hi(3),NVAR)

  integer :: i, j, k
  real(rt) :: kp, kr, nu, rho, temp, Ye
  logical, parameter :: comp_kp = .true.
  logical, parameter :: comp_kr = .true.

  ! scattering is assumed to be independent of nu.

  nu = nugroup(0)

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           rho = sta(i,j,k,URHO)
           temp = sta(i,j,k,UTEMP)
           if (naux > 0) then
              Ye = sta(i,j,k,UFX)
           else
              Ye = 0.e0_rt
           end if

           call get_opacities(kp, kr, rho, temp, Ye, nu, comp_kp, comp_kr)

           kps(i,j,k) = max(kr - kp, 0.e0_rt)
        end do
     end do
  end do

end subroutine ca_compute_scattering

subroutine ca_compute_scattering_2(lo, hi, &
                                   kps,kps_lo,kps_hi, &
                                   sta,sta_lo,sta_hi, &
                                   k0_p, m_p, n_p, &
                                   k0_r, m_r, n_r, &
                                   Tfloor, kfloor)

  use rad_params_module, only: ngroups, nugroup
  use meth_params_module, only: NVAR, URHO, UTEMP
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  integer, intent(in) :: lo(3), hi(3)
  integer, intent(in) :: kps_lo(3),kps_hi(3)
  integer, intent(in) :: sta_lo(3),sta_hi(3)
  real(rt), intent(inout) :: kps(kps_lo(1):kps_hi(1),kps_lo(2):kps_hi(2),kps_lo(3):kps_hi(3))
  real(rt), intent(in   ) :: sta(sta_lo(1):sta_hi(2),sta_lo(2):sta_hi(2),sta_lo(3):sta_hi(3),NVAR)
  real(rt), intent(in) :: k0_p, m_p, n_p
  real(rt), intent(in) :: k0_r, m_r, n_r
  real(rt), intent(in) :: Tfloor, kfloor

  integer :: i, j, k
  real(rt), parameter :: tiny = 1.0e-50_rt
  real(rt) :: Teff, k_p, k_r

  ! scattering is assumed to be independent of nu.

  if ( m_p.eq.0.e0_rt .and. n_p.eq.0.e0_rt .and. &
       m_r.eq.0.e0_rt .and. n_r.eq.0.e0_rt ) then
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              kps(i,j,k) = k0_r - k0_p
           end do
        end do
     end do
  else
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              Teff = max(sta(i,j,k,UTEMP), tiny)
              Teff = Teff + Tfloor * exp(-Teff / (Tfloor + tiny))
              k_p = k0_p * (sta(i,j,k,URHO) ** m_p) * (Teff ** (-n_p))
              k_r = k0_r * (sta(i,j,k,URHO) ** m_r) * (Teff ** (-n_r))
              kps(i,j,k) = max(k_r-k_p, kfloor)
           end do
        end do
     end do
  end if

end subroutine ca_compute_scattering_2
