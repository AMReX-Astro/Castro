  subroutine ca_react_state(lo,hi, &
                            state,s_l1,s_l2,s_h1,s_h2,&
                            time,dt_react)

      use eos_module
      use network           , only: nspec, naux
      use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                    small_dens, small_temp, allow_negative_energy, dual_energy_eta3
      use burner_module
      use bl_constants_module

      implicit none

      integer          :: lo(2),hi(2)
      integer          :: s_l1,s_h1,s_l2,s_h2
      double precision :: state(s_l1:s_h1,s_l2:s_h2,NVAR)
      double precision :: time,dt_react

      integer          :: i, j
      double precision :: rhoInv, rho_e_K, delta_rho_e

      type (eos_t_2D)  :: state_in
      type (eos_t_2D)  :: state_out

      state_in  = eos_t_2D(lo, hi)
      state_out = eos_t_2D(lo, hi)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            rhoInv = ONE / s_in(i,j,URHO)

            state_in % rho(i,j)   = s_in(i,j,URHO)
            state_in % T(i,j)     = s_in(i,j,UTEMP)

            rho_e_K = state(i,j,UEDEN) - HALF * rhoInv * (state(i,j,UMX)**2 + state(i,j,UMY)**2 + state(i,j,UMZ)**2)

            ! Dual energy formalism: switch between e and (E - K) depending on (E - K) / E.

            if ( rho_e_K / state(i,j,UEDEN) .lt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
               state_in % e(i,j) = rho_E_K * rhoInv
            else
               state_in % e(i,j) = state(i,j,UEINT) * rhoInv
            endif

            state_in % xn(i,j,:)  = s_in(i,j,UFS:UFS+nspec-1) * rhoInv
            state_in % aux(i,j,:) = s_in(i,j,UFX:UFX+naux-1) * rhoInv

         enddo
      enddo

      call burner(state_in, state_out, dt_react, time)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! Note that we want to update the total energy by taking the difference of the old
            ! rho*e and the new rho*e. If the user wants to ensure that rho * E = rho * e + rho * K,
            ! this reset should be enforced through an appropriate choice for the dual energy 
            ! formalism parameter dual_energy_eta2 in reset_internal_energy.

            delta_rho_e = state(i,j,URHO) * state_out % e(i,j) - state(i,j,UEINT)

            state(i,j,UEINT)           = state(i,j,UEINT) + delta_rho_e
            state(i,j,UEDEN)           = state(i,j,UEDEN) + delta_rho_e
            state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO) * state_out % xn(i,j,:)
            state(i,j,UFX:UFX+naux -1) = state(i,j,URHO) * state_out % aux(i,j,:)
            state(i,j,UTEMP)           = state_out % T(i,j)

         enddo
      enddo

  end subroutine ca_react_state

