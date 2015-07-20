  subroutine ca_react_state(lo,hi, &
                            state,s_l1,s_h1, &
                            time,dt_react)

      use eos_module
      use network           , only: nspec, naux
      use meth_params_module, only: NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                    small_temp, small_dens, allow_negative_energy, dual_energy_eta3
      use burner_module
      use bl_constants_module

      implicit none

      integer          :: lo(1),hi(1)
      integer          :: s_l1,s_h1
      double precision :: state(s_l1:s_h1,NVAR)
      double precision :: time,dt_react

      integer          :: i
      double precision :: rhoInv, rho_e_K, delta_rho_e

      type (eos_t_1D)  :: state_in
      type (eos_t_1D)  :: state_out

      state_in  = eos_t_1D(lo, hi)
      state_out = eos_t_1D(lo, hi)

      do i = lo(1), hi(1)

         rhoInv = ONE / s_in(i,URHO)

         state_in % rho(i)   = s_in(i,URHO)
         state_in % T(i)     = s_in(i,UTEMP)

         rho_e_K = state(i,UEDEN) - HALF * rhoInv * (state(i,UMX)**2 + state(i,UMY)**2 + state(i,UMZ)**2)

         ! Dual energy formalism: switch between e and (E - K) depending on (E - K) / E.

         if ( rho_e_K / state(i,UEDEN) .lt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
            state_in % e(i) = rho_E_K * rhoInv
         else
            state_in % e(i) = state(i,UEINT) * rhoInv
         endif

         state_in % xn(i,:)  = s_in(i,UFS:UFS+nspec-1) * rhoInv
         state_in % aux(i,:) = s_in(i,UFX:UFX+naux-1) * rhoInv

      enddo

      call burner(state_in, state_out, dt_react, time)

      do i = lo(1), hi(1)

         ! Note that we want to update the total energy by taking the difference of the old
         ! rho*e and the new rho*e. If the user wants to ensure that rho * E = rho * e + rho * K,
         ! this reset should be enforced through an appropriate choice for the dual energy 
         ! formalism parameter dual_energy_eta2 in reset_internal_energy.

         delta_rho_e = state(i,URHO) * state_out % e(i) - state(i,UEINT)

         state(i,UEINT)           = state(i,UEINT) + delta_rho_e
         state(i,UEDEN)           = state(i,UEDEN) + delta_rho_e
         state(i,UFS:UFS+nspec-1) = state(i,URHO) * state_out % xn(i,:)
         state(i,UFX:UFX+naux -1) = state(i,URHO) * state_out % aux(i,:)
         state(i,UTEMP)           = state_out % T(i)

      enddo
      
  end subroutine ca_react_state

