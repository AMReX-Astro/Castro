
  subroutine ca_react_state(lo,hi, &
                            state,s_lo,s_hi, &
                            time,dt_react)

      use eos_module
      use network           , only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, dual_energy_eta3, allow_negative_energy
      use burner_module
      use bl_constants_module

      implicit none

      integer          :: lo(3), hi(3)
      integer          :: s_lo(3), s_hi(3)
      double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
      double precision :: time, dt_react

      integer          :: i, j, k
      double precision :: rhoInv, rho_e_K, delta_rho_e

      type (eos_t_3D)  :: state_in
      type (eos_t_3D)  :: state_out

      state_in  = eos_t_3D(lo,hi)
      state_out = eos_t_3D(lo,hi)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               rhoInv = ONE / state(i,j,k,URHO)

               state_in % rho(i,j,k)   = state(i,j,k,URHO)
               state_in % T(i,j,k)     = state(i,j,k,UTEMP)
               
               rho_e_K = state(i,j,k,UEDEN) - HALF * rhoInv * (state(i,j,k,UMX)**2 + state(i,j,k,UMY)**2 + state(i,j,k,UMZ)**2)

               ! Dual energy formalism: switch between e and (E - K) depending on (E - K) / E.

               if ( rho_e_K / state(i,j,k,UEDEN) .lt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
                  state_in % e(i,j,k) = rho_E_K * rhoInv
               else
                  state_in % e(i,j,k) = state(i,j,k,UEINT) * rhoInv
               endif

               state_in % xn(i,j,k,:)  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
               state_in % aux(i,j,k,:) = state(i,j,k,UFX:UFX+naux-1) * rhoInv
               
            enddo
         enddo
      enddo

      if (allow_negative_energy .eq. 0) state_in % reset = .true.
      if (allow_negative_energy .eq. 0) state_out % reset = .true.
      
      call burner(state_in, state_out, dt_react, time)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Note that we want to update the total energy by taking the difference of the old
               ! rho*e and the new rho*e. If the user wants to ensure that rho * E = rho * e + rho * K,
               ! this reset should be enforced through an appropriate choice for the dual energy 
               ! formalism parameter dual_energy_eta2 in reset_internal_energy.

               delta_rho_e = state(i,j,k,URHO) * state_out % e(i,j,k) - state(i,j,k,UEINT)

               state(i,j,k,UEINT)           = state(i,j,k,UEINT) + delta_rho_e
               state(i,j,k,UEDEN)           = state(i,j,k,UEDEN) + delta_rho_e
               state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO) * state_out % xn(i,j,k,:)
               state(i,j,k,UFX:UFX+naux -1) = state(i,j,k,URHO) * state_out % aux(i,j,k,:)
               state(i,j,k,UTEMP)           = state_out % T(i,j,k)

            enddo
         enddo
      enddo

      call eos_deallocate(state_in)
      call eos_deallocate(state_out)      

  end subroutine ca_react_state
