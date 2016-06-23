module reactions_module

  implicit none

  public

contains

  subroutine ca_react_state(lo,hi, &
                            state,s_lo,s_hi, &
                            reactions,r_lo,r_hi, &
                            time,dt_react) bind(C, name="ca_react_state")

    use network           , only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
         UFS, UFX, dual_energy_eta3, allow_negative_energy, USHK
    use burner_module
    use burn_type_module
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: r_lo(3), r_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: reactions(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec+2)
    double precision :: time, dt_react

    integer          :: i, j, k, n
    double precision :: rhoInv, rho_e_K, delta_e, delta_rho_e

    type (burn_t) :: burn_state_in, burn_state_out
    type (eos_t) :: eos_state_in, eos_state_out

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)

             burn_state_in % rho = state(i,j,k,URHO)
             burn_state_in % T   = state(i,j,k,UTEMP)

             rho_e_K = state(i,j,k,UEDEN) - HALF * rhoInv * sum(state(i,j,k,UMX:UMZ)**2)

             ! Dual energy formalism: switch between e and (E - K) depending on (E - K) / E.

             if ( rho_e_K / state(i,j,k,UEDEN) .gt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
                burn_state_in % e = rho_E_K * rhoInv
             else
                burn_state_in % e = state(i,j,k,UEINT) * rhoInv
             endif

             do n = 1, nspec
                burn_state_in % xn(n) = state(i,j,k,UFS+n-1) * rhoInv
             enddo

             do n = 1, naux
                burn_state_in % aux(n) = state(i,j,k,UFX+n-1) * rhoInv
             enddo

             ! Ensure that the temperature going in is consistent with the internal energy.

             call burn_to_eos(burn_state_in, eos_state_in)
             call eos(eos_input_re, eos_state_in)
             call eos_to_burn(eos_state_in, burn_state_in)

             burn_state_in % i = i
             burn_state_in % j = j
             burn_state_in % k = k

             ! Now reset the internal energy to zero for the burn state.

             burn_state_in % e = ZERO

             burn_state_in % shock = .false.

#ifdef SHOCK_VAR
             if (state(i,j,k,USHK) > ZERO) then
                burn_state_in % shock = .true.
             endif
#endif

             call burner(burn_state_in, burn_state_out, dt_react, time)

             ! Note that we want to update the total energy by taking the difference of the old
             ! rho*e and the new rho*e. If the user wants to ensure that rho * E = rho * e + rho * K,
             ! this reset should be enforced through an appropriate choice for the dual energy 
             ! formalism parameter dual_energy_eta2 in reset_internal_energy.

             delta_e     = burn_state_out % e - burn_state_in % e
             delta_rho_e = burn_state_out % rho * delta_e

             state(i,j,k,UEINT) = state(i,j,k,UEINT) + delta_rho_e
             state(i,j,k,UEDEN) = state(i,j,k,UEDEN) + delta_rho_e

             do n = 1, nspec
                state(i,j,k,UFS+n-1) = state(i,j,k,URHO) * burn_state_out % xn(n)
             enddo

             do n = 1, naux
                state(i,j,k,UFX+n-1)  = state(i,j,k,URHO) * burn_state_out % aux(n)
             enddo

             ! Do an EOS call to get a temperature consistent with the new energy.
             ! Note that the burn state only contains the energy released during the burn,
             ! so we need to add to it the initial zone energy.

             call burn_to_eos(burn_state_out, eos_state_out)
             eos_state_out % e = eos_state_in % e + delta_e
             call eos(eos_input_re, eos_state_out)
             call eos_to_burn(eos_state_out, burn_state_out)
             
             state(i,j,k,UTEMP) = burn_state_out % T

             ! Add burning rates to reactions MultiFab, but be
             ! careful because the reactions and state MFs may
             ! not have the same number of ghost cells.

             if ( i .ge. r_lo(1) .and. i .le. r_hi(1) .and. &
                  j .ge. r_lo(2) .and. j .le. r_hi(2) .and. &
                  k .ge. r_lo(3) .and. k .le. r_hi(3) ) then

                do n = 1, nspec
                   reactions(i,j,k,n) = (burn_state_out % xn(n) - burn_state_in % xn(n)) / dt_react
                enddo
                reactions(i,j,k,nspec+1) = delta_e / dt_react
                reactions(i,j,k,nspec+2) = delta_rho_e / dt_react

             endif

          enddo
       enddo
    enddo

  end subroutine ca_react_state

end module reactions_module
