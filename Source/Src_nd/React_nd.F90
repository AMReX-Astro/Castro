module reactions_module

  implicit none

  public

contains

  subroutine ca_react_state(lo,hi, &
                            state,s_lo,s_hi, &
                            reactions,r_lo,r_hi, &
                            mask,m_lo,m_hi, &
                            time,dt_react) bind(C, name="ca_react_state")

    use network           , only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS, dual_energy_eta3
#if naux > 0
    use meth_params_module, only : UFX
#endif
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
#ifdef ACC
    use meth_params_module, only : do_acc
#endif
    use burner_module
    use burn_type_module
    use bl_constants_module

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: s_lo(3), s_hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: m_lo(3), m_hi(3)
    double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
    double precision :: reactions(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec+2)
    integer          :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    double precision :: time, dt_react

    integer          :: i, j, k, n
    double precision :: rhoInv, rho_e_K, delta_e, delta_rho_e

    type (burn_t) :: burn_state_in, burn_state_out
    type (eos_t) :: eos_state_in, eos_state_out

    !$acc data &
    !$acc copyin(lo, hi, r_lo, r_hi, s_lo, s_hi, m_lo, m_hi, dt_react, time) &
    !$acc copyin(mask) &
    !$acc copy(state, reactions) if(do_acc == 1)

    !$acc parallel if(do_acc == 1)

    !$acc loop gang vector collapse(3) &
    !$acc private(rhoInv, rho_e_K, delta_e, delta_rho_e) &
    !$acc private(eos_state_in, eos_state_out, burn_state_in, burn_state_out) &
    !$acc private(i,j,k)
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

#if naux > 0
             do n = 1, naux
                burn_state_in % aux(n) = state(i,j,k,UFX+n-1) * rhoInv
             enddo
#endif

             ! Ensure that the temperature going in is consistent with the internal energy.

             call burn_to_eos(burn_state_in, eos_state_in)
             call eos(eos_input_re, eos_state_in)
             call eos_to_burn(eos_state_in, burn_state_in)

             if (i >= lo(1) .and. i <= hi(1)) then
                burn_state_in % i = i
             else
                burn_state_in % i = -1
             endif

             if (j >= lo(2) .and. j <= hi(2)) then 
                burn_state_in % j = j
             else
                burn_state_in % j = -1
             endif

             if (k >= lo(3) .and. k <= hi(3)) then
                burn_state_in % k = k
             else
                burn_state_in % k = -1
             endif

             ! Now reset the internal energy to zero for the burn state.

             burn_state_in % e = ZERO

             burn_state_in % shock = .false.

#ifdef SHOCK_VAR
             if (state(i,j,k,USHK) > ZERO) then
                burn_state_in % shock = .true.
             endif
#endif

             ! Initialize the final state so that it has valid data in case this zone is masked out.

             burn_state_out = burn_state_in

             ! Don't burn on zones that we are intentionally masking out.

             if (mask(i,j,k) == 1) then
                call burner(burn_state_in, burn_state_out, dt_react, time)
             endif

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

#if naux > 0
             do n = 1, naux
                state(i,j,k,UFX+n-1)  = state(i,j,k,URHO) * burn_state_out % aux(n)
             enddo
#endif

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

    !$acc end parallel

    !$acc end data

  end subroutine ca_react_state

end module reactions_module
