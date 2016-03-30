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

    integer          :: i, j, k
    double precision :: rhoInv, rho_e_K, delta_x(nspec), delta_e, delta_rho_e

    type (burn_t) :: state_in
    type (burn_t) :: state_out

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             rhoInv = ONE / state(i,j,k,URHO)

             state_in % rho = state(i,j,k,URHO)
             state_in % T   = state(i,j,k,UTEMP)

             rho_e_K = state(i,j,k,UEDEN) - HALF * rhoInv * sum(state(i,j,k,UMX:UMZ)**2)

             ! Dual energy formalism: switch between e and (E - K) depending on (E - K) / E.

             if ( rho_e_K / state(i,j,k,UEDEN) .lt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
                state_in % e = rho_E_K * rhoInv
             else
                state_in % e = state(i,j,k,UEINT) * rhoInv
             endif

             state_in % xn  = state(i,j,k,UFS:UFS+nspec-1) * rhoInv
             state_in % aux = state(i,j,k,UFX:UFX+naux-1) * rhoInv

             state_in % shock = .false.

#ifdef SHOCK_VAR
             if (state(i,j,k,USHK) > ZERO) then
                state_in % shock = .true.
             endif
#endif

             call burner(state_in, state_out, dt_react, time)

             ! Note that we want to update the total energy by taking the difference of the old
             ! rho*e and the new rho*e. If the user wants to ensure that rho * E = rho * e + rho * K,
             ! this reset should be enforced through an appropriate choice for the dual energy 
             ! formalism parameter dual_energy_eta2 in reset_internal_energy.

             delta_x     = state_out % xn - state_in % xn
             delta_e     = state_out % e - state_in % e
             delta_rho_e = state_out % rho * delta_e

             state(i,j,k,UEINT)           = state(i,j,k,UEINT) + delta_rho_e
             state(i,j,k,UEDEN)           = state(i,j,k,UEDEN) + delta_rho_e
             state(i,j,k,UFS:UFS+nspec-1) = state(i,j,k,URHO) * state_out % xn
             state(i,j,k,UFX:UFX+naux-1)  = state(i,j,k,URHO) * state_out % aux
             state(i,j,k,UTEMP)           = state_out % T

             ! Add burning rates to reactions MultiFab, but be
             ! careful because the reactions and state MFs may
             ! not have the same number of ghost cells.

             if ( i .ge. r_lo(1) .and. i .le. r_hi(1) .and. &
                  j .ge. r_lo(2) .and. j .le. r_hi(2) .and. &
                  k .ge. r_lo(3) .and. k .le. r_hi(3) ) then

                reactions(i,j,k,1:nspec) = delta_x / dt_react
                reactions(i,j,k,nspec+1) = delta_e / dt_react
                reactions(i,j,k,nspec+2) = delta_rho_e / dt_react

             endif

          enddo
       enddo
    enddo

  end subroutine ca_react_state

end module reactions_module
