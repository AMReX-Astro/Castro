module reactions_module

  implicit none

  public

contains

#ifndef SDC

  subroutine ca_react_state(lo,hi, &
                            state,s_lo,s_hi, &
                            reactions,r_lo,r_hi, &
                            mask,m_lo,m_hi, &
                            time,dt_react) bind(C, name="ca_react_state")

    use network           , only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
         UFS, UFX, dual_energy_eta3
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use prob_params_module, only : dx_level, dim
    use amrinfo_module, only : amr_level
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
    double precision :: rhoInv, rho_e_K, delta_e, delta_rho_e, dx_min

    type (burn_t) :: burn_state_in, burn_state_out
    type (eos_t) :: eos_state_in, eos_state_out

    ! Minimum zone width

    dx_min = minval(dx_level(1:dim, amr_level))

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Don't burn on zones that we are intentionally masking out.

             if (mask(i,j,k) /= 1) cycle

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

             burn_state_in % dx = dx_min

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

#else

  subroutine ca_react_state(lo,hi, &
                            uold,uo_lo,uo_hi, &
                            unew,un_lo,un_hi, &
                            asrc,as_lo,as_hi, &
                            reactions,r_lo,r_hi, &
                            mask,m_lo,m_hi, &
                            time,dt_react) bind(C, name="ca_react_state")

    use network           , only : nspec, naux
    use meth_params_module, only : NVAR, URHO, UMX, UMZ, UEDEN, UEINT, UTEMP, &
                                   UFS, UFX, dual_energy_eta3, disable_shock_burning, &
                                   react_T_min, react_T_max, react_rho_min, react_rho_max
#ifdef SHOCK_VAR
    use meth_params_module, only : USHK
#endif
    use integrator_module, only : integrator
    use bl_constants_module, only : ZERO, HALF, ONE
    use sdc_type_module, only : sdc_t, SRHO, SMX, SMZ, SEDEN, SEINT, SFS

    implicit none

    integer          :: lo(3), hi(3)
    integer          :: uo_lo(3), uo_hi(3)
    integer          :: un_lo(3), un_hi(3)
    integer          :: as_lo(3), as_hi(3)
    integer          :: r_lo(3), r_hi(3)
    integer          :: m_lo(3), m_hi(3)
    double precision :: uold(uo_lo(1):uo_hi(1),uo_lo(2):uo_hi(2),uo_lo(3):uo_hi(3),NVAR)
    double precision :: unew(un_lo(1):un_hi(1),un_lo(2):un_hi(2),un_lo(3):un_hi(3),NVAR)
    double precision :: asrc(as_lo(1):as_hi(1),as_lo(2):as_hi(2),as_lo(3):as_hi(3),NVAR)
    double precision :: reactions(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec+2)
    integer          :: mask(m_lo(1):m_hi(1),m_lo(2):m_hi(2),m_lo(3):m_hi(3))
    double precision :: time, dt_react

    integer          :: i, j, k, n
    double precision :: rhooInv, rhonInv, rho_e_K, delta_e, delta_rho_e

    type (sdc_t) :: burn_state_in, burn_state_out

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

             ! Don't burn on zones that we are intentionally masking out.

             if (mask(i,j,k) /= 1) cycle

             ! Don't burn on zones inside shock regions, if the relevant option is set.

#ifdef SHOCK_VAR
             if (state(i,j,k,USHK) > ZERO .and. disable_shock_burning == 1) cycle
#endif

             ! Don't burn if we're outside of the relevant (rho, T) range.

             if (uold(i,j,k,UTEMP) < react_T_min .or. uold(i,j,k,UTEMP) > react_T_max .or. &
                 uold(i,j,k,URHO) < react_rho_min .or. uold(i,j,k,URHO) > react_rho_max) cycle

             ! Feed in the old-time state data.

             burn_state_in % y(SRHO)            = uold(i,j,k,URHO)
             burn_state_in % y(SMX:SMZ)         = uold(i,j,k,UMX:UMZ)
             burn_state_in % y(SEDEN)           = uold(i,j,k,UEDEN)
             burn_state_in % y(SEINT)           = uold(i,j,k,UEINT)
             burn_state_in % y(SFS:SFS+nspec-1) = uold(i,j,k,UFS:UFS+nspec-1)

             ! Tell the integrator about the non-reacting source terms.

             burn_state_in % ydot_a(SRHO)            = asrc(i,j,k,URHO)
             burn_state_in % ydot_a(SMX:SMZ)         = asrc(i,j,k,UMX:UMZ)
             burn_state_in % ydot_a(SEDEN)           = asrc(i,j,k,UEDEN)
             burn_state_in % ydot_a(SEINT)           = asrc(i,j,k,UEINT)
             burn_state_in % ydot_a(SFS:SFS+nspec-1) = asrc(i,j,k,UFS:UFS+nspec-1)

             ! Dual energy formalism: in doing EOS calls in the burn,
             ! switch between e and (E - K) depending on (E - K) / E.

             rhooInv = ONE / uold(i,j,k,URHO)

             rho_e_K = uold(i,j,k,UEDEN) - HALF * rhooInv * sum(uold(i,j,k,UMX:UMZ)**2)

             if ( rho_e_K / uold(i,j,k,UEDEN) .gt. dual_energy_eta3 .and. rho_e_K .gt. ZERO ) then
                burn_state_in % T_from_eden = .true.
             else
                burn_state_in % T_from_eden = .false.
             endif

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

             call integrator(burn_state_in, burn_state_out, dt_react, time)

             ! Update the state data.

             unew(i,j,k,URHO)            = burn_state_out % y(SRHO)
             unew(i,j,k,UMX:UMZ)         = burn_state_out % y(SMX:SMZ)
             unew(i,j,k,UEDEN)           = burn_state_out % y(SEDEN)
             unew(i,j,k,UEINT)           = burn_state_out % y(SEINT)
             unew(i,j,k,UFS:UFS+nspec-1) = burn_state_out % y(SFS:SFS+nspec-1)

             ! Add burning rates to reactions MultiFab, but be
             ! careful because the reactions and state MFs may
             ! not have the same number of ghost cells.

             if ( i .ge. r_lo(1) .and. i .le. r_hi(1) .and. &
                  j .ge. r_lo(2) .and. j .le. r_hi(2) .and. &
                  k .ge. r_lo(3) .and. k .le. r_hi(3) ) then

                rhonInv = ONE / unew(i,j,k,URHO)

                do n = 1, nspec
                   reactions(i,j,k,n) = (unew(i,j,k,UFS+n-1) * rhoninv - uold(i,j,k,UFS+n-1) * rhooinv) / dt_react
                enddo
                reactions(i,j,k,nspec+1) = (unew(i,j,k,UEINT) * rhonInv - uold(i,j,k,UEINT) * rhooInv) / dt_react
                reactions(i,j,k,nspec+2) = (unew(i,j,k,UEINT) - uold(i,j,k,UEINT)) / dt_react

             endif

          enddo
       enddo
    enddo

  end subroutine ca_react_state

#endif

end module reactions_module
