  subroutine ca_react_state(lo,hi, &
                            state,s_l1,s_l2,s_h1,s_h2,&
                            time,dt_react)

      use eos_module
      use network           , only: nspec, naux
      use meth_params_module, only: NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                    small_dens, small_temp, allow_negative_energy
      use burner_module
      use bl_constants_module

      implicit none

      integer          :: lo(2),hi(2)
      integer          :: s_l1,s_h1,s_l2,s_h2
      double precision :: state(s_l1:s_h1,s_l2:s_h2,NVAR)
      double precision :: time,dt_react

      integer          :: i, j
      double precision :: rhoInv

      type (eos_t_2D)  :: state_in
      type (eos_t_2D)  :: state_out

      state_in = eos_t_2D(lo, hi)
      state_out = eos_t_2D(lo, hi)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            rhoInv = ONE / s_in(i,j,URHO)

            state_in % rho(i,j)   = s_in(i,j,URHO)
            state_in % T(i,j)     = s_in(i,j,UTEMP)
            state_in % e(i,j)     = s_in(i,j,UEINT) * rhoInv
            state_in % xn(i,j,:)  = s_in(i,j,UFS:UFS+nspec-1) * rhoInv
            state_in % aux(i,j,:) = s_in(i,j,UFX:UFX+naux-1) * rhoInv

         enddo
      enddo

      call burner(state_in, state_out, dt_react, time)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! Note that we want to update the total energy by taking the difference of the old
            ! rho*e and the new rho*e rather than adding rho*e to KE because we don't necessarily
            ! want to force here that rho * E = rho * e + rho * K where K is defined in terms of
            ! the momenta. If this is desired that can always be forced through reset_internal_energy.

            state(i,j,UEINT)           = state(i,j,URHO) * state_out % e(i,j)
            state(i,j,UEDEN)           = state(i,j,UEDEN) + state(i,j,URHO) * (state_out % e(i,j) - state_in % e(i,j))
            state(i,j,UFS:UFS+nspec-1) = state(i,j,URHO) * state_out % xn(i,j,:)
            state(i,j,UFX:UFX+naux -1) = state(i,j,URHO) * state_out % aux(i,j,:)
            state(i,j,UTEMP)           = state_out % T(i,j)

         enddo
      enddo

  end subroutine ca_react_state

