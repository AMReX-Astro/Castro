
  subroutine ca_react_state(lo,hi,s_in,so_l1,so_h1,s_out,sn_l1,sn_h1, &
                            reaction_terms,r_l1,r_h1,time,dt_react)

      use eos_module
      use network           , only: nspec, naux
      use meth_params_module, only: NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                    small_temp, small_dens, allow_negative_energy
      use burner_module
      use bl_constants_module

      implicit none

      integer          :: lo(1),hi(1)
      integer          :: so_l1,so_h1
      integer          :: sn_l1,sn_h1
      double precision :: state(s_l1:s_h1,NVAR)
      double precision :: time,dt_react

      integer :: i
      double precision :: rhoInv

      type (eos_t_1D) :: state_in
      type (eos_t_1D) :: state_out

      state_in = eos_t_1D(lo, hi)
      state_out = eos_t_1D(lo, hi)

      do i = lo(1), hi(1)

         rhoInv = ONE / s_in(i,URHO)

         state_in % rho(i)   = s_in(i,URHO)
         state_in % T(i)     = s_in(i,UTEMP)
         state_in % e(i)     = s_in(i,UEINT) * rhoInv
         state_in % xn(i,:)  = s_in(i,UFS:UFS+nspec-1) * rhoInv
         state_in % aux(i,:) = s_in(i,UFX:UFX+naux-1) * rhoInv

      enddo

      call burner(state_in, state_out, dt_react, time)

      do i = lo(1), hi(1)

         ! Note that we want to update the total energy by taking the difference of the old
         ! rho*e and the new rho*e rather than adding rho*e to KE because we don't necessarily
         ! want to force here that rho * E = rho * e + rho * K where K is defined in terms of
         ! the momenta. If this is desired that can always be forced through reset_internal_energy.

         state(i,UEINT)           = state(i,URHO) * state_out % e(i)
         state(i,UEDEN)           = state(i,UEDEN) + state(i,URHO) * (state_out % e(i) - state_in % e(i))
         state(i,UFS:UFS+nspec-1) = state(i,URHO) * state_out % xn(i,:)
         state(i,UFX:UFX+naux -1) = state(i,URHO) * state_out % aux(i,:)
         state(i,UTEMP)           = state_out % T(i)

      enddo
      
  end subroutine ca_react_state

