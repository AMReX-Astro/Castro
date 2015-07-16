
  subroutine ca_react_state(lo,hi,s_in,so_l1,so_h1,s_out,sn_l1,sn_h1, &
                            reaction_terms,r_l1,r_h1,time,dt_react)

      use eos_module
      use network           , only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                     small_temp, small_dens, allow_negative_energy
      use burner_module
      use bl_constants_module

      implicit none

      integer :: lo(1),hi(1)
      integer :: so_l1,so_h1
      integer :: sn_l1,sn_h1
      integer ::  r_l1, r_h1
      double precision :: s_in (so_l1:so_h1,NVAR)
      double precision :: s_out(sn_l1:sn_h1,NVAR)
      double precision :: reaction_terms(r_l1:r_h1,nspec+2)
      double precision :: time,dt_react

      integer :: i
      double precision :: rhoInv

      type (eos_t_1D) :: state_in
      type (eos_t_1D) :: state_out

      state_in = eos_t_1D(lo, hi)
      state_out = eos_t_1D(lo, hi)

      do i = lo(1), hi(1)

         ! We do not want to integrate zones with rho < small_dens,
         ! but the hydro calls enforce_minimum_density, so we shouldn't
         ! be getting any zones like that. If so, throw an error.

         if (s_in(i,URHO) .lt. small_dens) then
            print *,'... rho < small_dens in react_state: ', i, s_in(i,URHO)
            call bl_error("Error:: React_1d.f90 :: ca_react_state")
         endif

         rhoInv            = ONE / s_in(i,URHO)

         state_in % rho(i)   = s_in(i,URHO)
         state_in % T(i)     = s_in(i,UTEMP)
         state_in % e(i)     = s_in(i,UEINT) * rhoInv
         state_in % xn(i,:)  = s_in(i,UFS:UFS+nspec-1) * rhoInv
         state_in % aux(i,:) = s_in(i,UFX:UFX+naux-1) * rhoInv

         ! The energy should never be negative coming into this call
         ! because of reset_internal_energy, so throw an error if it happens.

         if (allow_negative_energy .eq. 0 .and. state_in % e(i) .le. ZERO) then
            print *,'... e negative in react_state: ', i, state_in % e(i)
            call bl_error("Error:: React_1d.f90 :: ca_react_state")
         endif

      enddo

      call burner(state_in, state_out, dt_react, time)

      do i = lo(1), hi(1)

         ! Note that we want to update the total energy by taking the difference of the old
         ! rho*e and the new rho*e rather than adding rho*e to KE because we don't necessarily
         ! want to force here that rho * E = rho * e + rho * K where K is defined in terms of
         ! the momenta. If this is desired that can always be forced through reset_internal_energy.

         s_out(i,URHO)            = s_in(i,URHO)
         s_out(i,UEINT)           = s_out(i,URHO) * state_out % e(i)
         s_out(i,UEDEN)           = s_out(i,UEDEN) + s_out(i,URHO) * (state_out % e(i) - state_in % e(i))
         s_out(i,UFS:UFS+nspec-1) = s_out(i,URHO) * state_out % xn(i,:)
         s_out(i,UFX:UFX+naux -1) = state_in % aux(i,:)
         s_out(i,UTEMP)           = state_out % T(i)

         if (i.ge.r_l1 .and. i.le.r_h1) then
            reaction_terms(i,1:nspec) = reaction_terms(i,1:nspec) &
                 + (state_out % xn(i,:) - state_in % xn(i,:))
            reaction_terms(i,nspec+1) = reaction_terms(i,nspec+1) &
                 + (state_out % e(i) - state_in % e(i))
            reaction_terms(i,nspec+2) = reaction_terms(i,nspec+2) &
                 + s_out(i,URHO) * (state_out % e(i) - state_in % e(i))
         endif

      enddo
      
  end subroutine ca_react_state

