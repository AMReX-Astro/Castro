  subroutine ca_react_state(lo,hi, &
                            s_in ,so_l1,so_l2,so_h1,so_h2,&
                            s_out,sn_l1,sn_l2,sn_h1,sn_h2,&
                            reaction_terms,r_l1,r_l2,r_h1,r_h2, &
                            time,dt_react)

      use eos_module
      use network           , only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UTEMP, UFS, UFX, &
                                     small_dens, small_temp, allow_negative_energy
      use burner_module
      use bl_constants_module

      implicit none

      integer lo(2),hi(2)
      integer so_l1,so_h1,so_l2,so_h2
      integer sn_l1,sn_h1,sn_l2,sn_h2
      integer  r_l1, r_h1, r_l2, r_h2
      double precision s_in (so_l1:so_h1,so_l2:so_h2,NVAR)
      double precision s_out(sn_l1:sn_h1,sn_l2:sn_h2,NVAR)
      double precision reaction_terms(r_l1:r_h1,r_l2:r_h2,nspec+2)
      double precision time,dt_react

      integer          :: i,j
      double precision :: rhoInv

      type (eos_t) :: state_in(lo(1):hi(1),lo(2):hi(2))
      type (eos_t) :: state_out(lo(1):hi(1),lo(2):hi(2))

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! We do not want to integrate zones with rho < small_dens,
            ! but the hydro calls enforce_minimum_density, so we shouldn't
            ! be getting any zones like that. If so, throw an error.

            if (s_in(i,j,URHO) .lt. small_dens) then
               print *,'... rho < small_dens in react_state: ', i, j, s_in(i,j,URHO)
               call bl_error("Error:: React_2d.f90 :: ca_react_state")
            endif

            rhoInv              = ONE / s_in(i,j,URHO)

            state_in(i,j) % rho = s_in(i,j,URHO)
            state_in(i,j) % T   = s_in(i,j,UTEMP)
            state_in(i,j) % e   = s_in(i,j,UEINT) * rhoInv
            state_in(i,j) % xn  = s_in(i,j,UFS:UFS+nspec-1) * rhoInv
            state_in(i,j) % aux = s_in(i,j,UFX:UFX+naux-1) * rhoInv
            state_in(i,j) % loc = (/ i, j, -99 /)

            ! The energy should never be negative coming into this call
            ! because of reset_internal_energy, so throw an error if it happens.

            if (allow_negative_energy .eq. 0 .and. state_in(i,j) % e .le. ZERO) then
               print *,'... e negative in react_state: ', i, j, state_in(i,j) % e
               call bl_error("Error:: React_2d.f90 :: ca_react_state")
            endif
         enddo
      enddo

      call burner(state_in, state_out, dt_react, time)

      do j = lo(2), hi(2)
         do i = lo(1), hi(1)

            ! Note that we want to update the total energy by taking the difference of the old
            ! rho*e and the new rho*e rather than adding rho*e to KE because we don't necessarily
            ! want to force here that rho * E = rho * e + rho * K where K is defined in terms of
            ! the momenta. If this is desired that can always be forced through reset_internal_energy.

            s_out(i,j,URHO)            = s_in(i,j,URHO)
            s_out(i,j,UEINT)           = s_out(i,j,URHO) * state_out(i,j) % e
            s_out(i,j,UEDEN)           = s_out(i,j,UEDEN) + s_out(i,j,URHO) * (state_out(i,j) % e - state_out(i,j) % e)
            s_out(i,j,UFS:UFS+nspec-1) = s_out(i,j,URHO) * state_out(i,j) % xn(:)
            s_out(i,j,UFX:UFX+naux -1) = state_in(i,j) % aux(:)
            s_out(i,j,UTEMP)           = state_out(i,j) % T

            if (i.ge.r_l1 .and. i.le.r_h1 .and. j.ge.r_l2 .and. j.le.r_h2) then
               reaction_terms(i,j,1:nspec) = reaction_terms(i,j,1:nspec) &
                    + (state_out(i,j) % xn(:) - state_in(i,j) % xn(:))
               reaction_terms(i,j,nspec+1) = reaction_terms(i,j,nspec+1) &
                    + (state_out(i,j) % e - state_in(i,j) % e)
               reaction_terms(i,j,nspec+2) = reaction_terms(i,j,nspec+2) &
                    + s_out(i,j,URHO) * (state_out(i,j) % e - state_in(i,j) % e)
            endif

         enddo
      enddo

  end subroutine ca_react_state

