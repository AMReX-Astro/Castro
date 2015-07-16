
  subroutine ca_react_state(lo,hi, &
                            s_in ,so_l1,so_l2,so_l3,so_h1,so_h2,so_h3, &
                            s_out,sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                            reaction_terms,r_l1,r_l2,r_l3,r_h1,r_h2,r_h3, &
                            time,dt_react)

      use eos_module
      use network           , only : nspec, naux
      use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT, UTEMP, &
                                     UFS, UFX, small_dens, small_temp, allow_negative_energy
      use burner_module
      use bl_constants_module

      implicit none

      integer lo(3),hi(3)
      integer so_l1,so_h1,so_l2,so_h2,so_l3,so_h3
      integer sn_l1,sn_h1,sn_l2,sn_h2,sn_l3,sn_h3
      integer  r_l1, r_h1, r_l2, r_h2, r_l3, r_h3
      double precision s_in (so_l1:so_h1,so_l2:so_h2,so_l3:so_h3,NVAR)
      double precision s_out(sn_l1:sn_h1,sn_l2:sn_h2,sn_l3:sn_h3,NVAR)
      double precision reaction_terms(r_l1:r_h1,r_l2:r_h2,r_l3:r_h3,nspec+2)
      double precision time, dt_react

      integer          :: i,j,k
      double precision :: rhoInv

      type (eos_t_3D) :: state_in
      type (eos_t_3D) :: state_out

      state_in  = eos_t_3D(lo,hi)
      state_out = eos_t_3D(lo,hi)
      
      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! We do not want to integrate zones with rho < small_dens,
               ! but the hydro calls enforce_minimum_density, so we shouldn't
               ! be getting any zones like that. If so, throw an error.

               if (s_in(i,j,k,URHO) .lt. small_dens) then
                  print *,'... rho < small_dens in react_state: ', i, j, k, s_in(i,j,k,URHO)
                  call bl_error("Error:: React_3d.f90 :: ca_react_state")
               endif

               rhoInv                = ONE / s_in(i,j,k,URHO)

               state_in % rho(i,j,k)   = s_in(i,j,k,URHO)
               state_in % T(i,j,k)     = s_in(i,j,k,UTEMP)
               state_in % e(i,j,k)     = s_in(i,j,k,UEINT) * rhoInv
               state_in % xn(i,j,k,:)  = s_in(i,j,k,UFS:UFS+nspec-1) * rhoInv
               state_in % aux(i,j,k,:) = s_in(i,j,k,UFX:UFX+naux-1) * rhoInv

               ! The energy should never be negative coming into this call
               ! because of reset_internal_energy, so throw an error if it happens.
               
               if (allow_negative_energy .eq. 0 .and. state_in % e(i,j,k) .le. ZERO) then
                  print *,'... e negative in react_state: ', i, j, k, state_in % e(i,j,k)
                  call bl_error("Error:: React_3d.f90 :: ca_react_state")
               endif
            enddo
         enddo
      enddo

      call burner(state_in, state_out, dt_react, time)

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               ! Note that we want to update the total energy by taking the difference of the old
               ! rho*e and the new rho*e rather than adding rho*e to KE because we don't necessarily
               ! want to force here that rho * E = rho * e + rho * K where K is defined in terms of
               ! the momenta. If this is desired that can always be forced through reset_internal_energy.

               s_out(i,j,k,URHO)            = s_in(i,j,k,URHO)
               s_out(i,j,k,UEINT)           = s_out(i,j,k,URHO) * state_out % e(i,j,k)
               s_out(i,j,k,UEDEN)           = s_out(i,j,k,UEDEN) + s_out(i,j,k,URHO) * (state_out % e(i,j,k) - state_in % e(i,j,k))
               s_out(i,j,k,UFS:UFS+nspec-1) = s_out(i,j,k,URHO) * state_out % xn(i,j,k,:)
               s_out(i,j,k,UFX:UFX+naux -1) = state_in % aux(i,j,k,:)
               s_out(i,j,k,UTEMP)           = state_out % T(i,j,k)

               if (i.ge.r_l1 .and. i.le.r_h1 .and. j.ge.r_l2 .and. j.le.r_h2 .and. &
                   k.ge.r_l3 .and. k.le.r_h3) then
                  reaction_terms(i,j,k,1:nspec) = reaction_terms(i,j,k,1:nspec) &
                       + (state_out % xn(i,j,k,:) - state_in % xn(i,j,k,:))
                  reaction_terms(i,j,k,nspec+1) = reaction_terms(i,j,k,nspec+1) &
                       + (state_out % e(i,j,k) - state_in % e(i,j,k))
                  reaction_terms(i,j,k,nspec+2) = reaction_terms(i,j,k,nspec+2) &
                       + s_out(i,j,k,URHO) * (state_out % e(i,j,k) - state_in % e(i,j,k))
               endif

            enddo
         enddo
      enddo

  end subroutine ca_react_state
