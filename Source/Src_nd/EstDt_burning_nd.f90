
    subroutine ca_estdt_burning(u, u_lo, u_hi, &
                                reactions, r_lo, r_hi, &
                                lo, hi, dx, dt)

      use bl_constants_module, only: ZERO, TWO
      use network, only: nspec
      use meth_params_module, only : NVAR, URHO, UEINT, burning_timestep_factor

      implicit none

      integer          :: u_lo(3), u_hi(3)
      integer          :: r_lo(3), r_hi(3)
      integer          :: lo(3), hi(3)
      double precision :: u(u_lo(1):u_hi(1),u_lo(2):u_hi(2),u_lo(3):u_hi(3),NVAR)
      double precision :: reactions(r_lo(1):r_hi(1),r_lo(2):r_hi(2),r_lo(3):r_hi(3),nspec+2)
      double precision :: dx(3), dt

      double precision :: e, delta_e
      integer          :: i, j, k

      ! The reactions MultiFab contains the net changes in X (the
      ! first nspec values), e (the nspec+1 value), and rho*e (the
      ! nspec+2) value. 
      !
      ! What we want to do is limit so that the timestep is equal to
      ! burning_timestep_factor * (e / delta(e)).  If the timestep
      ! factor is equal to 1, this says that we don't want the
      ! internal energy to change by any more than its current
      ! magnitude in the next timestep. 
      !
      ! Note that since the MultiFab will only have been reacted for
      ! dt / 2 this is actually off by a factor of 2 from that
      ! estimate, so we correct for that. If the timestep factor is
      ! less than one, it functionally controls the fraction we will
      ! allow the internal energy to change in this timestep due to
      ! nuclear burning, provide that the last timestep's burning is a
      ! good estimate for the current timestep's burning.

      do k = lo(3), hi(3)
         do j = lo(2), hi(2)
            do i = lo(1), hi(1)

               e = u(i,j,k,UEINT) / u(i,j,k,URHO)
               delta_e = TWO * reactions(i,j,k,nspec+1)

               if (abs(delta_e) > 1.d-100) then
                  dt = min(dt, burning_timestep_factor * e / abs(delta_e))
               endif

            enddo
         enddo
      enddo

    end subroutine ca_estdt_burning
