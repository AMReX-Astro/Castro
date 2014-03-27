! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_compute_1d_gr_grav(var, r_l1, r_h1, grav, dx, problo)

      use probdata_module
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS, UFX
      use eos_module
      use network                     , only : nspec, naux
      use fundamental_constants_module, only : Gconst, c_light
      use bl_constants_module         

      implicit none

      integer         , intent(in   ) :: r_l1, r_h1
      double precision, intent(in   ) ::  var(r_l1:r_h1,NVAR)
      double precision, intent(  out) :: grav(r_l1:r_h1)
      double precision, intent(in   ) :: dx, problo(1)

      integer          :: i
      integer          :: pt_index(1)
      double precision :: rc,rlo,mass_encl,halfdx
      double precision :: ga, gb, gc

      double precision, parameter ::  fourpi       = FOUR * M_PI
      double precision, parameter ::  fourthirdspi = FOUR3RD * M_PI
      double precision, parameter ::  sqvc         = c_light**2

      type (eos_t) :: eos_state

      halfdx = HALF * dx

      do i = 0,r_h1
         rlo = problo(1) + dble(i) * dx
         rc  = rlo + halfdx
         if (i.gt.0) then
            mass_encl = mass_encl + fourthirdspi * halfdx * (rlo**2 + rlo*(rlo-halfdx) + (rlo-halfdx)**2) * var(i-1,URHO) + &
                                    fourthirdspi * halfdx * ( rc**2 +  rc* rlo         +  rlo**2        ) * var(i  ,URHO)
         else
            mass_encl = fourthirdspi * halfdx * (rc**2 + rc*rlo +  rlo**2) * var(i,URHO)
         end if

         grav(i) = -Gconst * mass_encl / rc**2

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This adds the post-Newtonian correction
!!       Added by Ken Chen, 6/9 2010
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!       Tolman-Oppenheimer-Volkoff(TOV) case

         if (var(i,URHO) .gt. ZERO) then

            eos_state % rho = var(i,URHO)
            eos_state % e   = var(i,UEINT) / var(i,URHO)
            eos_state % T   = var(i,UTEMP)
            eos_state % xn  = var(i,UFS:UFS+nspec-1) / var(i,URHO)
            eos_state % aux = var(i,UFX:UFX+naux-1) / var(i,URHO)

            pt_index(1) = i

            call eos(eos_input_re, eos_state, pt_index = pt_index)

            ga = (ONE + eos_state % p /(eos_state % rho * sqvc))
            gb = (ONE + fourpi * rc**3 * eos_state % p / (mass_encl*sqvc))
            gc =  ONE / (ONE - TWO * Gconst * mass_encl / (rc*sqvc))

            grav(i) = grav(i)*ga*gb*gc

         end if

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This ends the post-Newtonian correction
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

      if (problo(1) .eq. ZERO) then
         do i = r_l1,-1
             grav(i) = -grav(-i-1)
         end do
      end if

      end subroutine ca_compute_1d_gr_grav

