! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine ca_compute_1d_gr_grav(var, r_l1, r_h1, grav, dx, problo)

      use probdata_module
      use meth_params_module, only : NVAR, URHO, UEINT, UTEMP, UFS
      use eos_module
      use network                     , only : nspec
      use fundamental_constants_module, only : Gconst
      use bl_constants_module         , only : M_PI

      implicit none

      integer         , intent(in   ) :: r_l1, r_h1
      double precision, intent(in   ) ::  var(r_l1:r_h1,NVAR)
      double precision, intent(  out) :: grav(r_l1:r_h1)
      double precision, intent(in   ) :: dx, problo(1)

      integer          :: i,n
      integer          :: pt_index(1)
      double precision :: rc,rlo,mass_encl,halfdx
      double precision :: R,e,G, P, C, T,dpdr, dpde, X(nspec)
      double precision :: ga, gb, gc

      double precision, parameter ::  fourpi       = 4.d0 * M_PI
      double precision, parameter ::  fourthirdspi = 4.d0 * M_PI / 3.d0
      double precision, parameter ::  sqvc         = 29979245800.d0**2

      halfdx = 0.5d0 * dx

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

         if (R .gt. 0.d0) then

            R  =  var(i,URHO)
            e  =  var(i,UEINT)/R
            T  =  var(i,UTEMP)
   
            do n = 1, nspec
               X(n)= var(i,UFS+n-1)/R
            enddo

            pt_index(1) = i
            call eos_given_ReX(G, P, C, T, dpdr, dpde, R, e, X, pt_index=pt_index)

           ga = (1.d0 + P/(R*sqvc))
           gb = (1.d0 + fourpi * rc**3 * P / (mass_encl*sqvc))
           gc =  1.d0 / (1.d0 - 2.d0 * Gconst * mass_encl / (rc*sqvc))

           grav(i) = grav(i)*ga*gb*gc

         end if

!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!       This ends the post-Newtonian correction
!!       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo

      if (problo(1) .eq. 0.d0) then
         do i = r_l1,-1
             grav(i) = -grav(-i-1)
         end do
      end if

      end subroutine ca_compute_1d_gr_grav

