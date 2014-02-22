! Here we potentially reset the internal energy evolved separately
! from the total energy if we encounter a flow where we suspect that E
! - 0.5U^2 is not reliable.

subroutine reset_internal_e(u,u_l1,u_l2,u_h1,u_h2,lo,hi,verbose)

  use eos_module
  use network, only : nspec, naux
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT, UFS, UFX, &
                                 small_temp, allow_negative_energy
  use bl_constants_module

  implicit none

  integer          :: lo(2), hi(2), verbose
  integer          :: u_l1,u_l2,u_h1,u_h2
  double precision :: u(u_l1:u_h1,u_l2:u_h2,NVAR)
  
  ! Local variables
  integer          :: i,j
  integer          :: pt_index(2)
  double precision :: Up, Vp, ke
  double precision :: rho_eint, eint_new

  type (eos_t) :: eos_state
  
  ! Reset internal energy if negative.
  if (allow_negative_energy .eq. 0) then

     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           Up = u(i,j,UMX) / u(i,j,URHO)
           Vp = u(i,j,UMY) / u(i,j,URHO)

           ke = HALF * u(i,j,URHO) * (Up**2 + Vp**2)

           rho_eint = u(i,j,UEDEN) - ke

           if (rho_eint .gt. ZERO .and. rho_eint / u(i,j,UEDEN) .gt. 1.d-4) then
              ! Reset (rho e) if e is greater than 0.01% of E -- this
              ! is the conservative (and normal) method, i.e. we rely
              ! on the total E

               u(i,j,UEINT) = rho_eint

            else if (u(i,j,UEINT) .gt. ZERO) then
              ! If (e from E) < 0 or (e from E) < .0001*E but (e from e) > 0.
               
               u(i,j,UEDEN) = u(i,j,UEINT) + ke

           ! If not resetting and little e is negative ...
           else if (u(i,j,UEINT) .le. ZERO) then

              eos_state % rho = u(i,j,URHO)
              eos_state % T   = small_temp
              eos_state % xn  = u(i,j,UFS:UFS+nspec-1) / u(i,j,URHO)
              eos_state % aux = u(i,j,UFX:UFX+naux-1) / u(i,j,URHO)

              pt_index(1) = i
              pt_index(2) = j

              call eos(eos_input_rt, eos_state, pt_index = pt_index)

              eint_new = eos_state % e

              if (verbose .gt. 0) then
                 print *,'   '
                 print *,'>>> Warning: Castro_2d::reset_internal_energy  ',i,j 
                 print *,'>>> ... resetting neg. e from EOS using small_temp'
                 print *,'>>> ... from ',u(i,j,UEINT)/u(i,j,URHO),' to ', eint_new
                 print *,'    '
              end if

              u(i,j,UEINT) = u(i,j,URHO) * eint_new
              u(i,j,UEDEN) = u(i,j,URHO) * eint_new + ke

           end if

        enddo
     enddo

     ! If (allow_negative_energy .eq. 1) then just reset (rho e) from (rho E)
  else

     do j = lo(2),hi(2)
        do i = lo(1),hi(1)

           Up = u(i,j,UMX) / u(i,j,URHO)
           Vp = u(i,j,UMY) / u(i,j,URHO)
           ke = HALF * u(i,j,URHO) * (Up**2 + Vp**2)
           
           u(i,j,UEINT) = u(i,j,UEDEN) - ke
           
        enddo
     enddo
     
  endif

end subroutine reset_internal_e
