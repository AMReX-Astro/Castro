   subroutine ca_enforce_consistent_e(lo,hi,state,state_l1,state_h1)

     use meth_params_module, only : NVAR, URHO, UMX, UEDEN, UEINT

     implicit none

     integer          :: lo(1), hi(1)
     integer          :: state_l1,state_h1
     double precision :: state(state_l1:state_h1,NVAR)

     ! Local variables
     integer          :: i
     double precision :: u

     ! 
     ! Make sure to enforce (rho E) = (rho e) + 1/2 rho (u^2)
     !
     do i = lo(1), hi(1)

        u = state(i,UMX) / state(i,URHO)

        state(i,UEDEN) = state(i,UEINT) + &
            0.5d0 * state(i,URHO) * (u*u)

     end do

   end subroutine ca_enforce_consistent_e

