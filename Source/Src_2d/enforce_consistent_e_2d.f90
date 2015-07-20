   subroutine ca_enforce_consistent_e(lo,hi,state,state_l1,state_l2,state_h1,state_h2)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UEDEN, UEINT
     use bl_constants_module

     implicit none

     integer          :: lo(2), hi(2)
     integer          :: state_l1,state_l2,state_h1,state_h2
     double precision :: state(state_l1:state_h1,state_l2:state_h2,NVAR)

     ! Local variables
     integer          :: i,j
     double precision :: u,v

     ! 
     ! Make sure to enforce (rho E) = (rho e) + 1/2 rho (u^2 + v^2)
     !
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)

           u = state(i,j,UMX) / state(i,j,URHO)
           v = state(i,j,UMY) / state(i,j,URHO)

           state(i,j,UEDEN) = state(i,j,UEINT) + &
               HALF * state(i,j,URHO) * (u*u + v*v)

        end do
     end do

   end subroutine ca_enforce_consistent_e

