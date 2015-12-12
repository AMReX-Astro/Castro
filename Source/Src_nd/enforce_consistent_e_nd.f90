   subroutine ca_enforce_consistent_e(lo,hi,state,s_lo,s_hi,dx)

     use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ, UEDEN, UEINT
     use bl_constants_module

     implicit none

     integer          :: lo(3), hi(3)
     integer          :: s_lo(3), s_hi(3)
     double precision :: state(s_lo(1):s_hi(1),s_lo(2):s_hi(2),s_lo(3):s_hi(3),NVAR)
     double precision :: dx(3)
     
     ! Local variables
     integer          :: i,j,k
     double precision :: u, v, w, rhoInv
     
     ! 
     ! Enforces (rho E) = (rho e) + 1/2 rho (u^2 +_ v^2 + w^2)
     !
     do k = lo(3), hi(3)
        do j = lo(2), hi(2)
           do i = lo(1), hi(1)
              
              rhoInv = ONE/state(i,j,k,URHO)

              u = state(i,j,k,UMX) * rhoInv
              v = state(i,j,k,UMY) * rhoInv
              w = state(i,j,k,UMZ) * rhoInv

              state(i,j,k,UEDEN) = state(i,j,k,UEINT) + &
                     HALF * state(i,j,k,URHO) * (u*u + v*v + w*w)

           end do
        end do
     end do

   end subroutine ca_enforce_consistent_e

