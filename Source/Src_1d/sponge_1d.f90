     
      subroutine sponge( &
            uout, uout_l1 ,uout_h1, &
            lo,hi,dt)

      use meth_params_module, only : NVAR

      implicit none
      integer          :: lo(1), hi(1)
      integer          ::  uout_l1, uout_h1
      double precision :: uout(uout_l1:uout_h1,NVAR)
      double precision :: dt

      ! Nothing happens in this generic version of the sponge routine.

      end subroutine sponge

