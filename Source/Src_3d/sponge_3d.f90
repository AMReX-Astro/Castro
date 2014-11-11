module sponge_module

  implicit none

  private

  public sponge

contains

  subroutine sponge(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3,lo,hi, &
       time,dt,dx,dy,dz,domlo,domhi)

    use meth_params_module, only : NVAR
    
    implicit none
    integer          :: lo(3),hi(3),domlo(3),domhi(3)
    integer          :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
    double precision :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,NVAR)
    double precision :: time, dt
    double precision :: dx, dy, dz

    !
    ! Note that subroutine sponge is called inside OpenMP parallel region.
    ! So you must make sure static (i.e., save) variables, if there are any
    ! in your implementation of this subroutine, are threadprivate.
    !
    
    ! Nothing happens in this generic version of the sponge routine.
    
  end subroutine sponge

end module sponge_module
