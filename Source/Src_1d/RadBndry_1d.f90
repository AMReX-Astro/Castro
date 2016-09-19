! THIS IS ONLY A TEMPLATE!!!

! If radiation.[lo|hi]_bcflag is set to 1, this subroutine will be used to set
! boundary conditions for radiation diffusion.
!
! For LO_DIRICHLET, this should set Dirichlet value of radiation energy density.
! For LO_NEUMANN, this should set inward radiation flux.
! For LO_MARSHAK & LO_SANCHEZ_POMRANING, this should set incident radiation flux.

subroutine rbndry(  &
     bf, b_l1, b_h1, &
     &   d_l1, d_h1, &
     dx, xlo, t, dir, face) 

  use rad_params_module, only : ngroups

  implicit none

  integer, intent(in) :: b_l1, b_h1
  integer, intent(in) :: d_l1, d_h1  ! computational domain index
  integer, intent(in) :: dir         ! always 0 (i.e., x-direction) in 1D
  integer, intent(in) :: face        ! 0: low     1: high
  double precision, intent(out) :: bf(b_l1:b_h1,0:ngroups-1)
  double precision, intent(in) :: dx(1), xlo(1), t

  bf = 0.5d0

end subroutine rbndry
