
subroutine ca_test_type_flux_lab( lo, hi, &
     flab, fl_l1, fl_h1, &
     fcom, fc_l1, fc_h1, &
     lam,  lm_l1, lm_h1, &
     Er,   Er_l1, Er_h1, &
     S,     S_l1,  S_h1, &
     x,     x_l1,  x_h1, &
     nt, idim)
  use meth_params_module, only : NVAR, URHO, UMX
  use fluxlimiter_module, only : Edd_factor
  implicit none
  integer,intent(in) :: lo(1), hi(1)
  integer,intent(in) :: fl_l1, fl_h1
  integer,intent(in) :: fc_l1, fc_h1
  integer,intent(in) :: lm_l1, lm_h1
  integer,intent(in) :: Er_l1, Er_h1
  integer,intent(in) ::  S_l1,  S_h1  
  integer,intent(in) ::  x_l1,  x_h1
  integer,intent(in) :: nt, idim
  double precision           ::flab(fl_l1:fl_h1,0:nt-1)
  double precision,intent(in)::fcom(fc_l1:fc_h1) ! on face
  double precision,intent(in):: lam(lm_l1:lm_h1) ! on face
  double precision,intent(in)::  Er(Er_l1:Er_h1)
  double precision,intent(in)::   S( S_l1: S_h1,NVAR)
  double precision,intent(in)::   x( x_l1: x_h1)

  integer :: i
  double precision :: fr, lambda, eddf, vx

  do i=lo(1),hi(1)
     fr = 0.5d0 * (fcom(i)/(x(i)+1.d-50) + fcom(i+1)/x(i+1))
     lambda = 0.5d0 * (lam(i) + lam(i+1))
     eddf = Edd_factor(lambda)
     vx = S(i,UMX)/S(i,URHO)
     flab(i,idim) = fr + vx*(eddf+1.d0)*Er(i)
  end do
end subroutine ca_test_type_flux_lab
