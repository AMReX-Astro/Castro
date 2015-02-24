
subroutine ca_test_type_lambda( lo, hi, &
     test, test_l1, test_l2, test_h1, test_h2, ntest, n, &
     lamx, lamx_l1, lamx_l2, lamx_h1, lamx_h2, &
     lamy, lamy_l1, lamy_l2, lamy_h1, lamy_h2)

  implicit none

  integer, intent(in) :: lo(2), hi(2), test_l1, test_l2, test_h1, test_h2, ntest, n, &
       lamx_l1, lamx_l2, lamx_h1, lamx_h2, &
       lamy_l1, lamy_l2, lamy_h1, lamy_h2
  double precision, intent(in) :: lamx(lamx_l1:lamx_h1,lamx_l2:lamx_h2)
  double precision, intent(in) :: lamy(lamy_l1:lamy_h1,lamy_l2:lamy_h2)
  double precision             :: test(test_l1:test_h1,test_l2:test_h2,ntest)
  integer :: i,j

  do j = lo(2), hi(2)
     do i = lo(1), hi(1)
        test(i,j,n) = (lamx(i,j) + lamx(i+1,j) &
             & +       lamy(i,j) + lamy(i,j+1) ) / 4.d0
     end do
  end do

end subroutine ca_test_type_lambda


subroutine ca_test_type_flux_lab( lo, hi, &
     flab, fl_l1, fl_l2, fl_h1, fl_h2, &
     fcom, fc_l1, fc_l2, fc_h1, fc_h2, &
     lamx, lx_l1, lx_l2, lx_h1, lx_h2, &
     lamy, ly_l1, ly_l2, ly_h1, ly_h2, &
     Er,   Er_l1, Er_l2, Er_h1, Er_h2, &
     S,     S_l1,  S_l2,  S_h1,  S_h2, &
     x,x_l1,x_h1, nt, idim)
  use meth_params_module, only : NVAR, URHO, UMX, UMY
  use fluxlimiter_module, only : Edd_factor
  implicit none
  integer,intent(in) :: lo(2), hi(2)
  integer,intent(in) :: fl_l1, fl_l2, fl_h1, fl_h2
  integer,intent(in) :: fc_l1, fc_l2, fc_h1, fc_h2
  integer,intent(in) :: lx_l1, lx_l2, lx_h1, lx_h2
  integer,intent(in) :: ly_l1, ly_l2, ly_h1, ly_h2
  integer,intent(in) :: Er_l1, Er_l2, Er_h1, Er_h2
  integer,intent(in) ::  S_l1,  S_l2,  S_h1,  S_h2  
  integer,intent(in) ::  x_l1, x_h1, nt, idim
  double precision           ::flab(fl_l1:fl_h1,fl_l2:fl_h2,0:nt-1)
  double precision,intent(in)::fcom(fc_l1:fc_h1,fc_l2:fc_h2) ! on face
  double precision,intent(in)::lamx(lx_l1:lx_h1,lx_l2:lx_h2) ! on face
  double precision,intent(in)::lamy(ly_l1:ly_h1,ly_l2:ly_h2) ! on face
  double precision,intent(in)::  Er(Er_l1:Er_h1,Er_l2:Er_h2)
  double precision,intent(in)::   S( S_l1: S_h1, S_l2: S_h2,NVAR)
  double precision,intent(in):: x(x_l1:x_h1)

  call bl_error("ca_test_type_flux_lab not implemented in 2d")
  
end subroutine ca_test_type_flux_lab
