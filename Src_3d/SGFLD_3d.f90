
subroutine ca_test_type_lambda( lo, hi, &
     test, test_l1, test_l2, test_l3, test_h1, test_h2, test_h3, &
     ntest, n, &
     lamx, lamx_l1, lamx_l2, lamx_l3, lamx_h1, lamx_h2, lamx_h3, &
     lamy, lamy_l1, lamy_l2, lamy_l3, lamy_h1, lamy_h2, lamy_h3, &
     lamz, lamz_l1, lamz_l2, lamz_l3, lamz_h1, lamz_h2, lamz_h3 )

  implicit none

  integer, intent(in) :: lo(3), hi(3), ntest, n, &
       test_l1, test_l2, test_l3, test_h1, test_h2, test_h3, &
       lamx_l1, lamx_l2, lamx_l3, lamx_h1, lamx_h2, lamx_h3, &
       lamy_l1, lamy_l2, lamy_l3, lamy_h1, lamy_h2, lamy_h3, &
       lamz_l1, lamz_l2, lamz_l3, lamz_h1, lamz_h2, lamz_h3 
  double precision, intent(in) :: lamx(lamx_l1:lamx_h1,lamx_l2:lamx_h2,lamx_l3:lamx_h3)
  double precision, intent(in) :: lamy(lamy_l1:lamy_h1,lamy_l2:lamy_h2,lamy_l3:lamy_h3)
  double precision, intent(in) :: lamz(lamz_l1:lamz_h1,lamz_l2:lamz_h2,lamz_l3:lamz_h3)
  double precision             :: test(test_l1:test_h1,test_l2:test_h2,test_l3:test_h3,ntest)
  integer :: i,j,k

  do k = lo(3), hi(3)
     do j = lo(2), hi(2)
        do i = lo(1), hi(1)
           test(i,j,k,n) = (lamx(i,j,k) + lamx(i+1,j,k) &
                &   +       lamy(i,j,k) + lamy(i,j+1,k) &
                &   +       lamz(i,j,k) + lamz(i,j,k+1) ) / 6.d0
        end do
     end do
  end do

end subroutine ca_test_type_lambda


subroutine ca_test_type_flux_lab( lo, hi, &
     flab, fl_l1, fl_l2, fl_l3, fl_h1, fl_h2, fl_h3, &
     fcom, fc_l1, fc_l2, fc_l3, fc_h1, fc_h2, fc_h3, &
     lamx, lx_l1, lx_l2, lx_l3, lx_h1, lx_h2, lx_h3, &
     lamy, ly_l1, ly_l2, ly_l3, ly_h1, ly_h2, ly_h3, &
     lamz, lz_l1, lz_l2, lz_l3, lz_h1, lz_h2, lz_h3, &
     Er,   Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3, &
     S,     S_l1,  S_l2,  S_l3,  S_h1,  S_h2,  S_h3, &
     x,x_l1,x_h1, nt, idim)
  use meth_params_module, only : NVAR, URHO, UMX, UMY, UMZ
  use fluxlimiter_module, only : Edd_factor
  implicit none
  integer, intent(in) :: lo(3), hi(3)
  integer,intent(in) :: fl_l1, fl_l2, fl_l3, fl_h1, fl_h2, fl_h3
  integer,intent(in) :: fc_l1, fc_l2, fc_l3, fc_h1, fc_h2, fc_h3
  integer,intent(in) :: lx_l1, lx_l2, lx_l3, lx_h1, lx_h2, lx_h3
  integer,intent(in) :: ly_l1, ly_l2, ly_l3, ly_h1, ly_h2, ly_h3
  integer,intent(in) :: lz_l1, lz_l2, lz_l3, lz_h1, lz_h2, lz_h3
  integer,intent(in) :: Er_l1, Er_l2, Er_l3, Er_h1, Er_h2, Er_h3
  integer,intent(in) ::  S_l1,  S_l2,  S_l3,  S_h1,  S_h2,  S_h3  
  integer,intent(in) ::  x_l1, x_h1, nt, idim
  double precision           ::flab(fl_l1:fl_h1,fl_l2:fl_h2,fl_l3:fl_h3,0:nt-1)
  double precision,intent(in)::fcom(fc_l1:fc_h1,fc_l2:fc_h2,fc_l3:fc_h3) ! on face
  double precision,intent(in)::lamx(lx_l1:lx_h1,lx_l2:lx_h2,lx_l3:lx_h3) ! on face
  double precision,intent(in)::lamy(ly_l1:ly_h1,ly_l2:ly_h2,ly_l3:ly_h3) ! on face
  double precision,intent(in)::lamz(lz_l1:lz_h1,lz_l2:lz_h2,lz_l3:lz_h3) ! on face
  double precision,intent(in)::  Er(Er_l1:Er_h1,Er_l2:Er_h2,Er_l3:Er_h3)
  double precision,intent(in)::   S( S_l1: S_h1, S_l2: S_h2, S_l3: S_h3,NVAR)
  double precision,intent(in):: x(x_l1:x_h1)

  call bl_error("ca_test_type_flux_lab not implemented in 3d")
  
end subroutine ca_test_type_flux_lab
