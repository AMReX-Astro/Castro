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
