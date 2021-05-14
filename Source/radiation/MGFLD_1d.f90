
subroutine ca_flux_face2center( lo, hi, &
     t, t_l1, t_h1, &
     f, f_l1, f_h1, &
     x, x_l1, x_h1, &
     nt, idim, it) bind(C, name="ca_flux_face2center")

  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(1), hi(1)
  integer,intent(in)::t_l1,t_h1
  integer,intent(in)::f_l1,f_h1
  integer,intent(in)::x_l1,x_h1
  integer,intent(in) :: nt, idim, it
  real(rt)                   ::t(t_l1:t_h1,0:nt-1)
  real(rt)        ,intent(in)::f(f_l1:f_h1)
  real(rt)        ,intent(in)::x(x_l1:x_h1)

  integer i

  do i=lo(1),hi(1)
     t(i,it) = (f(i)/(x(i)+1.e-50_rt) + f(i+1)/x(i+1)) * 0.5e0_rt
  end do

end subroutine ca_flux_face2center


