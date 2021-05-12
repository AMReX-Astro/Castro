
subroutine ca_flux_face2center( lo, hi, &
     t, t_l1, t_l2, t_h1, t_h2, &
     f, f_l1, f_l2, f_h1, f_h2, &
     x, x_l1, x_h1, nt, idim, iflx) bind(C, name="ca_flux_face2center")

  use rad_params_module, only : ngroups
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(2), hi(2)
  integer,intent(in)::t_l1,t_h1,t_l2,t_h2
  integer,intent(in)::f_l1,f_h1,f_l2,f_h2
  integer,intent(in)::x_l1,x_h1
  integer,intent(in) :: nt, idim, iflx
  real(rt)                   ::t(t_l1:t_h1,t_l2:t_h2,0:nt-1)
  real(rt)        ,intent(in)::f(f_l1:f_h1,f_l2:f_h2)
  real(rt)        ,intent(in)::x(x_l1:x_h1)

  integer it, i, j

  it = idim*ngroups + iflx

  if (idim .eq. 0) then
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           t(i,j,it) = (f(i,j)/(x(i)+1.e-50_rt) + f(i+1,j)/x(i+1)) * 0.5e0_rt
        end do
     end do
  else 
     do j=lo(2), hi(2)
        do i=lo(1), hi(1)
           t(i,j,it) = (f(i,j)/x(i) + f(i,j+1)/x(i)) * 0.5e0_rt
        end do
     end do
  end if

end subroutine ca_flux_face2center
