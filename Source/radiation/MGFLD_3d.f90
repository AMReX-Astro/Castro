
subroutine ca_flux_face2center( lo, hi, &
     t, t_l1, t_l2, t_l3, t_h1, t_h2, t_h3, &
     f, f_l1, f_l2, f_l3, f_h1, f_h2, f_h3, &
     x, x_l1, x_h1, nt, idim, iflx) bind(C, name="ca_flux_face2center")

  use rad_params_module, only : ngroups
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  integer,intent(in):: lo(3), hi(3)
  integer,intent(in)::t_l1,t_h1,t_l2,t_h2,t_l3,t_h3
  integer,intent(in)::f_l1,f_h1,f_l2,f_h2,f_l3,f_h3
  integer,intent(in)::x_l1,x_h1
  integer,intent(in) :: nt, idim, iflx
  real(rt)                   ::t(t_l1:t_h1,t_l2:t_h2,t_l3:t_h3,0:nt-1)
  real(rt)        ,intent(in)::f(f_l1:f_h1,f_l2:f_h2,f_l3:f_h3)
  real(rt)        ,intent(in)::x(x_l1:x_h1)

  integer it, i, j, k

  it = idim*ngroups + iflx

  if (idim .eq. 0) then
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              t(i,j,k,it) = (f(i,j,k) + f(i+1,j,k)) * 0.5e0_rt
           end do
        end do
     end do
  else if (idim .eq. 1) then
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              t(i,j,k,it) = (f(i,j,k) + f(i,j+1,k)) * 0.5e0_rt
           end do
        end do
     end do
  else
     do k=lo(3), hi(3)
        do j=lo(2), hi(2)
           do i=lo(1), hi(1)
              t(i,j,k,it) = (f(i,j,k) + f(i,j,k+1)) * 0.5e0_rt
           end do
        end do
     end do
  end if

end subroutine ca_flux_face2center
