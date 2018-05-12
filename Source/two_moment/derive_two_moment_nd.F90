module derive_thornado_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  public

contains
  
! All subroutines in this file must be threadsafe because they are called
! inside OpenMP parallel regions.

  subroutine ca_der_J(J_avg,j_lo,j_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_J")

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: j_lo(3), j_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: J_avg(j_lo(1):j_hi(1),j_lo(2):j_hi(2),j_lo(3):j_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, ii

    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)
!            do is = 1, n_species
!            do im = 1, n_moments
!            do ie = 1, n_energy
!            do id = 1, n_rad_dof
!               U_R(i,j,k,ii) = uCR(id,ie,i,j,k,im,is)  ! This is just here to see the ordering
!               ii = (is-1)*(n_moments*n_energy*n_rad_dof) + &
!                    (im-1)*(n_energy*n_rad_dof) + &
!                    (ie-1)*n_rad_dof + (id-1)
!               J_avg(i,j,k,1) = some average over some number of the components
!            end do
!            end do
!            end do
!            end do
          end do
       end do
    end do

  end subroutine ca_der_J

end module derive_thornado_module
