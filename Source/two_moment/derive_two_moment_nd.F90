module derive_thornado_module

  use amrex_fort_module, only : rt => amrex_real
  use RadiationFieldsModule, only : nSpecies
  use ProgramHeaderModule, only : nE, nDOF

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
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icount, n_moments = 4

    print *,'NV IN CA_DER_J ',nv

    J_avg(:,:,:,:) = 0.0e0_rt  ! zero it out
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            icount = 0
            do is = 1, nSpecies
            do ie = 1, nE
            do id = 1, nDOF ! radiation degrees of freedom
              im = 1  ! J is first moment
              ii = (is-1)*(n_moments*nE*nDOF) + &
                   (im-1)*(nE*nDOF) + &
                   (ie-1)*nDOF + (id-1)
              icount = icount + 1  ! for averaging
              J_avg(i,j,k,:) = J_avg(i,j,k,:) + U_R(i,j,k,ii)
            end do
            end do
            end do
            J_avg(i,j,k,:) = J_avg(i,j,k,:) / icount  ! average

          end do
       end do
    end do

  end subroutine ca_der_J

  subroutine ca_der_J_per_E(J_avg,j_lo,j_hi,nv, &
                            U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                            domhi,delta,xlo,time,dt,bc,level,grid_no) &
                            bind(C, name="ca_der_J_per_E")

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: j_lo(3), j_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: J_avg(j_lo(1):j_hi(1),j_lo(2):j_hi(2),j_lo(3):j_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icount, n_moments = 4

    print *,'NV IN CA_DER_J_PER_E ',nv

    J_avg(:,:,:,:) = 0.0e0_rt  ! zero it out
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            icount = 0
            do is = 1, nSpecies
            do ie = 1, nE
            do id = 1, nDOF ! radiation degrees of freedom
              im = 1  ! J is first moment
              ii = (is-1)*(n_moments*nE*nDOF) + &
                   (im-1)*(nE*nDOF) + &
                   (ie-1)*nDOF + (id-1)
              icount = icount + 1  ! for averaging
              J_avg(i,j,k,ie) = J_avg(i,j,k,ie) + U_R(i,j,k,ii)
            end do
            end do
            end do
            J_avg(i,j,k,:) = J_avg(i,j,k,:) / icount  ! average

          end do
       end do
    end do

  end subroutine ca_der_J_per_E



  subroutine ca_der_Hx(Hx_avg,hx_lo,hx_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_Hx")

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: hx_lo(3), hx_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: Hx_avg(hx_lo(1):hx_hi(1),hx_lo(2):hx_hi(2),hx_lo(3):hx_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icount, n_moments = 4

    Hx_avg(:,:,:,:) = 0.0e0_rt  ! zero it out
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            icount = 0
            do is = 1, nSpecies
            do ie = 1, nE
            do id = 1, nDOF
              im = 2  ! Hx is second moment
              ii = (is-1)*(n_moments*nE*nDOF) + &
                   (im-1)*(nE*nDOF) + &
                   (ie-1)*nDOF + (id-1)
              icount = icount + 1  ! for averaging
              Hx_avg(i,j,k,:) = Hx_avg(i,j,k,:) + U_R(i,j,k,ii)
            end do
            end do
            end do
            Hx_avg(i,j,k,:) = Hx_avg(i,j,k,:) / icount  ! average

          end do
       end do
    end do

  end subroutine ca_der_Hx



  subroutine ca_der_Hy(Hy_avg,hy_lo,hy_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_Hy")

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: hy_lo(3), hy_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: Hy_avg(hy_lo(1):hy_hi(1),hy_lo(2):hy_hi(2),hy_lo(3):hy_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icount, n_moments = 4

    Hy_avg(:,:,:,:) = 0.0e0_rt  ! zero it out
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            icount = 0
            do is = 1, nSpecies
            do ie = 1, nE
            do id = 1, nDOF
              im = 3  ! Hy is third moment
              ii = (is-1)*(n_moments*nE*nDOF) + &
                   (im-1)*(nE*nDOF) + &
                   (ie-1)*nDOF + (id-1)
              icount = icount + 1  ! for averaging
              Hy_avg(i,j,k,:) = Hy_avg(i,j,k,:) + U_R(i,j,k,ii)
            end do
            end do
            end do
            Hy_avg(i,j,k,:) = Hy_avg(i,j,k,:) / icount  ! average

          end do
       end do
    end do

  end subroutine ca_der_Hy



  subroutine ca_der_Hz(Hz_avg,hz_lo,hz_hi,nv, &
                      U_R,d_lo,d_hi,nc,lo,hi,domlo, &
                      domhi,delta,xlo,time,dt,bc,level,grid_no) &
                      bind(C, name="ca_der_Hz")

    implicit none 

    integer, intent(in)     :: lo(3), hi(3)
    integer, intent(in)     :: hz_lo(3), hz_hi(3), nv
    integer, intent(in)     :: d_lo(3), d_hi(3), nc
    integer, intent(in)     :: domlo(3), domhi(3)
    integer, intent(in)     :: bc(3,2,nc)
    real(rt), intent(in)    :: delta(3), xlo(3), time, dt
    real(rt), intent(inout) :: Hz_avg(hz_lo(1):hz_hi(1),hz_lo(2):hz_hi(2),hz_lo(3):hz_hi(3),nv)
    real(rt), intent(in)    :: U_R(d_lo(1):d_hi(1),d_lo(2):d_hi(2),d_lo(3):d_hi(3),0:nc)
    integer, intent(in)     :: level, grid_no

    integer          :: i, j, k, is, ie, id, im, ii, icount, n_moments = 4

    Hz_avg(:,:,:,:) = 0.0e0_rt  ! zero it out
    do k = lo(3), hi(3)
       do j = lo(2), hi(2)
          do i = lo(1), hi(1)

            icount = 0
            do is = 1, nSpecies
            do ie = 1, nE
            do id = 1, nDOF
              im = 4  ! Hz is fourth moment
              ii = (is-1)*(n_moments*nE*nDOF) + &
                   (im-1)*(nE*nDOF) + &
                   (ie-1)*nDOF + (id-1)
              icount = icount + 1  ! for averaging
              Hz_avg(i,j,k,:) = Hz_avg(i,j,k,:) + U_R(i,j,k,ii)
            end do
            end do
            end do
            Hz_avg(i,j,k,:) = Hz_avg(i,j,k,:) / icount  ! average

          end do
       end do
    end do

  end subroutine ca_der_Hz

end module derive_thornado_module
