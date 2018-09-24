module MGutils_2D_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_apply_metric(lo, hi, &
       xlo, ylo, &
       xhi, yhi, &
       rhs, rl1, rl2, rh1, rh2,  &
       ecx, ecxl1, ecxl2, ecxh1, ecxh2, &
       ecy, ecyl1, ecyl2, ecyh1, ecyh2, dx, coord_type) &
       bind(C, name="ca_apply_metric")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3),hi(3),xlo(3),ylo(3),xhi(3),yhi(3)
    integer, intent(in) ::  rl1, rl2, rh1, rh2
    integer, intent(in) :: ecxl1, ecxl2, ecxh1, ecxh2
    integer, intent(in) :: ecyl1, ecyl2, ecyh1, ecyh2
    integer, intent(in) :: coord_type
    real(rt), intent(inout) :: rhs(rl1:rh1,rl2:rh2)
    real(rt)         ecx(ecxl1:ecxh1,ecxl2:ecxh2)
    real(rt)         ecy(ecyl1:ecyh1,ecyl2:ecyh2)
    real(rt)         dx(3)

    real(rt)         r
    integer i,j

    ! r-z
    if (coord_type .eq. 1) then

       ! At centers
       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          do j=lo(2),hi(2)
             rhs(i,j) = rhs(i,j) * r
          enddo
       enddo

       ! On x-edges
       do i=xlo(1),xhi(1)
          r = dble(i)*dx(1)
          do j=xlo(2),xhi(2)
             ecx(i,j) = ecx(i,j) * r
          enddo
       enddo

       ! On y-edges
       do i=ylo(1),yhi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          do j=ylo(2),yhi(2)
             ecy(i,j) = ecy(i,j) * r
          enddo
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in apply_metric ' ,coord_type
       call amrex_error("Error:: MGutils_2d.f90 :: ca_apply_metric")
#endif
    end if

  end subroutine ca_apply_metric



  subroutine ca_weight_cc(lo, hi, &
       cc, c_lo, c_hi,  &
       dx, coord_type) bind(C, name="ca_weight_cc")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) ::c_lo(3), c_hi(3)
    integer, value, intent(in) ::coord_type
    real(rt), intent(inout) ::cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_lo(3))
    real(rt), intent(in) :: dx(3)

    real(rt)         r
    integer i,j,k

    !$gpu

    ! r-z
    if (coord_type .eq. 1) then

       ! At centers
       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          do j=lo(2),hi(2)
             do k = lo(3), hi(3)
                cc(i,j,k) = cc(i,j,k) * r
             enddo
          enddo
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in weight_cc ' ,coord_type
       call amrex_error("Error:: MGutils_2d.f90 :: ca_weight_cc")
#endif
    end if

  end subroutine ca_weight_cc



  subroutine ca_unweight_cc(lo, hi, &
       cc, c_lo, c_hi,  &
       dx, coord_type) bind(C, name="ca_unweight_cc")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) ::c_lo(3), c_hi(3)
    integer, value, intent(in) ::coord_type
    real(rt), intent(inout) ::cc(c_lo(1):c_hi(1),c_lo(2):c_hi(2),c_lo(3):c_lo(3))
    real(rt), intent(in) :: dx(3)

    real(rt)         r
    integer i,j,k

    !$gpu

    ! r-z
    if (coord_type .eq. 1) then

       ! At centers
       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          do j=lo(2),hi(2)
             do k = lo(3),hi(3)
                cc(i,j,k) = cc(i,j,k) / r
             enddo
          enddo
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in unweight_cc ' ,coord_type
       call amrex_error("Error:: MGutils_2d.f90 :: ca_unweight_cc")
#endif
    end if

  end subroutine ca_unweight_cc



  subroutine ca_unweight_edges(lo, hi, &
       ec, ec_lo, ec_hi, dx, coord_type, idir) &
       bind(C, name="ca_unweight_edges")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) ::  lo(3),hi(3)
    integer, intent(in) :: ec_lo(3), ec_hi(3)
    integer, value, intent(in) :: coord_type, idir
    real(rt), intent(inout) :: ec(ec_lo(1):ec_hi(1),ec_lo(2):ec_hi(2),ec_lo(3):ec_hi(3))
    real(rt), intent(in) :: dx(3)

    real(rt)         :: r
    integer          :: i,j,k

    !$gpu

    ! r-z
    if (coord_type .eq. 1) then

       if (idir .eq. 0) then
          ! On x-edges
          do i = lo(1), hi(1)
             if (i .ne. 0) then
                r = dble(i)*dx(1)
                do j = lo(2),hi(2)
                   do k = lo(3),hi(3)
                      ec(i,j,k) = ec(i,j,k) / r
                   enddo
                enddo
             end if
          enddo
       else if (idir .eq. 1) then
          ! On y-edges
          do i = lo(1), hi(1)
             r = (dble(i)+0.5e0_rt) * dx(1)
             do j = lo(2),hi(2)
                do k = lo(3),hi(3)
                   ec(i,j,k) = ec(i,j,k) / r
                enddo
             enddo
          enddo
          ! On z-edges
          do i = lo(1), hi(1)
             r = (dble(i)+0.5e0_rt) * dx(1)
             do j = lo(2),hi(2)
                do k = lo(3),hi(3)
                   ec(i,j,k) = ec(i,j,k) / r
                enddo
             enddo
          enddo
       else

       end if

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in unweight_edges ' ,coord_type
       call amrex_error("Error:: MGutils_2d.f90 :: ca_unweight_edges")
#endif
    end if

  end subroutine ca_unweight_edges

end module MGutils_2D_module
