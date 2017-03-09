module MGutils_2D_module

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
    
    integer lo(2),hi(2),xlo(2),ylo(2),xhi(2),yhi(2)
    integer rl1, rl2, rh1, rh2
    integer ecxl1, ecxl2, ecxh1, ecxh2
    integer ecyl1, ecyl2, ecyh1, ecyh2
    integer coord_type
    real(rt)         rhs(rl1:rh1,rl2:rh2)
    real(rt)         ecx(ecxl1:ecxh1,ecxl2:ecxh2)
    real(rt)         ecy(ecyl1:ecyh1,ecyl2:ecyh2)
    real(rt)         dx(2)

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

    else 
       print *,'Bogus coord_type in apply_metric ' ,coord_type
       call bl_error("Error:: MGutils_2d.f90 :: ca_apply_metric")
    end if

  end subroutine ca_apply_metric



  subroutine ca_weight_cc(lo, hi, &
       cc, cl1, cl2, ch1, ch2,  &
       dx, coord_type) bind(C, name="ca_weight_cc")

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    integer lo(2),hi(2)
    integer cl1, cl2, ch1, ch2
    integer coord_type
    real(rt)         cc(cl1:ch1,cl2:ch2)
    real(rt)         dx(2)

    real(rt)         r
    integer i,j

    ! r-z
    if (coord_type .eq. 1) then

       ! At centers
       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          do j=lo(2),hi(2)
             cc(i,j) = cc(i,j) * r
          enddo
       enddo

    else 
       print *,'Bogus coord_type in weight_cc ' ,coord_type
       call bl_error("Error:: MGutils_2d.f90 :: ca_weight_cc")
    end if

  end subroutine ca_weight_cc



  subroutine ca_unweight_cc(lo, hi, &
       cc, cl1, cl2, ch1, ch2,  &
       dx, coord_type) bind(C, name="ca_unweight_cc")

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    integer lo(2),hi(2)
    integer cl1, cl2, ch1, ch2
    integer coord_type
    real(rt)         cc(cl1:ch1,cl2:ch2)
    real(rt)         dx(2)

    real(rt)         r
    integer i,j

    ! r-z
    if (coord_type .eq. 1) then

       ! At centers
       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          do j=lo(2),hi(2)
             cc(i,j) = cc(i,j) / r
          enddo
       enddo

    else 
       print *,'Bogus coord_type in unweight_cc ' ,coord_type
       call bl_error("Error:: MGutils_2d.f90 :: ca_unweight_cc")
    end if

  end subroutine ca_unweight_cc



  subroutine ca_unweight_edges(lo, hi, &
       ec, ecl1, ecl2, ech1, ech2, dx, coord_type, idir) &
       bind(C, name="ca_unweight_edges")

    use amrex_fort_module, only : rt => amrex_real
    implicit none
    
    integer lo(2),hi(2)
    integer ecl1, ecl2, ech1, ech2
    integer coord_type, idir
    real(rt)         ec(ecl1:ech1,ecl2:ech2)
    real(rt)         dx(2)

    real(rt)         :: r
    integer          :: i,j

    ! r-z
    if (coord_type .eq. 1) then

       if (idir .eq. 0) then
          ! On x-edges
          do i = lo(1), hi(1)
             if (i .ne. 0) then
                r = dble(i)*dx(1)
                do j = lo(2),hi(2)
                   ec(i,j) = ec(i,j) / r
                enddo
             end if
          enddo
       else
          ! On y-edges
          do i = lo(1), hi(1)
             r = (dble(i)+0.5e0_rt) * dx(1)
             do j = lo(2),hi(2)
                ec(i,j) = ec(i,j) / r
             enddo
          enddo
       end if

    else 
       print *,'Bogus coord_type in unweight_edges ' ,coord_type
       call bl_error("Error:: MGutils_2d.f90 :: ca_unweight_edges")
    end if

  end subroutine ca_unweight_edges

end module MGutils_2D_module
    
