module MGutils_1D_module

  use amrex_error_module
  use amrex_fort_module, only : rt => amrex_real
  implicit none

  public

contains

  subroutine ca_apply_metric(lo, hi, &
       xlo, xhi, &
       rhs, rl1, rh1,  &
       ecx, ecxl1, ecxh1, dx, coord_type) bind(C, name="ca_apply_metric")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(1),hi(1),xlo(1),xhi(1)
    integer, intent(in) :: rl1, rh1
    integer, intent(in) :: ecxl1, ecxh1
    integer, intent(in) :: coord_type
    real(rt), intent(inout) :: rhs(rl1:rh1)
    real(rt), intent(inout) :: ecx(ecxl1:ecxh1)
    real(rt), intent(in) :: dx(1)

    real(rt)         r,rlo,rhi
    integer i

    ! r-z
    if (coord_type .eq. 1) then

       ! At centers
       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          rhs(i) = rhs(i) * r
       enddo

       ! On edges
       do i=xlo(1),xhi(1)
          r = dble(i)*dx(1)
          ecx(i) = ecx(i) * r
       enddo

       ! spherical
    else if (coord_type .eq. 2) then

       ! At centers
       do i=lo(1),hi(1)
          !           r = (dble(i)+0.5e0_rt) * dx(1)
          !           rhs(i) = rhs(i) * r**2
          rlo = dble(i) * dx(1)
          rhi = rlo + dx(1)
          rhs(i) = rhs(i) * (rhi**3 - rlo**3) / (3.e0_rt * dx(1))
       enddo

       ! Note that (rhi**3 - rlo**3) / (3 dr) = ( (r + dr/2)**3 - (r - dr/2)**3 ) / (3 dr)
       !                                 = r^2 + dr^2 / 12

       ! On edges
       do i=xlo(1),xhi(1)
          r = dble(i)*dx(1)
          ecx(i) = ecx(i) * r**2
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in apply_metric ' ,coord_type
       call amrex_error("Error:: MGutils_1d.f90 :: ca_apply_metric")
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

    j = lo(2)
    k = lo(3)

    ! r-z
    if (coord_type .eq. 1) then

       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          cc(i,j,k) = cc(i,j,k) * r
       enddo

       ! spherical
    else if (coord_type .eq. 2) then

       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          cc(i,j,k) = cc(i,j,k) * r**2
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in weight_cc ' ,coord_type
       call amrex_error("Error:: MGutils_1d.f90 :: ca_weight_cc")
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

    j = lo(2)
    k = lo(3)

    ! r-z
    if (coord_type .eq. 1) then

       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          cc(i,j,k) = cc(i,j,k) / r
       enddo

       ! spherical
    else if (coord_type .eq. 2) then

       do i=lo(1),hi(1)
          r = (dble(i)+0.5e0_rt) * dx(1)
          cc(i,j,k) = cc(i,j,k) / r**2
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in unweight_cc ' ,coord_type
       call amrex_error("Error:: MGutils_1d.f90 :: ca_unweight_cc")
#endif
    end if

  end subroutine ca_unweight_cc



  subroutine ca_unweight_edges(lo, hi, &
       ecx, ecx_lo, ecx_hi, dx, coord_type, idir) &
       bind(C, name="ca_unweight_edges")

    use amrex_fort_module, only : rt => amrex_real
    implicit none

    integer, intent(in) :: lo(3),hi(3)
    integer, intent(in) :: ecx_lo(3), ecx_hi(3)
    integer, value, intent(in) :: coord_type, idir
    real(rt), intent(inout) :: ecx(ecx_lo(1):ecx_hi(1),ecx_lo(2):ecx_hi(2),ecx_lo(3):ecx_hi(3))
    real(rt), intent(in) :: dx(3)

    real(rt)         r
    integer i,j,k

    !$gpu

    j = lo(2)
    k = lo(3)

    ! r-z
    if (coord_type .eq. 1) then

       ! On edges
       do i=lo(1),hi(1)
          r = abs(dble(i))*dx(1)
          if (i.ne.0) ecx(i,j,k) = ecx(i,j,k) / r
       enddo

       ! spherical
    else if (coord_type .eq. 2) then

       ! On edges
       do i=lo(1),hi(1)
          r = dble(i)*dx(1)
          if (i.ne.0) ecx(i,j,k) = ecx(i,j,k) / r**2
       enddo

#ifndef AMREX_USE_CUDA
    else
       print *,'Bogus coord_type in unweight_edges ' ,coord_type
       call amrex_error("Error:: MGutils_1d.f90 :: ca_unweight_edges")
#endif
    end if

  end subroutine ca_unweight_edges

end module MGutils_1D_module
