module MGutils_2D_module

  use castro_error_module, only: castro_error
  use amrex_fort_module, only: rt => amrex_real

  implicit none

  public

contains

  subroutine ca_apply_metric(lo, hi, &
                             xlo, xhi, &
#if AMREX_SPACEDIM >= 2
                             ylo, yhi, &
#endif
                             rhs, rlo, rhi, &
                             ecx, ecxlo, ecxhi, &
#if AMREX_SPACEDIM >= 2
                             ecy, ecylo, ecyhi, &
#endif
                             dx, coord_type) &
                             bind(C, name="ca_apply_metric")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: xlo(3), xhi(3)
#if AMREX_SPACEDIM >= 2
    integer,  intent(in   ) :: ylo(3), yhi(3)
#endif
    integer,  intent(in   ) :: rlo(3), rhi(3)
    integer,  intent(in   ) :: ecxlo(3), ecxhi(3)
#if AMREX_SPACEDIM >= 2
    integer,  intent(in   ) :: ecylo(3), ecyhi(3)
#endif
    real(rt), intent(inout) :: rhs(rlo(1):rhi(1),rlo(2):rhi(2),rlo(3):rhi(3))
    real(rt), intent(inout) :: ecx(ecxlo(1):ecxhi(1),ecxlo(2):ecxhi(2),ecxlo(3):ecxhi(3))
#if AMREX_SPACEDIM >= 2
    real(rt), intent(inout) :: ecy(ecylo(1):ecyhi(1),ecylo(2):ecyhi(2),ecylo(3):ecyhi(3))
#endif
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: coord_type

    real(rt) :: r
    integer  :: i, j, k

    !$gpu

    ! r-z
    if (coord_type == 1) then

       ! At centers
       do i = lo(1), hi(1)
          r = (dble(i) + 0.5e0_rt) * dx(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                rhs(i,j,k) = rhs(i,j,k) * r
             end do
          end do
       end do

       ! On x-edges
       do i = xlo(1), xhi(1)
          r = dble(i) * dx(1)
          do k = xlo(3), xhi(3)
             do j = xlo(2), xhi(2)
                ecx(i,j,k) = ecx(i,j,k) * r
             end do
          end do
       end do

#if AMREX_SPACEDIM >= 2
       ! On y-edges
       do i = ylo(1), yhi(1)
          r = (dble(i) + 0.5e0_rt) * dx(1)
          do k = ylo(3), yhi(3)
             do j = ylo(2), yhi(2)
                ecy(i,j,k) = ecy(i,j,k) * r
             end do
          end do
       end do
#endif

#ifndef AMREX_USE_CUDA
    else

       print *,'Bogus coord_type in apply_metric ' ,coord_type
       call castro_error("Error:: MGutils_2d.f90 :: ca_apply_metric")
#endif

    end if

  end subroutine ca_apply_metric



  subroutine ca_weight_cc(lo, hi, &
                          cc, clo, chi, &
                          dx, coord_type) bind(C, name="ca_weight_cc")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: clo(3), chi(3)
    real(rt), intent(inout) :: cc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: coord_type

    real(rt) :: r
    integer  :: i, j, k

    !$gpu

    ! r-z
    if (coord_type == 1) then

       ! At centers
       do i = lo(1), hi(1)
          r = (dble(i) + 0.5e0_rt) * dx(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                cc(i,j,k) = cc(i,j,k) * r
             end do
          end do
       end do

#ifndef AMREX_USE_CUDA
    else

       print *,'Bogus coord_type in weight_cc ' ,coord_type
       call castro_error("Error:: MGutils_2d.f90 :: ca_weight_cc")
#endif

    end if

  end subroutine ca_weight_cc



  subroutine ca_unweight_cc(lo, hi, &
                            cc, clo, chi, &
                            dx, coord_type) bind(C, name="ca_unweight_cc")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: clo(3), chi(3)
    real(rt), intent(inout) :: cc(clo(1):chi(1),clo(2):chi(2),clo(3):chi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: coord_type

    real(rt) :: r
    integer  :: i, j, k

    !$gpu

    ! r-z
    if (coord_type == 1) then

       ! At centers
       do i = lo(1), hi(1)
          r = (dble(i) + 0.5e0_rt) * dx(1)
          do k = lo(3), hi(3)
             do j = lo(2), hi(2)
                cc(i,j,k) = cc(i,j,k) / r
             end do
          end do
       end do

#ifndef AMREX_USE_CUDA
    else

       print *,'Bogus coord_type in unweight_cc ' ,coord_type
       call castro_error("Error:: MGutils_2d.f90 :: ca_unweight_cc")
#endif

    end if

  end subroutine ca_unweight_cc



  subroutine ca_unweight_edges(lo, hi, &
                               ec, eclo, echi, &
                               dx, coord_type, idir) &
                               bind(C, name="ca_unweight_edges")

    use amrex_fort_module, only: rt => amrex_real

    implicit none

    integer,  intent(in   ) :: lo(3), hi(3)
    integer,  intent(in   ) :: eclo(3), echi(3)
    real(rt), intent(inout) :: ec(eclo(1):echi(1),eclo(2):echi(2),eclo(3):echi(3))
    real(rt), intent(in   ) :: dx(3)
    integer,  intent(in   ), value :: coord_type, idir

    real(rt) :: r
    integer  :: i, j, k

    !$gpu

    ! r-z
    if (coord_type == 1) then

       if (idir == 0) then

          ! On x-edges
          do i = lo(1), hi(1)
             if (i /= 0) then
                r = dble(i) * dx(1)
                do k = lo(3), hi(3)
                   do j = lo(2), hi(2)
                      ec(i,j,k) = ec(i,j,k) / r
                   end do
                end do
             end if
          end do

       else

          ! On y-edges
          do i = lo(1), hi(1)
             r = (dble(i) + 0.5e0_rt) * dx(1)
             do k = lo(3), hi(3)
                do j = lo(2), hi(2)
                   ec(i,j,k) = ec(i,j,k) / r
                end do
             end do
          end do
       end if

#ifndef AMREX_USE_CUDA
    else

       print *,'Bogus coord_type in unweight_edges ' ,coord_type
       call castro_error("Error:: MGutils_2d.f90 :: ca_unweight_edges")
#endif

    end if

  end subroutine ca_unweight_edges

end module MGutils_2D_module

