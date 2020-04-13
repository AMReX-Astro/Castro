void
Castro::apply_metric(const Box& bx,
                     const Box& xbx,
#if AMREX_SPACEDIM >= 2
                     const Box& ybx,
#endif
                     Array4<Real> const rhs, const Box& rbx,
                     Array4<Real> const ecx,
#if AMREX_SPACEDIM >= 2
                     Array4<Real> const ecy
#endif
                     )
{



  int coord_type = geom.Coord();

  // r-z
  if (coord_type == 1) {

    AMREX_PARALLEL_FOR_3D(bx, i, j, k,
    {

     IntVect idx(D_DECL(i, j, k));

     // at centers
     if (rbx.contains(idx)) {
       Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
       rhs(i,j,k) *= r;
     }

     // On x-edges
     if (xbx.contains(idx)) {
       Real r = static_cast<Real>(i) * dx[0];
       ecx(i,j,k) *= r;
     }

#if AMREX_SPACEDIM >= 2
     // On y-edges
     if (ybx.contains(idx)) {
       Real r = (static_cast<Real>(i) + 0.5_rt) * dx[0];
       ecy(i,j,k) *= r;
     }
#endif

    });

#ifndef AMREX_USE_CUDA
  } else {

    amrex::Print() << "Bogus coord_type in apply_metric " << coord_type << std::endl;

    amrex::Error("Error:: MGutils.cpp :: ca_apply_metric");
#endif

  }
}


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

